% clc
clear
rng(0)

ct = 200;
Pt = 10:5:30;
% Pt = 10;

pt = (10^-3)*db2pow(Pt);

height = 3;
R = 1.9;
M = 8;
N = 8;
noise = -90; %dBm
no = (10^-3)*db2pow(noise);

f = 28e9; % 28 GHz
c = 3e8; % Speed of light
lambda = c / f; % Free space wavelength
neff = 1.4; % Low-index materials like Teflon
lambda_g = lambda / neff;
eta = (c / (4 * pi * f))^2;
D2 = 10;
D_leng = 10;
P_max = 1; % Maximum power constraint
sigma2 = 1e-9; % Noise power

% Initialize power allocation
stepp = 2;
power_temp = flip(1:stepp:stepp*(M-1)+1);
power = power_temp/sum(power_temp);
power_eq  = (1/M)*ones(M,1);             % equal power (for conventional NOMA baseline)
% power_eq = power;

Rate_pro_avg = zeros(length(Pt), M);
Rate_fixpower_avg = zeros(length(Pt), M);
Rate_oma_avg = zeros(length(Pt), M);       % Oma with pinching
Rate_conv_avg = zeros(length(Pt), M);       % Oma with conventional
Rate_noma_avg = zeros(length(Pt), M);       % noma with conventional

for mi = 1:length(Pt) 
        mi;
        pt = (10^-3)*db2pow(Pt(mi));
        
        % Conventional antenna
        r_all_pin = zeros(M, ct);
        r_all_mul = zeros(M, ct);
        r_all_OMA = zeros(M, ct);
        r_all_mul_min = zeros(M, ct);
        r_oma_conv = zeros(M, ct);
        r_noma_conv  = zeros(M, ct);

       
       for i = 1:ct
           i
            %first general all five users' locations
        Mmax = M;
        locx = zeros(Mmax,2);
        locx(:,1) = D_leng*rand(Mmax,1)-D_leng/2; %length        
        locx(:,2) = D2*rand(Mmax,1)-D2/2; %width,  
        locx(Mmax,1) = -D2 + locx(Mmax,1); %best user, shift to the left 
        for m = Mmax-1:-1 : 1
            locx(m,:) = (Mmax-m)*D_leng + locx(m,:);%shift to the right and bottom
        end
        loc = [locx(1:M-1,:); locx(Mmax,:)];

    
        % ============== Conventional (no pinching) ================
            % Free-space pathloss to a single conventional antenna at origin
        dist_conv = max(1, sqrt(loc(:,1).^2 + loc(:,2).^2 + height^2));   % Mx1
        g_conv  = eta ./ (dist_conv.^2);                 % channel power gains
    
        % OMA baseline (conventional) ---
        % Each user gets 1/M of the time with full power pt
        r_oma_conv(:, i) = (1/M) * log2(1 + (g_conv * pt) / no);
    
        %  NOMA baseline (conventional, with SIC) ---
        % Order users by channel gain (weak -> strong)
        [g_sorted, order_ws] = sort(g_conv, 'ascend'); % weak first
          for m = 1:M
                rtemp = zeros(M-m+1,1);
                for mj = m:M
                    rtemp(mj-m+1) = log2(1 + g_sorted(mj) * power(mj) * pt / M  / ...
                        (no + g_sorted(mj) * sum(power(mj+1:M)) * pt / M ));
                end
                r_noma_conv(m,i) =  min(rtemp);
          end
        % Map back to original user indices (optional: you can keep sum only)
        % r_noma_eq_conv(:, i) = r_tmp;

            % Single Pinching antenna
            dist = max(1, sqrt(loc(:,2).^2 + height^2));
            r_all_pin(:,i) = log2(1 + eta * pt ./ dist.^2/no) / N;      
            
            %multiple pinching antennas         
        for mt = 1 : M %for each user, there are M distances to the M pinching antennas
            dist_temp = max(1,sqrt((loc(mt,1)-loc(:,1)).^2 + (loc(mt,2)-0).^2 +height^2));
            theta_temp = 2*pi* (loc(:,1)+100)/lambda_g; %just assume 100 as the location of the feed
            theta_channel = 2*pi*dist_temp/lambda;
            h(mt) =  abs(sum(1./dist_temp.*exp(-complex(0,1)*(theta_channel+theta_temp))))^2;
        end
    c = h .* eta * pt / N ;  % Effective channel gain (element-wise)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Total Power Minimization
       gamma_target = 0.01; % Example fixed SINR threshold (e.g., 3 dB = 10^(3/10) ≈ 2)

        cvx_begin gp quiet
            variables p(M)
            expressions interference(M, M)
             minimize(sum(p)) % Objective: minimize total power

            for m = 1:M
                for k = m:M
                    if k < M
                        interference(m, k) = sum(p(k+1:M)) + no / c(m);
                    else
                        interference(m, k) = no / c(m);
                    end

                    % SINR constraint for user m decoding user k's signal
                    interference(m, k) <= p(k) / gamma_target;
                end
            end

        cvx_end
       power_min = p./sum(p);

      %% %% Compute proposed form power allocation
            for m = 1:M
                rtemp = zeros(M-m+1,1);
                for mj = m:M
                    rtemp(mj-m+1) = log2(1 + h(mj) * power_min(mj) * eta * pt / N / ...
                        (no + h(mj) * sum(power_min(mj+1:M)) * eta * pt / N));
                end
                r_all_mul_min(m,i) = min(rtemp);
            end
      %% %% Compute fixed power allocation
            for m = 1:M
                rtemp = zeros(M-m+1,1);
                for mj = m:M
                    rtemp(mj-m+1) = log2(1 + h(mj) * power(mj) * eta * pt / N / ...
                        (no + h(mj) * sum(power(mj+1:M)) * eta * pt / N));
                end
                r_all_mul(m,i) = min(rtemp);
            end
      %% %% Compute OMA rate with PA
            for m = 1:M
                 r_all_OMA(m,i) = (1/M)*log2(1 + h(m) * eta * pt / N /no);
            end
         
      end
    
        %% Average Sum Rate results
        r_ave_sum(mi) = mean(sum(r_oma_conv,1));  % OMA conventional
        r_ave_sum_omaPA(mi) = mean(sum(r_all_OMA,1));  % OMA Pinching Antenna
        r_ave_pin_sum(mi) = mean(sum(r_all_pin,1));  % Single Pinching Antenna
        r_ave_mul_sum(mi) = mean(sum(r_all_mul,1));   % Fixed power PA NOMA
        r_ave_mul_sum_min(mi) = mean(sum(r_all_mul_min,1));  % Proposed Optimized PA NOMA
        r_ave_mul_sum_noma(mi) = mean(sum(r_noma_conv,1)); % NOMA conventional

        %% === Jain's Fairness Index ===
        Jain_oma_conv(mi)  = mean( sum(r_oma_conv,1).^2 ./ (M * sum(r_oma_conv.^2,1)) );
        Jain_oma_pin(mi)   = mean( sum(r_all_OMA,1).^2 ./ (M * sum(r_all_OMA.^2,1)) );
        Jain_single_pin(mi)   = mean( sum(r_all_pin,1).^2 ./ (M * sum(r_all_pin.^2,1)) );
        Jain_noma_fixed(mi) = mean( sum(r_all_mul,1).^2 ./ (M * sum(r_all_mul.^2,1)) );
        Jain_noma_prop(mi)  = mean( sum(r_all_mul_min,1).^2 ./ (M * sum(r_all_mul_min.^2,1)) );
        Jain_noma_conv(mi)   = mean( sum(r_noma_conv,1).^2 ./ (M * sum(r_noma_conv.^2,1)) );
end
Pt = 10:5:30;   % <-- replace with your actual x-axis variable if different

% === Colors and Markers for Methods ===
colors  = [0.47 0.67 0.19; ...        % OMA conventional  
           0.3 0.75 0.93; ...  % OMA-PA            
           1 0 1; ...      % Single PA         
           0 0 1; ...        % FPA-NOMA          
           1 0 0; ...        % Proposed OPA-NOMA 
           0.49 0.18 0.56];       % C-NOMA            

markers = {'p','s','d','^','o','x'};

% === Plot ===
figure(1); hold on; grid on; box on;

% M = 4 (solid lines)
plot(Pt, r_ave_sum,         '-', 'Color', colors(1,:), 'Marker', markers{1}, 'LineWidth', 1.8, 'MarkerSize', 7);
plot(Pt, r_ave_sum_omaPA,   '-', 'Color', colors(2,:), 'Marker', markers{2}, 'LineWidth', 1.8, 'MarkerSize', 7);
plot(Pt, r_ave_pin_sum,     '-', 'Color', colors(3,:), 'Marker', markers{3}, 'LineWidth', 1.8, 'MarkerSize', 7);
plot(Pt, r_ave_mul_sum,     '-', 'Color', colors(4,:), 'Marker', markers{4}, 'LineWidth', 1.8, 'MarkerSize', 7);
plot(Pt, r_ave_mul_sum_min, '-', 'Color', colors(5,:), 'Marker', markers{5}, 'LineWidth', 1.8, 'MarkerSize', 7);
plot(Pt, r_ave_mul_sum_noma,'-', 'Color', colors(6,:), 'Marker', markers{6}, 'LineWidth', 1.8, 'MarkerSize', 7);
% === Labels ===
xlabel('Transmit Power (dBm)', 'Interpreter','latex','FontSize', 14);
ylabel('Sum Rate (bps/Hz)', 'Interpreter','latex','FontSize', 14);

legend({'C-OMA', 'PA-OMA', 'SPA', 'FPA-NOMA', 'OPA-NOMA', 'C-NOMA', ...
     '$M=4$', '$M=8$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, ...
    'Box', 'off', 'Location','best', 'NumColumns', 2);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);




