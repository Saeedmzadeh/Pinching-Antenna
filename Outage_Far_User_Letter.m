% clc
clear
rng(0)

ct = 300;
% Pt = 10:5:30;
Pt = 10;

pt = (10^-3)*db2pow(Pt);

height = 3;
R = 1.9;
M = 5;
noise = -90; %dBm
no = (10^-3)*db2pow(noise);

f = 28e9; % 28 GHz
c = 3e8; % Speed of light
lambda = c / f; % Free space wavelength
neff = 1.4; % Low-index materials like Teflon
lambda_g = lambda / neff;
eta = (c / (4 * pi * f))^2;
% D = 10;
D2 = 10;
D1 = 10;
D_leng = D1;
P_max = 1; % Maximum power constraint
sigma2 = 1e-9; % Noise power

% Initialize power allocation
stepp = 2;
power_temp = flip(1:stepp:stepp*(M-1)+1);
power = power_temp/sum(power_temp);
% Define target rates (bps/Hz) for outage threshold
rate_thresh = 0.1:0.2:7; % example target rates from 0.1 to 3 bps/Hz

outage_prob_far_fix = zeros(length(rate_thresh),1);
outage_prob_far_Proposed = zeros(length(rate_thresh),1);
outage_prob_Conv = zeros(length(rate_thresh),1);
outage_prob_Single = zeros(length(rate_thresh),1);
outage_prob_NOMA = zeros(length(rate_thresh),1);
outage_prob_OMA = zeros(length(rate_thresh),1);
outage_prob_OMA_PA = zeros(length(rate_thresh),1);


for idx = 1:length(rate_thresh)
    idx;
    r_threshold = rate_thresh(idx);
    
    outage_count_fix = 0;
    outage_count_proposed = 0;
    outage_count_Conv = 0;
    outage_count_Single = 0;
    outage_count_NOMA = 0;
    outage_count_OMA = 0;
    outage_count_OMA_PA = 0;

        
       for i = 1:ct
               i;
           %first general all five users' locations
            Mmax = 5;
            locx = zeros(Mmax,2);
            locx(:,1) = D_leng*rand(Mmax,1)-D_leng/2; %length        
            locx(:,2) = D2*rand(Mmax,1)-D2/2; %width,  
            locx(Mmax,1) = -D2 + locx(Mmax,1); %best user, shift to the left 
            for m = Mmax-1:-1 : 1
                locx(m,:) = (Mmax-m)*D1 + locx(m,:);%shift to the right and bottom
            end
            loc = [locx(1:M-1,:); locx(Mmax,:)];
    
        
                % Conventional antenna
                r_all = zeros(M,1);
                dist = max(1, sqrt(loc(:,1).^2 + loc(:,2).^2 + height^2));
                r_all = log2(1 + eta * pt ./ dist.^2/no) / M;      
                
                % Pinching antenna
                r_all_pin = zeros(M,1);
                dist = max(1, sqrt(loc(:,2).^2 + height^2));
                r_all_pin = log2(1 + eta * pt ./ dist.^2/no) / M;     
               
                 % Free-space pathloss to a single conventional antenna at origin
                dist_conv = max(1, sqrt(loc(:,1).^2 + loc(:,2).^2 + height^2));   % Mx1
                g_conv    = eta ./ (dist_conv.^2);                                 % channel power gains
            
                % --- OMA baseline (conventional) ---
                % Each user gets 1/M of the time with full power pt
                r_oma_conv = zeros(M,1);
                r_oma_conv = (1/M) * log2(1 + (g_conv * pt) / no);
            
                % --- NOMA baseline (conventional, equal power with SIC) ---
                % Order users by channel gain (weak -> strong)
                [g_sorted, order_ws] = sort(g_conv, 'ascend'); % weak first

              
                %multiple pinching antennas 
                h = zeros(1, M);
            for mt = 1 : M %for each user, there are M distances to the M pinching antennas
                dist_temp = max(1,sqrt((loc(mt,1)-loc(:,1)).^2 + (loc(mt,2)-0).^2 +height^2));
                theta_temp = 2*pi* (loc(:,1)+100)/lambda_g; %just assume 100 as the location of the feed
                theta_channel = 2*pi*dist_temp/lambda;
                h(mt) =  abs(sum(1./dist_temp.*exp(-complex(0,1)*(theta_channel+theta_temp))))^2;
            end
        c = h .* eta * pt / M ;  % Effective channel gain (element-wise)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Total Power Minimization
           gamma_target = 0.1; % Example fixed SINR threshold (e.g., 3 dB = 10^(3/10) ≈ 2)
    
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
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %% %% Compute proposed form power allocation
                
            % Calculate achievable rates with fixed power allocation
            r_user = zeros(M,1);
            for m = 1:M
                rtemp = zeros(M-m+1,1);
                for mj = m:M
                    interference = no + h(mj) * sum(power(mj+1:end)) * eta * pt / M;
                    signal = h(mj) * power(mj) * eta * pt / M;
                    rtemp(mj - m + 1) = log2(1 + signal / interference);
                end
                r_user(m) = min(rtemp); % user m's rate
            end
            %  Proposed Method
            r_user_proposed = zeros(M,1);
                for m = 1:M
                    rtemp = zeros(M-m+1,1);
                    for mj = m:M
                        rtemp(mj-m+1) = log2(1 + h(mj) * power_min(mj) * eta * pt / M / ...
                            (no + h(mj) * sum(power_min(mj+1:M)) * eta * pt / M));
                    end
                    r_user_proposed(m) = min(rtemp);
                end
            %   Conventional NOMA
          r_noma = zeros(M,1);
          for m = 1:M
                rtemp = zeros(M-m+1,1);
                for mj = m:M
                    rtemp(mj-m+1) = log2(1 + g_sorted(mj) * power(mj) * pt / M / ...
                        (no + g_sorted(mj) * sum(power(mj+1:M)) * pt / M));
                end
                r_noma(m) =  min(rtemp);
          end
          % Compute OMA rate with PA
            r_oma_PA = zeros(M,1);
            for m = 1:M
                 r_oma_PA(m,1) = (1/M)*log2(1 + h(m) * eta * pt / M /no);
            end



%%%%%%% Check outage for the far user (assumed to be user M)  %%%%%
                if r_user(M) < r_threshold
                    outage_count_fix = outage_count_fix + 1;
                end
                if r_user_proposed(M) < r_threshold
                    outage_count_proposed = outage_count_proposed + 1;
                end
                if r_all(M) < r_threshold
                    outage_count_Conv = outage_count_Conv + 1;
                end
                if r_all_pin(M) < r_threshold
                    outage_count_Single = outage_count_Single + 1;
                end
                if r_noma(M) < r_threshold
                    outage_count_NOMA = outage_count_NOMA + 1;
                end
                if r_oma_conv(M) < r_threshold
                    outage_count_OMA = outage_count_OMA + 1;
                end
                if r_oma_PA(M) < r_threshold
                    outage_count_OMA_PA = outage_count_OMA_PA + 1;
                end
            
     
       end
    outage_prob_far_fix(idx) = outage_count_fix / ct;
    outage_prob_far_Proposed(idx) = outage_count_proposed / ct;  
    outage_prob_Conv(idx) = outage_count_Conv / ct;
    outage_prob_Single(idx) = outage_count_Single / ct;
    outage_prob_NOMA(idx) = outage_count_NOMA / ct;
    outage_prob_OMA(idx ) = outage_count_OMA / ct;
    outage_prob_OMA_PA(idx ) = outage_count_OMA_PA / ct;
    
end
filename = sprintf('results2');
save(filename);


figure(2)
p = plot(rate_thresh, outage_prob_far_Proposed, '-or', ...
         rate_thresh, outage_prob_far_fix, '-^b', ...
         rate_thresh, outage_prob_Single, '-*m', ...
         rate_thresh, outage_prob_NOMA, '-x', ...
         rate_thresh, outage_prob_OMA, '-p', ...
         rate_thresh, outage_prob_OMA_PA, '-s', ...
         'LineWidth', 2);

% Set marker sizes for each line
set(p(1), 'MarkerSize', 8);
set(p(2), 'MarkerSize', 8);
set(p(3), 'MarkerSize', 8);
set(p(4), 'MarkerSize', 8);
set(p(5), 'MarkerSize', 8);
set(p(6), 'MarkerSize', 8);

% Axis labels with font settings
xlabel('Target Rate (bits/s/Hz)', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Outage Probability (Far User)', 'Interpreter', 'latex', 'FontSize', 15);

% Legend with font settings
legend({'PA-OPA-NOMA', 'PA-FPA-NOMA', 'SPA', 'C-NOMA','C-OMA','PA-OMA'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

grid on

% Apply font to axis ticks
set(gca, 'Interpreter', 'latex', 'FontSize', 12);
