%% 1.---Simulate trajectories.
% % If you have done it already, comment for further analysis below
clear all
iM=1;
imin=1;
imax=1;%The number of simulations (After reaching the steady state)
w_cold=120;w_hot=240;
%T_c=120;
%n_c=1./(exp(w_c./T_c)-1);
n_h=10;T_h=w_hot/(log((n_h+1)/n_h));
n_c=1e-5;T_c=0;
sub_folder_name='Data';
mkdir(sub_folder_name)
for ur=0:0
    if ur==0
        i1=1;%don't touch, this is called in Factorisation;
        Factorisation;
        myVars = {"p1","p2","p3","p1_2","p2_2","p3_2","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
            'x_m_vec','p_m_vec','p1_vec','p2_vec','p1_2_vec','p2_2_vec','p3_2_vec','na_vec','t_vec_i1','w_hot','w_cold','w_cav','n_h','n_c','w_m','f','g'...
            ,'k','g_h','g_c','g_m','dt','J_h','J_m','J_cold','J_cav'};
        save([sub_folder_name,'/unconditional'],myVars{:});
        % plot(x_m_vec,1i*p_m_vec,'LineWidth',2);
        % xlim([-40 40])
        % ylim([-50 50])
    else
        for i1=imin:imax
            Factorisation;
            tvec_dN1=jump_times;
            myVars2={"tvec_dN1","p1","p2","p3","p1_2","p2_2","p3_2","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
                'x_m_vec','p_m_vec','p1_vec','p2_vec','p1_2_vec','p2_2_vec','p3_2_vec','na_vec','t_vec_i1','w_hot','w_cold','w_cav',...
                'n_h','n_c','w_m','f','g','k','g_h','g_c','g_m','dt','Q_h','Q_h_f',...
                'J_h','J_m','J_cold','J_cav'};
            save([sub_folder_name,'/conditional_traj',num2str(i1)],myVars2{:});
            [i1,imax]
        end
    end
end
%%%
iM=2;%This is the number of atoms, here it is two
% See multiple copies folder if you wanna change it. It shows up in some dependent codes
%%%
%% 2.---Entropy production and TUR analysis.
% Load unconditional data for heat current (needed for entropy production)
imax=1;
sub_folder_name='Data';
myVars_uncond = {'J_h','w_hot','w_cold','n_h','n_c','w_m'};
load([sub_folder_name,'/unconditional'],myVars_uncond{:});
J_hot = J_h;  % Store heat current from unconditional evolution
% Calculate temperatures
T_h = w_hot/(log((n_h+1)/n_h));
T_c = w_cold/(log((n_c+1)/n_c));
sigma = J_hot * (1/T_h - 1/T_c);  % Entropy production rate

% Now calculate filtered tick statistics
det_filt = 1;  % Enable detector filter
plot_filter = 0;
dtj = [];
muvec = zeros(1,imax);
varvec = zeros(1,imax);

for i1 = 1:imax
    myVars = {"tvec_dN1", 'w_m', 'w_hot', 'w_cold', 'w_cav', 'n_c'};
    load([sub_folder_name,'/conditional_traj',num2str(i1)], myVars{:})

    % Apply detector filter
    if det_filt == 1
        Detector_Filter_saturation;
        tvec_dN1 = tvec_dN1_I2(1:end);
    end

    % Normalize tick times
    tvec_dN1 = tvec_dN1 * w_m / pi;

    % Calculate inter-tick intervals
    dtjump = diff([0, tvec_dN1]);
    dtj = [dtj, dtjump];
    muvec(1,i1) = mean(dtjump(2:end));
    varvec(1,i1) = std(dtjump(2:end))^2;
end

% Remove first element and calculate statistics
dtj = dtj(2:end);
tau = mean(muvec, 'omitnan');  % Mean inter-tick time
var_ = mean(varvec, 'omitnan');  % Variance
N = tau^2 / var_;  % Accuracy

% Display results
fprintf('\n=== Results (with detector filter) ===\n');
fprintf('Mean inter-tick time (tau): %.4f (Omega_m/pi units)\n', tau);
fprintf('Accuracy (N): %.2f\n', N);
fprintf('Entropy production rate (sigma): %.6e\n', sigma);
fprintf('TUR quantity (2*N/(tau*sigma)): %.4f (should be >= 1)\n', 2*N/(tau*sigma));
fprintf('=====================================\n\n');