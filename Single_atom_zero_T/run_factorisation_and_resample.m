%% 1.---Simulate trajectories.
% % If you have done it already, comment for further analysis below
clear all
iM=1;
imin=1;
imax=10;%The number of simulations (After reaching the steady state)
w_cold=120;w_hot=240;
%T_c=120;
%n_c=1./(exp(w_c./T_c)-1);
n_h=10;T_h=w_hot/(log((n_h+1)/n_h));
n_c=0;T_c=0;
sub_folder_name='Data';
mkdir(sub_folder_name)
for ur=0:0
    if ur==0
        i1=1;%don't touch, this is called in Factorisation; 
        Factorisation;
        myVars = {"p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
            'x_m_vec','p_m_vec','p1_vec','p2_vec','na_vec','t_vec_i1','w_hot','w_cold','w_cav','n_h','n_c','w_m','f','g'...
            ,'k','g_h','g_c','g_m','dt'};
        save([sub_folder_name,'/unconditional'],myVars{:});
        % plot(x_m_vec,1i*p_m_vec,'LineWidth',2);
        % xlim([-40 40])
        % ylim([-50 50])
    else
        for i1=imin:imax
            Factorisation;
            tvec_dN1=jump_times;
            myVars2={"tvec_dN1","p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","x_m","p_m", ...
                'x_m_vec','p_m_vec','p1_vec','p2_vec','na_vec','t_vec_i1','w_hot','w_cold','w_cav',...
                'n_h','n_c','w_m','f','g','k','g_h','g_c','g_m','dt'};
            save([sub_folder_name,'/conditional_traj',num2str(i1)],myVars2{:});
            [i1,imax]
        end
    end
end
%%%
iM=1;%This is the number of atoms, here it is one, cannot be changed. 
% See multiple copies folder if you wanna change it. It shows up in some dependent codes
%%%
%% 2.---Analyze the data: conditional populations/phase space, and photon number

sub_folder_name='Data';
mkdir('Data/Pics')

% Load conditional trajectory
i1 = 5; % Use the First trajectory for this analysis
myVars2 = {'tvec_dN1','x_m_vec','p_m_vec','p1_vec','p2_vec','na_vec','t_vec_i1','w_m'};
load([sub_folder_name,'/conditional_traj',num2str(i1)],myVars2{:})

% Compute normalized time
t_norm = t_vec_i1 * w_m / pi;

% Select time window: 40.5 to 43.5
idx = (t_norm >= 40.5) & (t_norm <= 43.5);
t_plot = t_norm(idx);
x_m_plot = x_m_vec(idx);
p_m_plot = p_m_vec(idx);
p1_plot = p1_vec(idx);
p2_plot = p2_vec(idx);
p3_plot = 1 - p1_vec(idx) - p2_vec(idx);
na_plot = na_vec(idx);

% Filter tick times within the time window
tick_norm = tvec_dN1 * w_m / pi;
tick_idx = (tick_norm >= 40.5) & (tick_norm <= 43.5);
ticks_plot = tick_norm(tick_idx);

% Create figure with 3 subplots
figure('Position', [100, 100, 800, 600])

% (a) Mechanical oscillator quadratures
subplot(3,1,1)
hold on
plot(t_plot, x_m_plot, 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]) % Blue
plot(t_plot, p_m_plot, 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]) % Orange
% Add vertical lines for ticks
for i = 1:length(ticks_plot)
    xline(ticks_plot(i), 'k-', 'LineWidth', 1);
end
xlim([40.5, 43.5])
ylim([-30, 30])
xticks(41:43)
yticks(-30:10:30)
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('', 'FontSize', 24)
legend('$\langle \hat{x}_m \rangle$', '$\langle \hat{p}_m \rangle$', 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 24)
set(gca, 'LineWidth', 1, 'FontSize', 24, 'FontName', 'Helvetica')
set(gca, 'XGrid', 'off', 'YGrid', 'off')
box on
hold off

% (b) Populations
subplot(3,1,2)
hold on
plot(t_plot, p1_plot, 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]) % Blue
plot(t_plot, p2_plot, 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]) % Orange
plot(t_plot, p3_plot, 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]) % Purple
% Add vertical lines for ticks
for i = 1:length(ticks_plot)
    xline(ticks_plot(i), 'k-', 'LineWidth', 1);
end
xlim([40.5, 43.5])
ylim([-0.05, 1.05])
xticks(41:43)
yticks(0:0.25:1)
yticklabels({'0.00', '0.25', '0.50', '0.75', '1.00'})
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('', 'FontSize', 24)
legend('$\langle \hat{p}_1 \rangle$', '$\langle \hat{p}_2 \rangle$', '$\langle \hat{p}_3 \rangle$', 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 24)
set(gca, 'LineWidth', 1, 'FontSize', 24, 'FontName', 'Helvetica')
set(gca, 'XGrid', 'off', 'YGrid', 'off')
box on
hold off

% (c) Cavity photon number
subplot(3,1,3)
hold on
plot(t_plot, na_plot, 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]) % Blue
% Add vertical lines for ticks
for i = 1:length(ticks_plot)
    xline(ticks_plot(i), 'k-', 'LineWidth', 1);
end
xlim([40.5, 43.5])
ylim([-0.05, 1.05])
xticks(41:43)
yticks(0:0.25:1)
yticklabels({'0.00', '0.25', '0.50', '0.75', '1.00'})
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('', 'FontSize', 24)
legend('$\langle \hat{a}^\dagger \hat{a} \rangle$', 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 24)
set(gca, 'LineWidth', 1, 'FontSize', 24, 'FontName', 'Helvetica')
set(gca, 'XGrid', 'off', 'YGrid', 'off')
box on
hold off

% % Save figure
% saveas(gcf, [pwd '/Data/Pics/conditional_traj_ticks.png'])
% saveas(gcf, [pwd '/Data/Pics/conditional_traj_ticks.fig'])

%% 3.--- Tick Statistics: Inter-tick intervals (with and without filter)

sub_folder_name='Data';
imin = 1;
iM=1;
imax = 100; % Number of trajectories to analyze (increase as more data becomes available)

% Initialize arrays to store all inter-tick intervals
dtj_nofilter = [];  % Without filter
dtj_filter = [];    % With filter
muvec_nofilter = zeros(1,imax);
varvec_nofilter = zeros(1,imax);
muvec_filter = zeros(1,imax);
varvec_filter = zeros(1,imax);

fprintf('Analyzing tick statistics for %d trajectories...\n', imax);

for i1 = imin:imax
    % Load tick times and required parameters for filter
    myVars = {"tvec_dN1", 'w_m', 'w_hot', 'w_cold', 'w_cav', 'n_c'};
    load([sub_folder_name,'/conditional_traj',num2str(i1)], myVars{:})
    
    % --- WITHOUT FILTER ---
    % Normalize tick times for analysis
    tvec_dN1_norm = tvec_dN1 * w_m / pi;
    dtjump = diff([0, tvec_dN1_norm]);
    dtj_nofilter = [dtj_nofilter, dtjump];
    muvec_nofilter(i1) = mean(dtjump);
    varvec_nofilter(i1) = std(dtjump)^2;
    
    % --- WITH FILTER (apply detector dead time) ---
    % Reload original unnormalized tick times for filter
    plot_filter = 1;
    Detector_Filter_saturation;
    if exist('tvec_dN1_I2', 'var')
        tvec_dN1_filtered = tvec_dN1_I2*(w_m/pi); % Normalize filtered tick times
        dtjump_filter = diff([0, tvec_dN1_filtered]);
        %dtjump_filter = dtjump_filter(dtjump_filter < 2);
        dtj_filter = [dtj_filter, dtjump_filter];
        muvec_filter(i1) = mean(dtjump_filter);
        varvec_filter(i1) = std(dtjump_filter)^2;
    end
    
    fprintf('  Trajectory %d/%d processed\n', i1, imax);
end

% Calculate statistics
mu_nofilter = mean(muvec_nofilter);
var_nofilter = mean(varvec_nofilter);
N_nofilter = mu_nofilter^2 / var_nofilter;  % Accuracy
nu_nofilter = 1 / mu_nofilter;  % Tick frequency

mu_filter = mean(muvec_filter);
var_filter = mean(varvec_filter);
N_filter = mu_filter^2 / var_filter;  % Accuracy
nu_filter = 1 / mu_filter;  % Tick frequency

fprintf('\n--- Statistics Summary ---\n');
fprintf('WITHOUT FILTER: N = %.1f, ν = %.2f Ω_m/π\n', N_nofilter, nu_nofilter);
fprintf('WITH FILTER:    N = %.1f, ν = %.2f Ω_m/π\n', N_filter, nu_filter);

% Create histogram plots (separate figures)
bin = 200;

% Figure 1: Without filter
figure('Position', [100, 100, 700, 500])
histogram(dtj_nofilter(2:end), bin, 'Normalization', 'probability', 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'none')
xline([0 1 2 3 4], '-k', 'FontSize', 2, 'HandleVisibility', 'off');
xlim([-0.1, 3])
xlabel('$\Omega_{\rm m}\tau/\pi$', 'Interpreter', 'latex', 'FontSize', 24)
text(0.05, 0.95, sprintf('$\\mathcal{N} = %.1f$\n$\\nu = %.2f\\Omega_{\\rm m}/\\pi$', ...
    N_nofilter, nu_nofilter), 'Units', 'normalized', 'Interpreter', 'latex', ...
    'FontSize', 20, 'VerticalAlignment', 'top')
set(gca, 'LineWidth', 1, 'FontSize', 24, 'FontName', 'Helvetica')
set(gca, 'XGrid', 'off', 'YGrid', 'off')
box on
saveas(gcf, [pwd '/Data/Pics/tick_statistics_nofilter.png'])
saveas(gcf, [pwd '/Data/Pics/tick_statistics_nofilter.fig'])

% Figure 2: With filter
figure('Position', [100, 100, 700, 500])
histogram(dtj_filter(2:end), bin, 'Normalization', 'probability', 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'none')
xline([0 1 2 3 4], '-k', 'FontSize', 2, 'HandleVisibility', 'off');
xlim([-0.1, 3])
xlabel('$\Omega_{\rm m}\tau/\pi$', 'Interpreter', 'latex', 'FontSize', 24)
text(0.05, 0.95, sprintf('$\\mathcal{N} = %.1f$\n$\\nu = %.2f\\Omega_{\\rm m}/\\pi$', ...
    N_filter, nu_filter), 'Units', 'normalized', 'Interpreter', 'latex', ...
    'FontSize', 20, 'VerticalAlignment', 'top')
set(gca, 'LineWidth', 1, 'FontSize', 24, 'FontName', 'Helvetica')
set(gca, 'XGrid', 'off', 'YGrid', 'off')
box on
% saveas(gcf, [pwd '/Data/Pics/tick_statistics_filter.png'])
% saveas(gcf, [pwd '/Data/Pics/tick_statistics_filter.fig'])

fprintf('\nTick statistics plots saved!\n');
plot_filter = 0; % Reset filter flag

%% 4.--- Allan Variance (Filtered Data)

% Allan variance for filtered tick data
imax=100; % Number of trajectories to analyze (increase as more data becomes available)
det_filt = 1; % Enable detector filter
Allan

% Update labels to match the desired format
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 22)
ylabel('$\tilde{\sigma}^2_A(m)$', 'Interpreter', 'latex', 'FontSize', 22)
xlim([0, 525])
ylim([3e-3, 3e-1])
set(gca, 'LineWidth', 1, 'FontSize', 22, 'FontName', 'Helvetica')
grid on
box on
title('') % Remove title

fprintf('\nAllan variance plot saved!\n');

%% 5.---Autocorrelations
% please see another code with the same name!