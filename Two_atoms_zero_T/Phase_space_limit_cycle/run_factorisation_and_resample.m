%%%Warning, this is a slow code, only for phase space limit cycle!!
%% Simulate trajectories.
% % If you have done it already, comment for further analysis
clear all
iM=1;
imin=1;
imax=1;%The number of simulations (After reaching the steady state)
w_cold=120;w_hot=240;
%T_c=120;
%n_c=1./(exp(w_c./T_c)-1);
n_h=10;T_h=w_hot/(log((n_h+1)/n_h));
n_c=0;T_c=0;
sub_folder_name='Data';
mkdir(sub_folder_name)
for ur=0:0
    if ur==0
        i1=1;
        Factorisation;
        myVars = {'x_m_vec','p_m_vec','w_hot','w_cold','w_cav','n_h','n_c','w_m','f','g'...
            ,'k','g_h','g_c','g_m','dt',...
            "p1","p2","p3","na","re_ad_s12","im_ad_s12","na_p3","tvec"};
        save([sub_folder_name,'/unconditional'],myVars{:});
        % plot(x_m_vec,1i*p_m_vec,'LineWidth',2);
        % xlim([-40 40])
        % ylim([-50 50])
    else
        for i1=imin:imax
            [i1,imax]
            Factorisation;
            tvec_dN1=jump_times;
            myVars2={'x_m_vec','p_m_vec','w_hot','w_cold','w_cav',...
                'n_h','n_c','w_m','f','g','k','g_h','g_c','g_m','dt',...
                'p1_vec','p2_vec','p3_vec','na_vec','tvec_down','tvec_dN1'};
            save([sub_folder_name,'/conditional_traj',num2str(i1)],myVars2{:});
        end
    end
end
%% 2.--- Load and Plot the asymptotic limit cycles; and the populations
figure
sub_folder_name='Data';
myVars = {'x_m_vec','p_m_vec'};
load([sub_folder_name,'/unconditional'],myVars{:})
%subplot(1,4,iM);
% Plot full trajectory with thinner, darker blue line
hold on
plot(x_m_vec,p_m_vec,'LineWidth',1,'Color',[0 0.4470 0.7410])
% Plot last 2% as asymptotic limit with brighter, thicker line
cutoff = floor(0.98*length(x_m_vec));
plot(x_m_vec(cutoff:end),p_m_vec(cutoff:end),'LineWidth',2.5,'Color',[0.3010 0.7450 0.9330])
hold off
axis equal
set(gca,'linewidth',1.5,'FontSize',14)
ylabel('$\langle \widehat{p}_{\rm m} \rangle$','Interpreter','latex','FontSize',18);
xlabel('$\langle \widehat{x}_{\rm m} \rangle$','Interpreter','latex','FontSize',18);
box on
xlim([-15,20])
ylim([-20,20])

%% 2.1.---Plot time evolution of mechanical observables, populations, and cavity photon number
figure('Position', [100, 100, 800, 800])
myVars_time = {'x_m_vec','p_m_vec','p1','p2','p3','na','w_m','tvec'};
load([sub_folder_name,'/unconditional'],myVars_time{:})

% Create normalized time vector: Ω_m*t/π
t_normalized = w_m * tvec / pi;

% Select data range: 40.5 to 43.5
idx = (t_normalized >= 40.5) & (t_normalized <= 43.5);
t_plot = t_normalized(idx);
x_m_plot = x_m_vec(idx);
p_m_plot = p_m_vec(idx);
p1_plot = p1(idx);
p2_plot = p2(idx);
p3_plot = p3(idx);
na_plot = na(idx);

% Subplot (a): Mechanical position and momentum
subplot(3,1,1)
plot(t_plot, x_m_plot, 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
hold on
plot(t_plot, p_m_plot, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980])
hold off
set(gca, 'linewidth', 1.5, 'FontSize', 12)
ylabel('$\langle \widehat{x}_{\rm m} \rangle$, $\langle \widehat{p}_{\rm m} \rangle$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'$\langle \widehat{x}_{\rm m} \rangle$', '$\langle \widehat{p}_{\rm m} \rangle$'}, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
text(0.02, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top')
xlim([40.5 43.5])
box on
grid off

% Subplot (b): Atomic level populations
subplot(3,1,2)
plot(t_plot, p1_plot, 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
hold on
plot(t_plot, p2_plot, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980])
plot(t_plot, p3_plot, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
hold off
set(gca, 'linewidth', 1.5, 'FontSize', 12)
ylabel('$\langle \widehat{\rho}_1 \rangle$, $\langle \widehat{\rho}_2 \rangle$, $\langle \widehat{\rho}_3 \rangle$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'$\langle \widehat{\rho}_1 \rangle$', '$\langle \widehat{\rho}_2 \rangle$', '$\langle \widehat{\rho}_3 \rangle$'}, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
text(0.02, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top')
xlim([40.5 43.5])
box on
grid off

% Subplot (c): Cavity photon number
subplot(3,1,3)
plot(t_plot, na_plot, 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
set(gca, 'linewidth', 1.5, 'FontSize', 12)
ylabel('$\langle \widehat{a}^\dagger \widehat{a} \rangle$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'$\langle \widehat{a}^\dagger \widehat{a} \rangle$'}, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
text(0.02, 0.95, '(c)', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top')
xlim([40.5 43.5])
box on
grid off

%% 3. --- %% Load and Plot the asymptotic limit cycles conditional
figure('Position', [100, 100, 600, 600])
M=1;
sub_folder_name='Data';
myVars = {'x_m_vec','p_m_vec','w_m','tvec'};
load([sub_folder_name,'/conditional_traj1'],myVars{:})

% Calculate time in periods and select first 100 periods
t_normalized = w_m * tvec / (2*pi);  % time in periods
idx = t_normalized <= 100;
x_m_plot = x_m_vec(idx);
p_m_plot = p_m_vec(idx);

% Plot with thinner lines matching the style
plot(x_m_plot,p_m_plot,'LineWidth',0.5,'Color',[0 0.4470 0.7410])
axis equal
ylabel('$\langle \widehat{p}_{\rm m} \rangle$','Interpreter','latex','FontSize',24);
xlabel('$\langle \widehat{x}_{\rm m} \rangle$','Interpreter','latex','FontSize',24);
xlim([-50 50])
ylim([-50 50])
box on
grid on
set(gca,'GridColor',[0.85 0.85 0.85],'GridAlpha',1,'GridLineStyle','-','Layer','bottom','LineWidth',0.5,'FontSize',24,'TickLabelInterpreter','latex','linewidth',1)
%%%
%%%
%% 3.1.---Plot time evolution for conditional trajectory with jump times
figure('Position', [100, 100, 800, 800])
sub_folder_name='Data';
% Try to load conditional trajectory data
try
    myVars_cond = {'x_m_vec','p_m_vec','p1_vec','p2_vec','p3_vec','na_vec','w_m','tvec_down','tvec_dN1'};
    load([sub_folder_name,'/conditional_traj1'],myVars_cond{:})
catch ME
    error(['Conditional trajectory data incomplete. Required variables: p1_vec, p2_vec, p3_vec, na_vec, tvec_down, tvec_dN1.\n' ...
           'To generate this data:\n' ...
           '1. Change line 15 from "for ur=0:0" to "for ur=1:1"\n' ...
           '2. Run the script\n' ...
           '3. Change it back to "for ur=0:0" to avoid re-running simulations\n' ...
           'Original error: %s'], ME.message);
end

% Create normalized time vector: Ω_m*t/π
t_normalized = w_m * tvec_down / pi;
jump_times_norm = w_m * tvec_dN1 / pi;

% Select data range
t_range = [0 3.5];
idx = (t_normalized >= t_range(1)) & (t_normalized <= t_range(2));
t_plot = t_normalized(idx);
x_m_plot = x_m_vec(idx);
p_m_plot = p_m_vec(idx);
p1_plot = p1_vec(idx);
p2_plot = p2_vec(idx);
p3_plot = p3_vec(idx);
na_plot = na_vec(idx);
jump_idx = (jump_times_norm >= t_range(1)) & (jump_times_norm <= t_range(2));
jump_plot = jump_times_norm(jump_idx);

% Subplot (a): Mechanical position and momentum
subplot(3,1,1)
plot(t_plot, x_m_plot, 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
hold on
plot(t_plot, p_m_plot, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980])
% Add vertical lines for jump times
yl = ylim;
for ij = 1:length(jump_plot)
    plot([jump_plot(ij) jump_plot(ij)], yl, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off')
end
hold off
set(gca, 'linewidth', 1.5, 'FontSize', 12)
ylabel('$\langle \widehat{x}_{\rm m} \rangle$, $\langle \widehat{p}_{\rm m} \rangle$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'$\langle \widehat{x}_{\rm m} \rangle$', '$\langle \widehat{p}_{\rm m} \rangle$'}, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
text(0.02, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top')
xlim(t_range)
box on
grid off

% Subplot (b): Atomic level populations
subplot(3,1,2)
plot(t_plot, p1_plot, 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
hold on
plot(t_plot, p2_plot, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980])
plot(t_plot, p3_plot, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
% Add vertical lines for jump times
yl = ylim;
for ij = 1:length(jump_plot)
    plot([jump_plot(ij) jump_plot(ij)], yl, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off')
end
hold off
set(gca, 'linewidth', 1.5, 'FontSize', 12)
ylabel('$\langle \widehat{\rho}_1 \rangle$, $\langle \widehat{\rho}_2 \rangle$, $\langle \widehat{\rho}_3 \rangle$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'$\langle \widehat{\rho}_1 \rangle$', '$\langle \widehat{\rho}_2 \rangle$', '$\langle \widehat{\rho}_3 \rangle$'}, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
text(0.02, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top')
xlim(t_range)
box on
grid off

% Subplot (c): Cavity photon number
subplot(3,1,3)
plot(t_plot, na_plot, 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
hold on
% Add vertical lines for jump times
yl = ylim;
for ij = 1:length(jump_plot)
    plot([jump_plot(ij) jump_plot(ij)], yl, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off')
end
hold off
set(gca, 'linewidth', 1.5, 'FontSize', 12)
ylabel('$\langle \widehat{a}^\dagger \widehat{a} \rangle$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\Omega_{\rm m} t/\pi$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'$\langle \widehat{a}^\dagger \widehat{a} \rangle$'}, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
text(0.02, 0.95, '(c)', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top')
xlim(t_range)
box on
grid off

%% 4.--- The real deal (histogram etc)
figure('Position', [100, 100, 600, 600])
iM=1;
sub_folder_name='Data';
myVars = {'x_m_vec','p_m_vec'};
load([sub_folder_name,'/conditional_traj1'],myVars{:})

% Create 2D histogram
histogram2(x_m_vec', p_m_vec', [100, 100], ...
    'DisplayStyle', 'tile', 'ShowEmptyBins', 'on', 'EdgeColor', 'None')
axis equal
colormap('jet')
ylabel('$\langle \widehat{p}_{\rm m} \rangle$','Interpreter','latex','FontSize', 24);
xlabel('$\langle \widehat{x}_{\rm m} \rangle$','Interpreter','latex','FontSize', 24);
xlim([-30 30])
ylim([-30 30])
xticks(-30:10:30)
xticklabels({'-30','-20','','0','','20','30'})
yticks(-30:10:30)
box on
set(gca,'linewidth',1,'FontSize',24,'TickLabelInterpreter','latex')
colorbar off

hold on
% Plot the unconditional limit cycle on top
myVars_unconditional = {'x_m_vec','p_m_vec'};
load([sub_folder_name,'/unconditional'],myVars_unconditional{:})
% Use last 5% for the limit cycle
sx_m_vec = length(x_m_vec);
x_m_cycle = x_m_vec(1,floor(.95*sx_m_vec):end);
p_m_cycle = p_m_vec(1,floor(.95*sx_m_vec):end);
plot(x_m_cycle, p_m_cycle, 'LineWidth', 2, 'Color', 'r')
hold off

%% End of plotting sections
