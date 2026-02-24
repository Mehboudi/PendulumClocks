%% Plot Phase Space Limit Cycle from Conditional Trajectory
% This script loads x and p quadrature data from conditional_traj1.mat
% and plots the phase space trajectory and histogram

clear all
close all

%% Load data from conditional trajectory
sub_folder_name = 'Phase_space_limit_cycle/Data';
myVars = {'x_m_vec','p_m_vec'};
load([sub_folder_name,'/conditional_traj1'],myVars{:})

%% 1. Plot the full phase space trajectory
figure('Position', [100, 100, 600, 600])
plot(x_m_vec, p_m_vec, 'LineWidth', 0.5, 'Color', [0 0.4470 0.7410])
axis equal
ylabel('$\langle \widehat{p}_{\rm m} \rangle$', 'Interpreter', 'latex', 'FontSize', 24);
xlabel('$\langle \widehat{x}_{\rm m} \rangle$', 'Interpreter', 'latex', 'FontSize', 24);
xlim([-60 60])
ylim([-60 60])
box on
grid on
set(gca, 'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, 'GridLineStyle', '-', ...
    'Layer', 'bottom', 'LineWidth', 0.5, 'FontSize', 24, ...
    'TickLabelInterpreter', 'latex', 'linewidth', 1)
title('Phase Space Trajectory', 'Interpreter', 'latex', 'FontSize', 20)

%% 2. Plot 2D histogram of phase space
figure('Position', [100, 100, 600, 600])

% Add corner points to ensure full range coverage
x_extended = [x_m_vec(:); -60; -60; 60; 60];
p_extended = [p_m_vec(:); -60; 60; -60; 60];

% Create 2D histogram
histogram2(x_extended, p_extended, [100, 100], ...
    'DisplayStyle', 'tile', 'ShowEmptyBins', 'on', 'EdgeColor', 'None')
axis equal
colormap('jet')
ylabel('$\langle \widehat{p}_{\rm m} \rangle$','Interpreter','latex','FontSize', 24);
xlabel('$\langle \widehat{x}_{\rm m} \rangle$','Interpreter','latex','FontSize', 24);
xlim([-60 60])
ylim([-60 60])
xticks(-60:20:60)
yticks(-60:20:60)
box on
set(gca,'linewidth',1,'FontSize',24,'TickLabelInterpreter','latex')
colorbar off
title('Phase Space Histogram', 'Interpreter', 'latex', 'FontSize', 20)
