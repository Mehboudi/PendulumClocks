%% Factorisation.m - MATLAB Wrapper for Julia Simulation
% This script replaces the old simulation loop with a call to the robust Julia engine.
% It maintains the parameter definitions expected by the parent script.

% --- CONFIGURATION: JULIA PATH ---
% Set the specific path to the Julia executable found via 'which julia'
julia_path = '/opt/homebrew/bin/julia'; 
% ---------------------------------

% --- 1. Parameter Definitions (Preserved from original) ---
% These define the physics of the clock.
w_m = 1;                 % Mechanical frequency
f = 20 * w_m;            % Driving strength
g = 30 * w_m;            % Coupling
% Interaction terms
epsilon_1 = 0 * max(f,g);
epsilon_2 = 4 * max(f,g);
epsilon_3 = 8 * max(f,g);
Delta = 0 * max(f,g);
w_cav = w_hot - w_cold + Delta;
k = 10 * w_m;
g_h = 1 * w_m;
g_c = 100 * w_m;

% Thermals are set by parent script, but we ensure defaults if testing standalone
if ~exist('n_h', 'var'), n_h = 10; end
if ~exist('n_c', 'var'), n_c = 0; end

T_opt = T_c;
if T_opt == 0
    n_opt = 0;
else
    n_opt = 1./(exp(w_cav/T_opt)-1);
end

g_m = 0.01 * w_m;        % Mechanical damping

% Rescaling
r_alpha = .9;
w_m = r_alpha * w_m;

% Time settings
dt = 1e-4;               % Time step
if ~exist('ur', 'var')
    ur = 0; 
end

if ur == 0
    tmax = 10;           % Time for unconditional (Reduced by 10x)
    Np = 100;
else
    tmax = 50;           % Time for conditional (Reduced by 10x)
    Np = 1000;
end

% --- 2. Prepare Interface for Julia ---

% Structure parameters for export
julia_params.ur = ur;
julia_params.Np = Np;
julia_params.w_m = w_m;
julia_params.f = f;
julia_params.g = g;
julia_params.k = k;
julia_params.g_h = g_h;
julia_params.g_c = g_c;
julia_params.n_h = n_h;
julia_params.n_c = n_c;
julia_params.g_m = g_m;
julia_params.dt = dt;
julia_params.tmax = tmax;

% File names for communication
param_file = 'params_interop.mat';
result_file = 'results_interop.mat';

% Save parameters to .mat file (for both ur==0 and ur==1)
save(param_file, '-struct', 'julia_params');

if ur == 0
    result_op.rho0 = 0;
    result_op.xm0 = 0;
    result_op.pm0 = 0;
    result_op.exp_numc0 = 0;
    save(result_file, '-struct', 'result_op');
end

% --- 3. Call Julia Simulation ---
% We use 'system' to call the Julia script.
% Use the configured julia_path
% ADDED '-echo' to ensure progress messages are printed to MATLAB Command Window
cmd = sprintf('%s single_atom_simulation.jl "%s" "%s"', julia_path, param_file, result_file);
[status, cmdout] = system(cmd, '-echo');

if status ~= 0
    disp('--- Julia Output (Error) ---');
    disp(cmdout);
    error('Julia simulation failed. See output above. \nHint: If package installation fails repeatedly, try running this in your terminal:\n/opt/homebrew/bin/julia -e "using Pkg; Pkg.add([\"MAT\", \"DifferentialEquations\", \"QuantumOptics\"])"');
end

% --- 4. Load and Map Results ---
% Load the results produced by Julia
if exist(result_file, 'file')
    raw_results = load(result_file);
    % Debug: Show what fields are in the loaded file
    disp('Fields in loaded result file:');
    disp(fieldnames(raw_results));
else
    error('Results file not found. Julia simulation may have failed silently.');
end

if ur == 0
    % Unconditional Case: Expecting phase space vectors
    x_m_vec = raw_results.x_m_vec;
    p_m_vec = raw_results.p_m_vec;
    
    % Load population and cavity data
    p1 = raw_results.p1;
    p2 = raw_results.p2;
    p3 = raw_results.p3;
    na = raw_results.na;
    re_ad_s12 = raw_results.re_ad_s12;
    im_ad_s12 = raw_results.im_ad_s12;
    na_p3 = raw_results.na_p3;
    tvec = raw_results.tvec;
    
    % Calculate final values from vectors
    x_m = x_m_vec(end);
    p_m = p_m_vec(end);
    
    % Ensure correct orientation (row/col) if necessary
    if size(x_m_vec, 1) > 1, x_m_vec = x_m_vec'; end
    if size(p_m_vec, 1) > 1, p_m_vec = p_m_vec'; end
    if size(p1, 1) > 1, p1 = p1'; end
    if size(p2, 1) > 1, p2 = p2'; end
    if size(p3, 1) > 1, p3 = p3'; end
    if size(na, 1) > 1, na = na'; end
    if size(tvec, 1) > 1, tvec = tvec'; end
    
else
    % Conditional Case: Expecting downsampled time traces
    tvec_down = raw_results.tvec_down;
    x_m_vec = raw_results.t_x_m_down;
    p_m_vec = raw_results.t_p_m_down;
    p1_vec = raw_results.t_p1_down;
    p2_vec = raw_results.t_p2_down;
    p3_vec = raw_results.t_p3_down;
    na_vec = raw_results.t_na_down;
    re_ad_s12_vec = raw_results.t_re_ad_s12_down;
    im_ad_s12_vec = raw_results.t_im_ad_s12_down;
    
    % Jumps
    t_dN = raw_results.t_dN;
    tvec_dN1 = raw_results.tvec_dN1;
    
    % Ensure orientation (MATLAB often prefers rows for these legacy scripts)
    if size(tvec_down, 1) > 1, tvec_down = tvec_down'; end
    if size(tvec_dN1, 1) > 1, tvec_dN1 = tvec_dN1'; end
    
    % ALIAS: The parent script expects 'jump_times' to exist
    jump_times = tvec_dN1;
end

% Cleanup temporary files (Optional: keep them for debugging if needed)
%delete(param_file);
%delete(result_file);

fprintf('Julia simulation (ur=%d) completed successfully.\n', ur);