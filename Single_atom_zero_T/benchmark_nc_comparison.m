%% Benchmark script to compare Julia simulation performance for different Nc values
% Follows the pattern from run_factorisation_and_resample.m section 1

clear all
close all

%% Configuration
julia_path = '/opt/homebrew/bin/julia';
nc_values = [60, 20];  % Test both Nc values
results = struct();

%% Parameters (same as in Factorisation.m)
w_m = 1;                 
f = 20 * w_m;            
g = 30 * w_m;            
w_hot = 240;
w_cold = 120;
Delta = 0 * max(f,g);
w_cav = w_hot - w_cold + Delta;
k = 10 * w_m;
g_h = 1 * w_m;
g_c = 100 * w_m;
n_h = 10;
n_c = 0;
g_m = 0.01 * w_m;
r_alpha = 0.9;
w_m = r_alpha * w_m;
dt = 1e-4;
ur = 0;  % Run unconditional first
tmax = 10;
Np = 1;

%% Prepare Julia parameters
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

param_file = 'params_interop.mat';
result_file = 'results_interop.mat';

% Save parameters
save(param_file, '-struct', 'julia_params');

%% Run benchmark for each Nc value
fprintf('======================================\n');
fprintf('Julia Simulation Benchmark - Nc Comparison\n');
fprintf('======================================\n\n');

for idx = 1:length(nc_values)
    nc = nc_values(idx);
    
    fprintf('Running simulation with Nc = %d...\n', nc);
    fprintf('--------------------------------------\n');
    
    % Modify the Julia file to use the desired Nc
    % Read original file
    fid = fopen('single_atom_simulation.jl', 'r');
    content = fread(fid, '*char')';
    fclose(fid);
    
    % Create modified version
    modified_content = regexprep(content, 'Nc = \d+', sprintf('Nc = %d', nc));
    temp_file = sprintf('single_atom_simulation_nc%d.jl', nc);
    fid = fopen(temp_file, 'w');
    fprintf(fid, '%s', modified_content);
    fclose(fid);
    
    % Run Julia with timing
    cmd = sprintf('%s %s "%s" "%s"', julia_path, temp_file, param_file, result_file);
    
    tic;
    [status, cmdout] = system(cmd, '-echo');
    elapsed_time = toc;
    
    if status ~= 0
        fprintf('ERROR: Julia simulation failed for Nc = %d\n', nc);
        fprintf('%s\n', cmdout);
        continue;
    end
    
    % Store results
    results(idx).Nc = nc;
    results(idx).time = elapsed_time;
    results(idx).output = cmdout;
    
    fprintf('\n');
    fprintf('Completed Nc = %d in %.2f seconds\n', nc, elapsed_time);
    fprintf('======================================\n\n');
    
    % Clean up temporary file
    delete(temp_file);
end

%% Display comparison
fprintf('\n');
fprintf('======================================\n');
fprintf('BENCHMARK RESULTS\n');
fprintf('======================================\n');
for idx = 1:length(results)
    fprintf('Nc = %2d: %8.2f seconds\n', results(idx).Nc, results(idx).time);
end
fprintf('\n');

if length(results) == 2
    speedup = results(1).time / results(2).time;
    fprintf('Speedup (Nc=%d vs Nc=%d): %.2fx faster\n', ...
        results(2).Nc, results(1).Nc, speedup);
    
    if speedup > 1
        fprintf('Nc = %d is %.1f%% faster than Nc = %d\n', ...
            results(2).Nc, (speedup-1)*100, results(1).Nc);
    else
        fprintf('Nc = %d is %.1f%% slower than Nc = %d\n', ...
            results(2).Nc, (1-speedup)*100, results(1).Nc);
    end
end
fprintf('======================================\n');
