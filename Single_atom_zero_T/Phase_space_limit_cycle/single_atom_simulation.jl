# import Pkg; Pkg.add("QuantumOptics")
# import Pkg; Pkg.add("CairoMakie")
using MAT
using QuantumOptics
using DifferentialEquations
using LinearAlgebra
using SparseArrays

# Parameters (same as in paper)

# --- 1. Load Parameters ---
if length(ARGS) < 2
    error("Usage: julia single_atom_simulation.jl <param_file.mat> <result_file.mat>")
end

params = matread(ARGS[1])

# Extract Parameters (matching your naming convention)
ur = params["ur"]
println("Running simulation with ur = $ur")
Ω_m = params["w_m"]
f = params["f"]
g = params["g"]
κ = params["k"]
Γ_h = params["g_h"]
Γ_c = params["g_c"]
n_h = params["n_h"]
n_c = params["n_c"]
g_m = params["g_m"]
dt_sim = params["dt"]
tmax = params["tmax"]
Np = params["Np"]

α = 2; # ratio ω_{23}/ω_{13}, needed for entropy production

Γ = κ + Γ_h*n_h;

# --- 2. System Definition (QuantumOptics) ---

Nc = 20 # Cutoff for cavity
B_cav = FockBasis(Nc)
B_atom = NLevelBasis(3)
B = B_atom ⊗ B_cav

# Operators cavity
a = one(B_atom) ⊗ destroy(B_cav)
at = one(B_atom) ⊗ create(B_cav)
numc = at * a 

# Operators, atom

p1 = transition(B_atom,1,1) ⊗ one(B_cav);
p2 = transition(B_atom,2,2) ⊗ one(B_cav);
p3 = transition(B_atom,3,3) ⊗ one(B_cav);

σ = transition(B_atom, 1, 2) ⊗ one(B_cav)
σt = dagger(σ)

# Dissipator Operators
# Hot bath: 1 <-> 3
σ13 = transition(B_atom, 1, 3) ⊗ one(B_cav) # |1><3|
σ31 = dagger(σ13)                           # |3><1|

# Cold bath: 2 <-> 3
σ23 = transition(B_atom, 2, 3) ⊗ one(B_cav) # |2><3|
σ32 = dagger(σ23)                           # |3><2|

# Hamiltonian
function H(x)
    return f*(at*σ+a*σt)-sqrt(2)*g*numc*x;
end

# Jump Operators, different for zero or finite cold temperatures
if n_c == 0
    J = [sqrt(κ*(n_c+1))*a,sqrt(Γ_h*(n_h+1))*σ13,sqrt(Γ_h*n_h)*σ31,sqrt(Γ_c*(n_c+1))*σ23]
    Jnj = [sqrt(κ*(n_c+1))*a,sqrt(Γ_h*(n_h+1))*σ13,sqrt(Γ_h*n_h)*σ31] # jumps excluding the monitored one
else
    J = [sqrt(κ*(n_c+1))*a,sqrt(κ*n_c)*at,sqrt(Γ_h*(n_h+1))*σ13,sqrt(Γ_h*n_h)*σ31,sqrt(Γ_c*(n_c+1))*σ23,sqrt(Γ_c*n_c)*σ32]
    Jnj = [sqrt(κ*(n_c+1))*a,sqrt(κ*n_c)*at,sqrt(Γ_h*(n_h+1))*σ13,sqrt(Γ_h*n_h)*σ31,sqrt(Γ_c*n_c)*σ32] # jumps excluding the monitored one
end

Jc = sqrt(Γ_c*(n_c+1))*σ23;  #monitored jump operator

JdJ = dagger.(J).*J;

# --- 3. Simulation parameters ---

tsim = Np*2*pi/Ω_m; # simulation time
δt =0.001; # timestep for solving the equation. For the average evolution, we need δt*Ω_m << 1
tvec = [0:δt:tsim;]; # array with all discrete times
Δt = 20*δt; # timestep for plotting

# Initial state for ur==0. For ur==1 the last value of ur==0 is taken
if ur == 0
    ψ0 = nlevelstate(B_atom,1) ⊗ fockstate(B_cav,0);
    ρ0 = tensor(ψ0,dagger(ψ0));
    xm0 = 0; # initial mechanical position. 
    pm0 = 0;
    exp_numc0 = 0; # initial cavity photon number expectation value.
else
    # load the initial states from the file saved after ur==0
    if !isfile(ARGS[2])
        error("Cannot run ur==1: Result file $(ARGS[2]) does not exist. Please run with ur==0 first to generate initial conditions.")
    end
    params_ur0 = matread(ARGS[2])["params_ur0"]
    # Reconstruct the density operator from the saved matrix
    ρ0 = DenseOperator(B, params_ur0["rho0"]);
    xm0 = params_ur0["xm0"];
    pm0 = params_ur0["pm0"];
    exp_numc0 = params_ur0["exp_numc0"];
end

# --- 4. average evolution ur==0 evolution (i.e., for the unconditional)
if ur == 0
global j = 0; #used for saving data only every Δt

# setting initial state
global ρ = ρ0;
global xm = xm0;
global pm = pm0;
global exp_numc = exp_numc0;

# arrays that save data
global xm_avg = []; # mechanical position
global pm_avg = []; # mechanical momentum

global numc_avg = []; # cavity number

global p1_avg = []; # occupation state 1
global p2_avg = []; # occupation state 2
global p3_avg = []; # occupation state 3

global as_avg = []; # ⟨a^†σ⟩

global tplot_avg = []; # array of times for plotting

# Progress monitoring
start_time = time()
last_update = start_time
update_interval = 1.0  # Update every 1 second
total_steps = length(tvec)


for (step_idx, t) in enumerate(tvec)
    
    # Progress monitoring
    if time() - last_update > update_interval || step_idx == total_steps
        global last_update
        elapsed = time() - start_time
        progress = step_idx / total_steps
        eta = elapsed / progress - elapsed
        speed = step_idx / elapsed
        
        print("\rAverage evolution: $(round(Int, progress*100))% | "
              * "Elapsed: $(round(Int, elapsed))s | "
              * "ETA: $(round(Int, eta))s | "
              * "Speed: $(round(speed, digits=1)) steps/s")
        flush(stdout)
        last_update = time()
    end

    global ρ, exp_numc, xm, pm, j
    tdisc, ρ = timeevolution.master_nh([0,δt], ρ, H(xm)-1im/2*sum(JdJ), J); # using the masterequation solver to evolve the density matrix by a single time-step, tdisc is a vector of times that is not needed
    ρ = ρ[end]; # setting ρ equal to the final density matrix from the above simulation (i.e., the state at t+δt)
    exp_numc = real(expect(numc, ρ)) #average photon number
    
    # evolving the mechanics
    xm += δt*(Ω_m*pm-g_m/2*xm);
    pm += δt*(-Ω_m*xm-g_m/2*pm+sqrt(2)*g*exp_numc);
    
    # saving data only every Δt
    if t > j*Δt
        push!(xm_avg,xm)
        push!(pm_avg,pm)
        push!(numc_avg,exp_numc)
        push!(p1_avg,real(expect(p1, ρ)))
        push!(p2_avg,real(expect(p2, ρ)))
        push!(p3_avg,real(expect(p3, ρ)))
        push!(as_avg,expect(at*σ, ρ))
        push!(tplot_avg,t)
        j+=1;
    end
end
# Taking the final  state as the initial state for ur==1
ρ0 = ρ;
xm0 = xm;
pm0 = pm;
exp_numc0 = exp_numc;
# Save only the data matrix of the density operator, not the full Operator object
params_ur0 = Dict("rho0" => ρ0.data, "xm0" => xm0, "pm0" => pm0, "exp_numc0" => exp_numc0)
# Note: params_ur0 will be saved at the end along with other results

# Heat currents
ω13 = 1000; # Needed to plot the heat currents
Jh = Γ_h*((n_h+1)*p3_avg-n_h*p1_avg)/Ω_m; # Jh/(Ω_m*ω_{13}
Jm = -sqrt(2)*g*numc_avg.*pm_avg/ω13; # Note that only the time-average of this quantity really makes sense

println("\nAverage evolution completed in $(round(time() - start_time, digits=2))s")

end # end of average evolution ur==0 evolution

# --- 5. stochastic evolution ur==1 evolution (i.e., for the conditional)
if ur == 1
global j = 0; #used for saving data only every Δt

# setting initial state
global ρ = ρ0;
global xm = xm0;
global pm = pm0;
global exp_numc = exp_numc0;

# arrays that save data
global xm_st = []; # mechanical position
global pm_st = []; # mechanical momentum

global numc_st = []; # cavity number

global p1_st = []; # occupation state 1
global p2_st = []; # occupation state 2
global p3_st = []; # occupation state 3

global as_st = []; # ⟨a^†σ⟩

global tplot_st = []; # array of times for plotting

global Jtimes = []; # list of jumps

# Progress monitoring
start_time_st = time()
last_update_st = start_time_st
update_interval_st = 1.0  # Update every 1 second
total_steps_st = length(tvec)

for (step_idx_st, t) in enumerate(tvec)
    
    # Progress monitoring
    if time() - last_update_st > update_interval_st || step_idx_st == total_steps_st
        global last_update_st
        elapsed_st = time() - start_time_st
        progress_st = step_idx_st / total_steps_st
        eta_st = elapsed_st / progress_st - elapsed_st
        speed_st = step_idx_st / elapsed_st
        
        print("\rStochastic evolution: $(round(Int, progress_st*100))% | "
              * "Elapsed: $(round(Int, elapsed_st))s | "
              * "ETA: $(round(Int, eta_st))s | "
              * "Speed: $(round(speed_st, digits=1)) steps/s | "
              * "Jumps: $(length(Jtimes))")
        flush(stdout)
        last_update_st = time()
    end
    
    global ρ, exp_numc, xm, pm, j
    pj = Γ_c*δt*real(expect(p3,ρ));
    if pj < rand()
        #evolving the state, no jump
        tdisc, ρ = timeevolution.master_nh([0,δt], ρ, H(xm)-1im/2*sum(JdJ), Jnj);
        ρ = ρ[end]/tr(ρ[end]);
        exp_numc = real(expect(numc, ρ))

        # evolving the mechanics
        xm += δt*(Ω_m*pm-g_m/2*xm);
        pm += δt*(-Ω_m*xm-g_m/2*pm+sqrt(2)*g*exp_numc);
       
    else
        #jump happening
        ρ = δt*Jc*ρ*dagger(Jc)/pj;
        push!(Jtimes,t); # add jump time to list
    end

    # saving data only every Δt
    if t > j*Δt
        push!(xm_st,xm)
        push!(pm_st,pm)
        push!(numc_st,exp_numc)
        push!(p1_st,real(expect(p1, ρ)))
        push!(p2_st,real(expect(p2, ρ)))
        push!(p3_st,real(expect(p3, ρ)))
        push!(as_st,expect(at*σ, ρ))
        push!(tplot_st,t)
        j+=1;
    end
end

println("\nStochastic evolution completed in $(round(time() - start_time_st, digits=2))s")
println("Total jumps recorded: $(length(Jtimes))")

end # end of stochastic evolution ur==1 evolution

# --- 6. Data Processing and Saving ---
if ur == 0
    # Convert to proper numeric vectors for MATLAB
    x_m_vec = Float64.(xm_avg)
    p_m_vec = Float64.(pm_avg)
    p1 = Float64.(p1_avg)
    p2 = Float64.(p2_avg)
    p3 = Float64.(p3_avg)
    na = Float64.(numc_avg)
    re_ad_s12 = Float64.(real.(as_avg))
    im_ad_s12 = Float64.(imag.(as_avg))
    na_p3 = Float64.(numc_avg .* p3_avg)
    tvec = Float64.(tplot_avg)
    
    # Save results for ur==0 (including params_ur0 which was defined earlier in the ur==0 block)
    results = Dict(
        "params_ur0" => params_ur0,
        "x_m_vec" => x_m_vec,
        "p_m_vec" => p_m_vec,
        "p1" => p1,
        "p2" => p2,
        "p3" => p3,
        "na" => na,
        "re_ad_s12" => re_ad_s12,
        "im_ad_s12" => im_ad_s12,
        "na_p3" => na_p3,
        "tvec" => tvec
    )
    matwrite(ARGS[2], results)
end

if ur == 1
    # Process results for ur==1 - convert to proper numeric vectors
    t_p1_down = Float64.(p1_st)
    t_p2_down = Float64.(p2_st)
    t_p3_down = Float64.(p3_st)
    t_na_down = Float64.(numc_st)
    t_re_ad_s12_down = Float64.(real.(as_st))
    t_im_ad_s12_down = Float64.(imag.(as_st))
    
    # Process jumps
    if isempty(Jtimes)
        tvec_dN1 = Float64[]
        t_dN = zeros(Float64, length(tplot_st))
    else
        tvec_dN1 = Float64.(Jtimes)
        t_dN = zeros(Float64, length(tplot_st))
    end
    
    results = Dict(
        "tvec_down" => Float64.(tplot_st),
        "t_x_m_down" => Float64.(xm_st),
        "t_p_m_down" => Float64.(pm_st),
        "t_p1_down" => t_p1_down,
        "t_p2_down" => t_p2_down,
        "t_p3_down" => t_p3_down,
        "t_na_down" => t_na_down,
        "t_re_ad_s12_down" => t_re_ad_s12_down,
        "t_im_ad_s12_down" => t_im_ad_s12_down,
        "t_dN" => t_dN,
        "tvec_dN1" => tvec_dN1,
        # Add aliases for compatibility with MATLAB script
        "x_m_vec" => Float64.(xm_st),
        "p_m_vec" => Float64.(pm_st)
    )
    println("Saving ur==1 results to $(ARGS[2])")
    println("Result keys: ", keys(results))
    matwrite(ARGS[2], results)
    println("Successfully saved ur==1 results")
end