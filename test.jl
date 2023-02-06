include("dependencies.jl")
include("ghz_sim.jl")
include("ghz_exp.jl")

for x in ARGS; println(x); end


############## sweep nsites ###############

meas_fidelity = 0.99
T2 = 100.
n = 2.
meas_t = 10.
nsim = 10000
nsites_start = 3
nsites_end = 50
nsites_steps = 3
nsites_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l = sweepNsites(meas_fidelity,T2,n,meas_t,nsites_start,nsites_end,nsites_steps,nsim,true)

file_name="data/nsites_3_50_3.csv" 
saveData(nsites_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l,file_name)


############## sweep fidelity ###############

meas_fidelity_start = 1.0
meas_fidelity_end = 0.5
meas_fidelity_length = 50
T2 = 100.
n = 2.
meas_t = 10.
nsim = 10000
nsites = 20
meas_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l = sweepFidelity(meas_fidelity_start,meas_fidelity_end,meas_fidelity_length,T2,n,meas_t,nsites,nsim,true)

file_name="data/fidelity_1.0_0.5_50.csv" 
saveData(meas_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l,file_name)

############## sweep T2 ###############
meas_fidelity = 0.99
T2_start = 100.
T2_end = 10.
T2_length = 50
n = 2.
meas_t = 10.
nsim = 10000
nsites = 20
meas_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l = sweepT2(meas_fidelity,T2_start,T2_end,T2_length,n,meas_t,nsites,nsim,true)

file_name="data/T2_100_10_50_2.csv" 
saveData(meas_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l,file_name)

############### visualization ################

# file_name="data/T2_100_10_50_2.csv" 
nscan_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l = readData(file_name)
plot(nscan_l[1:30],population_ave_l[1:30],yaxis=:log)
plot(nscan_l,coherence_ave_l)


