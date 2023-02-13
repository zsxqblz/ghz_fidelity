include("dependencies.jl")
include("ghz_sim.jl")
include("ghz_exp.jl")

const meas_fidelity = 0.99
const T2_fidelity_start = parse(Float64,ARGS[1])
const T2_fidelity_end = parse(Float64,ARGS[2])
const T2_fidelity_length = parse(Int64,ARGS[3])
const n = 2.
const meas_t = 10.
const nsim = 10000
const nsites_start = parse(Int64,ARGS[4])
const nsites_end = parse(Int64,ARGS[5])
const nsites_steps = parse(Int64,ARGS[6])
const file_name = ARGS[7]

const T2_fidelity_l = range(T2_fidelity_start,stop=T2_fidelity_end,length=T2_fidelity_length)
const nsites_l = range(nsites_start,stop=nsites_end,step=nsites_steps)
const nsites_length = length(nsites_l)

population_ave_arr = Array{Float64,2}(undef,nsites_length,T2_fidelity_length)
population_std_arr = Array{Float64,2}(undef,nsites_length,T2_fidelity_length)
coherence_ave_arr = Array{Float64,2}(undef,nsites_length,T2_fidelity_length)
coherence_std_arr = Array{Float64,2}(undef,nsites_length,T2_fidelity_length)
@showprogress for (i,T2) in enumerate(T2_fidelity_l)
    ~,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l = sweepNsites(meas_fidelity,T2,n,meas_t,nsites_start,nsites_end,nsites_steps,nsim,false)
    population_ave_arr[:,i] = population_ave_l
    population_std_arr[:,i] = population_std_l 
    coherence_ave_arr[:,i] = coherence_ave_l 
    coherence_std_arr[:,i] = coherence_std_l  
end


save2DData(T2_fidelity_l,nsites_l,population_ave_arr,population_std_arr,coherence_ave_arr,coherence_std_arr,file_name)