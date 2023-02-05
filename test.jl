include("dependencies.jl")
include("ghz_sim.jl")
include("ghz_exp.jl")

# population_l, coherence_l = expGHZ(0.9,100.,2.,10.,10,10000)


# @printf "population | average: %.3f std dev: %.3f" mean(population_l) std(population_l)


# nsim_l = Int.(floor.(exp.(range(log(1),stop=log(10),length=10))))

population_ave_l,population_std_l,coherence_ave_l,coherence_std_l = nsimConv(0.9,100.,2.,10.,20,100,100000,10,true)

nsim_l = Int.(floor.(exp.(range(log(10),stop=log(100000),length=10))))
plot(nsim_l,population_ave_l,yerror=population_std_l,xaxis=:log)
plot(nsim_l,coherence_ave_l,yerror=coherence_std_l,xaxis=:log)
