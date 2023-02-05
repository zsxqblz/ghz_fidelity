function expGHZ(meas_fidelity::Float64,T2::Float64,n::Float64,meas_t::Float64,nsites::Int64,nsim::Int64,showPorg=false)
    population_l = Vector{Int64}(undef,nsim)
    coherence_l = Vector{Int64}(undef,nsim)

    if showPorg
        @showprogress for i in 1:nsim
            population,coherence = simGHZ(meas_fidelity,T2,n,meas_t,nsites)
            population_l[i] = population
            coherence_l[i] = coherence
        end
    else
        for i in 1:nsim
            population,coherence = simGHZ(meas_fidelity,T2,n,meas_t,nsites)
            population_l[i] = population
            coherence_l[i] = coherence
        end
    end

    return population_l, coherence_l
end

function nsimConv(meas_fidelity::Float64,T2::Float64,n::Float64,meas_t::Float64,nsites::Int64,nsim_start::Int64,nsim_end::Int64,nsim_length::Int64,showProg::Bool=false)
    nsim_l = Int.(floor.(exp.(range(log(nsim_start),stop=log(nsim_end),length=nsim_length))))
    population_ave_l = Vector{Float64}(undef,nsim_length)
    population_std_l = Vector{Float64}(undef,nsim_length)
    coherence_ave_l = Vector{Float64}(undef,nsim_length)
    coherence_std_l = Vector{Float64}(undef,nsim_length)
    
    if showProg
        for (i,nsim) in enumerate(nsim_l)
            @printf "current nsim: %d" nsim
            population_l, coherence_l = expGHZ(meas_fidelity,T2,n,meas_t,nsites,nsim)
            population_ave_l[i] = mean(population_l)
            population_std_l[i] = std(population_l)
            coherence_ave_l[i] = mean(coherence_l)
            coherence_std_l[i] = std(coherence_l)

        end
    else
        for (i,nsim) in enumerate(nsim_l)
            population_l, coherence_l = expGHZ(meas_fidelity,T2,n,meas_t,nsites,nsim)
            population_ave_l[i] = mean(population_l)
            population_std_l[i] = std(population_l)
            coherence_ave_l[i] = mean(coherence_l)
            coherence_std_l[i] = std(coherence_l)

        end
    end
    return population_ave_l,population_std_l,coherence_ave_l,coherence_std_l 
end