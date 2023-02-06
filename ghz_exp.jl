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
            @printf "current nsim: %d\n" nsim
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
    return nsim_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l 
end

function sweepNsites(meas_fidelity::Float64,T2::Float64,n::Float64,meas_t::Float64,nsites_start::Int64,nsites_end::Int64,nsites_steps::Int64,nsim::Int64,showProg::Bool=false)
    nsites_l = Int.(range(nsites_start,stop=nsites_end,step=nsites_steps))
    nsites_length = length(nsites_l)
    population_ave_l = Vector{Float64}(undef,nsites_length)
    population_std_l = Vector{Float64}(undef,nsites_length)
    coherence_ave_l = Vector{Float64}(undef,nsites_length)
    coherence_std_l = Vector{Float64}(undef,nsites_length)
    
    if showProg
        @showprogress for (i,nsites) in enumerate(nsites_l)
            population_l, coherence_l = expGHZ(meas_fidelity,T2,n,meas_t,nsites,nsim)
            population_ave_l[i] = mean(population_l)
            population_std_l[i] = std(population_l)
            coherence_ave_l[i] = mean(coherence_l)
            coherence_std_l[i] = std(coherence_l)

        end
    else
        for (i,nsites) in enumerate(nsites_l)
            population_l, coherence_l = expGHZ(meas_fidelity,T2,n,meas_t,nsites,nsim)
            population_ave_l[i] = mean(population_l)
            population_std_l[i] = std(population_l)
            coherence_ave_l[i] = mean(coherence_l)
            coherence_std_l[i] = std(coherence_l)

        end
    end
    return nsites_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l 
end

function sweepFidelity(meas_fidelity_start::Float64,meas_fidelity_end::Float64,meas_fidelity_length::Int64,T2::Float64,n::Float64,meas_t::Float64,nsites::Int64,nsim::Int64,showProg::Bool=false)
    meas_fidelity_l = range(meas_fidelity_start,stop=meas_fidelity_end,length=meas_fidelity_length)
    population_ave_l = Vector{Float64}(undef,meas_fidelity_length)
    population_std_l = Vector{Float64}(undef,meas_fidelity_length)
    coherence_ave_l = Vector{Float64}(undef,meas_fidelity_length)
    coherence_std_l = Vector{Float64}(undef,meas_fidelity_length)
    
    if showProg
        @showprogress for (i,meas_fidelity) in enumerate(meas_fidelity_l)
            population_l, coherence_l = expGHZ(meas_fidelity,T2,n,meas_t,nsites,nsim)
            population_ave_l[i] = mean(population_l)
            population_std_l[i] = std(population_l)
            coherence_ave_l[i] = mean(coherence_l)
            coherence_std_l[i] = std(coherence_l)

        end
    else
        for (i,meas_fidelity) in enumerate(meas_fidelity_l)
            population_l, coherence_l = expGHZ(meas_fidelity,T2,n,meas_t,nsites,nsim)
            population_ave_l[i] = mean(population_l)
            population_std_l[i] = std(population_l)
            coherence_ave_l[i] = mean(coherence_l)
            coherence_std_l[i] = std(coherence_l)

        end
    end
    return meas_fidelity_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l 
end

function sweepT2(meas_fidelity::Float64,T2_start::Float64,T2_end::Float64,T2_length::Int64,n::Float64,meas_t::Float64,nsites::Int64,nsim::Int64,showProg::Bool=false)
    T2_l = range(T2_start,stop=T2_end,length=T2_length)
    population_ave_l = Vector{Float64}(undef,meas_fidelity_length)
    population_std_l = Vector{Float64}(undef,meas_fidelity_length)
    coherence_ave_l = Vector{Float64}(undef,meas_fidelity_length)
    coherence_std_l = Vector{Float64}(undef,meas_fidelity_length)
    
    if showProg
        @showprogress for (i,T2) in enumerate(T2_l)
            population_l, coherence_l = expGHZ(meas_fidelity,T2,n,meas_t,nsites,nsim)
            population_ave_l[i] = mean(population_l)
            population_std_l[i] = std(population_l)
            coherence_ave_l[i] = mean(coherence_l)
            coherence_std_l[i] = std(coherence_l)

        end
    else
        for (i,T2) in enumerate(T2_length)
            population_l, coherence_l = expGHZ(meas_fidelity,T2,n,meas_t,nsites,nsim)
            population_ave_l[i] = mean(population_l)
            population_std_l[i] = std(population_l)
            coherence_ave_l[i] = mean(coherence_l)
            coherence_std_l[i] = std(coherence_l)

        end
    end
    return T2_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l 
end

function saveData(scan_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l,file_name)
    df = DataFrame()
    df.scan_l = scan_l
    df.population_ave_l = population_ave_l
    df.population_std_l = population_std_l
    df.coherence_ave_l = coherence_ave_l
    df.coherence_std_l = coherence_std_l
    CSV.write(file_name, df)
end

function readData(file_name)
    df = CSV.read(file_name, DataFrame)
    scan_l = df.scan_l
    population_ave_l = df.population_ave_l
    population_std_l = df.population_std_l
    coherence_ave_l = df.coherence_ave_l
    coherence_std_l = df.coherence_std_l
    return scan_l,population_ave_l,population_std_l,coherence_ave_l,coherence_std_l
end