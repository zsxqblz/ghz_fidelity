function genPlus(nsites::Int64)
    s = zero(Stabilizer, nsites)

    for i in 1:nsites
        s[i,i] = (true,false)
    end

    return s
end

function genGHZ(nsites::Int64)
    s = zero(Stabilizer, nsites)

    for i in 1:(nsites-1)
        s[i,i] = (false,true)
        s[i,i+1] = (false,true)
    end

    for i in 1:nsites
        s[nsites,i] = (true,false)
    end

    return s
end

function genRGHZ(zz_gnd_str::Vector{Int64},nsites::Int64)
    s = zero(Stabilizer, nsites)

    for i in 1:(nsites-1)
        s.tab.phases[i] = zz_gnd_str[i] > 0 ? 0x0 : 0x2
        s[i,i] = (false,true)
        s[i,i+1] = (false,true)
    end

    for i in 1:nsites
        s[nsites,i] = (true,false)
    end

    return s
end

function fidelityToErrProb(fidelity::Float64)
    return 1-fidelity
end

function T2ToErrProb(t::Float64,T2::Float64,n::Float64)
    return 1/2*(1-exp(-(t/T2)^n))
end

function genZZMeas(meas_fidelity::Float64,nsites::Int64)
    zz_gnd_truth_bit = bitrand(nsites-1)
    zz_gnd_truth = Vector{Int64}(undef,nsites-1)
    zz_meas = Vector{Int64}(undef,nsites-1)

    zz_err_prob = fidelityToErrProb(meas_fidelity)
    zz_err_rand = rand(nsites-1)
    
    for i in 1:(nsites-1)
        zz_gnd_truth[i] = zz_gnd_truth_bit[i] ? 1 : -1
        err = zz_err_rand[i]>zz_err_prob ? 1 : -1
        zz_meas[i] = zz_gnd_truth[i]*err
    end

    return zz_gnd_truth, zz_meas
end

function genZZErr(zz_str::Vector{Int64},meas_fidelity::Float64,nsites::Int64)
    zz_err_prob = fidelityToErrProb(meas_fidelity)
    zz_err_rand = rand(nsites-1)
    zz_meas = Vector{Int64}(undef,nsites-1)

    for i in 1:(nsites-1)
        err = zz_err_rand[i]>zz_err_prob ? 1 : -1
        zz_meas[i] = zz_str[i]*err
    end

    return zz_meas
end

function applyDephasing(s::Stabilizer,T2::Float64,n::Float64,meas_t::Float64,nsites::Int64)
    z_err_prob = T2ToErrProb(meas_t,T2,n)
    z_err_rand = rand(nsites)
    for i in 1:nsites
        if z_err_rand[i]<z_err_prob
            apply!(s,sZ(i))
        end
    end
end

function applyGHZCorrection(s::Stabilizer,zz_str::Vector{Int64},nsites::Int64)
    parity = 1
    for i in 2:nsites
        parity = parity * zz_str[i-1]
        if parity == -1
            apply!(s,sX(i))
        end
    end
end

function genZStr(nsites::Int64)
    p = P"Z"
    for i in 2:nsites
        p = p ⊗ Z
    end
    return p
end

function genXStr(nsites::Int64)
    p = P"X"
    for i in 2:nsites
        p = p ⊗ X
    end
    return p
end

function genYStr(nsites::Int64)
    p = P"Y"
    for i in 2:nsites
        p = p ⊗ Y
    end
    return p
end

function measureZZ(s::Stabilizer,nsites::Int64)
    zz_str = Vector{Int64}(undef,nsites-1)
    for i in 1:(nsites-1)
        zz_op = single_z(nsites,i)
        zz_op[i+1] = (false,true)
        new_s, parity = projectrand!(s,zz_op)
        zz_str[i] = parity == 0x0 ? 1 : -1
    end
    return zz_str
end

function measurePopulation(s::Stabilizer,nsites::Int64)
    boolStr = Vector{Bool}(undef,nsites)
    reg = Register(s, boolStr)

    for i in 1:nsites
        apply!(reg, sMZ(i,i)) 
    end

    for i in 2:nsites
        if reg.bits[1] ⊻ reg.bits[i]
            return 0
        end
    end

    return 1
end

function measureCoherence(s::Stabilizer,nsites::Int64)
    xstr = genXStr(nsites)
    project!(s, xstr)
    newstate, anticomindex, result = project!(s, xstr)
    if isnothing(result)
        return 0
    else
        return result == 0x0 ? 1 : -1
    end

end

function simGHZ_ZX(meas_fidelity::Float64,T2::Float64,n::Float64,meas_t::Float64,nsites::Int64)
    s = genPlus(nsites)
    applyDephasing(s,T2,n,meas_t,nsites)
    zz_str = measureZZ(s,nsites)
    zz_meas = genZZErr(zz_str,meas_fidelity,nsites)
    applyGHZCorrection(s,zz_meas,nsites)
    
    population = measurePopulation(copy(s),nsites)
    coherence = measureCoherence(s,nsites)

    return population,coherence
end

function simGHZ(meas_fidelity::Float64,T2::Float64,n::Float64,meas_t::Float64,nsites::Int64)
    zz_gnd_truth, zz_meas = genZZMeas(meas_fidelity,nsites)
    zz_state = genRGHZ(zz_gnd_truth,nsites)

    applyDephasing(zz_state,T2,n,meas_t,nsites)
    applyGHZCorrection(zz_state,zz_meas,nsites)
    
    population = measurePopulation(copy(zz_state),nsites)
    coherence = measureCoherence(zz_state,nsites)

    return population,coherence
end