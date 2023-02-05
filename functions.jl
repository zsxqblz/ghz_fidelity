function multPauli(lstr::PauliStr,rstr::PauliStr)
    out_str = copy(lstr.pauli_op)
    for i in rstr.idx_active
        out_str[i] = lstr.pauli_op[i]*rstr.pauli_op[i]
    end

    out_idx = Vector{Int64}(undef,0)
    for i in 1:lstr.nsites
        in_lsites = length(findall( x -> x == i, lstr.idx_active ))
        in_rsites = length(findall( x -> x == i, rstr.idx_active ))
        if in_lsites + in_rsites > 0
            append!(out_idx,i)
        end
    end
    return PauliStr(out_idx,out_str)
end

function evolPauli(pstr::PauliStr,ustr::PauliStr)
    for i in ustr.idx_active
        pstr.pauli_op[i] = ustr.pauli_op[i]*pstr.pauli_op[i]*ustr.pauli_op[i]'
    end
end

function kronPauli(str::PauliStr)
    kron_mat = ComplexF64[ 1.0+0.0im ]
        for si in str.pauli_op
            kron_mat = kron(kron_mat,si)
        end
    return kron_mat
end

function evolStabState(state::StabState,ustr::PauliStr)
    for stab in state.stablizers
        evolPauli(stab,ustr)
    end
end

function stabToVector(state::StabState)
    state_vec = rand(2^(state.nsites))
    for stab in state.stablizers
        proj = kronPauli(stab)
        
        for i in 1:2^(state.nsites)
            proj[i,i] = proj[i,i] + 1
        end

        state_vec = proj*state_vec
    end

    normalize!(state_vec)
    return state_vec
end

function expectVector(state_vec::Vector{ComplexF64},ostr::PauliStr)
    kron_ostr = kronPauli(ostr)
    return dot(state_vec,kron_ostr,state_vec)
end

function expectStabState(state::StabState,ostr::PauliStr)
    state_vec = stabToVector(state)
    return expectVector(state_vec,ostr)
end