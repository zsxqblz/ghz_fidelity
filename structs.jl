const s0 = ComplexF64[1 0 ; 0 1]
const sx = ComplexF64[0 1 ; 1 0]
const sy = ComplexF64[0 -1im ; 1im 0]
const sz = ComplexF64[1 0 ; 0 -1]

struct PauliStr
    idx_active::Vector{Int64}
    nsites::Int64
    pauli_op::Vector{Matrix{ComplexF64}}
    function PauliStr(idx_active::Vector{Int64},pauli_op::Vector{Matrix{ComplexF64}})
        nsites = length(pauli_op)
        new(idx_active,nsites,pauli_op)
    end
end

struct StabState
    stablizers::Vector{PauliStr}
    nstabl::Int64
    nsites::Int64
    function StabState(stablizers::Vector{PauliStr})
        nstabl = length(stablizers)
        nsites = stablizers[1].nsites
        new(stablizers,nstabl,nsites)
    end
end