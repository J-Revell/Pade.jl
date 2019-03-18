module Pade
    using SymPy
    using Main.Wynn

    # data structure to store the pade table
    # R[l/m] indexed by ints l & m, stored within dict.
    struct PadeTable{T}
        series::T
        terms::Vector{T}
        ptable::Dict{Tuple{Int,Int},T}
    end

    # filter function, when epsilon j-value is even
    flt(tuple::Tuple{Int,Int}) = iseven(tuple[2]) && (tuple[2] >= 0) &&
    (tuple[1]+tuple[2]/2 >= 0)

    function PadeTable(terms::Vector{T}; simplified::Bool = true) where T<:Union{Real,Sym}
        # compute epsilon algorithm coefficients
        series = sum(terms)
        etable = EpsilonTable(terms, simplified = simplified).etable

        # filter away auxiliary indexes
        idxs = Iterators.filter(flt, keys(etable))
        PadeTable(series, terms, Dict((Int(i+j/2),Int(j/2)) => etable[i,j] for (i,j) in idxs))
    end

    function pade(series::Vector{T}, l::Int, m::Int) where T<:Union{Real,Sym}
        if length(series) >= l+m+1
            return PadeTable(series[1:l+m+1]).ptable[l,m]
        else
            @error "A Pade approximation of order [$l/$m] requires at least $(l+m+1) series terms."
        end
    end

    export PadeTable, pade
end
