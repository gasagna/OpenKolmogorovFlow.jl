# DEFINE STUFF NECESSARY TO USE . NOTATION
import Base.Broadcast: _containertype, 
                       promote_containertype, 
                       broadcast_c!, 
                       broadcast_c

# when you are broadcasting over FTField and DiffOperator always go down to FTField
_containertype(::Type{<:FTField}) = FTField
_containertype(::Type{<:DiffOperator}) = FTField

# These are need for scalars and arrays
promote_containertype(::Type{FTField}, ::Type{Any})   = FTField
promote_containertype(::Type{Any}, ::Type{FTField})   = FTField
promote_containertype(::Type{FTField}, ::Type{Array}) = FTField
promote_containertype(::Type{Array}, ::Type{FTField}) = FTField

# extract underlying matrices and broadcast over that
@generated function broadcast_c!(f, ::Type{FTField}, ::Type{FTField}, dest, src...)
    rest = Expr(:tuple, ntuple(k -> :(unsafe_get(src[$k])), length(src))...)
    quote
        broadcast!(f, unsafe_get(dest), $rest...)
    end
end

# allocating version
@generated function broadcast_c(f, ::Type{FTField}, src...)
    rest = Expr(:tuple, ntuple(k -> :(unsafe_get(src[$k])), length(src))...)
    quote
        FTField(broadcast(f, $rest...))
    end
end

# extract the underlying matrix representation 
@inline Base.unsafe_get(U::FTField) = U.data
@inline Base.unsafe_get(op::DiffOperator) = op.data
