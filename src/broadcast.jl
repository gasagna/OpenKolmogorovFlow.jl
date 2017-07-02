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
@generated function broadcast_c!(f, ::Type{FTField}, ::Type{FTField}, dest, src::Vararg{Any, N}) where N
    args = [:(unsafe_get(src[$k])) for k = 1:N]
    quote
        broadcast!(f, unsafe_get(dest), $(args...))
    end
end

# allocating version
@generated function broadcast_c(f, ::Type{FTField}, src::Vararg{Any, N}) where N
    args = [:(unsafe_get(src[$k])) for k = 1:N]
    quote
        FTField(broadcast(f, $(args...)))
    end
end

# extract the underlying matrix representation 
@inline Base.unsafe_get(U::Union{FTField, DiffOperator}) = U.data
