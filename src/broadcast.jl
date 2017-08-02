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

# Extract underlying matrices in the fields and operators and then call broadcast for the arrays. 
@generated function Base.Broadcast.broadcast_c!(f, ::Type{FTField}, ::Type{FTField}, dest, args::Vararg{Any, M}) where {M}
    quote
        $(Expr(:meta, :inline))
        broadcast!(f, unsafe_get(dest), map(unsafe_get, args)...)
        return dest
    end
end

# Allocating version. This is called when you need to create a new object e.g.
# u = v .* 3. This creates a new object and returns it.
@generated function broadcast_c(f, ::Type{FTField}, src::Vararg{Any, N}) where N
    args = [:(unsafe_get(src[$k])) for k = 1:N]
    quote
        FTField(broadcast(f, $(args...)))
    end
end

# extract the underlying matrix representation 
@inline Base.unsafe_get(U::Union{FTField, DiffOperator}) = U.data
