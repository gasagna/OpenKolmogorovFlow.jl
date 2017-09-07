# DEFINE STUFF NECESSARY TO USE . NOTATION
import Base.Broadcast: _containertype,
                       promote_containertype,
                       broadcast_c!,
                       broadcast_c

# when you are broadcasting over FTField and DiffOperator always go down to FTField
_containertype(::Type{<:DiffOperator}) = FTField

# extract the underlying matrix representation
@inline Base.unsafe_get(U::Union{FTField, Field, DiffOperator}) = U.data

for T in (:FTField, :Field)
    @eval begin
        _containertype(::Type{<:$T}) = $T

        # These are need for scalars and arrays
        promote_containertype(::Type{$T}, ::Type{Any})   = $T
        promote_containertype(::Type{Any}, ::Type{$T})   = $T
        promote_containertype(::Type{$T}, ::Type{Array}) = $T
        promote_containertype(::Type{Array}, ::Type{$T}) = $T

        # Extract underlying matrices in the fields and operators and then call broadcast for the arrays.
        @generated function Base.Broadcast.broadcast_c!(f, ::Type{$T}, ::Type{$T}, dest, args::Vararg{Any, M}) where {M}
            quote
                $(Expr(:meta, :inline))
                broadcast!(f, unsafe_get(dest), map(unsafe_get, args)...)
                return dest
            end
        end

        # Allocating version. This is called when you need to create a new object e.g.
        # u = v .* 3. This creates a new object and returns it.
        @generated function broadcast_c(f, ::Type{$T}, src::Vararg{Any, N}) where N
            args = [:(unsafe_get(src[$k])) for k = 1:N]
            quote
                $($T)(broadcast(f, $(args...)))
            end
        end
    end
end