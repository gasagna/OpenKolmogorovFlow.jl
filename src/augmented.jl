import DualNumbers: Dual, value, epsilon

export AugmentedFTField, AugmentedField, state, prime

# ~~ VERSIONS OF FTFIELD AND FIELD WITH PERTURBATION OR ADJOINT VARIABLE ~~~
for (FT, numtype) in [(:FTField, Complex), (:Field, Real)]
    typename = Symbol(:Augmented, FT)
    abstractname = Symbol(:Abstract, FT)
    @eval begin
        # this is not the strictest typing but can't do otherwise
        # struct AugmentedFTField{n, T<:Dual{<:Complex}, F<:FTField{n, <:Complex}} <: AbstractFTField{n, T}
        # struct AugmentedField{n,   T<:Dual{<:Real},    F<:Field{n,   <:Real}}    <: AbstractField{n, T}
        struct $typename{n, T<:Dual{<:$numtype}, F<:$FT{n, <:$numtype}} <: $abstractname{n, T}
            x::F # state
            p::F # perturbation
        end

        # Constructors
        $typename(n::Int) = $typename(n, Dual{$(numtype == Complex ? Complex128 : Float64)})
        $typename(n::Int, ::Type{Dual{T}}) where {T} = $typename($FT(n, T), $FT(n, T))
        $typename(U::$FT{n, T}, u::$FT{n, T}) where {n, T} = $typename{n, Dual{T}, typeof(U)}(U, u)

        # accessors
        state(U::$typename) = U.x
        prime(U::$typename) = U.p
        state(U::$FT) = U # no op for FTField and Field

        # ~~ Array interface ~~
        Base.@propagate_inbounds @inline Base.getindex(U::$typename, i...) =
            Dual(U.x[i...], U.p[i...])
        Base.@propagate_inbounds @inline Base.setindex!(U::$typename, v::Dual, i...) = 
            (U.x[i...] = value(v); U.p[i...] = epsilon(v); v)
        Base.@propagate_inbounds @inline Base.setindex!(U::$typename, v::Real, i...) = 
            (U.x[i...] = v; v)

        Base.similar(U::$typename{n}, m::Int=n) where {n} = $typename(similar(U.x, m), similar(U.p, m))
    end
end