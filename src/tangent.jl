struct TangentEquation{n, T<:AbstractFloat64}
    tangImTerm::TangentEquationImTerm{n}
    tangExTerm::TangentEquationExTerm{n}
end

struct TangentEquationExTerm{n}
end

struct AugmentedState{T}
    x₁::T
    x₂::T
end

function (Eq::TangentEquationExTerm{n})(t::Real, 
                                        η::AugmentedState{FTField{n}}, 
                                        η̇::AugmentedState{FTField{n}}) where {n}
    
end

