import MacroTools: @capture, postwalk

# Extract the symbol used for indexing. Example:
#
#  _findindex(:(U[i]*V[i])) -> :i
#
# It is an error if two different indices are used
function _findindex(expr)
    indices = Set{Symbol}()
    postwalk(expr) do x
        if @capture(x, U_[i_])
            push!(indices, i)
        else
            x
        end
    end
    # index must be unique
    length(indices) == 0 &&
        throw(ArgumentError("no indexing found in expression"))
    length(indices) >  1 &&
        throw(ArgumentError("no unique index found in expression"))
    return pop!(indices)
end

# macro for sums in wave number space
macro Σ_jk(n, expr)
    # capture the running index
    i = _findindex(expr)
    return quote
        @inbounds begin
            $(esc(i)) = 1
            local d = $(esc(n)) >> 1
            # initialise sum variables
            Σ₁ = zero($(esc(expr)))
            Σ₂ = Σ₁
            # count column j = 0 with weight 1
            # skip the (0, 0) mode, assumed zero
            for $(esc(i)) = 2:$(esc(n))
                Σ₁ += $(esc(expr))
            end
            # count columns j=1:d with weight 2
            for $(esc(i)) = ($(esc(n))+1):($(esc(n))*d)
                Σ₂ += $(esc(expr))
            end
            # count column j=d+1 with weight 1
            for $(esc(i)) = ($(esc(n))*d+1):($(esc(n))*d + $(esc(n)))
                Σ₁ += $(esc(expr))
            end
        end
        # sum contributions
        bulk = 2*Σ₂ + Σ₁
        # must remove the contributions at the Nyquist
        # frequencies counted as exp^0 rather than cos^2
        for $(esc(i)) in (d+1, d*$(esc(n))+1, d*$(esc(n))+1+d)
            bulk -= 0.5*$(esc(expr))
        end
        bulk
    end
end