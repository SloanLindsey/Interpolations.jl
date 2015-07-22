type ReflectingExtrapolation{T,N,ITP,IT,GT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
end
ReflectingExtrapolation{T,ITP,IT,GT}(::Type{T}, N, itp::ITP, ::Type{IT}, ::Type{GT}) =
    ReflectingExtrapolation{T,N,ITP,IT,GT}(itp)

extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Reflect}) =
    ReflectingExtrapolation(T,N,itp,IT,GT)

function extrap_prep{T,N,ITP,IT,GT}(etp::Type{ReflectingExtrapolation{T,N,ITP,IT,GT}}, xs...)
    # First, translate x into the domain over [1,2N-1) (OnGrid) or [.5,2N+.5) (OnCell)
    # i.e. into twice the size of the domain.
    # Next, if x is now in the upper part of this "double-domain", reflect over the middle
    # to obtain a new value x' for which f(x') == f(x), but where x' is inside the domain
    if GT == OnGrid
        return quote
            @nexprs $N d->begin
                x_d = mod(xs[d] - one($xs[d]), 2 * size(etp.itp,d) - 2) + one($xs[d])
                x_d > size(etp.itp,d) && (x_d = convert($xs[d], 2) * size(etp.itp,d) - x_d)
            end
        end
    else # GT == OnCell
        return quote
            @nexprs $N d->begin
                x_d = mod(xs[d] - 1//2, 2*size(etp.itp,d) + 1//2)
                x_d > size(etp.itp,d) + 1//2 && (x = convert($xs[d],2) * size(etp.itp,d) + one($xs[d] - x_d))
            end
        end
    end
end
