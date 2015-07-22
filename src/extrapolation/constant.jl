type ConstantExtrapolation{T,N,ITP,IT,GT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
end
ConstantExtrapolation{T,ITP,IT,GT}(::Type{T}, N, itp::ITP, ::Type{IT}, ::Type{GT}) =
    ConstantExtrapolation{T,N,ITP,IT,GT}(itp)

extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Flat}) =
    ConstantExtrapolation(T,N,itp,IT,GT)

function extrap_prep{T,ITP,IT}(etp::Type{ConstantExtrapolation{T,1,ITP,IT,OnGrid}}, x)
    :(x = clamp(x, 1, size(etp,1)))
end
function extrap_prep{T,ITP,IT}(etp::Type{ConstantExtrapolation{T,1,ITP,IT,OnCell}}, x)
    :(x = clamp(x, .5, size(etp,1)+.5))
end
function extrap_prep{T,N,ITP,IT}(etp::Type{ConstantExtrapolation{T,N,ITP,IT,OnGrid}}, xs...)
    :(@nexprs $N d->(xs[d] = clamp(xs[d], 1, size(etp,d))))
end
function extrap_prep{T,N,ITP,IT}(etp::Type{ConstantExtrapolation{T,N,ITP,IT,OnCell}}, xs...)
    :(@nexprs $N d->(xs[d] = clamp(xs[d], .5, size(etp,d)+.5)))
end
