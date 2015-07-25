type ConstantExtrapolation{T,N,ITP,IT,GT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
end
ConstantExtrapolation{T,ITP,IT,GT}(::Type{T}, N, itp::ITP, ::Type{IT}, ::Type{GT}) =
    ConstantExtrapolation{T,N,ITP,IT,GT}(itp)

extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Flat}) =
    ConstantExtrapolation(T,N,itp,IT,GT)

function extrap_prep{T,N,ITP,IT}(etp::Type{ConstantExtrapolation{T,N,ITP,IT,OnGrid}}, xs...)
    :(@nexprs $N d->(x_d = clamp(xs[d], one($xs[d]), size(etp,d))))
end
function extrap_prep{T,N,ITP,IT}(etp::Type{ConstantExtrapolation{T,N,ITP,IT,OnCell}}, xs...)
    :(@nexprs $N d->(x_d = clamp(xs[d], 1//2, size(etp,d)+1//2)))
end
