type PeriodicExtrapolation{T,N,ITP,IT,GT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
end
PeriodicExtrapolation{T,ITP,IT,GT}(::Type{T}, N, itp::ITP, ::Type{IT}, ::Type{GT}) =
    PeriodicExtrapolation{T,N,ITP,IT,GT}(itp)

extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Periodic}) =
    PeriodicExtrapolation(T,N,itp,IT,GT)

function extrap_prep{T,N,ITP,IT,GT}(etp::Type{PeriodicExtrapolation{T,N,ITP,IT,GT}}, xs...)
    :(@nexprs $N d->(; x_d = mod(xs[d]-one($xs[d]), size(etp.itp,d)-one($xs[d]))+one($xs[d])))
end
