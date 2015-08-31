type LinearExtrapolation{T,N,ITP,IT,GT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
end
LinearExtrapolation{T,ITP,IT,GT}(::Type{T}, N, itp::ITP, ::Type{IT}, ::Type{GT}) =
    LinearExtrapolation{T,N,ITP,IT,GT}(itp)

extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Linear}) =
    LinearExtrapolation(T,N,itp,IT,GT)

function extrap_prep{T,N,ITP,IT}(etp::Type{LinearExtrapolation{T,N,ITP,IT,OnGrid}}, xs...)
    quote
        xc = collect(xs...)
        @nexprs $N d->begin
            if xc[d] < 1
                xc[d] = one(eltype(xc))
                k = gradient(etp.itp, xc...)[d]
                return etp.itp[xc...] + k * (xs[d] - one(typeof(xc[d])))
            elseif xc[d] > size(etp.itp, d)
                xc[d] = size(etp.itp, d)
                k = gradient(etp.itp, xc...)[d]
                return etp.itp[xc...] + k * (xs[d] - size(etp.itp, d))
            end
        end
    end
end

function extrap_prep{T,N,ITP,IT}(etp::Type{LinearExtrapolation{T,N,ITP,IT,OnCell}}, xs...)
    quote
        xc = collect(xs...)
        @nexprs $N d->begin
            if xc[d] < 1//2
                xc[d] = convert(typeof(xc[d], 1//2))
                k = gradient(etp.itp, xc...)[d]
                return etp.itp[xc...] + k * (xs[d] - convert(typeof(xc[d], 1//2)))
            elseif xc[d] > size(etp.itp, d) + 1//2
                xc[d] = size(etp.itp, d) + convert(typeof(xc[d], 1//2))
                k = gradient(etp.itp, xc...)[d]
                return etp.itp[xc...] + k * (xs[d] - size(etp.itp, d) - convert(typeof(xc[d], 1//2)))
            end
        end
        @nexprs $N d->(x_d = xs[d])
    end
end
