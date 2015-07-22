@generated function getindex{T,N,ITP,GT}(etp::AbstractExtrapolation{T,N,ITP,GT}, xs...)
    quote
        $(extrap_prep(etp, xs...))
        itp = etp.itp
        @nref $N itp x
    end
end
