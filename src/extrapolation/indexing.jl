@generated function getindex{T,N}(etp::AbstractExtrapolation{T,N}, xs...)
    quote
        $(extrap_prep(etp, xs...))
        itp = etp.itp
        @nref $N itp x
    end
end
