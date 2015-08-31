type ScaledInterpolation{T,N,ITPT,IT,GT,RT} <: AbstractInterpolationWrapper{T,N,ITPT,IT,GT}
    itp::ITPT
    ranges::RT
end
ScaledInterpolation{T,ITPT,IT,GT,RT}(::Type{T}, N, itp::ITPT, ::Type{IT}, ::Type{GT}, ranges::RT) =
    ScaledInterpolation{T,N,ITPT,IT,GT,RT}(itp, ranges)
function scale{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ranges::Range...)
    length(ranges) == N || error(string("Must scale $N-dimensional interpolation object with exactly $N ranges (you used $(length(ranges)))"))
    ScaledInterpolation(T,N,itp,IT,GT,ranges)
end

@generated function getindex{T,N}(sitp::ScaledInterpolation{T,N}, xs...)
    length(xs) == N || error(string("Must index into ScaledInterpolation{", N, "} with exactly ", N, " indices (you used ", length(xs), ")"))
    :(getindex(sitp.itp, coordlookup(sitp.ranges, tuple(xs...))...))
end

size(sitp::ScaledInterpolation, d) = size(sitp.itp, d)

coordlookup(r::FloatRange, x) = (r.divisor * x - r.start) ./ r.step + one(eltype(r))
coordlookup(r::StepRange, x) = (x - r.start) ./ r.step + one(eltype(r))
coordlookup(r::UnitRange, x) = x - r.start + one(eltype(r))
coordlookup{N}(r::NTuple{N}, x::NTuple{N}) = map(coordlookup,r,x)

