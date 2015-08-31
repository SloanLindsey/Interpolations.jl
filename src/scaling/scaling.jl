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

gradient{T,N}(sitp::ScaledInterpolation{T,N}, xs...) = gradient!(Array(T,N), sitp, xs...)
@generated function gradient!{T,N}(g, sitp::ScaledInterpolation{T,N}, xs...)
    ndims(g) == 1 || error(string("g must be a vector (but had ", ndims(g), " dimensions)"))
    length(xs) == N || error(string("Must index into ScaledInterpolation{", N, "} with exactly ", N, " indices (you used ", length(xs), ")"))

    quote
        length(g) == N || error(string("g must be a vector of length ", N, " (was ", length(g), ")"))
        gradient!(g, sitp.itp, coordlookup(sitp.ranges, tuple(xs...))...)
        for i in eachindex(g)
            g[i] = rescale_gradient(sitp.ranges[i], g[i])
        end
        g
    end
end

rescale_gradient(r::FloatRange, g) = g * r.divisor / r.step
rescale_gradient(r::StepRange, g) = g / r.step
rescale_gradient(r::UnitRange, g) = g