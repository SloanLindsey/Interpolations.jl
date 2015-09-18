type ScaledInterpolation{T,N,ITPT,IT,GT,RT} <: AbstractInterpolationWrapper{T,N,ITPT,IT,GT}
    itp::ITPT
    ranges::RT
end
ScaledInterpolation{T,ITPT,IT,GT,RT}(::Type{T}, N, itp::ITPT, ::Type{IT}, ::Type{GT}, ranges::RT) =
    ScaledInterpolation{T,N,ITPT,IT,GT,RT}(itp, ranges)
function scale{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ranges::Range...)
    length(ranges) == N || error("Must scale $N-dimensional interpolation object with exactly $N ranges (you used $(length(ranges)))")
    for d in 1:N
        if iextract(IT,d) != NoInterp 
            length(ranges[d]) == size(itp,d) || error("The length of the range in dimension $d ($(length(ranges[d]))) did not equal the size of the interpolation object in that direction ($(size(itp,d)))")
        elseif ranges[d] != 1:size(itp,d)
            error("NoInterp dimension $d must be scaled with unit range 1:$(size(itp,d))")
        end
    end

    ScaledInterpolation(T,N,itp,IT,GT,ranges)
end

@generated function getindex{T,N,ITPT,IT<:DimSpec}(sitp::ScaledInterpolation{T,N,ITPT,IT}, xs...)
    length(xs) == N || throw(ArgumentError("Must index into $N-dimensional scaled interpolation object with exactly $N indices (you used $(length(xs)))"))
    interp_types = length(IT.parameters) == N ? IT.parameters : tuple([IT.parameters[1] for _ in 1:N]...)
    interp_dimens = map(it -> interp_types[it] != NoInterp, 1:N)
    interp_indices = map(i -> interp_dimens[i] ? :(coordlookup(sitp.ranges[$i], xs[$i])) : :(xs[$i]), 1:N)
    return :(getindex(sitp.itp, $(interp_indices...)))
end

size(sitp::ScaledInterpolation, d) = size(sitp.itp, d)

coordlookup(r::FloatRange, x) = (r.divisor * x - r.start) / r.step + one(eltype(r))
coordlookup(r::StepRange, x) = (x - r.start) / r.step + one(eltype(r))
coordlookup(r::UnitRange, x) = x - r.start + one(eltype(r))
coordlookup(i::Bool, r::Range, x) = i ? coordlookup(r, x) : convert(typeof(coordlookup(r,x)), x)

gradient{T,N}(sitp::ScaledInterpolation{T,N}, xs...) = gradient!(Array(T,N), sitp, xs...)
@generated function gradient!{T,N,ITPT,IT}(g, sitp::ScaledInterpolation{T,N,ITPT,IT}, xs...)
    ndims(g) == 1 || throw(ArgumentError("g must be a vector (but had $(ndims(g)) dimensions)"))
    length(xs) == count_interp_dims(IT, N) || throw(ArgumentError("Must index into $N-dimensional scaled interpolation object with exactly $N indices (you used $(length(xs)))"))

    interp_types = length(IT.parameters) == N ? IT.parameters : tuple([IT.parameters[1] for _ in 1:N]...)
    interp_dimens = map(it -> interp_types[it] != NoInterp, 1:N)
    interp_indices = map(i -> interp_dimens[i] ? :(coordlookup(sitp.ranges[$i], xs[$i])) : :(xs[$i]), 1:N)

    quote
        length(g) == N || error(string("g must be a vector of length ", N, " (was ", length(g), ")"))
        gradient!(g, sitp.itp, $(interp_indices...))
        for i in eachindex(g)
            g[i] = rescale_gradient(sitp.ranges[i], g[i])
        end
        g
    end
end

rescale_gradient(r::FloatRange, g) = g * r.divisor / r.step
rescale_gradient(r::StepRange, g) = g / r.step
rescale_gradient(r::UnitRange, g) = g