module ScalingTests

using Interpolations
using Base.Test

# Model linear interpolation of y = -3 + .5x by interpolating y=x
# and then scaling to the new x range

itp = interpolate(1:1.0:10, BSpline(Linear), OnGrid)

sitp = scale(itp, -3:.5:1.5)

for (x,y) in zip(-3:.05:1.5, 1:.1:10)
    @test_approx_eq sitp[x] y
end

# Verify that it works in >1D, with different types of ranges

gauss(phi, mu, sigma) = exp(-(phi-mu)^2 / (2sigma)^2)
f(x,y) = gauss(x, 0.5, 4) * gauss(y, -.5, 2)

xs = -5:.5:5
ys = -4:.2:4
zs = Float64[f(x,y) for x in xs, y in ys]

itp2 = interpolate(zs, BSpline(Quadratic(Flat)), OnGrid)
sitp2 = scale(itp2, xs, ys)

for x in xs, y in ys
    @test_approx_eq f(x,y) sitp2[x,y]
end

# Test gradients of scaled grids
xs = -pi:.1:pi
ys = sin(xs)
itp = interpolate(ys, BSpline(Linear), OnGrid)
sitp = scale(itp, xs)

for x in -pi:.1:pi
	g = gradient(sitp, x)[1]
	@test_approx_eq_eps cos(x) g .05
end

end

module ScalingDimspecTests

using Interpolations, DualNumbers, Base.Test

xs = -pi:(2pi/10):pi
ys = -2:.1:2
f(x,y) = sin(x) * y^2

itp = interpolate(Float64[f(x,y) for x in xs, y in ys], Tuple{BSpline(Quadratic(Periodic)), BSpline(Linear)}, OnGrid)
sitp = scale(itp, xs, ys)

# Don't test too near the edges until issue #64 is resolved
for (ix,x0) in enumerate(xs[5:end-5]), (iy,y0) in enumerate(ys[2:end-1])
    x, y = x0 + 2pi/20, y0 + .05
    @test_approx_eq sitp[x0, y0] f(x0,y0)
    @test_approx_eq_eps sitp[x0, y0] f(x0,y0) 0.05

    g = gradient(sitp, x, y)
    fx = epsilon(f(dual(x,1), y))
    fy = (f(x, ys[iy+2]) - f(x, ys[iy+1])) / (ys[iy+2] - ys[iy+1])

    @test_approx_eq_eps g[1] fx 0.15
    @test_approx_eq_eps g[2] fy 0.05 # gradients for linear interpolation is "easy"
end

end