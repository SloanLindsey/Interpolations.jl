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

end
