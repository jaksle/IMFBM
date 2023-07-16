# Jakub Ślęzak
# Hugo Steinhaus Center, Wrocław University of Science and Technology
# 
# Ralf Metzler
# Institute of Physics & Astronomy, Potsdam Univeristy

# File estimation.jl is a part of IMFBM package
# IMFBM is free package: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.




#### estimation method for single trajectory based on variation and covariation

## fitting variation and covariation polygonal chains

using LsqFit, Optim

T = 10
m = 1000
dt = T/m

X = simulateIMFBM(0.3,0.4,2,1,5, dt, m, 1) # exemplary data

dX = diff(X,dims=1)

cv = dX[1:end-1] .* dX[2:end]
sq = dX[1:end-1] .^ 2

𝒞 = cumsum(cv)
𝒱 = cumsum(sq)



l2l(x,a,τ,c) = x < τ ? a*x : c*x + a*τ - c*τ # fitted polygonal chain

ts = dt:dt:T
l = m-2


a0 = 𝒞[end÷2]/l*2
b0 = l/2
c0 = (𝒞[end] - 𝒞[end÷2])/l*2
M(x,p) = l2l.(x,p[1],p[2],p[3])
fit = curve_fit(M,1:l,𝒞,[a0, b0, c0])



a1 = 𝒱[end÷2]/l*2
b1 = l/2
c1 = (𝒱[end] - 𝒱[end÷2])/l*2
fit2 = curve_fit(M,1:l,𝒱,[a1, b1, c1])

fopt(p) = sum( (𝒞[k]- l2l(k,p[1],p[3],p[2]))^2 +(𝒱[k]-l2l(k,p[4],p[3],p[5]))^2   for k in 1:l)
fitO = optimize(fopt,[a0,c0,b0,a1,c1])
par = Optim.minimizer(fitO)


## visualisation of the fitted variation and covariation
using Plots

plot(1:l,-𝒞,
    label = "estimated -covarition"
)
plot!(1:l,-l2l.(1:l,par[1],par[3],par[2]),
    label = "fitted - (polygonal chain)"
)

plot!(1:l,𝒱,
    label = "estimated covariation"
)
plot!(1:l,l2l.(1:l, par[4],par[3],par[5]),
    label = "fitted polygonal chain"
)



## final estimation

ae,ce,τe,ae2,ce2 = par

x = ae/ae2
estH1 = log2(2x+2)/2
estD1 = ae2/dt^(2estH1)
x2 = ce/ce2
estH2 = log2(2x2+2)/2
estD2 = ce2/dt^(2estH2)

println("H₁ ≈ $estH1")
println("H₂ ≈ $estH2")
println("D₁ ≈ $estD1")
println("D₂ ≈ $estD2")
println("τ ≈ $(T*τe/m)")



#### estimation method for sample of trajectories based on local msd

## generating exemplary data

T = 1
m, n = 1000, 100
dt = T/m
ts = dt:dt:T

X = simulateIMFBM(t->0.8 - 0.6t,t->1+t/2,dt,m,n)
dX = diff(X,dims=1)
m -= 1

## estimation of local msd

using Statistics
w = 10
Hs = Vector{Float64}(undef,size(1:m-w))
Ds = similar(Hs)
M(t,p) = p[1] .* (step(ts) .* t) .^ (2p[2])
for k in 1:m-w
     msd = mean(cumsum(dX[k:k+w-1,:],dims=1) .^2,dims=2)
     l = log.(msd)
     t = log.(1:w)
     X = [ones(10) t]
     b = (X'*X)\X'*l
     Hs[k] = b[2]/2
     Ds[k] = exp(b[1])
end

Ds ./= dt .^ (2Hs)

## visualisation of the results
using Plots

plot(dt:dt:length(Hs)*dt,Hs,
    label = "estimated H"
)
plot!(t->0.8 - 0.6t,0,1,
    label = "estimated local H"
)

plot(dt:dt:length(Ds)*dt,Ds,
    label = "estimated D",
    yscale = :log10
)
plot!(t->1+t/2,0,1,
    label = "estimated local D"
)