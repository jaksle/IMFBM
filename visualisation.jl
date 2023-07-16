# Jakub Ślęzak
# Hugo Steinhaus Center, Wrocław University of Science and Technology
# 
# Ralf Metzler
# Institute of Physics & Astronomy, Potsdam Univeristy

# File visualisation.jl is a part of IMFBM package
# IMFBM is free package: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

using Plots

## single trajectory, fig. 2

dt, m = 0.001, 1000
T = dt*m
ts = dt * (1:m)

X = simulateIMFBM(t->1 - t/T,t->1,dt,m,1)
n, j = 6, 6

dx = 0.1/n
dy = 8dx
x,y = dx,-dy
plot(
    xlim=(-dx/2,1),
    axisratio=1/8,
    framestyle = :box,
)
for x in 0:dx:1
    vals = X[x .< ts .< x+dx]
    isempty(vals) && continue 
    ys = -4:dy:4
    listYs = [y for y in ys if y+dy > minimum(vals) && y < maximum(vals)]
    for y in listYs
        plot!(Shape([x,x+dx,x+dx,x],[y, y, y + dy, y+ dy]),
            color = :white,
            linewidth = 1,
            linecolor = palette(:greys)[220],
            label = "",
        )
    end
end
plot!(ts, X,
    fontfamily = "Computer Modern",
    axisratio=1/8,
    linez  = LinRange(1,0,length(ts)),
    c = cgrad(palette(:turbo)[end-25:-1:25]),
    colorbar = true,
    linewidth = 1.5,
    colorbar_title = "Hurst exponent",
    label = "",
    xlabel = "time",
    ylabel = "displacement"
)

## 3 trajectories comparison, fig. 1

H1, H2 = 0.5, 0.3
D1, D2 = 1, 1
τ = 5 
ta, tb = 0.001, 10

C1 = covIMFBM(H1, H2, D1, D2, τ)

H1, H2 = 0.5, 0.7
D1, D2 = 1, 1
τ = 5 
ta, tb = 0.001, 10

C2 = covIMFBM(H1, H2, D1, D2, τ)

m = 10^4 # no of points
ts = LinRange(ta,tb,m)
S = [C1(s,t) for s in ts, t in ts]
A1 = cholesky(Symmetric(S)).U
S = [C2(s,t) for s in ts, t in ts]
A2 = cholesky(Symmetric(S)).U

n = 10 # sample

dX = randn(m, n) 

X = cumsum(dX,dims=1) .* sqrt(step(ts)) # Brownian motion
Y1 = A1'*dX
Y2 = A2'*dX

j = 5
cols = cgrad(palette(:turbo)[end-25:-1:25])

Plots.scalefontsizes()
plot(ts[:],X[:,j],
    xlabel = "time",
    ylabel = "displacement",
    fontfamily = "Computer Modern",
    xticks = 0:10,
    legend = :topleft,
    label = "H = 0.5",
    xlim = (0,10),
    framestyle = :box,
    linewidth = 0.5,
    c = :green, 
)

plot!(ts[end÷2:end],Y1[end÷2:end,j],
    c = palette(:default)[2],
    linewidth = 0.5,
    label = "H = 0.3",

)
plot!(ts[end÷2:end],Y2[end÷2:end,j],
    c = palette(:default)[1],
    linewidth = 0.5,
    label = "H = 0.7",
)

## spectrum plot, fig. 3 left

using FourierAnalysis, LaTeXStrings, Statistics

H1, H2 = 0.3, 0.697
D1, D2 = 1, 16
τ = 5
m, n = 1000, 1000
dt = 10/n # T = 10
X = simulateIMFBM(H1,H2,D1,D2,τ,dt,m,n)

sr, wl= 64, 256
s = 1 # choose part of trajectory 1 left
S = spectra(diff(X[s:s+500,:],dims=1),sr,wl)

c1,c2,c3 = palette(:tab10)[1],:black,palette(:tab10)[1]
w = (floor(Int,100*(s-1)/499))
plot(S.flabels,S,
    fontfamily = "Computer Modern",
    framestyle = :box,
    scale=:log10,
    title = "$(100-w)% left, $(w)% right",
    color = c1,
    alpha = 0.05,
    xlim = (0.25,8),
    xlabel = "frequency [Hz]",
    ylabel = "spectral density",
    yticks = ([2^-16,2^-12,2^-8],[L"2^{-16}",L"2^{-12}",L"2^{-8}"]),
    xticks = ([0.25,0.5,1,2,4,8],["1/4","1/2","1","2","4","8"]),
)
plot!(S.flabels,mean(S.y,dims=2),
    color=c2,
    linewidth = 3
)
q2 = [quantile(x,0.95) for x in eachrow(S.y)]
q1 = [quantile(x,0.05) for x in eachrow(S.y)]
plot!(S.flabels,q1,
    color = c3,
    linewidth = 2,
    linestyle = :dash,
)
plot!(S.flabels,q2,
    color = c3,
    linewidth = 2,
    linestyle = :dash,
)

plot!(f->1.2*2^-11*f^(1-2*0.3),S.flabels[1],S.flabels[end],
    color=:red,
    linestyle = :dot,
    linewidth = 2,
)


## spectrum plot, fig. 3 right

sr, wl= 64, 256
s = 500
S = spectra(diff(X[s:s+500,:],dims=1),sr,wl)

c1,c2,c3 = palette(:tab10)[1],:black,palette(:tab10)[1]
w = (floor(Int,100*(s-1)/499))
plot(S.flabels,S,
    fontfamily = "Computer Modern",
    framestyle = :box,
    scale=:log10,
    title = "$(100-w)% left, $(w)% right",
    color = c1,
    alpha = 0.05,
    xlim = (0.25,8),
    xlabel = "frequency [Hz]",
    ylabel = "spectral density",
    ylim = (2^-16,2^-6),
    yticks = ([2^-16,2^-12,2^-8],[L"2^{-16}",L"2^{-12}",L"2^{-8}"]),
    xticks = ([0.25,0.5,1,2,4,8],["1/4","1/2","1","2","4","8"]),
    label = "",
)
plot!([],[],
    color = c1,
    alpha = 0.3,
    label = "trejectory spectra"
)
plot!(S.flabels,mean(S.y,dims=2),
    color=c2,
    linewidth = 3,
    label = "mean spectrum"
)
q2 = [quantile(x,0.95) for x in eachrow(S.y)]
q1 = [quantile(x,0.05) for x in eachrow(S.y)]
plot!(S.flabels,q1,
    color = c3,
    linewidth = 2,
    linestyle = :dash,
    label = "",
)
plot!(S.flabels,q2,
    color = c3,
    linewidth = 2,
    linestyle = :dash,
    label = "90% quantiles"
)

plot!(f->4.5*2^-11*f^(1-2*0.7),S.flabels[1],S.flabels[end],
    color=:red,
    linestyle = :dot,
    linewidth = 2,
    label = L"power law $f^{1-2H}$",
    legend = :bottomleft
)



## spectrum, fig 3 center

sr, wl= 64, 256
s = 251
S = spectra(diff(X[s:s+500,:],dims=1),sr,wl)

c1,c2,c3 = palette(:tab10)[1],:black,palette(:tab10)[1]
w = (floor(Int,100*(s-1)/499))
plot(S.flabels,S,
    fontfamily = "Computer Modern",
    framestyle = :box,
    scale=:log10,
    title = "$(100-w)% left, $(w)% right",
    color = c1,
    alpha = 0.05,
    xlim = (0.25,8),
    xlabel = "frequency [Hz]",
    ylabel = "power spectrum",
    ylim = (2^-16,2^-6),
    yticks = ([2^-16,2^-12,2^-8],[L"2^{-16}",L"2^{-12}",L"2^{-8}"]),
    xticks = ([0.25,0.5,1,2,4,8],["1/4","1/2","1","2","4","8"]),
)
plot!(S.flabels,mean(S.y,dims=2),
    color=c2,
    linewidth = 3
)
q2 = [quantile(x,0.95) for x in eachrow(S.y)]
q1 = [quantile(x,0.05) for x in eachrow(S.y)]
plot!(S.flabels,q1,
    color = c3,
    linewidth = 2,
    linestyle = :dash,
)
plot!(S.flabels,q2,
    color = c3,
    linewidth = 2,
    linestyle = :dash,
)

plot!(f->1.7*2^-11*f^(1-2*0.7),S.flabels[1],1,
    color=:red,
    linestyle = :dot,
    linewidth = 2,
)
plot!(f->1.2*2^-11*f^(1-2*0.3),2,8,
    color=:red,
    linestyle = :dot,
    linewidth = 2,
)


