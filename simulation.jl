# Jakub Ślęzak
# Hugo Steinhaus Center, Wrocław University of Science and Technology
# 
# Ralf Metzler
# Institute of Physics & Astronomy, Potsdam Univeristy

# File simulation.jl is a part of IMFBM package
# IMFBM is free package: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.


using LinearAlgebra, SpecialFunctions

"""
    simulateIMFBM(H::Function, D::Function, dt, m, n)

For given local Hurst exponent H and local diffusivity D functions
returns `n` trajectories of IMFBM simulated at points dt, 2dt, ..., m*dt.

Requirements: 0 < H(t) < 1;    0 < D(t).

This is approximate IMFBM simulation scheme which converges to the true process as dt → 0.
"""
function simulateIMFBM(H::Function, D::Function, dt::Real, m::Int, n::Int) 

    ts = dt:dt:m*dt
    C = covIMFBM(H,D,dt)

    S = [C(s,t) for s in ts, t in ts]
    A = cholesky(Symmetric(S)).U
   
    ξ = randn(m, n)
    dX = A' * ξ
    return cumsum(dX, dims = 1)
end



"""
    simulateIMFBM(D₁, D₂, H₁, H₂, τ, ts, n)
    simulateIMFBM(H₁, H₂, D₁, D₂, τ, dt, m, n)

For given local Hurst exponents H and local diffusivity D which change from (D₁,H₁) to (D₂,H₂) at time τ
returns `n` trajectories of IMFBM simulated at points ts[1], ts[2], ...

Instead of ts timestep 0 < dt and number of points m::Int can be given.

Requirements: 0 < H₁, H₂ < 1;    0 < D₁, D₂;    0 < τ;    0 < ts[k].
Be wary of the condtion 0 < ts[k], 0 == ts[k] leads to covariance matrix not being positive definite.
Directly add value X(0) = 0 instead.

"""
function simulateIMFBM(H₁::Real, H₂::Real, D₁::Real, D₂::Real, τ::Real, ts::AbstractVector{<:Real}, n::Int) 

    m = length(ts)
    C = covIMFBM(H₁, H₂, D₁, D₂, τ)

    S = [C(s,t) for s in ts, t in ts]
    A = cholesky(Symmetric(S)).U
   
    ξ = randn(m, n)
    X = A' * ξ
    return X
end

simulateIMFBM(H₁::Real, H₂::Real, D₁::Real, D₂::Real, τ::Real, dt::Real, m::Int, n::Int) =
    simulateIMFBM(H₁, H₂, D₁, D₂, τ, dt:dt:m*dt, n::Int) 


_gammaFact(H) = H ≈ 1/2 ? -pi/2 : cospi(H)*gamma(-2H)
_crossTerm(H₁, H₂) = -_gammaFact((H₁+H₂)/2)/sqrt(_gammaFact(H₁)*_gammaFact(H₂))


"""
    covIMFBM(H::Function, D::Function, dt)

For given local Hurst exponent H and local diffusivity D functions
returns covariance function of the corresponding IMFBM.

Requirements: 0 < H(t) < 1;     0 < D(t).

This is covariance of the approximate IMFBM which converges to the true process as dt → 0.
"""
covIMFBM(H::Function, D::Function, dt::Real) = 
    function C(s,t)
        return _crossTerm(H(s),H(t))/2*sqrt(D(s)*D(t))*(abs(t-s+dt)^(H(s)+H(t))-2abs(t-s)^(H(s)+H(t))+abs(t-s-dt)^(H(s)+H(t)))
    end


"""
    msdIMFBM(H::Function, D::Function, dt)

For given local Hurst exponent H and local diffusivity D functions
returns mean square displacement (msd) function of the corresponding IMFBM.
Equal covIMFBM at arguments s=t, t=t.

Requirements: 0 < H(t) < 1;    0 < D(t).

This is msd of the approximate IMFBM which converges to the true process as dt → 0.
"""
msdIMFBM(H::Function, D::Function, dt::Real) = t -> covIMFBM(H,D,dt)(t,t)


"""
    covIMFBM(H₁, H₂, D₁, D₂, τ)

For given local Hurst exponents H and local diffusivity D which change from (D₁,H₁) to (D₂,H₂) at time τ
returns covariance function of the corresponding IMFBM.

Requirements: 0 < H₁, H₂ < 1;    0 < D₁, D₂;    0 < τ.
"""
function covIMFBM(H₁::Real, H₂::Real, D₁::Real, D₂::Real, τ::Real)
    c₁, c₂ = sqrt(D₁), sqrt(D₂)
    function C(s,t) 
        s > t && return C(t,s)
        t <= τ && return 1/2*( D₁*(s^(2H₁)+t^(2H₁)-(t-s)^(2H₁)) )
        s <= τ && return 1/2*( D₁*(s^(2H₁)+τ^(2H₁)-(τ-s)^(2H₁)) + 
            c₁*c₂*_crossTerm(H₁,H₂)*((τ-s)^(H₁+H₂)+t^(H₁+H₂)-τ^(H₁+H₂)-(t-s)^(H₁+H₂)) )
        return 1/2*( 2D₁*τ^(2H₁) + 
            c₁*c₂*_crossTerm(H₁,H₂)*(s^(H₁+H₂)+t^(H₁+H₂)-2τ^(H₁+H₂)-(t-τ)^(H₁+H₂)-(s-τ)^(H₁+H₂)) +
            D₂*((s-τ)^(2H₂)+(t-τ)^(2H₂)-(t-s)^(2H₂)) )
    end 
    return C
end

"""
    msdIMFBM(H₁, H₂, D₁, D₂, τ)

For given local Hurst exponents H and local diffusivity D which change from (D₁,H₁) to (D₂,H₂) at time τ
returns mean square displacement (msd) function of the corresponding IMFBM.
Equal covIMFBM at arguments s=t, t=t.

Requirements: 0 < H₁, H₂ < 1;    0 < D₁, D₂;    0 < τ.
"""
msdIMFBM(H₁::Real, H₂::Real, D₁::Real, D₂::Real, τ::Real) = t -> covIMFBM(H₁, H₂, D₁, D₂, τ)(t,t)
