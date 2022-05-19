
using FFTW

"""
    PhaseImprint(lb::LaserBeam)::PhaseImprint

Return the PhaseImprint type.

The PhaseImprint type is generated and returned. The `lb`
value has to be a LaserBeam type, and represents the
laser beam that will make the phase imprint.

# Example
```jldoctest
julia> lb = LaserBeam(Array(1:2), Array(1:2), Array(1:2), 1035e-9, 15e-6, 40e-6, 280e-15);

julia> PhaseImprint(lb)
PhaseImprint(LaserBeam([0 0; 0 0], [1, 2], [1, 2], [1, 2], 1.8199532051293268e15, 4.0e-5, 2.8e-13))
```

See also: [`Free`](@ref), [`Aperture`](@ref), [`Lense`](@ref), [`ElectronBeam`](@ref)
"""
struct PhaseImprint <: Component
    lb::LaserBeam
end


"""
Imprint the phase.

This function imprints a phase on the electron wave.
This function is only working if the coordinates of the wavefunction
are the same in x and y direction! This function only works if the
x and y components of the laser beam are the same as the x and y coordinates
of the wave object.
"""
function calculate!(wave::ElectronBeam, imprint::PhaseImprint)
    Δx = abs(imprint.lb.x[1] - imprint.lb.x[2])
    Δy = abs(imprint.lb.y[1] - imprint.lb.y[2])

    β = wave.v / c
    γ = 1 /sqrt(1-β^2)

    Ee = γ * m_e * c^2
    α = q^2 / (ħ * c) / (4 * π * ε_0)

    norm = sum(imprint.lb.I) * Δx * Δy

    prefactor = - α / (2*π * (1 + β)) * imprint.lb.E / Ee * imprint.lb.λ^2 / norm

    for j in eachindex(wave.y), i in eachindex(wave.x)
        wave.ψ[i, j] *= exp(1im * prefactor * interpolation(imprint.lb.I, imprint.lb.x, imprint.lb.y, wave.x[i], wave.y[j]))
    end
end


"""
    interpolation(A::Matrix{<:Real}, coords_x::Vector{<:Real},
                   coords_y::Vector{<:Real}, x::Real, y::Real)::Real

Return the bilinear interpolation.
Calculate and return the bilinear interpolation on the Matrix `A`
with the coordinates `coords_x` and `coords_y` at the point with
the position `x`, `y`.
"""
function interpolation(A::Matrix{<:Real}, coords_x::Vector{<:Real},
                       coords_y::Vector{<:Real}, x::Real, y::Real)::Real

    # if the value is not in the area of coords return 0
    if x < coords_x[1] || x > coords_x[end]
        return zero(eltype(A))
    elseif y < coords_y[1] || y > coords_y[end]
        return zero(eltype(A))
    end

    # calculate the Δx and Δy values
    Δx = abs(coords_x[1]-coords_x[2])
    Δy = abs(coords_y[1]-coords_y[2])

    # calculate the bins in which the x and y values are
    n = round(Int, (x-coords_x[1]) / Δx, RoundDown) + 1
    m = round(Int, (y-coords_y[1]) / Δy, RoundDown) + 1

    # interpolation in x direction on both y-points
    k1 = (A[n, m] - A[n+1, m]) / (coords_x[n] - coords_x[n+1])
    d1 = A[n, m] - k1 * coords_x[n]

    k2 = (A[n, m+1] - A[n+1, m+1]) / (coords_x[n] - coords_x[n+1])
    d2 = A[n, m+1] - k2 * coords_x[n]

    A1 = k1 * x + d1
    A2 = k2 * x + d2

    # use the interpolated values on the x-axis for the interpolation in the y axis
    k = (A1 - A2) / (coords_y[m] - coords_y[m+1])
    d = A1 - k * coords_y[m]

    # return the interpolated value
    return k * y + d
end
