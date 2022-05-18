"""
    ElectronBeam(ψ::Matrix{<:Complex}, λ::Real, x::Vector{<:Real}, y::Vector{<:Real})

Return the ElectronBeam type.

In the ElectronBeam type the wave function `ψ`, the de Broglie wavelength `λ` as well as
the coordinates `x` and `y` are saved.


# Example
```jldoctest
julia> ψ = ones(ComplexF64, 2, 2)
2×2 Matrix{ComplexF64}:
 1.0+0.0im  1.0+0.0im
 1.0+0.0im  1.0+0.0im

julia> ElectronBeam(ψ, 1e-12, [1, 2], [1, 2])
ElectronBeam(ComplexF64[1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im], 1.0e-12, [1, 2], [1, 2])

```
"""
mutable struct ElectronBeam <: Wave
    ψ::Matrix{<:Complex}
    λ::Real
    x::Vector{<:Real}
    y::Vector{<:Real}
end


"""
    deBroglieWavelength(U::Real)

Return the de Broglie wavelength of an electron.

This function calculates the relativistic de Broglie wavelength of
an electron using the acceleration Voltage `U` as an input.
The output unit will be in meters.

# Example
```jldoctest
julia> deBroglieWavelength(30e3)
6.979081574270336e-12
```

See also: [`ElectronBeam`](@ref)
"""
function deBroglieWavelength(U::Real)::Real
    return 2 * π * ħ / sqrt( 2 * U * q * m_e) / sqrt(1 + q * U / 2 / m_e / c^2)
end
