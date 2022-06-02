"""
    ElectronBeam(ψ::Matrix{<:Complex}, λ::Real, v::Real, x::Vector{<:Real}, y::Vector{<:Real})

Return the ElectronBeam type.

In the ElectronBeam type the wave function `ψ`, the de Broglie wavelength `λ`,
the electron velocity `v` as well as the coordinates `x` and `y` are saved.

# Example
```jldoctest
julia> ψ = ones(ComplexF64, 2, 2)
2×2 Matrix{ComplexF64}:
 1.0+0.0im  1.0+0.0im
 1.0+0.0im  1.0+0.0im

julia> ElectronBeam(ψ, 1e-12, 1e8, [1, 2], [1, 2])
ElectronBeam(ComplexF64[1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im], 1.0e-12, [1, 2], [1, 2])

```

See also: [`deBroglieWavelength`] [`LaserBeam`](@ref)
"""
mutable struct ElectronBeam
    ψ::Matrix{<:Complex}
    λ::Real
    v::Real
    x::Vector{<:Real}
    y::Vector{<:Real}
end


"""
    ElectronBeam(ψ::Matrix{<:Complex}, U::Real, x::Vector{<:Real}, y::Vector{<:Real})

Return the ElectronBeam type.

Calculate the de Broglie wavelength and the velocity, and returning the corresponding
ElectronBeam type. The wavefunction is `ψ`, the coordinates are `x` and `y`.
`U` is the electron acceleration voltage.

# Example
```jldoctest
julia> ψ = ones(ComplexF64, 2, 2)
2×2 Matrix{ComplexF64}:
 1.0+0.0im  1.0+0.0im
 1.0+0.0im  1.0+0.0im

julia> ElectronBeam(ψ, 30e3, [1, 2], [1, 2])
ElectronBeam(ComplexF64[1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im], 6.979081574270336e-12, 9.844470106002943e7, [1, 2], [1, 2])

```

See also: [`LaserBeam`](@ref)
"""
function ElectronBeam(ψ::Matrix{<:Complex}, U::Real, x::Vector{<:Real}, y::Vector{<:Real})
    λ = deBroglieWavelength(U)
    v = c * sqrt(1 - 1 / (1 + q*U/m_e/c^2)^2)
    return ElectronBeam(ψ, λ, v, x, y)
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
