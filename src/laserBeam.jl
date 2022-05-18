
"""
    LaserBeam(I::Matrix{<Real}, λ::Real, E::Real, x::Vector{<:Real}, y::Vector{<:Real})

Return the LaserBeam type.

Save the intensity distribution `I`, the wavelength `λ`, the pulse energy `E`
and the coordinates `x` and `y` in the LaserBeam struct that is returned.

# Example
```jldoctest
julia> I = rand(Float64, 2, 2)
2×2 Matrix{Float64}:
 0.956952  0.353413
 0.52275   0.874929

julia> LaserBeam(I, 1e-6, 1e-6, [1, 2], [1, 2])
LaserBeam([0.9569521001448673 0.35341262246269656; 0.5227499716351234 0.8749291118118734], 1.0e-6, 1.0e-6, [1, 2], [1, 2])

```

See also: [`ElectronBeam`](@ref)
"""
mutable struct LaserBeam <: Wave
    I::Matrix{<:Real}
    λ::Real
    E::Real
    x::Vector{<:Real}
    y::Vector{<:Real}
end
