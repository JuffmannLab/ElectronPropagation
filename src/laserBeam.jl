
using Images

"""
    LaserBeam(I::Matrix{<Real}, env::Vector{<:Real}, x::Vector{<:Real},
              y::Vector{<:Real}, t::Vector{<:Real}, λ::Real, E::Real)::LaserBeam

Return the LaserBeam type.

Take the 2D intensity pattern `I`, the temporal envelope function `env`, the coordinate system
`x`, `y` and `t`, the wavelength `λ` and the energy per laser pulse `E`. With these informations
a LaserBeam object is created and returned.

# Example
```jldoctest
julia> LaserBeam(ones(2, 2), Array(1:2), Array(1:2), Array(1:2), Array(1:2), 1035e-9, 13e-6)
LaserBeam([1.0 1.0; 1.0 1.0], [1, 2], [1, 2], [1, 2], [1, 2], 1.8199532051293268e15, 1.3e-5, 1.0)
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


"""
    loadintensity(s::String)::Matrix{<:Real}

Return the intensity matrix.

Load the intensity image from the file at `s`. The image will then
be normalized, such that the values range from 0 to 1, and then returned as
a Matrix with the datatype being Float64.

```jldoctest
julia> loadintensity("test.png")
2×2 Matrix{Float64}:
1.0  1.0
1.0  1.0
```
"""
function loadintensity(s::String)::Matrix{<:Real}
    # load the image, save it in a matrix
    input = Float64.(Gray.(load(s)))

    # normalize the image between 0 and 1
    input .-= minimum(input)
    input ./= maximum(input)

    # return the wanted image
    return input
end
