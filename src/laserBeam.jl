
using Images

# The Laser beam struct, it saves important information about the laser beam
mutable struct LaserBeam <: Wave
    I::Matrix{<:Real}   # Laser intensity distribution (assume uniform in z)
    env::Vector{<:Real} # temporal envelope
    x::Vector{<:Real}   # the x coordinates
    y::Vector{<:Real}   # the y coordinates
    t::Vector{<:Real}   # the t coordinates
    ω::Real             # the laser frequency
    E::Real             # energy per pulse
    norm::Real          # the normalization parameter
end

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
function LaserBeam(I::Matrix{<:Real}, env::Vector{<:Real}, x::Vector{<:Real},
                   y::Vector{<:Real}, t::Vector{<:Real}, λ::Real, E::Real)::LaserBeam
    # I    ...   laser beam intensity distribution
    # x    ...   x coordinates
    # y    ...   y coordinates, need to be the same as x!!!
    # λ    ...   wavelength
    # E    ...   energy that is in the light pulse

    # calculate the angular frequency
    ω = 2 * π * c / λ

    # set the normalization parameter to 1 (100% of the photons are still present)
    norm = 1.

    # return the LaserBeam struct
    return LaserBeam(I, env, x, y, t, ω, E, norm)
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
