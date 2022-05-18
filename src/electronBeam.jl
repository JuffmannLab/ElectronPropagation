

"""
    Create the electron beam type.
"""
mutable struct ElectronBeam <: Wave
    ψ::Array{<:Complex}     # The wave function
    λ::Real                 # de'Broglie wavelength
    x::Array{<:Real}        # x-axis
    y::Array{<:Real}        # y-axis
    norm::Real              # normalization parameter
end


"""
ElectronBeam(ψ::Matrix{<:Complex}, x::Vector{<:AbstractFloat},
             y::Vector{<:AbstractFloat}, U::Real)::ElectronBeam

Return the ElectronBeam type.

The electron beam type will be created and returned with this function.
The `x` value is the x-axis of the geometry, the `y` value is the y-axis
of the geometry, and `U` is the acceleration voltage of the Electrons.

The input wave `ψ` will be normalized and saved.

# Example
```jldoctest
julia> x = Array{Float64}(1:2);

julia> y = Array{Float64}(1:2);

julia> ψ = ones(ComplexF64, 2 , 2);

julia> U = 1;

julia> ElectronBeam(ψ, x, y, U)
ElectronBeam(ComplexF64[0.25 + 0.0im 0.25 + 0.0im; 0.25 + 0.0im 0.25 + 0.0im], 1.2264259654066947e-9, [1.0, 2.0], [1.0, 2.0], 2, 2, 1.0)


```

See also: [`LaserBeam`](@ref)
"""
function ElectronBeam(ψ::Matrix{<:Complex}, x::Vector{<:AbstractFloat},
                      y::Vector{<:AbstractFloat}, U::Real)::ElectronBeam
    # x   ...   the x axis of the geometry
    # y   ...   the y axis of the geometry
    # U   ...   the acceleration voltage of the electron beam

    # define the de'Broglie wavelength of the electrons
    λ = 2 * π * ħ / sqrt( 2 * U * q * m_e) / sqrt(1 + q * U / 2 / m_e / c^2)

    # normalize the input beam
    ψ ./= sqrt(sum(abs2.(ψ))*abs(x[1]-x[2])*abs(y[1]-y[2]))

    # set the normalization parameter to 1 (100% of the electrons are present)
    norm = 1.

    # return the electron beam struct
    return ElectronBeam(ψ, λ, x, y, norm)
end
