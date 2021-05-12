

"""
    Create the electron beam type.
"""
mutable struct ElectronBeam <: Wave
    ψ::Array{<:Complex}     # The wave function
    λ::Real                 # de'Broglie wavelength
    x::Array{<:Real}        # x-axis
    y::Array{<:Real}        # y-axis
    n_x::Integer            # Array size before zeropadding x
    n_y::Integer            # Array size before zeropadding y
    norm::Real              # normalization parameter
end


"""
    ElectronBeam(x::Vector{<:AbstractFloat}, y::Vector{<:AbstractFloat}, U::Real)::ElectronBeam

Return the ElectronBeam type.

The electron beam type will be created and returned with this function.
The `x` value is the x-axis of the geometry, the `y` value is the y-axis
of the geometry, and `U` is the acceleration voltage of the Electrons.

The ElectronBeam type will consist of a plane wave with a wavelength
that corresponds to the de'Broglie wavelength of the electrons at the
given acceleration Voltage.

# Example
```jldoctest
julia> x = Array{Float64}(1:2);

julia> y = Array{Float64}(1:2);

julia> U = 1;

julia> ElectronBeam(x, y, U)
ElectronBeam(Complex{Float64}[0.5 + 0.0im 0.5 + 0.0im; 0.5 + 0.0im 0.5 + 0.0im], 1.2264259654066947e-9, [1.0, 2.0], [1.0, 2.0], 2, 2)
```

See also: [`LaserBeam`](@ref)
"""
function ElectronBeam(x::Vector{<:AbstractFloat}, y::Vector{<:AbstractFloat}, U::Real)::ElectronBeam
    # x   ...   the x axis of the geometry
    # y   ...   the y axis of the geometry
    # U   ...   the acceleration voltage of the electron beam

    # define the de'Broglie wavelength of the electrons
    λ = 2 * π * ħ / sqrt( 2 * U * q / m_e ) / m_e

    # create a plane electron wave
    ψ = ones(complex(eltype(x)), size(x, 1), size(y, 1))

    # normalze the input electron beam
    ψ ./= sqrt(sum(abs2.(ψ)) * abs(x[1]-x[2]) * abs(y[1]-y[2]))

    # set the normalization parameter to 1 (100% of the electrons are present)
    norm = 1

    # return the electron beam struct
    return ElectronBeam(ψ, λ, x, y, size(x, 1), size(y, 1), norm)
end
