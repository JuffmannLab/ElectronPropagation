

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
end

"""
    Instantiate the electron beam type.

    The electron beam type will be created and returned with this function.
"""
function ElectronBeam(x, y, U)
    # x   ...   the x axis of the geometry
    # y   ...   the y axis of the geometry
    # U   ...   the acceleration voltage of the electron beam

    # define the de'Broglie wavelength of the electrons
    λ = 2 * π * ħ / sqrt( 2 * U * q / m_e ) / m_e

    # create a plane electron wave
    ψ = ones(complex(eltype(x)), size(x, 1), size(y, 1))

    # normalze the input electron beam
    ψ ./= sqrt(sum(abs2.(ψ)) * abs(x[1]-x[2]) * abs(y[1]-y[2]))

    # return the electron beam struct
    return ElectronBeam(ψ, λ, x, y, size(x, 1), size(y, 1))
end
