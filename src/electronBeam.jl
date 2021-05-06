

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

"""
    Add zeropadding.

    This function adds the zeropadding that is needed for fourier transform
    based propagation approaches.
"""
function zeropadding!(eb::ElectronBeam, d::Real)
    # eb   ...   the electron beam object
    # d    ...   the distance propagated

    # define the steopsize
    dx = abs(eb.x[1] - eb.x[2])
    dy = abs(eb.y[1] - eb.y[2])

    # calculate the critical sampling condition
    m_x = round(Int32, d * eb.λ / dx^2)
    m_y = round(Int32, d * eb.λ / dy^2)

    # error handling, if the wavefunction is oversampled
    if m_x < eb.n_x || m_y < eb.n_y
        throw(DomainError("The wavefunction is already oversampled"))
    end

    # calculate the shift of the original data in the new array
    shift_x = round(Int, (m_x-eb.n_x)/2)
    shift_y = round(Int, (m_y-eb.n_y)/2)

    # create the new coordinate system in x direction
    x_min = eb.x[1] - dx * shift_x
    x_max = x_min + dx * (m_x - 1)
    eb.x = Array{eltype(x)}(range(x_min, x_max, step=dx))

    # create the new coordinate system in y direction
    y_min = eb.y[1] - dy * shift_y
    y_max = y_min + dy * (m_y - 1)
    eb.y = Array{eltype(x)}(range(y_min, y_max, step=dy))

    # create the new zeropadded array
    ψ = zeros(eltype(eb.ψ), m_x, m_y)

    # fill the new array
    ψ[shift_x:shift_x+eb.n_x-1, shift_y:shift_y+eb.n_y-1] = eb.ψ

    # save the zeropadded array
    eb.ψ = ψ
end


"""
    Remove zeropadding.

    This function removes the previously added zeropadding again, such that
    only the incident fov is shown.
"""
function removezeropadding!(eb::ElectronBeam)
    # eb   ...   the electron beam that needs the zeropadding removed

    # check if removing the zeropadding is needed
    if size(eb.x) == eb.n_x && size(eb.y) == eb.n_y
        return
    end

    # calculate the shift, such that the correct data is retrieved
    shift_x = round(Int, (size(eb.x, 1)-eb.n_x)/2)
    shift_y = round(Int, (size(eb.y, 1)-eb.n_y)/2)

    # remove the zeropadding from the coordinate systems
    eb.x = eb.x[shift_x+1:shift_x+eb.n_x]
    eb.y = eb.y[shift_y+1:shift_y+eb.n_y]

    # remove the zeropadding from the wavefunction
    eb.ψ = eb.ψ[shift_x:shift_x+eb.n_x-1, shift_y:shift_y+eb.n_y-1]
end
