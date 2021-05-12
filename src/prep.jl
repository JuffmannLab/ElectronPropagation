"""
    zeropadding!(eb::ElectronBeam, d::Real)    

Apply zeropaddinig to `eb`.

This function adds the zeropadding that is needed for fourier transform
based propagation approaches. It calculates the critical sampling for
the zeropadding from the information that is stored in the ElectronBeam
type `eb` and the propagation distance `d`. The new, zeropadded wavefunction
is saved in the ElectronBeam struct.

See also: [`removezeropadding!`](@ref)
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
    eb.x = Array{eltype(eb.x)}(range(x_min, x_max, step=dx))

    # create the new coordinate system in y direction
    y_min = eb.y[1] - dy * shift_y
    y_max = y_min + dy * (m_y - 1)
    eb.y = Array{eltype(eb.x)}(range(y_min, y_max, step=dy))

    # create the new zeropadded array
    ψ = zeros(eltype(eb.ψ), m_x, m_y)

    # fill the new array
    ψ[shift_x:shift_x+eb.n_x-1, shift_y:shift_y+eb.n_y-1] = eb.ψ

    # save the zeropadded array
    eb.ψ = ψ
end


"""
    removezeropadding!(eb::ElectronBeam)

Remove the zeropadding from `eb`.

Removes the zeropadding from the ElectronBeam struct that is
given with it. It is compatible with the zeropadding! function
and does the oposite of it.

See also [`zeropadding!`](@ref)
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
