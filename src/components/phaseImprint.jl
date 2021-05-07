
using FFTW

"""
    Struct for the phase imprint.

    This struct is the basis vor the Phase imprint that can be done.
"""
struct PhaseImprint <: Component
    lb::LaserBeam     # the laser beam type
end


"""
    Imprint the phase.

    This function imprints a phase on the electron wave.
    This function is only working if the coordinates of the wavefunction
    are the same in x and y direction! This function only works if the
    x and y components of the laser beam are the same as the x and y coordinates
    of the wave object.
"""
function calculate!(wave::Wave, imprint::PhaseImprint)
    # wave      ...   the electron beam
    # imprint   ...   the imprint object

    @info "Calculate the phase imprint..."

    # extract the laserbeam type out of the imprint type
    lb = imprint.lb
    x = lb.x
    y = lb.y
    z = lb.z

    # define the dx dy dz values
    dx = abs(x[2]-x[1])
    dy = abs(y[2]-y[1])
    dz = abs(z[2]-z[1])

    # define the 3D envelope function
    envelope = Array{complex(eltype(lb.I))}(undef, size(x, 1), size(y, 1), size(z, 1))

    # create the 3D Intensity function
    intensity = similar(envelope)

    # fill those two functions
    for k = 1:size(z, 1)

        # the value has to be calculated only once
        val = exp(-z[k]^2 / (c*lb.Δt)^2)

        for j = 1:size(y, 1)
            for i = 1:size(x, 1)
                envelope[i, j, k] = val
                intensity[i, j, k] = lb.I[i, j]
            end
        end
    end

    # calculate the normalization constant
    I0 = lb.E * c / (sum(intensity .* envelope) * dx * dy * dz)

    # calculate the convolution of the intensity and the envelope function
    fft!(envelope, 3)
    fft!(intensity, 3)
    intensity .*= envelope
    ifft!(intensity, 3)
    intensity .*= dz

    # calculate constants
    α = q^2 * I0 / (2 * m_e * ε_0 * c^2 * lb.ω^2 * ħ)

    # apply the phase to the electron beam
    # first get out the size of the intenisty array
    n = size(intensity, 1)

    # calculate the position of the left upper array in
    # the electron wave
    posx = round(Int, (size(wave.ψ, 1)-n)/2)
    posy = round(Int, (size(wave.ψ, 1)-n)/2)


    # apply the phase to the electron beam on the
    # calculated positions
    wave.ψ[posx+1:posx+n, posy+1:posy+n] .*=
           exp.(1im * α * real(intensity[:, :, end]))
end
