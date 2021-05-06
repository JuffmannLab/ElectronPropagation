
using ZernikePolynomials

# The Laser beam struct, it saves important information about the laser beam
mutable struct LaserBeam
    I::Matrix{<:Real}   # Laser intensity distribution (assume uniform in z)
    x::Vector{<:Real}   # the x coordinates
    y::Vector{<:Real}   # the y coordinates
    z::Vector{<:Real}   # the z coordinates
    ω::Real             # the laser frequency
    E::Real             # energy per pulse
    Δt::Real            # pulse length
end

"""
    LaserBeam(x::Vector{<:Real}, y::Vector{<:Real}, z::Vector{<:Real}, λ::Real
              d::Real, E::Real, Δt::Real)::LaserBeam

Construct and return the LaserBeam type.

This function creates a LaserBeam type and creates the Zernike Intensity pattern
that is used to calculate the phase imprint onto the electron beam.
"""
function LaserBeam(x::Vector{<:Real}, y::Vector{<:Real}, z::Vector{<:Real},
                   λ::Real, d::Real, E::Real, Δt::Real)::LaserBeam
    # x    ...   x coordinates
    # y    ...   y coordinates, need to be the same as x!!!
    # z    ...   z coordinates
    # λ    ...   wavelength
    # d    ...   diameter of the beam
    # E    ...   energy that is in the light pulse
    # Δt   ...   pulse length of the laser

    # calculate the angular frequency
    ω = 2 * π * c / λ

    # calculate the zernike intensity pattern
    I = _zernikeAmplitude(x, d)

    # return the LaserBeam struct
    return LaserBeam(I, x, y, z, ω, E, Δt)
end


"""
    _zernikeAmplitude(x::Vector{<:Real}, d::Real)::Matrix{<:Real}

Return a Zernike Beam.

This function takes the coordinates `x`, and the diameter `d` to construct
a Amplitude of the spherical aberration zernike polynomial, with the given
diameter. The polynimial will be calculated in 2D with the sides being the `x`
coordinate. Due to the nature of the ZernikePolynomials package, only quadratic
shapes can be returned.
"""
function _zernikeAmplitude(x::Vector{<:Real}, d::Real)::Matrix{<:Real}
    # x   ...   coordinates for the amplitude
    # d   ...   the diameter of the beam

    # Create the empty laser intensity matrix
    I = Matrix{eltype(x)}(undef, size(x, 1), size(x, 1))

    # change the diameter to the radius that is needed four our calculation
    r = d / 2

    # calculate the dx value (dy MUST be the same... because the zernike
    # package only supplies with symmetric outputs, and i am lazy)
    dx = abs(x[1] - x[2])

    # calculate the amount of pixels that are the zernike phase diameter
    # +2 because the evaluateZernike function leaves some room that is not used
    m = round(Int, d / dx) + 3

    # shift the matrix elements to be positive
    ϕ = evaluateZernike(m, 4, 1., index=:Noll)
    ϕ .+= abs(minimum(ϕ))

    # calculate the relative shift to apply the zernike polynomials to the
    # wavefunction
    ishift = round(Int, (size(I, 1)-size(ϕ, 1))/2)
    jshift = round(Int, (size(I, 2)-size(ϕ, 2))/2)

    # fill the intensity array with either 0 or the zernike amplitude
    for j = 1:size(I, 2)
        for i = 1:size(I, 1)
            if r^2 >= x[i]^2+y[j]^2
                I[i, j] = ϕ[i-ishift, j-jshift]
            else
                I[i, j] = 0
            end
        end
    end

    # return the intesity matrix
    return I
end
