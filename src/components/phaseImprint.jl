
using FFTW

"""
    PhaseImprint(lb::LaserBeam)::PhaseImprint

Return the PhaseImprint type.

The PhaseImprint type is generated and returned. The `lb`
value has to be a LaserBeam type, and represents the
laser beam that will make the phase imprint.

# Example
```jldoctest
julia> lb = LaserBeam(Array(1:2), Array(1:2), Array(1:2), 1035e-9, 15e-6, 40e-6, 280e-15);

julia> PhaseImprint(lb)
PhaseImprint(LaserBeam([0 0; 0 0], [1, 2], [1, 2], [1, 2], 1.8199532051293268e15, 4.0e-5, 2.8e-13))
```

See also: [`PropTf`](@ref), [`Aperture`](@ref), [`Lense`](@ref),
[`PropDirect`](@ref), [`Edge`](@ref)
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

    # extract the laserbeam type out of the imprint type
    lb = imprint.lb
    x = lb.x
    y = lb.y
    t = lb.t

    # define the dx dy dz values
    dx = abs(x[2]-x[1])
    dy = abs(y[2]-y[1])
    dt = abs(t[1]-t[2])

    # create the phase function
    phase = similar(lb.I)

    # calculate the normalization constant
    env_int = sum(lb.env) * dt
    I0 = lb.E / (sum(lb.I .* env_int) * dx * dy)

    # calculate constants
    α = q^2 * I0 / (2 * m_e * ε_0 * c * lb.ω^2 * ħ)

    # calculate the phase imprint.
    for j = 1:size(lb.I, 2)
        for i = 1:size(lb.I, 1)
            phase[i, j] = I0 * α * lb.I[i, j] * env_int
        end
    end

    # apply the phase to the electron beam
    @. wave.ψ *= exp(1im * phase)
end
