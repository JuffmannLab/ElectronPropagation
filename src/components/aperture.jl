

"""
    Aperture(d::Real)::Aperture

Return the Aperture type.

This struct simulates an aperture, where `d` is the
diameter of the aperture.

# Example
```jldoctest
julia> Aperture(1)
Aperture(1)
```

See also: [`Free`](@ref), [`PhaseImprint`](@ref), [`Lens`](@ref)
"""
struct Aperture <: Component
    d::Real
end

"""
    Calculate the aperture.

    Apply the aperture to a given wavefunction.
"""
function calculate!(wave::ElectronBeam, aperture::Aperture)
    for i = 1:size(wave.x, 1), j = 1:size(wave.y, 1)
        if wave.x[i]^2+wave.y[j]^2 <= (aperture.d/2)^2
        else
            wave.ψ[i, j] = 0
        end
    end
end
