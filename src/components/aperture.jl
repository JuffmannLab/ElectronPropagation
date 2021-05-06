

"""
    The aperture struct.
"""
struct Aperture <: Component
    d::Real
end

"""
    Calculate the aperture.

    Apply the aperture to a given wavefunction.
"""
function calculate!(wave::Wave, aperture::Aperture)
    # wave       ...   the wave type
    # aperture   ...   the aperture type

    @info "Calculate the aperture..."

    # apply the aperute to the wavefunction
    for i = 1:size(wave.x, 1)
        for j = 1:size(wave.y, 1)
            if wave.x[i]^2+wave.y[j]^2 <= aperture.d^2
            else
                wave.ψ[i, j] = 0
            end
        end
    end
end
