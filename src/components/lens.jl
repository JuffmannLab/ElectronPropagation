

"""
    Lens(f::Real)::Lens

Return a Lens type.

Here a lens type is returned. The `f` value is the focal length of the lens.

# Example
```jldoctest
julia> Lens(1)
Lens(1)
```

See also: [`PropTf`](@ref), [`PhaseImprint`](@ref), [`Aperture`](@ref),
[`PropDirect`](@ref), [`Edge`](@ref)
"""
struct Lens <: Component
    f::Real   # focal length of the lense
end

"""
    calculate!(wave::Wave, lens::Lens)

Calculate the lense.

This function simulates the effect of a perfect lens on the wavefunction.
"""
function calculate!(wave::ElectronBeam, lens::Lens)
    if lens.f == 0
        return
    else
        wave.ψ .*= @. exp(-1im * π / wave.λ / lens.f * (wave.x^2 + wave.y'^2))
    end
end
