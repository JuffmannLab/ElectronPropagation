

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
function calculate!(wave::Wave, lens::Lens)
    # wave    ...   the wave type
    # lense   ...   the lense type

    # define the wavevector
    k = 2 * π / wave.λ

    # add the lense phase onto the wavefunction
    for j = 1:size(wave.ψ, 2)
        for i = 1:size(wave.ψ, 1)
            wave.ψ[i, j] *= exp(-1im * k / 2 / lens.f * (wave.x[i]^2 + wave.y[j]^2))
        end
    end
end
