

"""
    Lense(f::Real)::Lense

Return a Lense type.

Here a lense type is returned. The `f` value is the focal length of the lense.

# Example
```jldoctest
julia> Lense(1)
Lense(1)
```

See also: [`PropTf`](@ref), [`PhaseImprint`](@ref), [`Aperture`](@ref),
[`PropDirect`](@ref), [`Edge`](@ref)
"""
struct Lense <: Component
    f::Real   # focal length of the lense
end

"""
    calculate!(wave::Wave, lense::Lense)

Calculate the lense.

This function simulates the effect of a perfect lense on the wavefunction.
"""
function calculate!(wave::Wave, lense::Lense)
    # wave    ...   the wave type
    # lense   ...   the lense type

    @info "Calculate the lense..."

    # define the wavevectoro
    k = 2 * π / wave.λ

    # add the lense phase onto the wavefunction
    for i = 1:size(wave.ψ, 1)
        for j = 1:size(wave.ψ, 2)
            wave.ψ[i, j] *= exp(-1im * k / 2 / lense.f * (wave.x[i]^2 + wave.y[j]^2))
        end
    end
end
