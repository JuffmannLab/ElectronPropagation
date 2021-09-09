

# create the abstract type of the component
abstract type Component end

# include all the components
include("./components/aperture.jl")
include("./components/freePropagation.jl")
include("./components/knifeEdge.jl")
include("./components/lens.jl")
include("./components/mcp.jl")
include("./components/phaseImprint.jl")

# The setup struct
struct Setup
    setup::Array{<:Component}   # an array of components
end

"""
    Setup(comps::Component...)::Setup

Return the setup.

Create and return the setup type. This type will save the current
components in the order that they should be calculated. The components
are in the varargs `comps`.

See also: [`propagation!`](@ref)
"""
function Setup(comps::Component...)::Setup
    # comps...  varargs for the input of the components
    return Setup(collect(comps))
end

"""
    propagation!(wave::Wave, setup::Setup)

Propagate the setup.

This function propagates the given wave `wave` through the given setup `setup`.
The wave will be altered according to the propagation.

See also: [`Setup`](@ref)
"""
function propagation!(wave::Wave, setup::Setup)
    # wave    ...   The wave that should be propagated
    # setup   ...   The setup through which the propagation should take place

    # iterate through the setup
    for comp in setup.setup
        calculate!(wave, comp)
    end

    # normalization
    _normalization!(wave)
end

"""
    _normalization!(wave::Wave)

Normalize the wave.

This function will normalize the wave that is given in `wave`.
For that it calculates the absoltue square and matches it to the normalizaion
paramter `wave.norm` that is given in the wave type.
"""
function _normalization!(wave::Wave)
    # wave   ...   either an electron or light wave
 
    # Normalize the wave function
    wave.ψ .*= sqrt(wave.norm / sum(abs2.(wave.ψ)) /
                    abs(wave.x[1] - wave.x[2]) /
                    abs(wave.y[1] - wave.y[2]))
end
