

# create the abstract type of the component
abstract type Component end

# include all the components
include("./components/aperture.jl")
include("./components/freePropagation.jl")
include("./components/knifeEdge.jl")
include("./components/lense.jl")
include("./components/mcp.jl")
include("./components/phaseImprint.jl")

# The setup struct
struct Setup
    setup::Array{<:Component}   # an array of components
end

"""
    Create the setup.

    In this struct the whole setup is saved as an array of components.
    The setup will be calculated in the order of appearence in the function
    definition.
"""
function Setup(comps::Component...)
    # comps...  varargs for the input of the components
    return Setup(collect(comps))
end

"""
    Propagate the setup.

    This function propagates the given wave through the given setup.
    The wave will be altered according to the propagation.
"""
function propagation!(wave::Wave, setup::Setup)
    # wave    ...   The wave that should be propagated
    # setup   ...   The setup through which the propagation should take place

    # iterate through the setup
    for comp in setup.setup
        calculate!(wave, comp)
    end

    # renormalize the wave function
    wave.ψ ./= sqrt(sum(abs2.(wave.ψ)) * abs(wave.x[1]-wave.x[2])
                    * abs(wave.y[1]-wave.y[2]))
end

"""

"""
function _normalization(wave)
end
