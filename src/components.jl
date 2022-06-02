

# create the abstract type of the component
abstract type Component end

# include all the components
include("./components/aperture.jl")
include("./components/freePropagation.jl")
include("./components/lens.jl")
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
Setup(comps::Component...) = Setup(collect(comps))

"""
    propagation!(wave::ElectronBeam, setup::Setup)

Propagate the setup.

This function propagates the given wave `wave` through the given setup `setup`.
The wave will be altered according to the propagation.

See also: [`Setup`](@ref)
"""
function propagation!(wave::ElectronBeam, setup::Setup)
    for comp in setup.setup
        calculate!(wave, comp)
    end
end
