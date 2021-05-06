
struct Edge <: Component
    offset::Int    # number of pixel offset
end

"""
    calculate!(wave.:Wave, edge::Edge)

Calculate the effect of the knife edge.

This function calculates the impact on the wavefunction of `wave`
from a knife edge that is define by the `edge` object. In the edge
object there is an offset saved, that can be added, and that will
offset the knife edge from the middle of the coordinate system with
`edge.offset` number of pixels.
"""
function calculate!(wave::Wave, edge::Edge)
    # wave     ...   the wave type
    # edge     ...   the edge type

    @info "Calculate the knife edge..."

    for i = 1:round(Int, size(wave.ψ, 2)/2+1)+edge.offset
        wave.ψ[:, i] .*= 0
    end
end
