
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

    # define the length of the x axis, and the dy slice
    l = abs(wave.x[1] - wave.x[end])
    dy = abs(wave.y[1] - wave.y[2])

    # iterate over the columns with no transmission
    for i = 1:round(Int, size(wave.ψ, 2)/2+1)+edge.offset
        
        # subtract the electron probability that gets lost in this column
        wave.norm -= sqrt(sum(abs2.(wave.ψ[:, i])) * l * dy)

        # apply the edge
        wave.ψ[:, i] .*= 0
    end
end
