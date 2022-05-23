
# import statements
using FFTW

struct Free <: Component
    transferfunction::Matrix{<:Complex}
end

"""
    Free(wave::Wave, distance::Real)::PropTf

Return the Free type.

Create and return the PropTf struct that is needed to calculate the
wavefunction a given distance `distance` away. It also needs the
`wave`, which is an ElectronBeam type. Be aware that
for this function to calculate something meaningfull, critical sampling
must be mainitained!

# Example
```jldoctest
julia> eb = ElectronBeam(Array{Float64}(1:2), Array{Float64}(1:2), 1);

julia> Free(eb, 1)
Free(Complex{Float64}[1.0 - 1.926465401546721e-9im 1.0 - 1.926465401546721e-9im;
1.0 - 1.926465401546721e-9im 1.0 - 1.926465401546721e-9im])
```

See also: [`Aperture`](@ref), [`PhaseImprint`](@ref), [`Lens`](@ref), [`Edge`](@ref)
"""
function Free(wave::ElectronBeam, distance::Real)::Free

    # get the x axis out of the wave struct
    x = wave.x
    y = wave.y
    λ = wave.λ

    # get some information out of the transverse coordinates
    dx = abs(x[1] - x[2])
    dy = abs(y[1] - y[2])

    # define the frequency coordinates
    fx = Vector{Float64}(range(-1/(2*dx), 1/(2*dx), length=size(x, 1)))
    fy = Vector{Float64}(range(-1/(2*dy), 1/(2*dy), length=size(y, 1)))

    # fill the transferfunction arrays with the proper values
    transferfunction = @. exp(-1im * π * λ * distance * (fx^2 + fy'^2))

    # shift the transferfunction and return the Free object
    return Free(fftshift(transferfunction))
end


"""
    calculate!(wave::Wavem free::Free)

Calculate the Free.

This function calculates the light field in a given direction. After the
calculation the changed LightField object is returned.
"""
function calculate!(wave::ElectronBeam, free::Free)

    # Fouriertransform the input wave
    Ψ = fft(fftshift(wave.ψ))

    # apply the transferfunction
    Ψ .*= free.transferfunction

    # calculate the inverse fouriertransform
    wave.ψ = ifftshift(ifft(Ψ))
end


