
# import statements
using FFTW
using ProgressMeter


struct PropTf <: Component
    transferfunction::Array{<:Complex}
end

"""
    PropTf(wave::Wave, distance::Real)::PropTf

Return the PropTf type.

Create and return the PropTf struct that is needed to calculate the
wavefunction a given distance `distance` away. It also needs the
`wave` type, which is either a LightBeam or ElectronBeam. Be aware that
for this function to calculate something meaningfull, critical sampling
must be mainitained!

# Example
```jldoctest
julia> eb = ElectronBeam(Array{Float64}(1:2), Array{Float64}(1:2), 1);

julia> PropTf(eb, 1)
PropTf(Complex{Float64}[1.0 - 1.926465401546721e-9im 1.0 - 1.926465401546721e-9im;
1.0 - 1.926465401546721e-9im 1.0 - 1.926465401546721e-9im])
```

See also: [`Aperture`](@ref), [`PhaseImprint`](@ref), [`Lense`](@ref),
[`PropDirect`](@ref), [`Edge`](@ref)
"""
function PropTf(wave::Wave, distance::Real)::PropTf
    # wave       ...   some kind of wave struct
    # distance   ...   propagation distance

    # get the x axis out of the wave struct
    x = wave.x
    y = wave.y
    λ = wave.λ

    # get some information out of the transverse coordinates
    dx = abs(x[1] - x[2])
    dy = abs(y[1] - y[2])
    n = size(x, 1)
    m = size(y, 1)

    # define the frequency coordinates
    fx = Array(range(-1/(2*dx), 1/(2*dx), length=n))
    fy = Array(range(-1/(2*dy), 1/(2*dy), length=m))

    # create the empty transferfunction arrays
    transferfunction = similar(wave.ψ)

    # fill the transferfunction arrays with the proper values
    @. transferfunction = exp(-1im * π * λ * distance * (fx^2 + fy'^2))

    # shift the transferfunction and return thhe PropTf object
    return PropTf(fftshift(transferfunction))
end


"""
    calculate!(wave::Wavem proptff::PropTf)

Calculate the PropTf.

This function calculates the light field in a given direction. After the
calculation the changed LightField object is returned.
"""
function calculate!(wave::Wave, proptf::PropTf)
    # wave        ...   some kind of wave struct
    # proptf      ...   propagation object

    @info "Calculate transfer function propagation..."

    # Fouriertransform the input wave
    Ψ = fft(fftshift(wave.ψ))

    # apply the transferfunction
    Ψ .*= proptf.transferfunction

    # calculate the inverse fouriertransform
    wave.ψ = ifftshift(ifft(Ψ))
end


