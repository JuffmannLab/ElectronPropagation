
# import statements
using FFTW
using ProgressMeter
using SharedArrays
using Distributed


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


"""
    PropDirect(x::Vector{<:Real}, y::Vector{<:Real}, d::Real)::PropDirect

Return the PropDirect type.

Create the PropDirect type, that can be used to calculate the wavefunction
propagated a distance `d` in z direction, with the target coordinates
`x` and `y`. This calculation uses the direct method, and is thus really
slow. It also can happen that there are artefacts that appear periodically
in calculated image. This happens when the source and target coordinates do
not match to a for me unknown criterium (i have to investigate!!).

# Example
```jldoctest
julia> PropDirect(Array{Float64}(1:2), Array{Float64}(1:2), 1)
PropDirect([1.0, 2.0], [1.0, 2.0], 1)
```

See also: [`PropTf`](@ref), [`PhaseImprint`](@ref), [`Lense`](@ref),
[`Aperture`](@ref), [`Edge`](@ref)
"""
struct PropDirect <: Component
    x::Vector{<:Real}       # Target x coordinates
    y::Vector{<:Real}       # Target y coordinates
    d::Real                 # propagation distance
end


"""
    calculate!(wave::Wave, prop::PropDirect)

Calculate the propagation.

This function calculates the propagation, using direct integration
of the Fresnel integral. This method has the advantage that it is
flexible (you can change the coordinate system to change arbitrarly),
but the disadvantage that it doesn't use the Fourier transform, and such
is a lot slower.
"""
function calculate!(wave::Wave, prop::PropDirect)
    # wave   ...   wave object
    # prop   ...   PropDirect struct

    @info "Calculate direct Propagation..."

    # get data from the wave and prop struct
    u = wave.x
    v = wave.y
    du = abs(u[1]-u[2])
    dv = abs(v[1]-v[2])
    x = prop.x
    y = prop.y
    z = prop.d

    # calculate the wavenumber
    k = 2 * π / wave.λ

    # define the ouput array
    ψ = SharedArray{eltype(wave.ψ)}(size(x, 1), size(y, 1))

    # calculate the constant factors
    α = 1im * k / 2 / z
    β = exp(1im * k * z) / (1im * wave.λ * z)

    # integrate over the given wavefunction with the calculation
    @showprogress pmap(i->parallel_calculate(ψ,wave.ψ,x,y,u,v,du,dv,z,k,α,β,i),1:size(ψ,2))

    # set the wave struct to the new conditions
    wave.x = x
    wave.y = y
    wave.ψ = convert(Array, ψ)
end


"""
    Parallelized function.

    This function is parallelized to accelerate the direct calculation of
    the light field.
"""
function parallel_calculate(ψ::SharedArray{<:Complex}, ψ0::Array{<:Complex},
                            x::Array{<:Real}, y::Array{<:Real},
                            u::Array{<:Real}, v::Array{<:Real},
                            du::Real, dv::Real, z::Real, k::Real,
                            α::Complex, β::Complex, i::Int)

    # iterate over the new array
    for j = 1:size(ψ, 1)
        # set the ψ value to 0
        ψ[j, i] = 0

        for l = 1:size(ψ0, 2)
            for m = 1:size(ψ0, 1)

                # calculate the new wavefunction fore each pixel
                # @inbounds removes the expression checking if the index is
                # inside of the array... if something goes wrong, things might
                # crash or so... don't let it go wrong!!!
                @inbounds ψ[j, i] += ψ0[m,l]*exp(α*((x[j]-u[m])^2+(y[i]-v[l])^2))
            end
        end

        # multiply the constant β and the pixelsize of u and v
        ψ[j, i] *= β * du * dv
    end
end
