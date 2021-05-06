
# import statements
using FFTW
using ProgressMeter
using SharedArrays
using Distributed

"""
    Create the PropTf struct.

    This struct saves the transfer functions that are needed to calculate the
    lightfield at a given propagation distance with the transfer function
    approach. There is one transfer function saved for the positive propagation
    direction and one for the negative direction.
"""
struct PropTf <: Component
    transferfunction::Array{<:Complex}
end


"""
    Construct the PropTf struct.

    This method calculates the transfer functions that are saved in the struct
    and then creates the struct. For this it needs the propagation distance,
    the wavelnegth of the light field, and the geometry on which the lightfield
    is operating.
"""
function PropTf(wave::Wave, distance::Real)
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
    Struct for slow, but flexible propagation.

    This struct takes the new transverse coordinates as an input, aswell
    as the propagation distance that is desired.
"""
struct PropDirect <: Component
    x::Array{<:Real}       # Target x coordinates
    y::Array{<:Real}       # Target y coordinates
    d::Real                # propagation distance
end


"""
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
