
module ElectronPropagation

# export the different kind of waves
export ElectronBeam, LaserBeam, loadintensity

# export helping functions
export zeropadding!, removezeropadding!

# export the different kinds of setup possibilities
export PropTf, PhaseImprint, Aperture, Lens
export PropDirect, Edge, Setup, propagation!

# create the abstract wave type
abstract type Wave end

# define some physical constants
global const c = 299792458           # the speed of light
global const ħ = 1.054571817e-34     # the reduced planck constant
global const m_e = 9.1093837015e-31  # the electron mass
global const q = 1.602176634e-19     # electron charge
global const ε_0 = 8.8541878128e-12  # vacuum permitivity

# include all the needed code
include("./electronBeam.jl")
include("./laserBeam.jl")
include("./components.jl")
include("./prep.jl")

# end module
end
