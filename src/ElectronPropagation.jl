
module ElectronPropagation

# export the different kind of waves
export ElectronBeam, LaserBeam, deBroglieWavelength

# export the different kinds of setup possibilities
export Free, PhaseImprint, Aperture, Lens
export Setup, propagation!

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

# end module
end
