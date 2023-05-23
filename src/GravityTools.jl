module GravityTools

# Required libraries
using Distributions
using Measurements
using Distributed
using ProgressMeter

# Utility methods and functions
include("Utils.jl")

# Abstract Framework definitions
include("Frameworks.jl")

# Physics related methods and definitions
include("Physics.jl")

# Astrophysical Objects related methods and definitions
include("AstrophysicalObjects.jl")

# Astrophysical Framework methods and definitions
include("AstrophysicalFramework.jl")

# Export important types, functions, and abstract types
export DEF, GR, SimpleEOS, TabularKernel, Physics, StellarObject, AstrophysicalFramework, simulate!, calculate!, export_abstract_types

# Export abstract types separately
export_abstract_types() = begin
    export AbstractStellarObject, AbstractStellarObjectQuantities
    export AbstractAstrophysicalObject, AbstractAstrophysicalObjectQuantities, AbstractDoubleStellarObjectQuantities
end

# Mathematical constants and variables
export lvl_1σ, lvl_2σ, lvl_3σ, lvl_4σ, lvl_5σ, lvl_6σ, lvl_7σ, lvl_68CL, lvl_90CL, lvl_95CL, lvl_99CL
export LinRule, LogRule, RangeVariable, ValueVariable, Variable, Var

end # module
