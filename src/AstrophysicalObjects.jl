# AstrophysicalObjects.jl is a Julia module that defines various types and structures related to astrophysical objects. This module provides a set of abstract types and concrete implementations to represent and model different astrophysical entities, including stellar objects, stellar systems, and associated quantities.

# The module includes the following key components:

# 1. Stellar Objects: Abstract types, such as `AbstractStellarObject`, represent the general concept of a stellar object, while concrete subtypes like `StellarObject` provide specific implementations. Stellar objects are characterized by properties such as type, mass, and quantities.

# 2. Stellar Object Quantities: Abstract types like `AbstractStellarObjectQuantities` define the general concept of quantities associated with a stellar object. Concrete subtypes, such as `GRStellarObjectQuantities` and `DEFStellarObjectQuantities`, represent specific quantity implementations for different gravity models, such as General Relativity and DEFGravity.

# 3. Double Stellar Object Quantities: Abstract types like `AbstractDoubleStellarObjectQuantities` capture the concept of quantities specific to double stellar objects, such as binary star systems. Concrete implementations, like `DEFDoubleStellarObjectQuantities`, provide details and properties related to DEFGravity.

# 4. Stellar Systems: The module includes an abstract type `AbstractStellarSystem` that represents the general concept of a stellar system. Concrete subtypes, such as `TripleSystem`, provide specific implementations for different types of stellar systems.

# These components allow users to create and manipulate astrophysical objects, associate them with relevant quantities, and model various stellar systems. The module provides a foundation for simulating and analyzing astrophysical phenomena, enabling researchers to study the behavior and characteristics of different objects and systems in a structured and flexible manner.

# AstrophysicalObjects.jl serves as a crucial building block for constructing higher-level frameworks and simulations in the field of astrophysics. It provides a convenient and modular approach to representing, organizing, and working with astrophysical objects, contributing to the advancement of research and understanding in the realm of astrophysics.

# Define abstract type for astrophysical objects
abstract type AbstractAstrophysicalObject end

# Define abstract type for stellar objects, which is a subtype of AbstractAstrophysicalObject
abstract type AbstractStellarObject <: AbstractAstrophysicalObject end

# Define abstract type for stellar object quantities
abstract type AbstractStellarObjectQuantities end

# Define specific quantities type for General Relativity
mutable struct GRStellarObjectQuantities <: AbstractStellarObjectQuantities
    I::Float64
    GRStellarObjectQuantities() = new()
end

const GRQuantities = GRStellarObjectQuantities

# Define specific quantities type for DEFGravity
mutable struct DEFStellarObjectQuantities <: AbstractStellarObjectQuantities
    alphaA::Float64
    betaA::Float64
    kA::Float64
    IA::Float64
    DEFStellarObjectQuantities() = new()
end

const DEFQuantities = DEFStellarObjectQuantities

# Function to get the corresponding quantities type for a given gravity type
StellarQuantities(gravity_type::Type{GeneralRelativity}) = GRStellarObjectQuantities()
StellarQuantities(gravity_type::Type{DEFGravity}) = DEFStellarObjectQuantities()

# Define abstract type for double stellar object quantities
abstract type AbstractDoubleStellarObjectQuantities end

# Define specific double quantities type for DEFGravity
mutable struct DEFDoubleStellarObjectQuantities <: AbstractDoubleStellarObjectQuantities
    alphaA::Float64
    betaA::Float64
    kA::Float64
    IA::Float64
    alphaB::Float64
    betaB::Float64
    kB::Float64
    IB::Float64
end

# Define the StellarObject struct with a type, mass, and quantities
mutable struct StellarObject{T <: Union{AbstractStellarObjectQuantities, Nothing}}
    type::Symbol
    mass::Float64
    quantities::T  # Quantities field can hold T subtype or be empty (Nothing)
end

# Constructor for StellarObject without quantities
StellarObject(type::Symbol, mass::Float64) = StellarObject{Nothing}(type, mass, nothing)

# Constructor for StellarObject with gravity-specific quantities
function StellarObject(gravity_type::Type{<:AbstractGravity}, type::Symbol, mass::Float64)
    quantities = StellarQuantities(gravity_type)  # Get the quantities based on gravity type
    StellarObject(type, mass, quantities)
end

# Function to initialize the StellarObject based on the gravity type
function StellarObject(gravity_type::Type{<:AbstractGravity}, object::StellarObject)
    StellarObject(gravity_type, object.type, object.mass)
end

function update_StellarObject!(gravity_type::Type{<:AbstractGravity}, stellarObject::StellarObject)
    stellarObject.quantities = StellarQuantities(gravity_type)  # Get the quantities based on gravity type
end

# Define abstract type for stellar systems, which is a subtype of AbstractAstrophysicalObject
abstract type AbstractStellarSystem <: AbstractAstrophysicalObject end

# Define specific stellar system type TripleSystem
struct TripleSystem <: AbstractStellarSystem
end

# Define abstract type for PSR binary system parameters
abstract type PSRBinarySystemParameters end

# Define specific parameters type for Keplerian parameters
mutable struct KeplerianParameters{T <: Union{Float64, Measurement{Float64}}} <: PSRBinarySystemParameters
    P0::T
    T0::T
    e0::T
    omega0::T
    x0::T
end

# Define specific parameters type for Post-Keplerian parameters
mutable struct PostKeplerianParameters{T <: Union{Float64, Measurement{Float64}}} <: PSRBinarySystemParameters
    k::T
    gamma::T
    Pbdot_GW::T
    r::T
    s::T
    dtheta::T
    dr::T
end

# Define specific parameters type for extra binary parameters
mutable struct ExtraBinaryParameters{T <: Union{Float64, Measurement{Float64}}} <: PSRBinarySystemParameters
    mc::T
    q::T
    omegadot::T
    h3::T
    varsigma::T
    Delta_N::T
end

# Define the PSRBinarySystem struct with various parameters and objects
mutable struct PSRBinarySystem{T1 <: AbstractStellarObject, T2 <: AbstractDoubleStellarObjectQuantities}
    comp_type::Symbol
    mp::Float64
    mc::Float64
    K_params::KeplerianParameters
    pulsar::T1
    companion::T1
    Obj_quantities::T2
    PK_params::PostKeplerianParameters
    X_params::ExtraBinaryParameters
end
