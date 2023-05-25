# AstrophysicalObjects.jl is a Julia module that defines various types and structures related to astrophysical objects. This module provides a set of abstract types and concrete implementations to represent and model different astrophysical entities, including stellar objects, stellar systems, and associated quantities.

# The module includes the following key components:

# 1. Stellar Objects: Abstract types, such as `AbstractStellarObject`, represent the general concept of a stellar object, while concrete subtypes like `StellarObject` provide specific implementations. Stellar objects are characterized by properties such as type, mass, and quantities.

# 2. Stellar Object Quantities: Abstract types like `AbstractStellarObjectQuantities` define the general concept of quantities associated with a stellar object. Concrete subtypes, such as `GRStellarObjectQuantities` and `DEFStellarObjectQuantities`, represent specific quantity implementations for different gravity models, such as General Relativity and DEFGravity.

# 3. Double Stellar Object Quantities: Abstract types like `AbstractDoubleStellarObjectQuantities` capture the concept of quantities specific to double stellar objects, such as binary star systems. Concrete implementations, like `DEFDoubleStellarObjectQuantities`, provide details and properties related to DEFGravity.

# 4. Stellar Systems: The module includes an abstract type `AbstractStellarSystem` that represents the general concept of a stellar system. Concrete subtypes, such as `TripleSystem`, provide specific implementations for different types of stellar systems.

# These components allow users to create and manipulate astrophysical objects, associate them with relevant quantities, and model various stellar systems. The module provides a foundation for simulating and analyzing astrophysical phenomena, enabling researchers to study the behavior and characteristics of different objects and systems in a structured and flexible manner.

# AstrophysicalObjects.jl serves as a crucial building block for constructing higher-level frameworks and simulations in the field of astrophysics. It provides a convenient and modular approach to representing, organizing, and working with astrophysical objects, contributing to the advancement of research and understanding in the realm of astrophysics.

# Define abstract type for astrophysical objects
abstract type AbstractAstrophysicalObject <: AbstractGravityToolsType end

abstract type AbstractAstrophysicalObjectQuantities <: AbstractGravityToolsType end

# Define abstract type for stellar objects, which is a subtype of AbstractAstrophysicalObject
abstract type AbstractStellarObject <: AbstractAstrophysicalObject end

# Define abstract type for stellar object quantities
abstract type AbstractStellarObjectQuantities <: AbstractAstrophysicalObjectQuantities end

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

# Define the StellarObject struct with a type, mass, and quantities
mutable struct StellarObject{T <: Union{AbstractStellarObjectQuantities, Nothing}} <: AbstractStellarObject
    type::Symbol
    mass::Float64
    quantities::T  # Quantities field can hold T subtype or be empty (Nothing)
end

input_parameters(stellarObject::StellarObject) = (:type, :mass)

# Constructor for StellarObject 
StellarObject(type::Symbol=Symbol(""), mass::Float64=0.0; quantities::Union{AbstractStellarObjectQuantities, Nothing}=nothing) = StellarObject(type, mass, quantities)

# Constructor for StellarObject with gravity-specific quantities
function StellarObject(gravity_type::Type{<:AbstractGravity})
    quantities = StellarQuantities(gravity_type)  # Get the quantities based on gravity type
    StellarObject(quantities = quantities)
end

# Function to initialize the StellarObject based on the gravity type
function StellarObject(gravity_type::Type{<:AbstractGravity}, object::StellarObject)
    quantities = StellarQuantities(gravity_type)
    StellarObject(object.type, object.mass, quantities=quantities)
end

function update_StellarObject!(gravity_type::Type{<:AbstractGravity}, stellarObject::StellarObject)
    stellarObject.quantities = StellarQuantities(gravity_type)  # Get the quantities based on gravity type
end

# Define abstract type for stellar systems, which is a subtype of AbstractAstrophysicalObject
abstract type AbstractStellarSystem <: AbstractAstrophysicalObject end

# Define abstract type for stellar object quantities
abstract type AbstractStellarSystemQuantities <: AbstractAstrophysicalObjectQuantities end

# Define specific stellar system type TripleSystem
struct TripleSystem <: AbstractStellarSystem
end

# Define abstract type for PSR binary system parameters
abstract type PSRBinarySystemParameters <: AbstractAstrophysicalObjectQuantities end

# Define specific parameters type for Keplerian parameters
mutable struct KeplerianParameters{T <: Union{Float64, Measurement{Float64}}} <: PSRBinarySystemParameters
    P0::T
    T0::T
    e0::T
    omega0::T
    x0::T
end

KeplerianParameters() = KeplerianParameters(zeros(Float64,5)...)

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

PostKeplerianParameters() = PostKeplerianParameters(zeros(Float64,7)...)

# Define specific parameters type for extra binary parameters
mutable struct ExtraBinaryParameters{T <: Union{Float64, Measurement{Float64}}} <: PSRBinarySystemParameters
    mc::T
    q::T
    omegadot::T
    h3::T
    varsigma::T
    Delta_N::T
end

ExtraBinaryParameters() = ExtraBinaryParameters(zeros(Float64, 6)...)

# Define abstract type for double stellar object quantities
abstract type AbstractDoubleStellarObjectQuantities <: AbstractAstrophysicalObjectQuantities end 

# # Define specific double quantities type for DEFGravity
# mutable struct DEFDoubleStellarSystemQuantities <: AbstractDoubleStellarObjectQuantities
#     alphaA::Float64
#     betaA::Float64
#     kA::Float64
#     IA::Float64
#     alphaB::Float64
#     betaB::Float64
#     kB::Float64
#     IB::Float64
# end




# # Define the PSRBinarySystem struct with various parameters and objects
# mutable struct PSRBinarySystem{T1 <: AbstractStellarObject, T2 <: Union{AbstractDoubleStellarObjectQuantities, Nothing}} <: AbstractStellarSystem
#     comp_type::Symbol
#     mp::Float64
#     mc::Float64
#     K_params::KeplerianParameters{Float64}
#     pulsar::T1
#     companion::T1
#     Obj_quantities::T2
#     PK_params::PostKeplerianParameters{Float64}
#     X_params::ExtraBinaryParameters{Float64}
# end

# input_parameters(object::PSRBinarySystem) = (:comp_type, :mp, :mc, :K_params)

# # Constructor for StellarObject without quantities
# function PSRBinarySystem(comp_type::Symbol=Symbol(""), mp::Float64=0.0, mc::Float64=0.0)
#     K_params = KeplerianParameters()
#     pulsar = StellarObject(:NS, mp)
#     companion = StellarObject(comp_type, mc)
#     Obj_quantities = nothing
#     PK_params = PostKeplerianParameters()
#     X_params = ExtraBinaryParameters()
#     return PSRBinarySystem(comp_type, mp, mc, K_params, pulsar, companion, Obj_quantities, PK_params, X_params)
# end

mutable struct BinarySystemQuantities{T <: Union{Nothing, AbstractStellarObjectQuantities}} <: AbstractStellarSystemQuantities
    K_params::KeplerianParameters{Float64}
    PK_params::PostKeplerianParameters{Float64}
    X_params::ExtraBinaryParameters{Float64}
    SO_params::T
end

abstract type BinaryStellarObjectQuantities <: AbstractStellarObjectQuantities end

mutable struct GRBinaryStellarObjectQuantities
    IA::Float64
    IB::Float64
end

mutable struct DEFBinaryStellarObjectQuantities
    alphaA::Float64
    betaA::Float64
    kA::Float64
    IA::Float64
    alphaB::Float64
    betaB::Float64
    kB::Float64
    IB::Float64
end

BinaryStellarObjectQuantities(gravity_type::Type{GeneralRelativity}) = GRBinaryStellarObjectQuantities(0.0, 0.0)
BinaryStellarObjectQuantities(gravity_type::Type{DEFGravity}) = DEFBinaryStellarObjectQuantities(zeros(Float64,8)...)


# Define the BinarySystem struct with various parameters and objects
mutable struct BinarySystem{T1 <: AbstractStellarObject, T2 <: Union{Nothing, AbstractStellarSystemQuantities}} <: AbstractStellarSystem
    primary::T1
    companion::T1
    quantities::T2
end

property_mappings(obj::BinarySystem) = Dict{Symbol, Vector{Symbol}}(
    :prim_type => [:primary, :type],
    :mp => [:primary, :mass],
    :comp_type => [:companion, :type],
    :mc => [:companion, :mass]
)

input_parameters(binarySystem::BinarySystem) = (:prim_type, :mp, :comp_type, :mc)

# Constructor for BinarySystem 
BinarySystem(primary::T1=StellarObject(), companion::T1=StellarObject(); quantities::Union{AbstractStellarObjectQuantities, Nothing}=nothing) where {T1 <: AbstractStellarObject} = BinarySystem(primary, companion, quantities)

# Constructor for BinarySystem with gravity-specific quantities
function BinarySystem(gravity_type::Type{<:AbstractGravity})
    primary = StellarObject(gravity_type)
    companion = StellarObject(gravity_type)
    quantities = nothing  # Get the quantities based on gravity type
    BinarySystem(primary, companion, quantities)
end

# Function to initialize the BinarySystem based on the gravity type
function BinarySystem(gravity_type::Type{<:AbstractGravity}, binarySystem::BinarySystem)
    primary = StellarObject(gravity_type, binarySystem.primary)
    companion = StellarObject(gravity_type, binarySystem.companion)
    quantities = nothing
    BinarySystem(primary, companion, quantities)
end



