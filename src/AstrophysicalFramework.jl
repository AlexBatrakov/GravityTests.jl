# AstrophysicalFramework.jl is a Julia module that defines the `AstrophysicalFramework` struct, which serves as a central component for organizing and managing astrophysical simulations and analyses. This module provides a framework for combining physics models, astrophysical objects, and computational tools into a unified structure.

# The `AstrophysicalFramework` struct encapsulates the physics model and the astrophysical object within a single entity. It allows researchers and scientists to conveniently associate a specific physics model, such as a gravity model, with a particular astrophysical object, such as a stellar object or stellar system. This combination enables the study of physical phenomena and the simulation of astrophysical systems under specific physical conditions.

# The module includes a constructor for creating instances of the `AstrophysicalFramework` struct. The constructor takes the physics model and astrophysical object as arguments and initializes an `AstrophysicalFramework` object that encapsulates both. This allows for a straightforward and organized approach to setting up simulations and analyses.

# By providing a unified structure, AstrophysicalFramework.jl facilitates the management and manipulation of astrophysical systems. It enables researchers to easily access and utilize the underlying physics model and the associated astrophysical object for computations, simulations, and analysis tasks.

# Overall, AstrophysicalFramework.jl provides a convenient and modular framework for organizing and working with astrophysical simulations and analyses. It helps researchers leverage the power of physics models and astrophysical objects in a cohesive and structured manner, enabling more efficient and insightful studies in the field of astrophysics.


# Define the AstrophysicalFramework struct with a physics and object
struct AstrophysicalFramework{T1 <: Physics, T2 <: AbstractAstrophysicalObject} <: AbstractFramework
    physics::T1
    object::T2
    function AstrophysicalFramework{T1, T2}(physics::T1, astrophysicalObject::T2) where {T1 <: Physics, T2 <: AbstractAstrophysicalObject}
        return new{T1, T2}(physics, astrophysicalObject)
    end
end

# Outer constructor for AstrophysicalFramework
function AstrophysicalFramework(physics::T1, astrophysicalObject::T2) where {T1 <: Physics, T2 <: AbstractAstrophysicalObject}
    # Initialize the astrophysical object based on the provided physics object
    initialized_object = initialize_AstrophysicalObject(physics, astrophysicalObject)

     # Create and return an AstrophysicalFramework object
    return AstrophysicalFramework{T1, typeof(initialized_object)}(physics, initialized_object)
end


# Function to initialize the astrophysical object based on the physics object
function initialize_AstrophysicalObject(gravity_type::Type{<:AbstractGravity}, object::AbstractAstrophysicalObject)
    if object isa StellarObject
        StellarObject(gravity_type, object)
    elseif object isa AbstractStellarSystem
        StellarSystem(gravity_type, object)
    else
        error("Unsupported astrophysical object type: $(typeof(object))")
    end
end

