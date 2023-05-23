# Physics.jl is a Julia module that provides a framework for representing physical systems and their properties in the field of astrophysics. The module defines several key components, including gravity models, equations of state, and physical kernels. It also includes structures to represent various astrophysical objects such as stellar objects and stellar systems.

# The Physics module is designed to facilitate the study and analysis of astrophysical phenomena by providing a flexible and extensible framework. It allows users to define and manipulate physical systems, incorporate different gravity models and equations of state, and compute various properties and quantities associated with astrophysical objects.

# By separating the different components into individual types and providing modular functions for their initialization and manipulation, Physics.jl offers a convenient and structured approach to working with astrophysical systems. This modular design allows for easy extension and customization, making it suitable for a wide range of astrophysical simulations, calculations, and analyses.

# Overall, Physics.jl provides a powerful toolkit for researchers, scientists, and enthusiasts working in the field of astrophysics, enabling them to model and explore complex physical systems with ease and flexibility.

# Define abstract type for gravity
abstract type AbstractGravity <: AbstractGravityToolsType end

# Define specific gravity type GeneralRelativity
struct GeneralRelativity <: AbstractGravity
end

const GR = GeneralRelativity

# Define specific gravity type DEFGravity
mutable struct DEFGravity <: AbstractGravity
    alpha0::Float64
    beta0::Float64
    function DEFGravity()
        return new()
    end
    function DEFGravity(alpha0, beta0)
        return new(alpha0, beta0)
    end
end

const DEF = DEFGravity


# Define abstract type for EOS
abstract type AbstractEOS <: AbstractGravityToolsType end

# Define specific EOS type SimpleEOS
mutable struct SimpleEOS <: AbstractEOS
    eos_name::Symbol
    function SimpleEOS(eos_name::Symbol)
        return new(eos_name)
    end
    function SimpleEOS()
        return new()
    end
end

# Define abstract type for physical kernel with gravity parameter
abstract type AbstractPhysicalKernel{G <: AbstractGravity} <: AbstractGravityToolsType end

# Define specific type for structure solver kernel
abstract type StructureSolverKernel{G} <: AbstractPhysicalKernel{G} end

# Define abstract type for tabular kernel with gravity parameter
abstract type AbstractTabularKernel{G} <: AbstractPhysicalKernel{G} end

# Define a parametric struct for tabular kernel with path_to_grids parameter
struct TabularKernel{PathToGrids}
    path_to_grids::PathToGrids
end

# Constructor for TabularKernel
TabularKernel(path_to_grid::String) = TabularKernel{Symbol(path_to_grid)}

# Function to retrieve path_to_grids from a TabularKernel type
path_to_grid(kernel::Type{T}) where {T <: TabularKernel}  = String(kernel.parameters[1])

# Define specific tabular kernel type for GeneralRelativity
mutable struct GRTabularKernel <: AbstractTabularKernel{GR}
    path_to_grids::String
    mass_grids
    function GRTabularKernel(path_to_grids)
        return new(path_to_grids)
    end
    # Fields specific to GRTabularKernel
end

# Define specific tabular kernel type for DEFGravity
mutable struct DEFTabularKernel <: AbstractPhysicalKernel{DEF}
    path_to_grids::String
    alpha0::Float64
    beta0::Float64
    eos_name::Symbol
    full_grids
    mass_grids
    function DEFTabularKernel(path_to_grids)
        return new(path_to_grids)
    end
end

# Define method to create GRTabularKernel from GeneralRelativity and TabularKernel type
TabularKernel(::Type{GeneralRelativity}, kernel::Type{T}) where {T <: TabularKernel} = GRTabularKernel(path_to_grid(kernel))

# Define method to create DEFTabularKernel from DEFGravity and TabularKernel type
TabularKernel(::Type{DEFGravity}, kernel::Type{T}) where {T <: TabularKernel} = DEFTabularKernel(path_to_grid(kernel))

# Define the Physics struct with gravity, EOS, and kernel parameters
struct Physics{G <: AbstractGravity, E <: AbstractEOS, K <: AbstractPhysicalKernel} <: AbstractGravityToolsType
    gravity::G
    eos::E
    kernel::K
end

# Constructor for Physics with gravity, EOS, and kernel type
function Physics(gravity::AbstractGravity, eos::AbstractEOS, kernel::Type{T}) where {T <: TabularKernel}
    return Physics(gravity, eos, TabularKernel(typeof(gravity), kernel))
end



