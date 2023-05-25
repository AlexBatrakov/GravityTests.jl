#5--------------------------------------------------
#4--------------------------------------------------
#3----------------------------------------------------------------
#2------------------------------------------------------------------------------
#1--------------------------------------------------------------------------------------------

abstract type AbstractGravityToolsType end

abstract type AbstractFramework <: AbstractGravityToolsType end 

calculate!(inputpool) = nothing


# function input_parameters(inputpool, fr::AbstractFramework)
#     return map(param_name -> getfield(inputpool, param_name), input_parameters(typeof(fr)))
# end


function print_input_parameters(object, level::Int = 0)
    indent = "    "  # Number of spaces for each level of indentation
    
    # Print the structure name and type
    name = string(repeat(indent, level), typeof(object).name.name)
    println(name)
    # Print the fields and their values
    for field in fieldnames(typeof(object))
        value = isdefined(object, field) ? getfield(object, field) : nothing
        if fieldnames(typeof(value)) != ()
            print_input_parameters(value, level + 1) # Recursive call for nested structures
        else
            fieldname = string(repeat(indent, level + 1), field)
            fieldvalue = value isa Symbol ? string(value) : repr(value)
            if field in input_parameters(object)
                println("$fieldname = $fieldvalue")
            end
        end
    end
end

function input_parameters(object::T) where {T <: AbstractGravityToolsType}
    paramspool = ()
    for field in fieldnames(typeof(object))
        value = isdefined(object, field) ? getfield(object, field) : nothing

        if fieldnames(typeof(value)) != ()
            paramspool = (paramspool..., input_parameters(value)...) # Recursive call for nested structures
        else

        end
    end
    return paramspool 
end


input_parameters(object) = ()

function update_object!(object, inputpool::Dict)
#    println("object = $object")
    for field in fieldnames(typeof(object))
#        println("field = $field")
        if haskey(inputpool, field) && field in input_parameters(object)
#            println("setfield!")
            setfield!(object, field, getindex(inputpool, field))
        elseif isdefined(object, field) && input_parameters(getfield(object, field)) != ()
#            println("update_object!")
            update_object!(getfield(object, field), inputpool)
        end
    end
    return object
end

function update_framework!(framework::AbstractFramework, inputpool::Dict)
    return update_object!(framework, inputpool)
end

function show_objecttree(io::IO, object, level::Int = 0)
    indent = "    "  # Number of spaces for each level of indentation
    
    # Print the structure name and type
    name = string(repeat(indent, level), typeof(object).name.name)
    println(io, name)
    
    # Print the fields and their values
    for field in fieldnames(typeof(object))
        value = isdefined(object, field) ? getfield(object, field) : nothing
        if fieldnames(typeof(value)) != ()
            show_objecttree(io, value, level + 1)  # Recursive call for nested structures
        else
            fieldname = string(repeat(indent, level + 1), field)
            fieldvalue = value isa Symbol ? string(value) : repr(value)
            println(io, "$fieldname = $fieldvalue")
        end
    end
end

function Base.show(io::IO, object::T) where {T <: AbstractGravityToolsType}
    show_objecttree(io, object)
end




abstract type PropertyMappings <: AbstractGravityToolsType end

property_mappings(obj::T) where {T<:AbstractGravityToolsType} = Dict{Symbol, Vector{Symbol}}()

function Base.getproperty(obj::T, prop::Symbol) where {T<:AbstractGravityToolsType}
    mappings = property_mappings(obj)
    if prop in keys(mappings)
        field_path = mappings[prop]
        value = obj
        for field in field_path
            value = getfield(value, field)
        end
        return value
    else
        return getfield(obj, prop)
    end
end

function Base.setproperty!(obj::T, prop::Symbol, value) where {T<:AbstractGravityToolsType}
    mappings = property_mappings(obj)
    if prop in keys(mappings)
        field_path = mappings[prop]
        parent = obj
        for field in field_path[1:end-1]
            parent = getfield(parent, field)
        end
        setfield!(parent, field_path[end], value)
    else
        setfield!(obj, prop, value)
    end
end