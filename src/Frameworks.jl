#5--------------------------------------------------
#4--------------------------------------------------
#3----------------------------------------------------------------
#2------------------------------------------------------------------------------
#1--------------------------------------------------------------------------------------------

abstract type AbstractGravityToolsType end

abstract type AbstractFramework <: AbstractGravityToolsType end 

calculate!(inputpool) = nothing


function input_parameters(inputpool, fr::AbstractFramework)
    return map(param_name -> getfield(inputpool, param_name), input_parameters(typeof(fr)))
end



input_parameters(object) = fieldnames(typeof(object))

function update_object!(object, inputpool::Dict)
    for field in input_parameters(typeof(object))
        if haskey(inputpool, field)
            setfield!(object, field, getindex(inputpool, field))
        elseif input_parameters(getfield(object, field)) != ()
            update_input!(inputpool, getfield(object, field))
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




