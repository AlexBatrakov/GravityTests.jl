#5--------------------------------------------------
#4--------------------------------------------------
#3----------------------------------------------------------------
#2------------------------------------------------------------------------------
#1--------------------------------------------------------------------------------------------

abstract type AbstractFramework end

calculate!(inputpool) = nothing

#function input_parameters(fr::AbstractFramework)
#    return input_parameters(typeof(fr))
#end

function initialize(fr_type::Type{T}) where {T <: AbstractFramework}
    return fr_type(a, b, NaN)
end

function input_parameters(inputpool, fr::AbstractFramework)
    return map(param_name -> getfield(inputpool, param_name), input_parameters(typeof(fr)))
end



#function initialize(inputpool, fr_type::Type{T}) where {T <: AbstractFramework}
#    return fr_type(a, b, NaN)
#end

#struct GeneralInputPool
#
#end


input_parameters(object) = fieldnames(typeof(object))

function update_input!(inputpool, object)
    for field in input_parameters(object)
        if haskey(inputpool, field)
            setfield!(object, field, getindex(inputpool, field))
#        elseif input_parameters(getfield(object, field)) != ()
#            update_input!(inputpool, getfield(object, field))
        end
    end
    return object
end




