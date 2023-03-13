#--------------------------------------------------------------------------------------------------------------
# Test parameters

abstract type AbstractTestParameters end

struct TestParameters{T1 <: AbstractRangeRule, T2 <: AbstractRangeRule, T3 <: ValueVariable, T4 <: RangeVariable} <: AbstractTestParameters
    x::RangeVariable{T1}
    y::RangeVariable{T2}
    int::Vector{T3}
    ext::Vector{T4}
end

function Base.show(io::IO, params::TestParameters)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Test parameters:")
    println(io, ' '^(indent + 4), "X range parameter:")
    println(IOContext(io, :indent => indent+8), params.x)
    println(io, ' '^(indent + 4), "Y range parameter:")
    print(IOContext(io, :indent => indent+8), params.y)
    if !isempty(params.int)
        println(io, "\n", ' '^(indent + 4), "Internal value parameters:")
        print_array(IOContext(io, :indent => indent+8), params.int)
    end
    if !isempty(params.ext)
        println(io, "\n", ' '^(indent + 4), "External range parameters:")
        print_array(IOContext(io, :indent => indent+8), params.ext)
    end
	return nothing
end


#--------------------------------------------------------------------------------------------------------------
# Test Framework

struct TestFramework{T1 <: TestParameters, T2 <: AbstractKernel, T3}
    params::T1
    kernels::Vector{T2}
    data::T3
end

function TestFramework(params::T1, kernels::Vector{T2}, ref_sets::RefinementSettings) where {T1 <: TestParameters, T2 <: AbstractKernel}
    if isempty(params.ext)
        grid = Refinement2DGrid(params.x, params.y, ref_sets)
        return TestFramework(params, kernels, grid)
    end
end

function Base.show(io::IO, test::TestFramework)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Test framework:")
    println(IOContext(io, :indent => indent+4), test.params)
    println(io, ' '^(indent+4), "Kernels:")
    print_array(IOContext(io, :indent => indent+8), test.kernels)
    println(io, ' '^(indent+4), "Data:")
    println(IOContext(io, :indent => indent+8), test.data)
	return nothing
end


function calculate!(test::TestFramework)

    function target_function(x,y)
        
    end

    params_function!(grid) = nothing

    calculate_2DGrid!(tf.grid, target_function, params_function!)
    return tf
end


#--------------------------------------------------------------------------------------------------------------
# Tests

abstract type AbstractTest end

struct GeneralTest <: AbstractTest
    rparams::Vector{RangeVariable}
    vparams::Vector{ValueVariable}
end

function GeneralTest(args...)
    params = collect(args)
    rparams = Vector{RangeVariable}()
    vparams = Vector{ValueVariable}()
    for param in params
        if typeof(param) <: RangeVariable
            push!(rparams, param)
        elseif typeof(param) <: ValueVariable
            push!(vparams, param)
        end
    end
    return GeneralTest(rparams, vparams)
end

function Base.show(io::IO, test::GeneralTest)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "General test:")
    println(io, ' '^(indent + 4), "Ranged parameters:")
    for rparam in test.rparams
        println(IOContext(io, :indent => indent+8), rparam)
    end
    println(io, ' '^(indent + 4), "Valued parameters:")
    for i in 1:length(test.vparams)-1
        println(IOContext(io, :indent => indent+8), test.vparams[i])
    end
    print(IOContext(io, :indent => indent+8), test.vparams[end])
	return nothing
end



