#--------------------------------------------------------------------------------------------------------------
#Abstract kernel

abstract type AbstractKernel end

struct GeneralInput{T1, T2 <: ValueVariable}
    raw::Vector{T1}
    in::Vector{T2}
end

function Base.show(io::IO, input::GeneralInput)
    indent = get(io, :indent, 0)
    println(io, ' '^(indent), "Raw information:")
    print_array(IOContext(io, :indent => indent+4), input.raw)
    println(io, "\n", ' '^(indent), "Desired information:")
    print_array(IOContext(io, :indent => indent+4), input.in)
	return nothing
end

struct GeneralOutput{T1, T2 <: ValueVariable}
    raw::Vector{T1}
    out::Vector{T2}
end

function Base.show(io::IO, output::GeneralOutput)
    indent = get(io, :indent, 0)
    println(io, ' '^(indent), "Raw information:")
    print_array(IOContext(io, :indent => indent+4), output.raw)
    println(io, "\n", ' '^(indent), "Desired information:")
    print_array(IOContext(io, :indent => indent+4), output.out)
	return nothing
end

#--------------------------------------------------------------------------------------------------------------
#Simple kernel

struct SimpleKernel <: AbstractKernel
    sets
    input
    output
end

function Base.show(io::IO, kernel::SimpleKernel)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Simple kernel:")
    println(io, ' '^(indent+4), "Settings:")
    println(io, ' '^(indent+8), kernel.sets)
    println(io, ' '^(indent+4), "Input:")
    println(IOContext(io, :indent => indent+8), kernel.input)
    println(io, ' '^(indent+4), "Output:")
    println(IOContext(io, :indent => indent+8), kernel.output)
	return nothing
end

SimpleKernel() = SimpleKernel(nothing,[],[])

function target_function(kernel::SimpleKernel)
    var1, var2 = kernel.input.raw
    kernel.output["sum"] = var1 + var2
    kernel.output["mul"] = var1 * var2
end
