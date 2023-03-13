#--------------------------------------------------------------------------------------------------------------
#Abstract kernel

abstract type AbstractKernel end

function Base.show(io::IO, kernel::T) where {T <: AbstractKernel}
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Simple kernel:")
    println(io, ' '^(indent+4), "Settings:")
    println(IOContext(io, :indent => indent+8), kernel.sets)
    println(io, ' '^(indent+4), "Input:")
    println(IOContext(io, :indent => indent+8), kernel.input)
    println(io, ' '^(indent+4), "Output:")
    println(IOContext(io, :indent => indent+8), kernel.output)
	return nothing
end

function clear_input!(kernel::T) where {T <: AbstractKernel}
    kernel.input = typeof(kernel.input)()
end
function clear_output!(kernel::T) where {T <: AbstractKernel}
    kernel.output = typeof(kernel.output)()
end

mutable struct GeneralInput{T1, T2 <: ValueVariable}
    int::Vector{T1}
    ext::Vector{T2}
end

GeneralInput{T1, T2}() where {T1, T2 <: ValueVariable} = GeneralInput{T1, T2}(Vector{T1}(undef,0), Vector{T2}(undef,0))

function Base.show(io::IO, input::GeneralInput)
    indent = get(io, :indent, 0)
    println(io, ' '^(indent), "Internal flow:")
    print_array(IOContext(io, :indent => indent+4), input.int)
    println(io, "\n", ' '^(indent), "External flow:")
    print_array(IOContext(io, :indent => indent+4), input.ext)
	return nothing
end

mutable struct GeneralOutput{T1, T2 <: ValueVariable}
    int::Vector{T1}
    ext::Vector{T2}
end

GeneralOutput{T1, T2}() where {T1, T2 <: ValueVariable} = GeneralOutput{T1, T2}(Vector{T1}(undef,0), Vector{T2}(undef,0))

function Base.show(io::IO, output::GeneralOutput)
    indent = get(io, :indent, 0)
    println(io, ' '^(indent), "Internal flow:")
    print_array(IOContext(io, :indent => indent+4), output.int)
    println(io, "\n", ' '^(indent), "External flow:")
    print_array(IOContext(io, :indent => indent+4), output.ext)
	return nothing
end

#--------------------------------------------------------------------------------------------------------------
#Simple kernel

mutable struct SimpleKernel <: AbstractKernel
    sets::Nothing
    input::GeneralInput{Float64, ValueVariable{Float64}}
    output::GeneralOutput{Float64, ValueVariable{Float64}}
end

function Base.show(io::IO, kernel::SimpleKernel)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Simple kernel:")
    println(io, ' '^(indent+4), "Input:")
    println(IOContext(io, :indent => indent+8), kernel.input)
    println(io, ' '^(indent+4), "Output:")
    println(IOContext(io, :indent => indent+8), kernel.output)
	return nothing
end


function SimpleKernel()
    input = GeneralInput{Float64, ValueVariable{Float64}}()
    output = GeneralOutput{Float64, ValueVariable{Float64}}()
    return SimpleKernel(nothing, input, output)
end

function target_function!(kernel::SimpleKernel)
    if !isempty(kernel.input.ext)
        var1, var2 = kernel.input.ext[1].value, kernel.input.ext[2].value
    elseif  !isempty(kernel.input.int)
        var1, var2 = kernel.input.int
    end
    kernel.output.int = [sin(var1^2*var2), var1 + var2, var1 * var2]
    names = ["out", "sum", "mul"]
    kernel.output.ext = [Var(name=names[i], value = kernel.output.int[i]) for i in eachindex(names)]
end