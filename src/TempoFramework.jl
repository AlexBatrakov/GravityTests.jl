
#--------------------------------------------------------------------------------------------------------------
# tempo parameters

mutable struct TempoParameter{T1, T2, T3}
    name::String
    name_symbol::Symbol
    line::String
    value::T1
    flag::T2
    uncertainty::T3
end

TP = TempoParameter

function TempoParameter(line::String)
    line_split = split(line)
    n = length(line_split)
    line_parsed = parse_tparam_field.(line_split)
    line_parsed_types = typeof.(line_parsed)
    if n >= 3 && line_parsed_types[1:3] == [String, String, String]
        n_name = 3
        name =        String(line_split[1] * " " * line_split[2] * " " * line_split[3])
        name_symbol = Symbol(name)
    else
        n_name = 1
        name =        String(line_split[1])
        name_symbol = Symbol(name)
    end

    value = n_name < n ? parse_tparam_field(line_split[n_name + 1]) : nothing

    flag = n_name + 1 < n ? parse_tparam_field(line_split[n_name + 2]) : nothing

    uncertainty = n_name + 2 < n ? parse_tparam_field(line_split[n_name + 3]) : nothing

    tparam = TempoParameter(name, name_symbol, line, value, flag, uncertainty)
    
    update_tparam_line!(tparam)

    return tparam
end

function TempoParameter(name, value, flag=nothing, uncertainty=nothing)
    name_symbol = Symbol(name)
    tparam = TempoParameter(name, name_symbol, "", value, flag, uncertainty)
    update_tparam_line!(tparam)
    return tparam
end

TempoParameter(var::ValueVariable) = TempoParameter(var.name, var.value)

function Base.show(io::IO, tparam::TempoParameter)
    indent = get(io, :indent, 0)
    print(io, " "^indent, tparam.name)
    print(io, " ", tparam.value)
    if tparam.flag !== nothing
        print(io, " ", tparam.flag)
    end
    if tparam.uncertainty !== nothing
        print(io, " ", tparam.uncertainty)
    end
	return nothing
end

function update_tparam_line!(tparam::TempoParameter)
    n_name = 20
    n_value = 27
    n_flag = 6
    n_uncertainty = 27
    line = tparam.name * " "^(n_name - length(tparam.name))
    value_str = string(tparam.value)
    line *=  value_str * " "^(n_value - length(value_str))
    flag_str = string(tparam.flag)
    line *= tparam.flag !== nothing ? flag_str * " "^(n_flag - length(flag_str)) : ""
    uncertainty_str = string(tparam.uncertainty)
    line *= tparam.uncertainty !== nothing ? uncertainty_str * " "^(n_uncertainty - length(uncertainty_str)) : ""
    tparam.line = line
    return line
end

function parse_tparam_field(value_str)
    n = length(value_str)
    value_int64 = tryparse(Int64, value_str)
    if value_int64 !== nothing
        return value_int64::Int64
    end
    value_float64 = tryparse(Float64, value_str)
    if (value_float64 !== nothing) && ((value_float64 > 0.0 && n <= 20) || (value_float64 < 0.0 && n <= 21))
        return value_float64::Float64
    end
    return String(value_str)::String
end

function update_tparam!(tparam::TempoParameter; value=tparam.value, flag=tparam.flag, uncertainty=tparam.uncertainty)
    tparam.value = value
    tparam.flag = flag
    tparam.uncertainty = uncertainty
    tparam.line = update_tparam_line!(tparam)
    return tparam
end

#TempoParameter(name::String, value::T, flag::Int64=-1, uncertainty::Float64=0.0) where {T} = TempoParameter{T}(name, Symbol(name), value, flag, uncertainty)
#TempoParameter(name_symbol::Symbol, value::T, flag::Int64=-1, uncertainty::Float64=0.0) where {T} = TempoParameter{T}(String(name), name_symbol, value, flag, uncertainty)

#--------------------------------------------------------------------------------------------------------------
# tempo par files

mutable struct TempoParFile
    name::String
    tparams::Dict{Symbol,TempoParameter}
    order::Vector{Symbol}
end

TempoParFile(name::String) = TempoParFile(name, Dict{Symbol,TempoParameter}(), Vector{Symbol}())

function Base.show(io::IO, par_file::TempoParFile)
    println(io, "Tempo parameter file $(par_file.name): ")
    for (i, name_symbol) in enumerate(par_file.order)
        #print(IOContext(io, :indent => indent+4), par_file.tparams[name_symbol])
        print("    ", par_file.tparams[name_symbol].line)
        if i < length(par_file.order)
            print("\n")
        end
    end
	return nothing
end

function read_par_file(par_file::TempoParFile)
    par_file.order = Vector{Symbol}()
    open(par_file.name, "r") do file_in
        for line in eachline(file_in)
            if startswith(line, "C ") || startswith(line, "c ")
                continue
            end
            tparam = TempoParameter(line)
            par_file.tparams[tparam.name_symbol] = tparam
            push!(par_file.order, tparam.name_symbol)
        end
    end
    return par_file
end

function write_par_file(par_file::TempoParFile, name_out=par_file.name)
    open(name_out, "w") do file_out
        for name_symbol in par_file.order
            tparam = par_file.tparams[name_symbol]
            println(file_out, tparam.line)
        end
    end
    return par_file
end

function update_par_file(par_file::TempoParFile, tparam::TempoParameter)
    return par_file
end

function update_par_file()

end
#=
struct TempoIteration
    nits::Int64
    gain::Float64
    tparams::Vector{TempoParameter}
end

function Base.show(io::IO, iter::TempoIteration)
    indent = get(io, :indent, 0)
    println(io, " "^indent, "Number of iterations: ", iter.nits)
    println(io, " "^indent, "GAIN value for convergence stage: ", iter.gain)
    println(io, " "^indent, "Parameters: ", iter.tparams)
    if tparam.flag != -1
        print(io, " ", tparam.flag)
    end
    if tparam.uncertainty != 0.0
        print(io, " ", tparam.uncertainty)
    end
	return nothing
end
=#

#--------------------------------------------------------------------------------------------------------------
# tempo settings

struct TempoSettings
    version::String
    par_file_init::String
    par_file_work::String
    tim_file::String
    add_flag::String
    fit_XPBDOT::Bool
    iters::Vector{Vector{TempoParameter}}
end

function Base.show(io::IO, tsets::TempoSettings)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Tempo settings:")
    println(io, ' '^(indent + 4), "Version: ", tsets.version)
    println(io, ' '^(indent + 4), "Initial par file: ", tsets.par_file_init)
    println(io, ' '^(indent + 4), "Working par file: ", tsets.par_file_work)
    println(io, ' '^(indent + 4), "Working tim file: ", tsets.tim_file)
    println(io, ' '^(indent + 4), "Selected additional flags: ", tsets.add_flag)
    println(io, ' '^(indent + 4), "Fit PBDOT to GR value: ", tsets.fit_XPBDOT)
    for (i,iter) in enumerate(tsets.iters)
        println(io, ' '^(indent + 4), "Tempo parameters in iteration #$i:")
        for (j, tparam) in enumerate(iter)
            print(io, ' '^(indent + 8), "$tparam")
            if i != length(tsets.iters) || j != length(iter)
                print("\n")
            end
        end
    end
	return nothing
end

TempoSettings(;version, par_file_init, par_file_work=par_file_init[1:end-4]*"_work.par", tim_file, add_flag, fit_XPBDOT, iters=Vector{Vector{TempoParameter}}()) = TempoSettings(version, par_file_init, par_file_work, tim_file, add_flag, fit_XPBDOT, iters)

TempoSettings(args... ;version, par_file_init, par_file_work=par_file_init[1:end-4]*"_work.par", tim_file, add_flag, fit_XPBDOT) = TempoSettings(version, par_file_init, par_file_work, tim_file, add_flag, fit_XPBDOT, collect(args))

