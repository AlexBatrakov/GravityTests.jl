#--------------------------------------------------------------------------------------------------------------
# Some statistics stuff

lvl_1σ = quantile(Chisq(2), 0.682689492137086)
lvl_2σ = quantile(Chisq(2), 0.954499736103642)
lvl_3σ = quantile(Chisq(2), 0.997300203936740)
lvl_4σ = quantile(Chisq(2), 0.999936657516334)
lvl_5σ = quantile(Chisq(2), 0.999999426696856)
lvl_6σ = quantile(Chisq(2), 0.999999998026825)
lvl_7σ = quantile(Chisq(2), 0.999999999997440)

lvl_68CL = quantile(Chisq(2), 0.68)
lvl_90CL = quantile(Chisq(2), 0.90)
lvl_95CL = quantile(Chisq(2), 0.95)
lvl_99CL = quantile(Chisq(2), 0.99)
lvl_997CL = quantile(Chisq(2), 0.997)

#--------------------------------------------------------------------------------------------------------------
# Variables

abstract type AbstractVariable end

struct ValueVariable{T} <: AbstractVariable
    name::String
    name_symbol::Symbol
    value::T
end

ValueVariable(name::String, value::T) where {T} = ValueVariable{T}(name, Symbol(name), value)

ValueVariable(;name::String, value::T) where {T} = ValueVariable{T}(name, Symbol(name), value)

function Base.show(io::IO, var::ValueVariable)
    indent = get(io, :indent, 0)
    print(io, " "^indent, var.name)
    print(io, " ", var.value)
	return nothing
end

abstract type AbstractRangeRule end

function transformation(range_rule::T, min, max, N::Int64) where {T <: AbstractRangeRule}
end

struct LinRule <: AbstractRangeRule end

function transformation(range_rule::LinRule, min, max, N::Int64)
    lin_values = collect(LinRange(min, max, N))
    values = lin_values
    return lin_values, values
end

struct LogRule <: AbstractRangeRule
    sign::Int64
end

function transformation(range_rule::LogRule, min, max, N::Int64)
    min_log = log10(abs(min))
    max_log = log10(abs(max))
    lin_values = collect(LinRange(min_log, max_log, N))
    values = range_rule.sign .* 10.0 .^ lin_values
    return lin_values, values
end

struct RangeVariable{T <: AbstractRangeRule}
    name::String
    name_symbol::Symbol
    min::Float64
    max::Float64
    N::Int64
    range_rule::T
    lin_values::Vector{Float64}
    values::Vector{Float64}
    function RangeVariable{T}(name::String, min, max, N::Int64, range_rule::T) where {T <: AbstractRangeRule}
        name_symbol = Symbol(name)
        lin_values, values = transformation(range_rule, min, max, N)
        return new(name, name_symbol, min, max, N, range_rule, lin_values, values)
    end
end

RangeVariable(name::String, min, max, N::Int64, range_rule::T) where {T <: AbstractRangeRule} = RangeVariable{T}(name, min, max, N, range_rule)

function RangeVariable(;name::String, min=nothing, max=nothing, N=nothing, range_rule=nothing)
    if range_rule == :lin
        return RangeVariable(name, min, max, N, LinRule())
    elseif range_rule == :log
        return RangeVariable(name, min, max, N, LogRule(sign(min)))
    end
end

function Base.show(io::IO, var::RangeVariable)
    indent = get(io, :indent, 0)
    print(io, " "^indent, var.name)
    print(io, " [$(var.min), $(var.max)]")
    print(io, " $(var.N) $(typeof(var).parameters[1]) entries")
	return nothing
end

function Variable(;name::String, value=nothing, min=nothing, max=nothing, N=nothing, range_rule=nothing)
    if value === nothing
        return RangeVariable(name=name, min=min, max=max, N=N, range_rule=range_rule)
    else
        return ValueVariable(name=name, value=value)
    end
end

Var = Variable

#--------------------------------------------------------------------------------------------------------------
# Print utils
function print_array(io, arr)
    if length(arr) > 0
        for i in 1:length(arr)-1
            println(io, arr[i])
        end
        print(io, arr[end])
    else
        print(io)
    end
end

function print_array(io, arr::Vector{Float64})
    indent = get(io, :indent, 0)
    print(io, ' '^indent, arr)
end


#--------------------------------------------------------------------------------------------------------------
#=
function read_param!(param, alpha0, log10alpha0, beta0)
    if param.name == "alpha0"
        alpha0 = (min = param.min, max = param.max, N = param.N)
        log10alpha0 = (min = log10(abs(param.min)), max = log10(abs(param.max)), N = param.N)
    elseif param.name == "log10alpha0"
        log10alpha0 = (min = param.min, max = param.max, N = param.N)
        alpha0 = (min = -exp10(param.min), max = -exp10(param.max), N = param.N)
    elseif param.name == "beta0"
        beta0 = (min = param.min, max = param.max, N = param.N)
    end
    return alpha0, log10alpha0, beta0
end

function GeneralTest(;psrname::String, eosname::String="", alpha0=NaN, log10alpha0=NaN, beta0=NaN, param1, param2)
    if !isnan(alpha0)
        log10alpha0 = log10(abs(alpha0))
    elseif !isnan(log10alpha0)
        log10alpha0 = Float64(log10alpha0)
        alpha0 = -exp10(log10alpha0)
    end
    alpha0, log10alpha0, beta0 = read_param!(param1, alpha0, log10alpha0, beta0)
    alpha0, log10alpha0, beta0 = read_param!(param2, alpha0, log10alpha0, beta0)
    return GeneralTest(psrname, eosname, alpha0, log10alpha0, beta0, param1, param2)
end

struct GridSetttings
    N_refinement::Int64
    CL::Vector{Float64}
    contours::Vector{Float64}
    refinement_type::String
    delta_chisqr_max::Float64
    delta_chisqr_diff::Float64
    gr_in_chisqr::Bool
    function GridSetttings(N_refinement, CL, contours, refinement_type, delta_chisqr_max, delta_chisqr_diff, gr_in_chisqr)
        CL = typeof(CL) == Float64 ? [CL] : CL
        return new(N_refinement, isempty(CL) ? cdf(Chisq(2), contours) : CL, isempty(contours) ? quantile.(Chisq(2), CL) : contours, refinement_type, delta_chisqr_max, delta_chisqr_diff, gr_in_chisqr)
    end
end

function Base.show(io::IO, gsets::GridSetttings)
    println(io, "Grid settings:")
	println(io, "   Desired refinement level: ", gsets.N_refinement)
    println(io, "   Desired confidence levels: ", gsets.CL)
    println(io, "   Desired Δχ2 contours: ", gsets.contours)
    println(io, "   Type of the refinements: ", gsets.refinement_type)
    println(io, "   Maximum Δχ2 value: ", gsets.delta_chisqr_max)
    println(io, "   Maximum difference in Δχ2 value: ", gsets.delta_chisqr_diff)
    print(io,   "   Include GR in Δχ2 value: ", gsets.gr_in_chisqr)
	return nothing
end


GridSetttings(;N_refinement=1, CL=Float64[], contours=Float64[], refinement_type="nice", delta_chisqr_max=10.0, delta_chisqr_diff=1.0, gr_in_chisqr=false) = GridSetttings(N_refinement, CL, contours, refinement_type, delta_chisqr_max, delta_chisqr_diff, gr_in_chisqr)
=#

#=
function get_label(name)
    if name == "PBDOT"
        return L"\dot{P}_\mathrm{b}\, (10^{-12})"
    elseif name == "M2"
        return L"m_{\mathrm{c}}\,(M_\odot)"
    elseif name == "MTOT"
        return L"m_{\mathrm{tot}}\,(M_\odot)"
    elseif name == "GAMMA"
        return L"\gamma"
    elseif name == "XDOT"
        return L"\dot{x}\,(10^{-12} s/s)"
    elseif name == "OMDOT"
        return L"\dot{\omega}"
    elseif name == "COSI"
        return L"\cos i"
    elseif name == "DTHETA"
        return L"\delta_\theta"
    elseif name == "log10alpha0"
        return L"\log_{10}|\alpha_0|"
    elseif name == "beta0"
        return L"\beta_0"
    elseif name == "H3"
        return L"h_3"
    elseif name == "VARSIGMA"
        return L"\varsigma"
    end
end
=#

