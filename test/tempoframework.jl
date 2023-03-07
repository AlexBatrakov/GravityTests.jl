using Revise
using PyPlot

using Distributed
addprocs(4)
@everywhere using GravityTests

test = GeneralTest(
    Var(name="EOS", value="MPA1"),
    Var(name="COMP_TYPE", value="WD"),
    Var(name = "ALPHA0", min = -1e-4, max = -1e-1, N = 9, range_rule=:log),
    Var(name = "BETA0", min = -6.0, max = 6.0, N = 9, range_rule=:lin)
    )

tsets = TempoSettings(
    version = "tempo",
    par_file_init = "J2222-0137_T1_DDSTG_DMX.par",
    tim_file = "J2222-0137_T1.tim",
    add_flag = "-c -j -ni npulses.dat",
    fit_XPBDOT = false,
    iters = [
        [TP("GAIN", 1), TP("NITS", 3)],
        [TP("GAIN", 1), TP("NITS", 3)]]
    )

ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = true,
    DiffContourUnit(:chisqr, max = 10.0, diff = 1.0, contours = [4.0])
    )

tf = TempoFramework(test, tsets, ref_sets)

calculate!(tf)