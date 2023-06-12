using Revise
using GravityTools




#--------------------------------------------------------------------------------------------------------------
test_params = TestParameters(
    Var(name = "var1", min = 0.0, max = 4.0, N = 5, range_rule=:lin),
    Var(name = "var2", min = 0.0, max = 3.0, N = 5, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

kernels = [SimpleKernel()]

ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = true,
#    DiffUnit(:val1, diff = 0.1),
#    ContourUnit(:val1, contours = [0.5])
    DiffContourUnit(:val1, diff = 0.05, contours = [0.5])
    )

test = TestFramework(
    test_params,
    kernels,
    ref_sets
)





#--------------------------------------------------------------------------------------------------------------
test_params = TestParameters(
    Var(name = "BETA0", min = -6.0, max = 6.0, N = 9, range_rule=:lin),
    Var(name = "ALPHA0", min = -1e-4, max = -1e-1, N = 9, range_rule=:log),
    [Var(name="EOS", value="MPA1"),
    Var(name="COMP_TYPE", value="WD")],
    RangeVariable[]
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

    tempo_kernel = TempoKernel(tsets)

ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = true,
    DiffContourUnit(:chisqr, max = 10.0, diff = 1.0, contours = [4.0])
    )

test = TestFramework(test_params, [tempo_kernel], ref_sets)



#--------------------------------------------------------------------------------------------------------------

tempo_sets = TempoSettings(
    version = "tempo",
    par_file_init = "J2222-0137_T1_DDSTG_DMX.par",
    tim_file = "J2222-0137_T1.tim",
    add_flag = "-c -j -ni npulses.dat",
    fit_XPBDOT = false,
    iters = [
        [TP("GAIN", 1), TP("NITS", 3)],
        [TP("GAIN", 1), TP("NITS", 3)]]
    )

tempo_framework

test_params = TestParameters(
    Var(name = "BETA0", min = -6.0, max = 6.0, N = 9, range_rule=:lin),
    Var(name = "ALPHA0", min = -1e-4, max = -1e-1, N = 9, range_rule=:log),
    [Var(name="EOS", value="MPA1"),
    Var(name="COMP_TYPE", value="WD")],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = true,
    DiffContourUnit(:chisqr, max = 10.0, diff = 1.0, contours = [4.0])
    )

test_framework

calculate!(test_framework, tempo_framework)


