using Revise
using GravityTools


test_params = TestParameters(
    Var(name = "BETA0", min = -6.0, max = 6.0, N = 9, range_rule=:lin),
    Var(name = "ALPHA0", min = -1e-4, max = -1e-1, N = 9, range_rule=:log),
    [Var(name="EOS", value="MPA1"),
    Var(name="COMP_TYPE", value="WD")],
    RangeVariable[]
)



#--------------------------------------------------------------------------------------------------------------
test_params = TestParameters(
    Var(name = "var1", min = 0.0, max = 4.0, N = 5, range_rule=:lin),
    Var(name = "var2", min = 0.0, max = 3.0, N = 5, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

test = TestFramework(
    test_params,
    SimpleKernel(),
    Float64[]
)