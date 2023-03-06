lin_var = RangeVariable(name = "BETA0", min = -6.0, max = 6.0, N = 9, range_rule=:lin)
log_var = RangeVariable(name = "ALPHA0", min = -1e-4, max = -1e-1, N = 9, range_rule=:log)

test = GeneralTest(
    Var(name="EOS", value="MPA1"),
    Var(name="COMP_TYPE", value="WD"),
    Var(name = "ALPHA0", min = -1e-4, max = -1e-1, N = 9, range_rule=:log),
    Var(name = "BETA0", min = -6.0, max = 6.0, N = 9, range_rule=:lin)
    )