using Revise
using GravityTests
using PyPlot

test = GeneralTest(
    psrname = "J1952+2630",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    )

tsets = TempoSettings(
    version = "tempo2",
    par_file_init = "J1952+2630_DDSTG.par",
    tim_file = "TOAs_pdev_puppi_fast_T2_23-26_efac_gauss.tim",
    add_flag = "-c -j -ni npulses.dat",
    fit_XPBDOT = false,
    iters = [
        [TP("GAIN", 1), TP("NITS", 3)],
        [TP("GAIN", 1), TP("NITS", 3)]]
    )

ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = false,
    DiffContourUnit(:chisqr, max = 10.0, diff = 1.0, contours = [4.0])
    )