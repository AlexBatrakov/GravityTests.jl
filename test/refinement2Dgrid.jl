using GravityTests


ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = true,
    DiffUnit(:chisqr, 0.0, 10.0, 0.1),
    ContourUnit(:chisqr, 0.0, 10.0, [0.5, 0.7])
    )

x = [0.0, 1.0, 2.0]
y = [0.0, 1.0, 2.0, 3.0]
func(x,y) = sin(x^2*y)
target_function(x,y) = (val1 = sin(x^2*y),)
params_function!(grid) = nothing

grid = Refinement2DGrid(x, y, ref_sets)
precalculate_2DGrid(grid, target_function, params_function!)