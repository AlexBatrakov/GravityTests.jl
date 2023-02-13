#--------------------------------------------------------------------------------------------------------------
# Single_core routines
using Revise
using GravityTests
using PyPlot


ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = false,
#    DiffUnit(:val1, diff = 0.1),
#    ContourUnit(:val1, contours = [0.5])
    DiffContourUnit(:val1, diff = 0.05, contours = [0.5])
    )

x = [0.0, 1.0, 2.0]
y = [0.0, 1.0, 2.0, 3.0]
target_function(x,y) = (val1 = sin(x^2*y),)
params_function!(grid) = nothing

grid_init = Refinement2DGrid(x, y, ref_sets)

@time precalculate_2DGrid(grid_init, target_function, params_function!)
grid_init.value[:val1]

@time grid_refined = refine_2DGrid(grid_init, target_function, params_function!)
grid_refined.value[:val1]

grid_array = Array{Refinement2DGrid}(undef, 9)
grid_array[1] = grid_init
for i in 2:9
    grid_array[i] = refine_2DGrid(grid_array[i-1], target_function, params_function!)
end

#--------------------------------------------------------------------------------------------------------------
# Parallel routines
using Revise
using GravityTests
using PyPlot
using Distributed

addprocs(8)

ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = true,
#    DiffUnit(:val1, diff = 0.1),
#    ContourUnit(:val1, contours = [0.5])
    DiffContourUnit(:val1, diff = 0.05, contours = [0.5])
    )

x = [0.0, 1.0, 2.0]
y = [0.0, 1.0, 2.0, 3.0]
@everywhere target_function(x,y) = (val1 = sin(x^2*y),)
@everywhere params_function!(grid) = nothing

grid_init = Refinement2DGrid(x, y, ref_sets)

@time precalculate_2DGrid(grid_init, target_function, params_function!)
grid_init.value[:val1]

@time grid_refined = refine_2DGrid(grid_init, target_function, params_function!)
grid_refined.value[:val1]

grid_array = Array{Refinement2DGrid}(undef, 9)
grid_array[1] = grid_init
for i in 2:9
    grid_array[i] = refine_2DGrid(grid_array[i-1], target_function, params_function!)
end