using Revise
using GravityTools

input = GravityTools.GeneralInput([5.0, 6.0], [Var(name="var1", value=5.0), Var(name="var2", value=6.0)])

output = GravityTools.GeneralOutput([11.0, 30.0], [Var(name="sum", value=11.0), Var(name="mul", value=30.0)])

kernel = SimpleKernel(
    nothing,
    input,
    output
)