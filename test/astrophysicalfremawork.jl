# Create an instance of the gravity type DEFGravity
gravity = DEF()

# Create an instance of the EOS type SimpleEOS
eos = SimpleEOS()

# Create an instance of the kernel type TabularKernel with the specified path to grids
kernel = TabularKernel("path/to/grids")

# Create an instance of the Physics struct with the gravity, EOS, and kernel
physics = Physics(gravity, eos, kernel)


# Create an instance of the StellarObject struct with the specified type and mass
stellarObject = StellarObject()

# Create an instance of the AstrophysicalFramework struct with the physics and stellarObject
framework = AstrophysicalFramework(physics, stellarObject)

inputpool = Dict(:type => :WD, :mass => 1.4, :alpha0 => -1e-4, :beta0 => -4.0, :eos_name => :MPA1)

update_framework!(framework, inputpool)


binarySystem = BinarySystem(StellarObject(), StellarObject(), nothing)

framework = AstrophysicalFramework(physics, binarySystem)