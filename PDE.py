import numpy as np
import matplotlib.pyplot as plt

from fenics import *

from TubularTuber import *
import SolutionStorage as SS

# only display messages of things that 
# may go boom later
set_log_level(30)

# defines mesh, boundaries

# inlet forcing
np.random.seed(1)
noise_level = 0.1
inlet_forcing = lambda t, freq, freq2: np.cos(freq * t)  #+ 0.5 * np.cos(freq2 * t)

# frequency of forcing
frequency = 10 # Hz
high_freq = 2000 # Hz
# maximum time, long enough to stabilize
Tmax = 0.1 # s

# averaging of end of chronologic signal
signal_averaging_length = 1000
# nondimensionalize
speed_of_sound = 343 # m / s
c = speed_of_sound * 39.37 # convert to in / s

# time scale
t0 = x0 / c

# nondimensional angular frequency
w_low = 2 * np.pi * frequency * t0
w_high = 2 * np.pi * high_freq * t0
# nondimensional maximum time
tmax = Tmax / t0

print("Angular frequency: ", w_low)
print("Oscillation period: ", 2 * np.pi / w_low)
print("Maximum time: ", tmax)

# time steps, 
nt_steps = 1000
# saving frequency
save_freq = 10
count_freq = 100

dt = tmax / nt_steps

print("Time step: ", dt)

print("Time steps per period: ", int(2 * np.pi / w_low / dt))

input("Press any key to continue...")

# functional analysis definitions
# build problem on mixed function space
V = FunctionSpace(mesh, "CG", 2)

# trial and test functions
# u - pressure, v = du/dt 
u = TrialFunction(V)
u_ = TestFunction(V)

# temporary storage
u_curr = Function(V)
u_next = Function(V)
v_curr = Function(V)
v_next = Function(V)

# boundary conditions
# tube bc is zero Neumann, no contribution
tube_bc = DirichletBC(V, Constant(0), tube_boundary)
inlet_bc = lambda t: DirichletBC(V, inlet_forcing(t, w_low, w_high), inlet_boundary)

# variational PDE formulation, newmark-beta method
# g, b - time stepping parameters
g = 0.5
b = 0.25
# equivalent to midpoint method, unconditionally stable


# IC
u_curr.assign(Constant(0.0))
v_curr.assign(Constant(0.0))

# storage for chronologic solutions
u_sols = []

exte = ExhaustExtractor()
ext = interpolate(exte, V)

# variational form, newmark-beta method, written as A[u, u_] = F
Au = u * u_ * dx + b * dt**2 * dot(grad(u_), grad(u)) * dx 
F1uv = u_curr * u_ * dx + dt * v_curr * u_ * dx - (1 - 2 * b) / 2 * dt ** 2 * dot(grad(u_curr), grad(u_)) * dx  
#        ext * v_curr * u_ * ds

Av = u * u_ * dx
F2uv = v_curr * u_ * dx - dt * g * dot(grad(u_), grad(u_next)) * dx - dt * (1 - g) * dot(grad(u_), grad(u_curr)) * dx

# iterate
for j in range(nt_steps):
    # update bcs
    #bcs = [tube_bc, inlet_bc(j * dt)]
    bcs = [inlet_bc(j * dt)]

    # bilinear form remains the same
    # forcing term updates automatically via c++ pointers


    # solve for next u variational problem and store in u_next
    solve(Au == F1uv, u_next, bcs=bcs)

    # solve for next v, given next u
    # solve(Av == F2uv, v_next)
    # update v
    v_bcs = [DirichletBC(V, -u_next.dx(1), outlet_boundary)]
    solve(Av == F2uv, v_next, bcs=v_bcs)

    # updates
    # append solution to list
    u_sols.append(Function(V, name="Pressure"))
    u_sols[j].assign(u_next)

    # increment solutions
    u_curr.assign(u_next)
    v_curr.assign(project(v_next, V))

    if j % count_freq == 0:
        print("Completed No. ", j)


# save data to file
xdmf_fname = "soln.xdmf"
xdmf = SS.build_xdmf_file(xdmf_fname)
xdmf = SS.write_soln_to_file(xdmf, u_sols, dt * np.arange(nt_steps), 
        save_freq=save_freq)

# final snapshot
u_square = list(map(lambda x: x * x, u_sols))
u_final = 2 * sum(u_square[-signal_averaging_length::2]) / (dt * signal_averaging_length)

# plotting
plt.figure()
p = plot(u_final)
plt.colorbar(p)
plt.savefig(folder + "final-snapshot-fs-{:d}.pdf".format(nfins))
#plt.show()

# extract output signal
output_signal = 2 * np.array(list(map(lambda u: assemble(u * ext * dx), 
    u_sols[-signal_averaging_length::2])))

time = np.arange(nt_steps)[-signal_averaging_length:] * dt

# save to file
np.savetxt(folder + "signal-output-fs-{:d}.txt".format(nfins), output_signal)

# plot
#plt.figure()
#plt.plot(time, output_signal)
#plt.xlabel("t")
#plt.ylabel("amplitude")
#plt.show()

