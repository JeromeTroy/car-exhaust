from fenics import *
from mshr import *
import matplotlib.pyplot as plt

from WaveBreaker import *

use_breaker = True
display_mesh = True

if use_breaker: 
    folder = "figures/fins/"
else:
    folder = "figures/one-plugged/"

# mesh resolution
res = 50

# boundary tolerance
btol = 1e-5

# parameters, lengths in inches
length = 18
diameter = 2.5
radius = 0.25
bracket_length = 1.5
engine_offshoot_length = 2.5

fin_length = 1.5
fin_depth = 0.125

# nondimensionalize
x0 = diameter

l = length / x0
d = diameter / x0
r = radius / x0
bl = bracket_length / x0
eol = engine_offshoot_length / x0

fl = fin_length / x0
fd = fin_depth / x0

# construct domain
eo_corner1 = Point(0, d / 2)
eo_corner2 = Point(eol + d, -d / 2)
eo_box = Rectangle(eo_corner1, eo_corner2)

ring_center = Point(eol + d + r, 0)
ring_splitter = Circle(ring_center, r + d) - Circle(ring_center, r)
outer_corner1 = Point(eol + d + r, r + d)
outer_corner2 = Point(eol + 2 * (d + r), -r - d)
outer_rectangle = Rectangle(outer_corner1, outer_corner2)
ring_splitter = ring_splitter - outer_rectangle

tube1_corner1 = Point(eol + d, -r - d)
tube1_corner2 = Point(eol + d + l, -r)
tube1 = Rectangle(tube1_corner1, tube1_corner2)
tube2_corner1 = Point(eol + d, r + d)
tube2_corner2 = Point(eol + d + l, r)
tube2 = Rectangle(tube2_corner1, tube2_corner2)

exhaust = eo_box + ring_splitter + tube1 + tube2

# fin construction parameters
sup_start = [eol + d + l, -r - 0.5 * d + 0.5 * fd]
sup_stop = [eol + d, -r - 0.5 * d - 0.5 * fd]
fin_start = [eol + 1.5 * d, -r - 0.5 * d]
length = l - d
nfins = 3

# construct fins
breaker = build_wave_breaker(fin_start, length, nfins, fl, fd, sup_start, sup_stop)

if use_breaker:
    exhaust = exhaust - breaker

# generate mesh
mesh = generate_mesh(exhaust, res)


# write details to file
fname = folder + "details.txt"
fout = open(fname, "w")

fout.write("mesh resolution, ")
fout.write(str(res))
fout.write("\n")

fout.write("No. fins, ")
fout.write(str(nfins))
fout.write("\n")

fout.write("Fin spacing, ")
fout.write(str(round(length / nfins, 3)))
fout.write("\n")

fout.write("dx, ")
fout.write(str(mesh.hmin()))
fout.write("\n")

fout.close()

# boundaries
def inlet_boundary(x, on_boundary):
    return on_boundary and x[0] < btol

def outlet_boundary(x, on_boundary):
    return on_boundary and abs(x[0] - (eol + d + l)) < btol and x[1] > 0
#def outlet_boundary(x, on_boundary):
#    return on_boundary and abs(x[0] - (eol + d + l)) < btol

def tube_boundary(x, on_boundary):
    return on_boundary and (not inlet_boundary(x, on_boundary)) and (not outlet_boundary(x, on_boundary))

# expression to extract solution on outbound boundary
class ExhaustExtractor(UserExpression):
    """
    Zero everywhere except at output of exhaust
    """

    def eval(self, value, x):
        """
        Evaluate function

        Input:
            value - output variable for result (replaces return statement)
            x - input for function

        Function evaluates to 1 on outbound boundary, 0 elsewhere
        """

        if abs(x[0] - eol - d - l) < btol and x[1] > 0:
            value[0] = 1.0
        else:
            value[0] = 0.0

    def value_shape(self):
        """
        What is the shape of the outpu
        In this case - scalar output, no shape
        """
        return ()

# test plot
plt.figure()
plot(mesh)
plt.savefig(folder + "mesh.pdf")
if display_mesh:
    plt.show()
