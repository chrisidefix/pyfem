# Include PyFEM libraries
import os,sys
sys.path.append(os.getcwd()+"/..")

from math import *
from pyfem.mesh import *
from pyfem.fem import *
from pyfem.equilib import *

# Generate mesh
block = Block3D()
block.set_coords([
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 6.0, 0.0,
    0.0, 6.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 6.0, 1.0,
    0.0, 6.0, 1.0,
    0.5, 0.0, 0.0,
    1.0, 4.0, 0.0,
    0.5, 6.0, 0.0,
    0.0, 4.0, 0.0,
    0.5, 0.0, 1.0,
    1.0, 4.0, 1.0,
    0.5, 6.0, 1.0,
    0.0, 4.0, 1.0,
    0.0, 0.0, 0.5,
    1.0, 0.0, 0.5,
    1.0, 6.0, 0.5,
    0.0, 6.0, 0.5 ])

block.set_divisions (1,20,1)

iblock = BlockInset  ()
iblock.set_coords    ([0.5, 2.05, 0.2,  0.5, 6.0, 0.8])

mesh = Mesh()
mesh.blocks.append  (block)
mesh.blocks.append  (iblock)

alp = atan((0.8-0.2)/(6.0-2.05))

# Mesh generation
mesh.generate  ()
mesh.write_mesh("tmesh.vtk")


#################################################################################################


# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)

# Setting element types and parameters
Dm  = 0.025
Pm  = Dm*3.14159
A   = 3.14159*Dm**2/4.0
#A   = 0.00
C   = 10.0
phi = 30*3.14159/180.0

domain.elems.solids.set_elem_model(EqElasticSolid({"E": 1.E4, "nu": 0.0}))
domain.elems.solids.set_state({"sxx": -100.0, "syy": -100.0, "szz": -100.0 })

domain.elems.lines .set_elem_model(EqPerfectPlasticBar({"E": 210.0E6, "A": A, "sig_max":500000}))
#domain.elems.lines .set_elem_model(EqElasticBar({"E": 210.0E6, "A": 0.005}))
domain.elems.line_joints.set_elem_model(EqMohrCoulombJoint({"Ks": 100.0E3, "Dm": Dm, "C": C, "phi": phi}))

# Nodes
bar_nodes = domain.elems.lines.nodes
bar_nodes.sort_in_y()
hook_node = bar_nodes[-1]

# Calculating loads
load_levels = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 0.999]
load_incs   = [ i-j for i,j in zip(load_levels[1:],load_levels)]

nstages     = len(load_incs)
tload       = (3.14159*Dm*4.0)*(C + 100.0*tan(phi))

# Setting solver 
domain.set_solver(SolverEq())
domain.solver.set_scheme("MNR")
domain.solver.set_scheme("NR")
domain.solver.set_precision(1.0E-3)
domain.solver.track(bar_nodes, 'inter')
#inter_data = domain.solver.track(bar_nodes, 'inter')
domain.solver.set_track_per_inc(False)

# Loop along load stages
for i in range(nstages):
    domain.nodes.clear_brys()
    domain.elems.solids.nodes.set_brys({ "ux": 0.0, "uy": 0.0, "uz": 0.0 })
    hook_node.set_brys({"fz": sin(alp)*tload*load_incs[i]})
    hook_node.set_brys({"fy": cos(alp)*tload*load_incs[i]})
    domain.solver.set_incs(4)
    domain.solver.solve()
    domain.solver.write_output()

print " total load = ", tload 
print "End of analysis"


#################################################################################################

#line_data = LineData([], [], ndiv=10)
#domain.solver.track(line_data)
#
#from pyfem.plot import *
#
#inter_data.plot()
#inter_data.plot(["tau"],[])
#inter_data[0, "ux"]
#
#bar_nodes.plot
#
#plot(node, ["",""], ["lab1","lab2"], ["xlab", "ylab"])
#plot_stages(bar_nodes, range(0,8), ["tau"]
#

from plot import *
import pylab

tau_max = C + 100.*tan(phi)

def normalize(l):
    for i in range(0,len(l)):
        l[i] = -l[i]/tau_max
    return l

tables = read_tables('inter.dat')

pylab.subplot(111)
l0, = pylab.plot(tables[0]["dist"], normalize(tables[0]["tau"]), 'k*-', linewidth=1.5, markersize=8)
l2, = pylab.plot(tables[1]["dist"], normalize(tables[1]["tau"]), 'kx-', linewidth=1.5, markersize=8)
l3, = pylab.plot(tables[2]["dist"], normalize(tables[2]["tau"]), 'kv-', linewidth=1.5, markersize=8)
l4, = pylab.plot(tables[3]["dist"], normalize(tables[3]["tau"]), 'k^-', linewidth=1.5, markersize=8)
l5, = pylab.plot(tables[4]["dist"], normalize(tables[4]["tau"]), 'ks-', linewidth=1.5, markersize=8)
l6, = pylab.plot(tables[5]["dist"], normalize(tables[5]["tau"]), 'kp-', linewidth=1.5, markersize=8)
l7, = pylab.plot(tables[6]["dist"], normalize(tables[6]["tau"]), 'ko-', linewidth=1.5, markersize=8)
pylab.xlabel('d', fontsize='large')
pylab.ylabel('tau')
pylab.tick_params(labelsize='large')
pylab.legend((l0,l2,l3,l4,l5,l6,l7),("0.20", "0.40", "0.60", "0.80", "0.90", "0.98", "1.00", ), 'lower right', shadow=True)
pylab.grid(True)
pylab.show()
