# Include PyFEM libraries

import os,sys; sys.path.insert(0, os.getcwd()+"/../..")

from pyfem import *

# Mesh generation
# -----------------------------------------------------------------------------------

# Generate mesh
block = Block3D()
block.make_box((0,0,0),(1,6,1))

block.set_divisions (1,20,1)
#block.set_quadratic()

iblock = BlockInset  ()
punctual_model  = False
iblock.punctual = punctual_model
iblock.set_coords ([ (0.5, 2.045256, 0.2),  (0.5, 6.0 , 0.8)]) #sloped
#iblock.set_quadratic(False)

mesh = Mesh()
mesh.add_blocks(block, iblock)

alp = atan((0.8-0.2)/(6.0-2.045256))

# Mesh generation
mesh.generate  ()

# Setting domain
# -----------------------------------------------------------------------------------

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)

# Setting element types and parameters
Dm  = 0.4/pi
#As  = pi*Dm**2/4.0
As  = 0.005
C   = 10.0
phi = 30*pi/180.0
Es  = 210.0E3
tau_max = C + 100.*tan(phi)

domain.elems.solids.set_elem_model(EqElasticSolid(E=1.E4, nu=0.0))
domain.elems.solids.set_state     (sxx=-100.0, syy=-100.0, szz=-100.0)
domain.elems.lines .set_elem_model(EqElasticTruss(E=Es, A=As))

if punctual_model:
    domain.elems.joints.set_elem_model(EqMCPunctualJoint(Ks=100.0E3, Kn=100.0E3, Dm=Dm, C=C, phi=phi))
else:
    domain.elems.joints.set_elem_model(EqMohrCoulombJoint(Ks=100.0E3, Kn=100.0E3, Dm=Dm, C=C, phi=phi))

# Reinforcement nodes
bar_nodes = domain.elems.lines.nodes
bar_nodes.sort_in_y()

# Node to apply external force (last node)
hook_node = bar_nodes[-1]

# Defining load levels
load_levels = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 0.9999]
load_incs   = [ i-j for i,j in zip(load_levels[1:],load_levels)]

# Number of stages
nstages     = len(load_incs)
# Total load
tload       = (3.14159*Dm*4.0)*(C + 100.0*tan(phi))

# Setting solver 
domain.set_solver(SolverEq())
domain.solver.set_scheme("NR")
domain.solver.set_precision(1.0E-3)
domain.solver.track(bar_nodes, 'inter')
domain.solver.set_track_per_inc(False)

# Loop along load stages
for i in range(nstages):
    domain.elems.solids.nodes.set_bc(ux=0.0, uy=0.0, uz=0.0 )
    hook_node.set_bc(fz = sin(alp)*tload*load_incs[i])
    hook_node.set_bc(fy = cos(alp)*tload*load_incs[i])
    domain.solver.set_incs(4)
    domain.solver.solve()
    domain.solver.write_output()

print " total load = ", tload

# Plot results
# -----------------------------------------------------------------------------------
bar_nodes.plot("tau", coef=-1.0/tau_max, legend=load_levels[1:])

