# Include libraries
import os,sys; sys.path.insert(0, os.getcwd()+'/../..')

from pyfem import *

# Generate mesh
mesh = Mesh()
bl0 = Block3D()
bl0.make_box((0,0,0),(1,1,10))
bl0.set_divisions(1,1,10)
#bl0 = Block3D(box=[(0,0,0),(1,1,10)], div=(1,1,10))
mesh.add_blocks(bl0)
mesh.generate()
mesh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)

# Setting element types and parameters
# input information
load   = -10
hd     = 10
E      = 5000.0
nu     = 0.25
k      = 1.0e-5
gammaw = 10.0
mv     = (1+nu)*(1-2*nu)/(E*(1-nu))
cv     = k/(mv*gammaw)
print "mv = " , mv, "cv = " , cv

emodel = HydromecLin(E=5000, nu=0.25, k=1.0e-5, gammaw=10.)
domain.elems.set_elem_model(emodel)

outnodes = domain.nodes.sub(x=0, y=0)
outnodes.sort_in_z()

#print domain.elems[0].elem_model.calcQh()
#exit()

#Setting solver 
domain.set_solver( SolverHydromec() )
domain.solver.set_scheme("MNR")
domain.solver.set_scheme("FE")
domain.solver.track(outnodes)

# Stage 1 (to get steady state)
domain.nodes.set_bc(ux=0, uy=0)
domain.nodes.sub(z=0 ).set_bc(uz=0)
domain.nodes.sub(z=10).set_bc(wp=0)
domain.solver.set_incs(3)
domain.solver.set_scheme("FE")

domain.solver.solve(1000000.0)
domain.solver.write_output()

# Stage 2 (load application)
domain.nodes.set_bc(ux=0, uy=0)
domain.nodes.sub(z=0 ).set_bc(uz=0)
domain.nodes.sub(z=10).set_bc(wp=0)
domain.nodes.sub(z=10).set_bc(fz=-load)
domain.solver.set_incs(1)
domain.solver.solve(10.0)
domain.solver.write_output()

#generating times increments
T  = [0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0]  # time
t  = [Ti*hd**2/cv for Ti in T]
Dt = [t[i] - t[i-1] for i in range(1,len(t))]

# Next stage (consolidation)
for dt in Dt:
    OUT('dt')
    domain.nodes.set_bc(ux=0, uy=0)
    domain.nodes.sub(z=0).set_bc(uz=0)
    domain.nodes.sub(z=10).set_bc(wp=0)
    domain.solver.set_incs(2)
    domain.solver.solve(dt)
    domain.solver.write_output()

print outnodes

# Plot results
outnodes.plot("d","wp")

