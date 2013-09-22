# Include PyFEM libraries
from pyfem import *

# Generate mesh
mesh = Mesh()
mesh.set_ndim(3)

block = Block3D()
block.make_box([0,0,0], [1,1,0.5])
block.set_divisions(2,2,2)
#block.set_quadratic(True)

mesh.blocks.append(block)
mesh.set_verbose(True)
mesh.generate()

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)

deg = pi/180.

# Setting element types and parameters
#elem_model = EqMohrCoulombSolid(c=0.01, phi=1*deg, E=100., nu=0.25)
elem_model = EqMohrCoulombSolid(c=0.084, phi=0.125, E=100., nu=0.25)
domain.elems.set_elem_model(elem_model)
domain.elems.set_state(sxx=0.0, syy=0.0, szz=0.0)

# Setting boudary conditions
top_nodes  = domain.nodes.sub(z=0.5)

domain.nodes.sub(z=0.0).set_bc(ux=0, uy=0, uz=0)
#domain.nodes.sub(z=0.5).set_bc(uz=-0.033)
domain.nodes.sub(x=1.0).sub(y=1.0).set_bc(fz=-0.4)
domain.nodes.sub(x=[0.0, 1.0]).set_bc(ux=0, uy=0)
domain.nodes.sub(y=[0.0, 1.0]).set_bc(ux=0, uy=0)

#Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.set_verbose(True)
domain.solver.set_incs(10)
#domain.elems[0].set_save_history(True)
domain.solver.track(domain.elems[0], 'elem_0')
domain.solver.set_track_per_inc(True)
domain.solver.set_scheme("FE")

domain.solver.solve()

domain.solver.set_incs(10)
top_nodes.set_bc({"uz": 0.008})
domain.solver.solve()

print "Residue: ", domain.solver.residue

domain.solver.write_output()

domain.elems[0].plot("J1", "srJ2D")
domain.elems[0].plot("ezz", "szz")
