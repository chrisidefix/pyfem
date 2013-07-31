try: import sys; sys.path.insert(0, "c:\\Dropbox\\codes\\pyfem")
except: pass

# Include PyFEM libraries
from pyfem.mesh    import *
from pyfem.fem     import *
from pyfem.equilib import *

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

# Setting element types and parameters
elem_model = EqDruckerPragerSolid({"alpha": 0.05, "kappa": 0.1, "E": 100.0, "nu": 0.25})
domain.elems.set_elem_model(elem_model)
domain.elems.set_state(sxx=0.0, syy=0.0, szz=0.0)

# Setting boudary conditions
base_nodes = domain.nodes.filter(lambda n: n.z==0)
top_nodes  = domain.nodes.filter(lambda n: n.z==0.5)
lat_nodes  = domain.nodes.filter(lambda n: n.x==0 or n.x==1 or n.y==0 or n.y==1)

domain.nodes.with_z(0.0).set_brys(ux=0, uy=0, uz=0)
domain.nodes.with_z(0.5).set_brys(uz=-0.033)
domain.nodes.with_x_in(0.0, 1.0).set_brys(ux=0, uy=0)
domain.nodes.with_y_in(0.0, 1.0).set_brys(ux=0, uy=0)

#Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.set_verbose(True)
domain.solver.set_incs(20)
#domain.elems[0].set_save_history(True)
domain.solver.track(domain.elems[0], 'elem_0')
domain.solver.set_track_per_inc(True)
domain.solver.set_scheme("FE")

domain.solver.solve()

domain.solver.set_incs(20)
top_nodes.set_brys({"uz": 0.008})
domain.solver.solve()

print "Residue: ", domain.solver.residue

domain.solver.write_output()

domain.elems[0].plot("J1", "srJ2D")
domain.elems[0].plot("ezz", "szz")
