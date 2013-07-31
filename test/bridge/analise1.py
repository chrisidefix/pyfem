from pyfem  import *

# Generate mesh
malha = Mesh()
#malha.load_file("C:/Users/Renato/Documents/UnB/PJTPC/Modelos de pontes/mesh.vtk")
malha.load_file("mesh.vtk")

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(malha)
#domain.set_analysis_type("plane_strain")

# Setting element types and parameters
mat1 = EqPlasticBar(E=200000, A=0.03, sig_max=200.0E3)
mat2 = EqPlasticBar(E=210000, A=0.03, sig_max=200.0E3)
mat3 = EqPlasticBar(E=220000, A=0.03, sig_max=200.0E3)

# Selections
vert_elems = domain.elems.with_dx(0.0)
hori_elems = domain.elems.with_dy(0.0)
slop_elems = domain.elems - vert_elems - hori_elems

#barras_tipo2 = domain.elems.with_tag("tipo2")

# Material application
vert_elems.set_elem_model(mat1)
hori_elems.set_elem_model(mat2)
slop_elems.set_elem_model(mat3)

#Setting initial conditions
#domain.elems.set_state(sa=0.1)

#Setting boundary conditions
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -50.0
domain.nodes.with_y(7.5).set_brys(fy=load)

#Setting solver and solving
solver = SolverEq()
solver.set_domain(domain)
solver.set_scheme("NR")
solver.set_incs(10)

#for i in range(5):
solver.solve()
solver.reset_displacements()

solver.write_output()

# Stage 2

domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).set_brys(fy=load)


solver.solve()
solver.write_output()

#    solver.trim_elems(min_sa=0.03)

