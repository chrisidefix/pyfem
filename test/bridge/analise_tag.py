from pyfem import *

# Generate mesh
malha = Mesh ()
malha.load_file("tag.vtk")

# setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(malha)
#domain.set_analysis_type("plane_strain")

#Setting element types and parameters
mat1 = EqPlasticBar (E=213800000, A=0.0086  , sig_max=294.0E6)
mat2 = EqPlasticBar (E=213800000, A=0.0204  , sig_max=294.0E6)
mat3 = EqPlasticBar (E=213800000, A=0.0066  , sig_max=294.0E6)
mat4 = EqPlasticBar (E=213800000, A=0.0076  , sig_max=294.0E6)
mat5 = EqPlasticBar (E=213800000, A=0.0146  , sig_max=294.0E6)
mat6 = EqPlasticBar (E=213800000, A=0.0122  , sig_max=294.0E6)

# Selections

barras_S1 = domain.elems.with_tag("S-1")
barras_S2 = domain.elems.with_tag("S-2")
barras_S3 = domain.elems.with_tag("S-3")
barras_S4 = domain.elems.with_tag("S-4")
barras_S5 = domain.elems.with_tag("S-5")
barras_S6 = domain.elems.with_tag("S-6")

#Material application
barras_S1.set_elem_model(mat1)
barras_S2.set_elem_model(mat2)
barras_S3.set_elem_model(mat3)
barras_S4.set_elem_model(mat4)
barras_S5.set_elem_model(mat5)
barras_S6.set_elem_model(mat6)

#Setting boundary conditions
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
#load = -10.0
#domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
#domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
#domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
#domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
#domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

domain.elems.set_body_force(-0.1)

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

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 3

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 4

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 5

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 6

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 7

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 8

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 9

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 10

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 11

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 12

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 13

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 14

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 15

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 16

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 17

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 18

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 19

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 20

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 21

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()

# Stage 22

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()


# Stage 23

#domain.nodes.clear_brys()
domain.nodes.with_y(0.0).with_x(0.0).set_brys(ux=0, uy=0.0)
domain.nodes.with_y(0.0).with_x(40.0).set_brys(uy=0.0)
load = -10.0
domain.nodes.with_y(7.5).with_x(6.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(13.333).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(20.000).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(26.667).set_brys(fy=load, fx=load)
domain.nodes.with_y(7.5).with_x(33.333).set_brys(fy=load, fx=load)

solver.solve()
solver.write_output()
