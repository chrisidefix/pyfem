# Include PyFEM libraries
from pyfem.mesh import *
from pyfem.fem import *
from pyfem.equilib import *

# Define mesh

# ... Solid elements
block = Block2D()
block.make_box([0,0], [40, 2])
block.set_divisions(4,3)
block.set_quadratic(True)

# ... Generating the mesh
mesh = Mesh([block])
mesh.generate()

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)
#domain.nodes.sort_in_y()

# Setting element model
domain.elems.solids.set_elem_model(EqElasticSolid({"E": 10.0E6, "nu": 0.0}))

# Setting boundary conditions
#domain.faces.with_x(0.0).set_brys({"ux": 0})
#domain.nodes.with_x(0.0).with_y(1.0).set_brys({"ux": 0, "uy": 0})
domain.faces.with_x(0.0).set_brys({"ux": 0, "uy": 0})
domain.faces.with_y(2.0).set_brys({"ty": -0.24})

# Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.set_incs(1)
domain.solver.set_scheme("NR")
domain.solver.set_precision(1E-3)
domain.solver.set_plane_stress(True)
domain.solver.set_precision(1E-3)
tnodes = domain.elems.nodes.with_y(1.0)
tnodes.sort_in_x()
domain.solver.track(tnodes, "tnodes")
domain.solver.solve()

# Write output file
domain.solver.write_output()


#for ip in domain.elems[0].ips:
#    print ip

data = []
for i, e in enumerate(domain.elems):
    emodel = e.elem_model
    nn = 12
    dd = [-1+j*2.0/nn for j in range(nn+1)]

    DD = []
    tau = []

    for d in dd:
        R = [d, 0.0]
        D = emodel.ips[0].mat_model.stiff()
        B, detJ = emodel.calcB(R, emodel.coords())
        U = emodel.U()
        sig = dot(D, dot(B, U))
        tau.append(-sig[3]/1.4142) # sxy
        DD.append(i*10 + (d+1)/2.*10)
    data.append((DD, tau))

an_data = []
DD = [0,10,20,30,40]
for D in DD:
    tau = 1.5*0.24*(40-D)/2.0
    an_data.append(tau)


from plot import *
import pylab

pylab.subplot(111)
pylab.title("Graph")

l0, = pylab.plot(data[0][0], data[0][1], 'bo-', markeredgecolor='blue', markerfacecolor='white', markersize=6, color='b')
pylab.plot(data[1][0], data[1][1], 'bo-', markeredgecolor='blue', markerfacecolor='white', markersize=6, color='b')
pylab.plot(data[2][0], data[2][1], 'bo-', markeredgecolor='blue', markerfacecolor='white', markersize=6, color='b')
pylab.plot(data[3][0], data[3][1], 'bo-', markeredgecolor='blue', markerfacecolor='white', markersize=6, color='b')

l1, = pylab.plot(DD, an_data, '-', color='k')

l2, = pylab.plot(tnodes.data_book[-1].get_col("dist"), tnodes.data_book[-1].get_col("sxy",-1.0),'k^', markersize=7)

pylab.legend((l0,l1,l2),("Unsmoothed shear stress", "Beam theory", "Proposed smoothing"), 'upper right', shadow=True)

pylab.grid(True)
pylab.show()
