# Include PyFEM libraries
from pyfem import *

# Define mesh

# ... Solid elements
block = Block3D()
block.make_box([0,0,0], [1.0, 40, 2])
block.set_divisions(1,4,3)
block.set_quadratic()

# Define domain and conditions

# ... Generating the mesh
mesh = Mesh([block])
mesh.generate()

# ... Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)

# ... Setting element model
elem_model = EqElasticSolid(E=20.0E6, nu=0.0)
domain.elems.solids.set_elem_model(elem_model)

# ... Setting boundary conditions
domain.faces.sub(x=0.0).set_bc(ux=0)
domain.faces.sub(x=1.0).set_bc(ux=0)
domain.faces.sub(y=0.0).set_bc(ux=0, uy=0, uz=0)
domain.faces.sub(z=2.0).set_bc(tz=-0.24)

# Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.set_incs(1)
domain.solver.set_scheme("NR")
domain.solver.set_precision(1E-3)
tnodes = domain.elems.nodes.sub(x=0.0).sub(z=1.0)
tnodes.sort_in_y()
domain.solver.track(tnodes, "tnodes")
domain.solver.solve()

#dl = domain.get_data_line([],[])
#dl.plot("syz")

# Write output file
domain.solver.write_output()

# Getting data from integration points
data = []
for i, e in enumerate(domain.elems):
    emodel = e.elem_model
    nn = 12
    dd = [-1+j*2.0/nn for j in range(nn+1)]

    DD = []
    sigl = []

    for d in dd:
        R = [-1.0, d, 0.0]
        D = emodel.ips[0].mat_model.stiff()
        B, detJ = emodel.calcB(R, emodel.coords())
        U = emodel.U()
        sig = dot(D, dot(B, U))
        sigl.append(-sig[4]/1.4142) # syz
        DD.append(i*10 + (d+1)/2.*10)
    data.append((DD, sigl))

#print sigl

import pylab

pylab.subplot(111)

pylab.plot(data[0][0], data[0][1], 'ko-', markerfacecolor='white', label='unsmoothed shear stress')
pylab.plot(data[1][0], data[1][1], 'ko-', markerfacecolor='white')
pylab.plot(data[2][0], data[2][1], 'ko-', markerfacecolor='white')
pylab.plot(data[3][0], data[3][1], 'ko-', markerfacecolor='white')
#pylab.plot(DD, sigl, 'ko-', markerfacecolor='white')
pylab.plot(tnodes.data_book[-1].get_col("dist"), tnodes.data_book[-1].get_col("syz",-1.0),'k^', label='proposed extrapolation')
pylab.plot(tnodes.data_book[-1].get_col("dist"), tnodes.data_book[-1].get_col("syz",+1.0),'k-', label='beam theory')

from pylab import rc
rc('text', usetex=True)
rc('text', usetex=False)
rc('font', family='serif')
pylab.legend(loc=1)
pylab.xlabel(r'Distance')
pylab.ylabel(r'Shear stress $\tau$')

pylab.grid(True)
pylab.show()



