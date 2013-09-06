# -*- coding: utf-8 -*-
# Teste análise de trelica usando Python e elementos finitos

# Importar as bibliotecas para trabalhar com elem. finitos
from pyfem import *

# Geracaoo da malha

malha = Mesh()

bloco1 = BlockLine()
bloco1.set_coords([[0,0], [1,0]])

bloco2 = BlockLine()
bloco2.set_coords([[0,0], [0,1]])

bloco3 = BlockLine()
bloco3.set_coords([[0,1], [1,0]])

malha.add_blocks(bloco1, bloco2, bloco3)
malha.generate()

# Análise via elementos finitos
dom = Domain()  # Criação do domínio de elementos finitos
dom.load_mesh(malha)

# definição dos tipos de elementos e suas propriedades
tipo_elem = EqElasticBar(E=1000.0, A=1.0)   # E em kN/m2=kPa    e     A em m2
dom.elems.set_elem_model(tipo_elem)

# condições de contorno de deslocamento
dom.nodes.sub(x=0).sub(y=0).set_bc(ux=0, uy=0)
dom.nodes.sub(x=1).sub(y=0).set_bc(uy=0)

# cond de contorno de força
dom.nodes.sub(y=1).set_bc(fx=10.0)

# Solução
dom.set_solver( SolverEq() ) # definição do solver só para equilíbrio
dom.solver.solve()

# Check
N, E = dom.elems[0].elem_model.get_nodal_and_elem_vals()

print "Checking stresses:"
compare(N['sa'], [10.0, 10.0])

print "Checking strains:"
compare(N['ea'],[0.01,0.01])

#dom.solver.write_output()






