# -*- coding: utf-8 -*-
# Teste análise de trelica usando Python e elementos finitos

# Importar as bibliotecas para trabalhar com elem. finitos
from pyfem import *

# Geracaoo da malha

malha = Mesh()

malha.set_ndim(2)  # duas dimensoees 

bloco1 = BlockLine()
bloco1.set_coords([[0,0], [1,0]])

bloco2 = BlockLine()
bloco2.set_coords([[0,0], [0,1]])

bloco3 = BlockLine()
bloco3.set_coords([[0,1], [1,0]])


malha.blocks.append(bloco1)
malha.blocks.append(bloco2)
malha.blocks.append(bloco3)

malha.generate()


# Análise via elementos finitos

dom = Domain()  # Criacaoo do domínio de elementos finitos
dom.load_mesh(malha)


# definição dos tipos de elementos e suas propriedades
tipo_elem = EqElasticBar(E=100000.0, A=0.01)   # E em kN/m2=kPa    e     A em m2

dom.elems.set_elem_model(tipo_elem)

# cond de contorno de deslocamento
dom.nodes.sub(x=0).sub(y=0).set_brys(ux=0, uy=0)
dom.nodes.sub(x=1).sub(y=0).set_brys(uy=0)

# cond de contorno de força
dom.nodes.with_y(1.0).set_brys(fx=10.0)

# Solução

dom.set_solver( SolverEq() ) # definicaoo do solver so para equiliibrio
dom.solver.solve()
dom.solver.write_output()






