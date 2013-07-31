#File:mesh.py
from utils import *
from Tkinter import *

def generateMesh(soilPos, rigidPos, nx, ny):
	"""Generates system mesh"""
	meshX= []
	meshY= []
	geoDensFunc = []
	airDensFunc = []
	rigDensFunc = []
	wallSwitch = []
	n = nx*ny
	#gera a malha em si
	for i in range(nx):
		for j in range(ny):
			meshX.append(i+1)
			meshY.append(j+1)
	#gera a funcao de densidade do solo
	for i in range(n):
		if (meshX[i] >= soilPos[0]) and (meshX[i] <= soilPos[1]) and (meshY[i] >= soilPos[2]) and (meshY[i] <= soilPos[3]):
			geoDensFunc.append(1.0)
		else:
			geoDensFunc.append(0.0)	
	#gera a funcao de densidade do corpo rigido
	for i in range(n):
		if (meshX[i] >= rigidPos[0]) and (meshX[i] <= rigidPos[1]) and (meshY[i] >= rigidPos[2]) and (meshY[i] <= rigidPos[3]):
			rigDensFunc.append(1.0)
		else:
			rigDensFunc.append(0.0)	
	#gera a funcao de densidade do ar
	for i in range(n):
		if (geoDensFunc[i] == 0.0) and (rigDensFunc == 0.0):
			airDensFunc.append(1.0)
		else:
			airDensFunc.append(0.0)	
	#gera o controle de parede	
	for i in range(n):	
		wallSwitch.append(0)
		
	return meshX, meshY, geoDensFunc, rigDensFunc, airDensFunc, wallSwitch

           
def drawMesh(nx, ny, meshX, meshY, geoDensFunc, rigDensFunc, largura, altura):
    """Draws system mesh"""   
    
    n= nx*ny #numero total de dados
    #largura = 10 #largura dos retangulos da malha em x
    #altura = 10 #largura dos reatangulos da malha em y
    w = largura*nx #largura da janela
    h = largura*ny #altura da janela
    canvas = Canvas(width=w, height=h, bg = 'white') #fundo do desenho
    canvas.pack(expand=YES, fill=BOTH)
    
    ##desenha o solo    
    for i in range(n):
    	xCoord = largura*nx - meshX[i]*largura
    	yCoord = altura*ny - meshY[i]*altura
    	
    	if (geoDensFunc[i] <=0):
    	    canvas.create_rectangle(xCoord,yCoord,xCoord+largura,yCoord+altura,
    	    width=1,fill='green')
        else:
            canvas.create_rectangle(xCoord,yCoord,xCoord+largura,yCoord+altura,
    	    width=1,fill='red')
    
    ##desenha o corpo rigido
    for i in range(n):
    	xCoord = largura*nx - meshX[i]*largura
    	yCoord = altura*ny - meshY[i]*altura
    	
    	if (geoDensFunc[i] <=0) and (rigDensFunc[i] <=0):
    	    canvas.create_rectangle(xCoord,yCoord,xCoord+largura,yCoord+altura,
    	    width=1,fill='green')
        elif (rigDensFunc[i] >=1.0):
            canvas.create_rectangle(xCoord,yCoord,xCoord+largura,yCoord+altura,
    	    width=1,fill='grey')	    
    
    mainloop()
    
    	
	

