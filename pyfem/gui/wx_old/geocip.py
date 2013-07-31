#!/usr/bin/python
# File: geocip.py


from Tkinter import *
from tkMessageBox import *
import input_window as iw
from os import *
from utils import *
from mesh import *
from time import *
from easygui import *



class GeoCip( Frame ):
    """Main Window of GEOCIP program"""

    def __init__( self ):
        """Create four buttons, pack them and bind events"""

        Frame.__init__( self )
        self.pack( expand = NO, fill = BOTH )
        self.master.title( "Programa GEOCIP (Simulacao de Grandes Deformacoes)" )
        
        showinfo( "Bem-Vindo ao GEOCIP", 
        'Caso queira simplesmente testar o programa, clique no botao "Escrever Dados/Gerar Malha"' + 
        'e depois em "Processar Dados". Sera rodado um exemplo.' )
        
 

        #variaveis de instancia
        self._NX = 50 #: Number of meshes (x-direction)
        self._NY = 20#: Number of meshes (y-direction)
        self._CFL = 5.0e-2 #: Courant Number (Courant-Friedrich-Levy)
        self._DTMAX = 1.0e-3 #: Maximum value of time increment
        self._DX = 0.01#: Mesh size (x-direction)
        self._DY = 0.01#: Mesh size (y-direction)
        self._TOTALTIME = 3.0#: Total calculation time
        self._GOUT = 1.0e-1#: Time interval for output
        self._ISCREEN = 20#: Time interval for output on screen
        self._GX = 0.0#: Gravity acceleration (x-direction)
        self._GY = 9.81#: Gravity acceleration (y-direction)
        self._ROS = 2000.0#: Density of geomaterial
        self._ROA = 1.25#: Density of air
        self._ROR = 100.0#: Density of rigid body
        self._VIS = 1.0#: Viscosity of geomaterial
        self._VIA = 1.0e-5#: Viscosity of air
        self._VIR = 1.0e-3#: Viscosity of rigid body
        self._IBISWT = 1#: Switch for viscosity (0:Newton 1:Bingham)
        self._BINA = 500.0#: Cohesion
        self._BINB = 0.0#: Internal friction angle
        self._ETAMAX = 1.0e10#: Maximum value of viscosity
        self._ETAMIN = 1.0#: Minimum value of viscosity
        self._GMIN = 1.0e-10#: Minimum value of shear strain rate
        self._YURIN = 0.0#: Initial velocity of rigid body (x-direction)
        self._YVRIN = 2.0#: Initial velocity of rigid body (y-direction)
        self._YRRIN = 0.0#: Initial angular velocity of rigid body
        self._MATRIXSW = 2#: Switch for matrix solver (1:PCGS 2:BICGS-TAB)

        self._SoilNXI = 1
        self._SoilNXF = 50
        self._SoilNYI = 1
        self._SoilNYF = 10

        self._RBNXI = 23
        self._RBNXF = 28
        self._RBNYI = 11
        self._RBNYF = 13

        # create input data button
        self.genDataButton = Button( self, text = "Entrar Dados Gerais",
           command = self.enterGenData )
        self.genDataButton.bind( "<Enter>", self.rolloverEnter )
        self.genDataButton.bind( "<Leave>", self.rolloverLeave )
        #self.genDataButton.pack( side = LEFT, padx = 5, pady = 5 )        
        self.genDataButton.grid( row = 0, column = 0 )        

        # create input data button
        self.matDataButton = Button( self, text = "Entrar Dados dos Materiais",
           command = self.enterMatData )
        self.matDataButton.bind( "<Enter>", self.rolloverEnter )
        self.matDataButton.bind( "<Leave>", self.rolloverLeave )
        #self.matDataButton.pack( side = LEFT, padx = 5, pady = 5 )
        self.matDataButton.grid( row = 0, column = 1 )       

        # create input data button
        self.meshDataButton = Button( self, text = "Entrar Dados da Malha",
           command = self.enterMeshData )
        self.meshDataButton.bind( "<Enter>", self.rolloverEnter )
        self.meshDataButton.bind( "<Leave>", self.rolloverLeave )
        #self.meshDataButton.pack( side = LEFT, padx = 5, pady = 5 )
        self.meshDataButton.grid( row = 1, column = 0 )        

        # create processing button
        self.writeDataButton = Button( self, text = "Escrever Dados/Gerar Malha",
           command = self.writeData )
        self.writeDataButton.bind( "<Enter>", self.rolloverEnter )
        self.writeDataButton.bind( "<Leave>", self.rolloverLeave )
        #self.writeDataButton.pack( side = LEFT, padx = 5, pady = 5 )
        self.writeDataButton.grid( row = 1, column = 1 )        
        # create processing button
        self.procDataButton = Button( self, text = "Processar Dados",
           command = self.procData )
        self.procDataButton.bind( "<Enter>", self.rolloverEnter )
        self.procDataButton.bind( "<Leave>", self.rolloverLeave )
        #self.procDataButton.pack( side = LEFT, padx = 5, pady = 5 )
        self.procDataButton.grid( row = 2, column = 0 )        

        #create post-processing button
        self.visResButton = Button( self, text = "Visualizar Resultado",
           command = self.visRes )
        self.visResButton.bind( "<Enter>", self.rolloverEnter )
        self.visResButton.bind( "<Leave>", self.rolloverLeave )
        #self.visResButton.pack( side = LEFT, padx = 5, pady = 5 )
        self.visResButton.grid( row = 2, column = 1 )        


    def enterGenData( self ):
        #showinfo( "Entrada de dados", "Insira os dados nos campos correspondentes" )
        values = []
        values = iw.genPropsWindow()
        self._NX = long(values[0])
        self._NY = long(values[1])
        self._CFL = float(values[2])
        self._DTMAX = float(values[3])
        self._DX = float(values[4])
        self._DY = float(values[5])
        self._TOTALTIME = float(values[6])
        self._GOUT = float(values[7])
        self._ISCREEN = float(values[8])
        self._GX = float(values[9])
        self._GY = float(values[10])
        self._MATRIXSW = long(values[11])
        
    def enterMatData( self ):
        #showinfo( "Entrada de dados", "Insira os dados nos campos correspondentes" )
        values = []
        values = iw.matPropsWindow()
        self._ROS = float(values[0])
        self._ROA = float(values[1])
        self._ROR = float(values[2])
        self._VIS = float(values[3])
        self._VIA= float(values[4])
        self._VIR = float(values[5])
        self._BINA = float(values[6])
        self._BINB = float(values[7])
        self._ETAMAX = float(values[8])
        self._ETAMIN = float(values[9])
        self._GMIN = float(values[10])
        self._YURIN = float(values[11])
        self._YVRIN = float(values[12])
        self._YRRIN = float(values[13])

    def enterMeshData( self ):
        #showinfo( "Entrada de dados", "Insira os dados nos campos correspondentes" )
        values = []
        values = iw.meshDataWindow()

        self._SoilNXI = long(values[0])
        self._SoilNXF = long(values[2])
        self._SoilNYI = long(values[1])
        self._SoilNYF = long(values[3])

        self._RBNXI = long(values[4])
        self._RBNXF = long(values[6])
        self._RBNYI = long(values[5])
        self._RBNYF = long(values[7])
        
        
        
    def writeData( self ):
        showinfo( "Escrevendo dados em arquivo.", "Gerar arquivos data.txt e mesh.txt" )
        #gera o arquivo data.txt
        inputData = []
        inputData = ['%10d%10d'%(self._NX,self._NY),
                   '%10.2e%10.2e%10.2e%10.2e'%(self._CFL,self._DTMAX,self._DX,self._DY),
                   '%10.2e%10.2e%10d%10.2e%10.2e'%(self._TOTALTIME,self._GOUT,self._ISCREEN,self._GX,self._GY),
                   '%10.2e%10.2e%10.2e'%(self._ROS,self._ROA,self._ROR),
                   '%10.2e%10.2e%10.2e'%(self._VIS,self._VIA,self._VIR),
                   '%10d%10.2e%10.2e%10.2e%10.2e%10.2e'%(self._IBISWT,self._BINA,self._BINB,self._ETAMAX,self._ETAMIN,self._GMIN),
                   '%10.2e%10.2e%10.2e'%(self._YURIN,self._YVRIN,self._YRRIN),
                   '%10d'%(self._MATRIXSW)]
        #escreve em arquivo
        writef("data.txt",inputData)
        
        #gera a malha - arquivo malha.txt
        meshX, meshY, geoDensFunc, rigDensFunc, airDensFunc, wallSwitch = generateMesh([self._SoilNXI,self._SoilNXF,self._SoilNYI,self._SoilNYF],
        [self._RBNXI,self._RBNXF,self._RBNYI,self._RBNYF],self._NX,self._NY)

        n = self._NX*self._NY

        data = []

        for i in range(n):
            data.append('%5d%5d%10.3f%10.3f%10.3f%5d'%(meshX[i],meshY[i],geoDensFunc[i],
            airDensFunc[i],rigDensFunc[i],wallSwitch[i]))
               
        writef("mesh.txt",data)
        
        showinfo( "Visualizar a malha", "Visualizar a malha gerada" )

        #desenha a malha
        largura = long(1024/self._NX)
        altura = largura 
        drawMesh(self._NX,self._NY, meshX, meshY, geoDensFunc, rigDensFunc, largura, altura)
       
        
    def procData( self ):
        showinfo( "Processamento...", "Clique em OK e aguarde o fim do processamento. " )
        ###################################################
        ### Comandos dependentes do sistema operacional ###
        ################################################### 
        system('./cip_proc')
                       
        ##gera nome dos diretorios onde serao armazenados os arquivos de saida
        t = ctime(time()) 
        tp = t.split()
        t = ctime(time()) 
        tp = t.split()
        dir1 = 'dados'+tp[4]+tp[1]+tp[2]+"-"+tp[3].replace(':','_')
        dir2 = 'imgs'+tp[4]+tp[1]+tp[2]+"-"+tp[3].replace(':','_')
        dir3 = 'datain'+tp[4]+tp[1]+tp[2]+"-"+tp[3].replace(':','_')
        ##cria os diretorios e copia os arquivos
        system('rm -rf anim')
        system('mkdir '+dir1)
        system('mkdir '+dir2)
        system('mkdir '+dir3)
        system('mkdir anim')
        system('cp *.bmp anim')
        system('mv *.bmp '+dir2)
        system('mv *.DAT '+dir1)
        system('mv *.txt '+dir3)
        
        showinfo( "Processamento encerrado!", 
        "Os arquivos de entrada de dados estao armazenados no diretorio" + dir3
        +" e os de saida estao armazenados nos diretorios " + dir1 +"," + dir2)

                
    def visRes( self):
        #showinfo( "Processamento...", "Clique em OK e aguarde o fim do processamento. " )
        ###################################################
        ### Comandos dependentes do sistema operacional ###
        ################################################### 
        fps = integerbox(message='Entre o intervalo de tempo entre os quadros da animacao.', title='Numero de fps para a animacao', argDefault=100, argLowerBound=5, argUpperBound=300)

        system('animate -delay '+ str(fps) + ' anim/*.bmp')
        

    def rolloverEnter( self, event ):
        event.widget.config( relief = GROOVE )

    def rolloverLeave( self, event ):
        event.widget.config( relief = RAISED )




def main():
    GeoCip().mainloop()

if __name__ == "__main__":
    main()

############################################################################
### (C) Copyright 2008 by Carlos Eduardo Veras Neves                       #
### Programa de Pos-Graduacao em Geotecnia da Universidade de Brasilia     #
###                                                                        #
### DISCLAIMER: The authors and publisher of this book have used their     #
### best efforts in preparing the book. These efforts include the          #
### development, research, and testing of the theories and programs        #
### to determine their effectiveness. The authors and publisher make       #
### no warranty of any kind, expressed or implied, with regard to these    #
### programs or to the documentation contained in these books. The authors #
### and publisher shall not be liable in any event for incidental or       #
### consequential damages in connection with, or arising out of, the       #
### furnishing, performance, or use of these programs.                     #
############################################################################



