from easygui import *
from Tkinter import *

def genPropsWindow():
    msg = "Insira os dados gerais da simulacao em unidades coerentes"
    title = "Parametros gerais da simulacao"
    fieldNames = ["Numero de celulas na direcao x:",
              "Numero de celulas na direcao y:",
              "Valor do numero de Courant-Friedrich-Levy:",
              "Valor maximo do incremento de tempo:",
              "Tamanho da celula na direcao x:",
              "Tamanho da celula na direcao y:",
              "Tempo total de processamento:",
              "Intervalo de tempo para a saida de dados:",
              "Intervalo de tempo para a impressao dos dados na tela:",
              "Aceleracao da gravidade na direcao x:",
              "Aceleracao da gravidade na direcao y:",
              "Algoritmo de solucao do sistema linear (1:PCGS/2:BICGS-TAB):"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)

    # make sure that none of the fields was left blank
    while 1:
        if fieldValues == None: break
        errmsg = ""
        for i in range(len(fieldNames)):
            if fieldValues[i].strip() == "":
                errmsg = errmsg + ('"%s" nao pode ficar vazio.' % fieldNames[i])
        if errmsg == "": break # no problems found
        fieldValues = multenterbox(errmsg, title, fieldNames, fieldValues)

    #print "Reply was:", fieldValues

    return fieldValues

def matPropsWindow():
    msg = "Insira as propriedades dos materiais em unidades coerentes"
    title = "Propriedades dos materiais"
    fieldNames = ["Densidade do geomaterial:",
              "Densidade do ar:",
              "Densidade do corpo rigido:",
              "Viscosidade do geomaterial:",
              "Viscosidade do ar:",
              "Viscosidade do corpo rigido:",
              "Tipo de viscosidade (0:Newton/1:Bingham):",
              "Coesao:",
              "Angulo de Atrito interno (graus):",
              "Valor maximo de viscosidade:",
              "Valor minimo de viscosidade:",
              "Valor minimo de deformacao cisalhante:",
              "Velocidade inicial de corpo rigido na direcao x:",
              "Velocidade inicial de corpo rigido na direcao y:",
              "Velocidade inicial de corpo rigido na direcao z:",
              "Velocidade inicial angular de corpo rigido:"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)

    # make sure that none of the fields was left blank
    while 1:
        if fieldValues == None: break
        errmsg = ""
        for i in range(len(fieldNames)):
            if fieldValues[i].strip() == "":
                errmsg = errmsg + ('"%s" nao pode ficar vazio.' % fieldNames[i])
        if errmsg == "": break # no problems found
        fieldValues = multenterbox(errmsg, title, fieldNames, fieldValues)

    #print "Reply was:", fieldValues

    return fieldValues

def meshDataWindow():
    msg = "Insira a geometria da malha. Lembre-se que se faz referencia ao numero das celulas correspondentes"
    title = "Dados da Malha"
    fieldNames = ["Celula inicial para o solo em x:",
              "Celula inicial para o solo em y:",
              "Celula final para o solo em x:",
              "Celula final para o solo em y:",
              "Celula inicial para o corpo rigido em x:",
              "Celula inicial para o corpo rigido em y:",
              "Celula final para o corpo rigido em x:",
              "Celula final para o corpo rigido em y:"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)

    # make sure that none of the fields was left blank
    while 1:
        if fieldValues == None: break
        errmsg = ""
        for i in range(len(fieldNames)):
            if fieldValues[i].strip() == "":
                errmsg = errmsg + ('"%s" nao pode ficar vazio.' % fieldNames[i])
        if errmsg == "": break # no problems found
        fieldValues = multenterbox(errmsg, title, fieldNames, fieldValues)

    #print "Reply was:", fieldValues

    return fieldValues
