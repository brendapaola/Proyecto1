# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 04:21:21 2019

@author: brend
"""

import Bio
from Bio import SeqIO
from Bio import Entrez
import re
from Bio.Seq import Seq
import math
import pandas as pd
import scipy.stats as ss
import numpy as np 
import matplotlib.pyplot as plt
from timeit import timeit

#Actualizamos la matrix de estados acultos dados un conjunto de genes
def ActMatOcu(genes):
    r = len(genes)
    oculta = list(range(len(genes)))
    for i in range(len(genes)):
        oculta[i] =[]
        for j in range(len(genes)):
            oculta[i].append(2/3*1/(r))
        oculta[i].append(1/3)
        
    oculta.append([1/len(genes)]*(r+1))
    
    return(oculta)

def ActMatObs(genes): #Funcion que atualiza la matriz observada
    dist = list(range(len(genes)))

    for i in range(len(genes)):
        temp = ''
        for j in range(len(genes[i])):
            temp = temp + genes[i][j]
        aux = Seq(temp)        
        dist[i] = [aux.count('A'), aux.count('C'), aux.count('T'), aux.count('G')]
        
    #Ditribucion de los intrones
    AI = Seq(cadena).count('A')
    CI = Seq(cadena).count('C')
    TI = Seq(cadena).count('T')
    GI = Seq(cadena).count('G')

    totales = []

    for j in range(4):
        conta = 0
        for i in range(len(dist)):
            conta = conta + dist[i][j]
        totales.append(conta)
    
    I = [AI - totales[0], CI -totales[1], TI -totales[2], GI -totales[2]]
    dist.append(I)
    
    matrizTrans = list(range(len(dist))) #np.zeros((len(dist),4))
    for i in range(len(dist)):
        total = sum(dist[i])
        matrizTrans[i] = []
        for j in range(4):
            matrizTrans[i].append(dist[i][j]*(1/total)) 
    
    return(matrizTrans)

def ActCadOcu(genes): #Funcion qeuc crea la nueva cadena oculta de acuerdo a los nuevos genes estimados
    gen = list(range(len(genes)))
    cadOculta = ''

    for i in range(len(genes)):
        for j in range(len(genes[i])):
            laLista = re.split('ATG'+ genes[i][j], cadena)
            cadOculta = cadOculta  + len(laLista[0])*'I' + (len(genes[i][j]) + 3)*(str(gen[i])) 
            cade = laLista[1]
    cadOculta = cadOculta + 'I'*len(cade)

    return(cadOculta)

def ObtenerGenes(cadena):
    posible = re.split('ATG', cadena)
    posible2 = re.split('ATG', cadena)
    posible.pop(0)
    posible2.pop(0)
    for i in range(len(posible)):
        posible[i] = 'ATG' + posible[i]
    patrone = '((TAA)|(TGA)|(TAG))'

    #Separamos por los STOPS encontramos los posibles genes
    posgen = list(range(len(posible2)))
    posgen2 = list(range(len(posible2)))

    for i in range(len(posible)):
        if posible2[i] == re.split(patrone, posible2[i])[0]:
            posgen[i] = 'rep'
            posgen2[i] = 'rep'
        else:
            posgen[i] = 'ATG' + re.split(patrone, posible2[i])[0] + 'TAA'
            posgen2[i] =re.split(patrone, posible2[i])[0]

    cuenta = posgen.count('rep')
    for i in range(cuenta):
        posgen.remove('rep')
        posgen2.remove('rep')
    
    #Estos son nuestros posbles genes candidatos
    #posgen pero eliminiemos alas secuencias que no pueden ser frgmentos como los de longitud 1

    ayuda = []
    for i in range(len(posgen2)):
        if len(posgen2[i]) < 3:
            ayuda.append(i)

            #posgen.remove(i)

    cuen = 0
    for i in range(len(ayuda)):
        posgen2.pop(ayuda[i]-cuen)
        posgen.pop(ayuda[i]-cuen)
        cuen = cuen + 1
        
    #Localizamos el dinucleotido CpG para poder agrupar genes
    info = []
    for i in range(len(posgen2)):
        info.append(Seq(posgen2[i]).count('CG')/len(posgen2[i]))
        
    #Ya agrupamos los genes con sus fragmentos asociados
    genes = []
    aux = []
    for i in range(len(info)):
        if info[i] == 0:
            aux.append(posgen2[i])
        else:
            genes.append(aux)
            aux = [posgen2[i]]

    genes2 = []
    aux2 = []
    for i in range(len(info)):
        if info[i] == 0:
            aux2.append(posgen[i])
        else:
            genes2.append(aux2)
            aux = [posgen[i]]
            
    return(genes)
#Distancia Eclidania
def dist2(x, y):
    cuento = 0
    for i in range(len(x)):
        cuento = cuento + (x[i]-y[i])**2
    return math.sqrt(cuento)

#Si un exon tiene una distribucion similar a otra los asociaremos si queremos disminuir el número estimado de r

def varia(matrizTrans):
    var = list(range(len(matrizTrans)))
    
    for i in range(len(matrizTrans)):
        var[i] = []
        for j in range(len(matrizTrans)):
            var[i].append(dist2(matrizTrans[i], matrizTrans[j]))
        
    return(var)

def SimR(r):
    rnueva = ss.poisson.rvs(r)
    return(rnueva)

def rep(num, n):
    x = list(range(n))
    for i in range(n):
        x[i] = num
    return(x)
    
    
def verosimilitud(cadena, cadOculta, matrizTrans, oculta):
    gen = list(range(len(oculta)-1))
    
    for i in range(len(gen)):
        gen[i] = str(i)
    gen = gen + ['I']
    
    matObs =  pd.DataFrame(matrizTrans)
    matObs.columns = ['A', 'C', 'T', 'G']
    matObs.index = gen
    
    matOcu = pd.DataFrame(oculta)
    matOcu.columns = gen
    matOcu.index = gen
    
    #print(matObs)
    #print(matOcu)
    
    cuentaOcu = 2 + 1/len(gen) #Probabilidad de iniciar en cualquier estado oculto
    edoOcu = re.split('', cadOculta)
    edoObs = re.split('', cadena)
    cuentaObs = 2 + matObs[edoObs[1]][edoOcu[1]]
    
    for i in range(len(oculta)-1):
        cuentaOcu = np.log(matOcu[edoOcu[i+2]][edoOcu[i+1]]+2)*cuentaOcu
        cuentaObs = cuentaObs*np.log(matObs[edoObs[i+1]][edoOcu[i+1]]+2)
        #print([cuentaOcu, cuentaObs])
    
    total = np.log(cuentaOcu*cuentaObs)
    return(total)
    
#Optimizando
#Elegimos un parametro y lo volvemos a siomular
def SimPar(r, cadOcu, matrizTrans, oculta, genes):
    unif = np.random.rand() 
    var = varia(matrizTrans)
    var2 = varia(matrizTrans)
    
    gen = list(range(len(genes)))
    for i in range(len(gen)):
        gen[i] = str(i)
    gen = gen + ['I']
    
    if unif < 1/3:
        r2 = SimR(r) #Tambien tenemos que volver a simular los otras variables pues 
        if r2 < r:
            elegidos = list(set(list(np.random.randint(0, r, r-r2))))  #Borrrar a colapasar uno al azar
 
            var = varia(matrizTrans)
            var2 = varia(matrizTrans)

            for i in range(len(elegidos)):
                var2[elegidos[i]].remove(min(var2[elegidos[i]]))
                masCercana = min(var2[elegidos[i]])
                indice = var[elegidos[i]].index(masCercana)
                
                aux = list(range(len(genes)))
                if (indice in aux) and (elegidos[i] in aux):
                    genes[indice] = genes[indice] + list(genes[elegidos[i]])
                    genes.pop(elegidos[i])
                
            matrizTrans = ActMatObs(genes)
            oculta = ActMatOcu(genes)
            cadOcu = ActCadOcu(genes)
            r = len(genes)
            
        elif (r < r2 and r >r2): #if (r2 < (len(genesglob)-1)) and (r2>r): #Creamos otro renglon
            div = list(np.random.randint(0, r, r2-r)) 
            div = list(set(div))
            #print(div)
            for i in range(len(div)):
                #print(genes[div[i]])
                nuevoGen = int(np.random.randint(0, len(genes[div[i]]), 1))
                #print(genes[div[i]][nuevoGen])
                aux = [genes[div[i]][nuevoGen]]
                genes.append(aux)

            matrizTrans = ActMatObs(genes)
            oculta = ActMatOcu(genes)
            cadOcu = ActCadOcu(genes)
            r = len(genes) 
    
        else:
            matrizTrans = matrizTrans
        #print('simulamos r')

    elif (1/3 <= unif) and (unif < 2/3):
        #Sim Ocul
        nueva = []
        for i in range(len(oculta)):
            alphai = list(oculta[i])
            if 0 in alphai:
                alphai = tuple(rep(1/len(oculta[i]), len(oculta[i])))
                nueva.append(list(np.random.dirichlet(alphai, 1)[0]))
            else:
                alphai = tuple(oculta[i])
                nueva.append(list(np.random.dirichlet(alphai, 1)[0]))
        oculta = nueva
        #print('simulamos MatOcu')
    else:
        nueva = []
                               
        for i in range(len(matrizTrans)):
            alphai = list(matrizTrans[i])
            if 0 in alphai:
                alphai = tuple(rep(1/4, 4))
                nueva.append(list(np.random.dirichlet(alphai, 1)[0]))
            else:
                alphai = tuple(matrizTrans[i])
                nueva.append(list(np.random.dirichlet(alphai, 1)[0]))
                
        matrizTrans = nueva
        #print('simulamos MatObs')
        
    return([r, cadOcu, matrizTrans, oculta, genes])
    
#Muestreo de Gibbs
def MuestreoGibs(cadena, cadOculta, matrizTrans, oculta, genes, num):
    compara = verosimilitud(cadena, cadOculta, matrizTrans, oculta) #número de la verosimilitud :v
    graf = [compara]
    r = len(genes)
    numGen = [r]
    
    for i in range(num):
        theta = SimPar(r, cadOculta, matrizTrans, oculta, genes) #lista de esta forma: [r, cadOcu, matrizTrans, oculta, genes]
        compara2 = verosimilitud(cadena, theta[1], theta[2], theta[3]) 
        
        if compara > compara2: #Si no mejora le damos chance tal vez ponerna

            proba = compara2/compara
            proba = abs(proba) - abs(int(proba))
            #if proba > 1:
            #    proba = 1/proba
                
            desci = int(np.random.binomial(1, 1-proba, 1))
            
            if desci == 1:
                r = theta[0]
                cadOculta = theta[1]
                matrizTrans = theta[2]
                oculta = theta[3]
                genes = theta[4]
                compara = compara2
        else: #Si mejora nos quedamos con la nueva
            r = theta[0]
            cadOculta = theta[1]
            matrizTrans = theta[2]
            oculta = theta[3]
            genes = theta[4]
            compara = compara2
        graf.append(compara)
        numGen.append(len(genes))

        #print(compara)
    return([cadena, cadOculta, matrizTrans, oculta, graf, numGen])
    
genesglob = ObtenerGenes(cadena)
ocultaglob = ActMatOcu(genesglob)
matrizTransglob = ActMatObs(genesglob)
cadOcuglob = ActCadOcu(genesglob)
#pd.DataFrame(matrizTransglobal)
verosimilitud(cadena, cadOcuglob, matrizTransglob, ocultaglob)
MuestreoGibs(cadena, cadOcuglob, matrizTransglob, ocultaglob, genesglob, 80)