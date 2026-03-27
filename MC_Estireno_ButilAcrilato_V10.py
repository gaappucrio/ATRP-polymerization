'''
Programa de distribuição de Monte Carlo para copolimerização de estireno n-butil acrilato
'''
#################################################################################################
#               Importação de bibliotecas                                                       #
#################################################################################################

import time
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO
import lib_MC

#################################################################################################
#               Constantes                                                                      #
#################################################################################################

MW1,MW2,MW_sol = 104.14,128.17,138.16
rho1,rho2,rho_sol = 0.906,0.894,1.053
Rcte, Temp, N = 1.987,273.15+110, 6.022e23
V = 5.e-19
start = time.time() 
#################################################################################################
#               Definição das constantes cinéticas                                              #
#################################################################################################

#Mathematical Modeling of Atom-Transfer Radical Copolymerization
#Dynamic Monte Carlo Simulation of ATRP in a Batch Reactor
#General Method for Determination of the Activation, Deactivation, and Initiation Rate Constants in Transition Metal-Catalyzed Atom Transfer Radical Processes

#Iniciação
kai,kdi = 26e-3,29e6

#Equilíbrio e propagação
ka1,ka2 = 0.45,0.055
kp11, kp22 = 1.06e7*np.exp(-7067/(Rcte*Temp)), 7.37e5*np.exp(-2299/(Rcte*Temp))
r1,r2 = 0.79,0.26
kp12,kp21 = kp11/r1,kp22/r2
kd1,kd2 = 1.15e7,8.e7

#Transferência para monômero
ktr11,ktr22,ktr12 = kp11*2.198e-1*np.exp(-2820/Temp),kp22*(1.3e-4),2.997e4*np.exp(-7835.8/(Rcte*Temp))
ktr21 = ktr12

#Terminação por combinação
ktc12,ktc11 = 7.681e9*np.exp(-2690.42/(Rcte*Temp)),kp11**2*1.1e-5*np.exp(12452.2/(Rcte*Temp))
ktc21 = ktc12
ktc22 = kp22/(2.5e-4)

#################################################################################################
#               Definição das constantes cinéticas de Monte Carlo                               #
#################################################################################################

#Iniciação
kaiMC,kdiMC = kai/(N*V),kdi/(N*V)

#Equilíbrio e propagação
ka1MC,ka2MC,kp11MC,kp22MC,kp12MC,kp21MC,kd1MC,kd2MC = ka1/(N*V),ka2/(N*V),kp11/(N*V),kp22/(N*V),kp12/(N*V),kp21/(N*V),kd1/(N*V),kd2/(N*V)

#Transferência para monômero
ktr11MC,ktr22MC,ktr12MC,ktr21MC = ktr11/(N*V),ktr22/(N*V),ktr12/(N*V),ktr21/(N*V)

#Terminação por combinação
ktc11MC,ktc12MC,ktc21MC,ktc22MC = 2*ktc11/(N*V),ktc12/(N*V),ktc21/(N*V),2*ktc22/(N*V)

#################################################################################################
#               Condições iniciais                                                              #
#################################################################################################

#Atom Transfer Radical Copolymerization of Styrene and n-Butyl Acrylate
f = 0.510
n1,n2 = 0.04*f,0.04*(1-f)
Vol = 1e-3*(n1*MW1/rho1+n2*MW2/rho2+1e-3*MW_sol/rho_sol)
n = 11
X = np.zeros(n,dtype=int)
X[0] = 0.4e-3/Vol*N*V                           #Catalisador
X[1] = 0.e0                                     #Catalisador ativado
X[2] = n1/Vol*N*V                               #Monômero A
X[3] = n2/Vol*N*V                               #Monômero B
X[4] = 0.e0                                     #Polímero vivo A
X[5] = 0.e0                                     #Polímero vivo B
X[6] = 0.e0                                     #Polímero dormente A
X[7] = 0.e0                                     #Polímero dormente B
X[8] = 0.e0                                     #Polímero morto
X[9] = 0.e0                                     #Polímero dormente iniciado
X[10] = 0.4e-3/Vol*N*V                          #Iniciador
n_max = int(0)
PA,PB,L,DA,DB = np.zeros(n_max,dtype=int),np.zeros(n_max,dtype=int),np.zeros(1000,dtype=int),np.zeros(n_max,dtype=int),np.zeros(n_max,dtype=int)

x = np.zeros(len(L),dtype=int)
for i in range(1,len(L)):
    x[i] = i


ti,tf = 0.0,4.e4
lista_t,lista_p,lista_Mn,lista_Mw,lista_PDI = [],[],[],[],[]
M0 = X[2]+X[3]
#################################################################################################
#               Taxas das reações                                                               #
#################################################################################################

n_reacoes = 20
Reacao = np.zeros(n_reacoes,dtype=float)

#Equilíbrio e propagação
Reacao[0] = ka1MC*X[6]*X[0]
Reacao[1] = ka2MC*X[7]*X[0]
Reacao[2] = kp11MC*X[4]*X[2]
Reacao[3] = kp22MC*X[5]*X[3]
Reacao[4] = kp12MC*X[4]*X[3]
Reacao[5] = kp21MC*X[5]*X[2]
Reacao[6] = kd1MC*X[4]*X[1]
Reacao[7] = kd2MC*X[5]*X[1]

#Transferência para monômeros
Reacao[8] = ktr11MC*X[4]*X[2]
Reacao[9] = ktr12MC*X[4]*X[3]
Reacao[10] = ktr22MC*X[5]*X[3]
Reacao[11] = ktr21MC*X[5]*X[2]

#Terminação por combinação
Reacao[12] = ktc11MC*X[4]*(X[4]-1)/4
Reacao[13] = ktc12MC*X[4]*X[5]
Reacao[14] = ktc22MC*X[5]*(X[5]-1)/4
Reacao[15] = ktc21MC*X[5]*X[4]

#Iniciação
Reacao[16] = kaiMC*X[0]*X[10]
Reacao[17] = kdiMC*X[9]*X[1]
Reacao[18] = kp11MC*X[9]*X[2]
Reacao[19] = kp22MC*X[9]*X[3]

#################################################################################################
#               Simulação de MonteCarlo                                                         #
#################################################################################################
t = ti
conv = 1-(X[2]+X[3])/(2*0.02/Vol*N*V)
lista_Mn.append(0)
lista_Mw.append(0)
lista_PDI.append(float("NaN"))
lista_t.append(t)
lista_p.append(conv)
while t<tf:
    Reacao_soma = np.sum(Reacao)
    r = np.random.exponential()
    dt = r/Reacao_soma
    r = np.random.rand()
    R = r*Reacao_soma
    if R <= Reacao[0]:                          #Reação 1
        i = int(np.random.rand()*X[6])
        PA = np.append(PA, DA[i])
        DA = np.delete(DA,i)
        X[0]-=1
        X[6]-=1
        X[1]+=1
        X[4]+=1
    elif R <= sum(Reacao[0:2]):                 #Reação 2
        i = int(np.random.rand()*X[7])
        PB = np.append(PB, DB[i])
        DB = np.delete(DB,i)
        X[0]-=1
        X[7]-=1
        X[1]+=1
        X[5]+=1
    elif R <= sum(Reacao[0:3]):                 #Reação 3      
        i = int(np.random.rand()*X[4])
        PA[i]+=1
        X[2]-=1  
    elif R <= sum(Reacao[0:4]):                 #Reação 4
        i = int(np.random.rand()*X[5])
        X[3]-=1
        PB[i]+=1
    elif R <= sum(Reacao[0:5]):                 #Reação 5
        i = int(np.random.rand()*X[4])
        PB = np.append(PB, PA[i]+1)
        PA = np.delete(PA,i)
        X[3]-=1
        X[4]-=1
        X[5]+=1
    elif R <= sum(Reacao[0:6]):                 #Reação 6
        i = int(np.random.rand()*X[5])
        PA = np.append(PA, PB[i]+1)
        PB = np.delete(PB,i)
        X[2]-=1
        X[5]-=1
        X[4]+=1
    elif R <= sum(Reacao[0:7]):                 #Reação 7
        i = int(np.random.rand()*X[4])
        DA = np.append(DA, PA[i])
        PA = np.delete(PA,i)
        X[0]+=1
        X[6]+=1
        X[1]-=1
        X[4]-=1
    elif R <= sum(Reacao[0:8]):                 #Reação 8
        i = int(np.random.rand()*X[5])
        DB = np.append(DB, PB[i])
        PB = np.delete(PB,i)
        X[0]+=1
        X[1]-=1
        X[5]-=1
        X[7]+=1
    elif R <= sum(Reacao[0:9]):                 #Reação 9
        i = int(np.random.rand()*X[4])
        L[PA[i]]+=1
        PA[i] = 1
        X[8]+=1
        X[2]-=1
        Mn,Mw = lib_MC.Mn(DA,DB,x,L)
        lista_Mn.append(Mn)
        lista_Mw.append(Mw)
        lista_PDI.append(Mw/Mn)
        conv = 1-(X[2]+X[3])/(2*0.02/Vol*N*V)
        lista_t.append(t)
        lista_p.append(conv)
    elif R <= sum(Reacao[0:10]):                #Reação 10
        i = int(np.random.rand()*X[4])
        L[PA[i]]+=1
        PB = np.append(PB, 1)
        PA = np.delete(PA,i)
        X[4]-=1
        X[5]+=1
        X[8]+=1
        X[3]-=1
        Mn,Mw = lib_MC.Mn(DA,DB,x,L)
        lista_Mn.append(Mn)
        lista_Mw.append(Mw)
        lista_PDI.append(Mw/Mn)
        conv = 1-(X[2]+X[3])/(2*0.02/Vol*N*V)
        lista_t.append(t)
        lista_p.append(conv)
    elif R <= sum(Reacao[0:11]):                #Reação 11
        i = int(np.random.rand()*X[5])
        L[PB[i]]+=1
        PB[i] = 1
        X[8]+=1
        X[3]-=1
        Mn,Mw = lib_MC.Mn(DA,DB,x,L)
        lista_Mn.append(Mn)
        lista_Mw.append(Mw)
        lista_PDI.append(Mw/Mn)
        conv = 1-(X[2]+X[3])/(2*0.02/Vol*N*V)
        lista_t.append(t)
        lista_p.append(conv)        
    elif R <= sum(Reacao[0:12]):                #Reação 12
        i = int(np.random.rand()*X[5])
        L[PB[i]]+=1
        PA = np.append(PA,1)
        PB = np.delete(PB,i)
        X[8]+=1
        X[2]-=1
        X[5]-=1
        X[4]+=1
        Mn,Mw = lib_MC.Mn(DA,DB,x,L)
        lista_Mn.append(Mn)
        lista_Mw.append(Mw)
        lista_PDI.append(Mw/Mn)
        conv = 1-(X[2]+X[3])/(2*0.02/Vol*N*V)
        lista_t.append(t)
        lista_p.append(conv)
    elif R <= sum(Reacao[0:13]):                #Reação 13
        i,j = int(np.random.rand()*X[4]),int(np.random.rand()*X[4])
        if i == j:
            if i ==0:
                j+=1
            else:
                j-=1
        L[PA[i]+PA[j]]+=1
        PA = np.delete(PA,[i,j])
        X[8]+=1
        X[4]-=2
    elif R <= sum(Reacao[0:14]):                #Reação 14
        i,j = int(np.random.rand()*X[4]),int(np.random.rand()*X[5])
        L[PA[i]+PB[j]]+=1
        PA = np.delete(PA,i)
        PB = np.delete(PB,j)
        X[8]+=1
        X[4]-=1
        X[5]-=1
    elif R <= sum(Reacao[0:15]):                #Reação 15
        i,j = int(np.random.rand()*X[5]),int(np.random.rand()*X[5]) 
        if i == j:
            if i ==0:
                j+=1
            else:
                j-=1
        L[PB[i]+PB[j]]+=1
        PB = np.delete(PB,[i,j])
        X[8]+=1
        X[5]-=2
    elif R <= sum(Reacao[0:16]):                #Reação 16
        i,j = int(np.random.rand()*X[5]),int(np.random.rand()*X[4])
        L[PB[i]+PA[j]]+=1
        PB = np.delete(PB,i) 
        PA = np.delete(PA,j)
        X[8]+=1
        X[4]-=1
        X[5]-=1
    elif R<= sum(Reacao[0:17]):                 #Reação 17                 
        X[10]-=1
        X[0]-=1
        X[1]+=1
        X[9]+=1
    elif R<= sum(Reacao[0:18]):                 #Reação 18
        X[10]+=1
        X[0]+=1
        X[1]-=1
        X[9]-=1  
    elif R<= sum(Reacao[0:19]):                 #Reação 19
        PA = np.append(PA,1)
        X[9]-=1
        X[2]-=1
        X[4]+=1
    else:                                       #Reação 20
        PB = np.append(PB,1)
        X[9]-=1
        X[3]-=1
        X[5]+=1
    
    #Equilíbrio e propagação
    Reacao[0] = ka1MC*X[6]*X[0]
    Reacao[1] = ka2MC*X[7]*X[0]
    Reacao[2] = kp11MC*X[4]*X[2]
    Reacao[3] = kp22MC*X[5]*X[3]
    Reacao[4] = kp12MC*X[4]*X[3]
    Reacao[5] = kp21MC*X[5]*X[2]
    Reacao[6] = kd1MC*X[4]*X[1]
    Reacao[7] = kd2MC*X[5]*X[1]

    #Transferência para monômeros
    Reacao[8] = ktr11MC*X[4]*X[2]
    Reacao[9] = ktr12MC*X[4]*X[3]
    Reacao[10] = ktr22MC*X[5]*X[3]
    Reacao[11] = ktr21MC*X[5]*X[2]

    #Terminação por combinação
    Reacao[12] = ktc11MC*X[4]*(X[4]-1)/4
    Reacao[13] = ktc12MC*X[4]*X[5]
    Reacao[14] = ktc22MC*X[5]*(X[5]-1)/4
    Reacao[15] = ktc21MC*X[5]*X[4]

    #Iniciação
    Reacao[16] = kaiMC*X[0]*X[10]
    Reacao[17] = kdiMC*X[9]*X[1]
    Reacao[18] = kp11MC*X[9]*X[2]
    Reacao[19] = kp22MC*X[9]*X[3]
    
    t+=dt

conv = 1-(X[2]+X[3])/(2*0.02/Vol*N*V)
lista_t.append(t)
lista_p.append(conv)
Mn,Mw = lib_MC.Mn(DA,DB,x,L)
lista_Mn.append(Mn)
lista_Mw.append(Mw)
lista_PDI.append(Mw/Mn)
tempo_total = time.time()-start
print("Tempo total de simulação: %f" %(time.time()-start))

#################################################################################################
#               Tratamento dos dados                                                            #
#################################################################################################
L_soma = np.zeros(501,dtype=int)
M_soma = np.zeros(len(L_soma),dtype=float)
W_soma = np.zeros(len(L_soma),dtype=float)

#Distribuição de tamanho 200
x_soma = np.zeros(len(L_soma),dtype=int)
for i in range(len(L_soma)):
    x_soma[i] = i
    L_soma[i] = L[i]
    M_soma[i] = i*(MW1+MW2)/2

for i in range(len(L_soma)):
    W_soma[i] = L_soma[i]*M_soma[i]/sum(L_soma*M_soma)

#################################################################################################
#               Arquivos                                                                        #
#################################################################################################

file1 = open('distribuicao.txt','w')
for i in range(len(L_soma)):
    file1.write("%i,%i,%.2f,%.2f\n" %(x_soma[i],L_soma[i],M_soma[i],W_soma[i]))
file1.close

file2 = open('log.txt','w')
file2.write("Tempo (s): %.6f\nVolume de controle (L): %.2e" %(tempo_total,V))
file2.close

file3 = open('conversao.txt','w')
for i in range(len(lista_t)):
    file3.write('%.6f,%.6f,%.2f,%.2f,%.3f\n' %(lista_t[i],lista_p[i],lista_Mn[i],lista_Mw[i],lista_PDI[i]))
file3.close

file4 = open('distribuicao_total.txt','w')
for i in range(len(L)):
    file4.write('%i;%i\n' %(x[i],L[i]))
file4.close