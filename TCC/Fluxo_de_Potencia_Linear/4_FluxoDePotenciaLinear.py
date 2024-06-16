import requests
from amplpy import AMPL, tools
import pandas as pd
#from requests import PandaRequest
import numpy as np
import matplotlib.pyplot as plt

tools.modules.load()
ampl = AMPL()
#print(ampl.get_option('version'))

ampl.read('C:\\TCC_1\\virt\\4_FluxoDePotenciaNaoLinear_PQ.mod')
ampl.readData('C:\\TCC_1\\virt\\4_Sistema_69_Barras_PQ.dat')

ampl.setOption('solver', 'cplex')
ampl.solve()

Ob = ampl.getSet("Ob").getValues().toPandas().index.astype(int)

Vnom = ampl.getParameter("Vnom").getValues().toPandas().iloc[0,0]
Vsqr = ampl.getVariable("Vsqr").getValues().toPandas()
Vpu = {}

print("\n|-------------------------------------|")
print("|  Solucao do Fluxo de potencia para  |")
print("| Sistemas de distribuicao de energia |")
print("|-------------------------------------|\n")

for i in Ob:
    Vpu[i] = Vsqr.loc[i, 'Vsqr.val'] / Vnom**2

print("|-------------------------------------|")
print("|              Tensões                |")
print("|-----------|-----------|-------------|")
print("|   Barra   |  V [kV]   |  V [p.u.]   |")
print("|-----------|-----------|-------------|")
for i in Ob:
    print("|%10d |%10.4f |%12.4f |" % (i, Vsqr.loc[i, 'Vsqr.val'], Vpu[i]))
print("|-----------|-----------|-------------|\n")
#print("PROXIMOOOO\n")


Ol1, Ol2 = zip(*ampl.getSet("Ol").getValues().toPandas().index)

Ot = ampl.getSet("Ot").getValues().toPandas().index.astype(int)


Pperdas = {}
Qperdas = {}
R = ampl.getParameter("R").getValues().toPandas()
X = ampl.getParameter("X").getValues().toPandas()
#Q = ampl.getParameter("Q").getValues().toPandas()

Isqr = ampl.getVariable("Isqr").getValues().toPandas()
P = ampl.getVariable("P").getValues().toPandas()


for (i,j,t) in zip(Ol1, Ol2, Ot):
    Pperdas[(i,j)] = R.at[(i,j), 'R'] * Isqr.at[(i,j,t),'Isqr.val']
    Qperdas[(i,j)] = X.at[(i,j), 'X'] * Isqr.at[(i,j,t), 'Isqr.val']


print ('Reultado Linhas\n')
 
print('  i    j    I[amp]    P[kw]    P[kwar]    Pperdas[kw]    Qperdas[kwar]  \n')
  
for (i,j) in zip(Ol1, Ol2):
 print(i,j,Isqr.at[(i,j), 'Isqr.val'],P.at[(i,j), 'P.val'],P.at[(i,j), 'P.val']*1000,Pperdas[i,j],Qperdas[i,j])

 
 
 
print('\n\n')
eixo = list(range(0,len(Ob)))

plt.plot(eixo,Vsqr)
plt.xlabel('Barras')
plt.ylabel('Tensão Kv')
plt.grid()
plt.show()

 
ampl.close