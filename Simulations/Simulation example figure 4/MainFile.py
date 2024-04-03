import numpy as np
import matplotlib.pyplot as plt
from shutil import copy
import random
import math
import pandas as pd
import sys

#Classes: vertices-cells-tissue
exec(open("classes_and_functions.py").read())

# Load initial configuration data
path = "./Initial-condition-fig4/"
celdas_data = pd.read_csv(str(path)+'/cell_vertex_id.csv')
celda = list(celdas_data["cell id"])
v_celda = list(celdas_data["vertex id"])

junctions_data = pd.read_csv(str(path)+'/junctions_data.csv')
lista_i = list(junctions_data["vertex i"])
lista_j = list(junctions_data["vertex j"])
lij_save = list(junctions_data["lij"])
lij0_save = list(junctions_data["l0ij"])
tensiones_distribution = list(junctions_data["Tij"])

vertex_data = pd.read_csv(str(path)+'/vertex_data.csv')
x = list(vertex_data["x"])
y = list(vertex_data["y"])
z = np.zeros(len(x))

top = pd.read_csv(str(path)+'/adjacents_vertices_cells.csv')
v_top = list(top["vertex i"])
adj1 = list(top["adj 1"])
adj2 = list(top["adj 2"])
adj3 = list(top["adj 3"])
cell1 = list(top["cell 1"])
cell2 = list(top["cell 2"])
cell3 = list(top["cell 3"])
adjs_tipo1 = list(zip(zip(adj1, adj2), zip(adj3, adj1), zip(adj2, adj3)))
cells_tipo1 =list(zip(cell1, cell2, cell3))
n_cells = int(celda[len(celda)-1]+1)
n_vertices = int(len(x))



# Basic parameters
A0=1
mu  = 0.636
KA = 1
Kl = 0.7
l_T1 = 0.01
T00 = 0.1
e_critic = 0.05
Gamma_act = 0.5

# Strain thresholds for active tension toggling
strainOff = e_critic * 3 * 1e10
strainOn = 0.1

# Turnover time for active tension to toggle off
# Cooldown time is minimum interval between turning off and on again
turnoverTime = 1.6
cooldownTime = turnoverTime

# box size
factor = np.sqrt(2/(3*np.sqrt(3)))
L_x = 13 * np.sqrt(3)*factor
L_y = 16 * np.sqrt(3)*np.sqrt(3)/2*factor
L_z = 0

# Set step size, steps per save, and total saves
saveResolutionFactor = 1
dt = 0.01
d_pasos = 25 // saveResolutionFactor
total_txt = 4 * saveResolutionFactor * 15


#activity information pulses
data_act = np.loadtxt('data_active_pulse.txt', ndmin=2)
i_act = data_act[:,0].tolist()
j_act = data_act[:,1].tolist()
ti_act = data_act[:,2].tolist()
delta_act = data_act[:,3].tolist()

#Creation of the tissue
exec(open("tissue_creation.py").read())


#Main program:
T1.cal_polygons()
T1.cal_areas()
T1.data()
copy('./results/data_vertices.txt', './results/0_vertices.txt')
copy('./results/data_cells.txt', './results/0_cells.txt')
copy('./results/data_celda.txt', './results/0_celda.txt')

lista_i = []
lista_j = []
lij_save = []
lij0_save = []
Tij_save = []
Tijact_save = []
flats = []
t_act_ij = []
state_ij=[]

for i in range(n_vertices):
    flat = [T1.vs[i].adj[0][0],T1.vs[i].adj[0][1],T1.vs[i].adj[1][0]]
    flats.append(flat)
for i in range(n_vertices):
    for j in flats[i]:
        if j>i and j in np.array(T1.vs[i].adj).flatten():
            lista_i.append(i)
            lista_j.append(j)
            lij = T1.modulo(T1.resta_vectores(T1.vs[i].r,T1.vs[j].r))
            lij_save.append(lij)
            lij0_save.append(LargosNat[i][j])
            Tijact_save.append(Tact[i][j])
            if junctionActive[i][j]==True:
                t_act_ij.append(0-timeActivated[i][j])
            else:
                t_act_ij.append(-1)
            state_ij.append(JunctionState[i][j])

np.savetxt('./lengths_tensions/0.txt', np.c_[lista_i, lista_j, lij_save, lij0_save, Tijact_save,t_act_ij,state_ij], fmt='%1.5f')

for i in range(int(d_pasos * total_txt)-1):
    N = (i+1)*dt
    T1.evol_vertex((i+1)*dt)
    T1.evol_Tact((i+1)*dt)
    c = (i + 1) % d_pasos
    d = int((i + 1) / d_pasos)
    if c == 0 :
        print((i+1)*dt)

        T1.data()
        copy('./results/data_vertices.txt', './results/'+ str(d) + '_vertices.txt')
        copy('./results/data_cells.txt', './results/'+str(d) + '_cells.txt')
        copy('./results/data_celda.txt', './results/'+str(d) + '_celda.txt')

        n_i3 = []
        lista_i = []
        lista_j = []
        lij_save = []
        lij0_save = []
        Tij_save = []
        Tijact_save = []
        flats3 = []
        t_act_ij = []
        state_ij=[]
        for i in range(n_vertices):
            n_i3.append(i)
            flat = [T1.vs[i].adj[0][0],T1.vs[i].adj[0][1],T1.vs[i].adj[1][0]]
            flats3.append(flat)
        for ii in range(len(n_i3)):
            i = n_i3[ii]
            for j in flats3[ii]:
                if j>i:
                    lista_i.append(i)
                    lista_j.append(j)
                    lij = T1.modulo(T1.resta_vectores(T1.vs[i].r,T1.vs[j].r))
                    lij_save.append(lij)
                    lij0_save.append(LargosNat[i][j])
                    Tijact_save.append(Tact[i][j])
                    if junctionActive[i][j]==True:
                        t_act_ij.append(N-timeActivated[i][j])

                    else:
                        t_act_ij.append(-1)
                    state_ij.append(JunctionState[i][j])

        np.savetxt('./lengths_tensions/'+str(d)+'.txt', np.c_[lista_i, lista_j, lij_save, lij0_save, Tijact_save, t_act_ij,state_ij], fmt='%1.5f')
