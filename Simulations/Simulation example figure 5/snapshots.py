import numpy as np
import matplotlib.pyplot as plt
from shutil import copy
import math

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})

class Vertex:
    def __init__(self, position):
        self.r = position

    @property
    def x(self):
        return self.r[0]

    @property
    def y(self):
        return self.r[1]

class Tissue:
    def __init__(self, vertices,vertex_i,vertex_j, Stateij):
        self.vs = vertices
        self.vi = vertex_i
        self.vj = vertex_j
        self.state_ij = Stateij

    def draw_tissue(self):
        for i in range(len(self.state_ij)):
            fac = self.state_ij[i]
            colores = ["maroon","dodgerblue","#C1CDCD"]
            if int(fac)==1:
                color_ij = colores[0]
                z = 40
                lw=2.5
            elif int(fac)==2:
                color_ij = colores[1]
                z = 50
                lw=2.5
            elif int(fac)==3:
                color_ij = colores[2]
                z = 30
                lw=2.5
            v1=int(self.vi[i])
            v2=int(self.vj[i])

            x1,y1=self.vs[v1].x,self.vs[v1].y
            x2,y2=self.vs[v2].x,self.vs[v2].y

            if np.abs(x2-x1)>L_x/2 or np.abs(y2-y1)>L_y/2:
                if np.abs(y2-y1)<L_y/2:
                    if x2-x1>L_x/2:
                        x2left = x2-L_x
                        x1right = x1+L_x
                        plt.plot([x1right,x2],
                                 [y1,y2], '-', linewidth=lw, color=color_ij,zorder=z)
                        plt.plot([x1,x2left],
                                 [y1,y2], '-',linewidth=lw, color=color_ij,zorder=z)
                    else:
                        x2left = x2+L_x
                        x1right = x1-L_x
                        plt.plot([x1right,x2],
                                 [y1,y2], '-', linewidth=lw,color=color_ij,zorder=z)
                        plt.plot([x1,x2left],
                                 [y1,y2], '-', linewidth=lw,color=color_ij,zorder=z)

                if np.abs(x2-x1)<L_x/2:
                    if y1-y2>L_y/2:
                        y2left = y2+L_y
                        y1right = y1-L_y
                        plt.plot([x1,x2],
                                 [y1right,y2], '-', linewidth=lw,color=color_ij,zorder=z)
                        plt.plot([x1,x2],
                                 [y1,y2left], '-', linewidth=lw,color=color_ij,zorder=z)
                    #else y2-y1>L_y/2:
                    else:
                        y2left = y2-L_y
                        y1right = y1+L_y
                        plt.plot([x1,x2],
                                 [y1right,y2], '-',linewidth=lw, color=color_ij,zorder=z)
                        plt.plot([x1,x2],
                                 [y1,y2left], '-',linewidth=lw, color=color_ij,zorder=z)

                else:
                    if x2-x1>L_x/2 and y2-y1>L_y/2:
                        x2left = x2-L_x
                        x1right = x1+L_x
                        y2left = y2-L_y
                        y1right = y1+L_y
                        plt.plot([x1right,x2],
                                 [y1right,y2], '-', linewidth=lw,color=color_ij,zorder=z)
                        plt.plot([x1,x2left],
                                 [y1,y2left], '-', linewidth=lw,color=color_ij,zorder=z)

                    elif x2-x1>L_x/2 and y1-y2>L_y/2:
                        x2left = x2-L_x
                        x1right = x1+L_x
                        y2left = y2+L_y
                        y1right = y1-L_y
                        plt.plot([x1right,x2],
                                 [y1right,y2], '-',linewidth=lw, color=color_ij,zorder=z)
                        plt.plot([x1,x2left],
                                 [y1,y2left], '-',linewidth=lw, color=color_ij,zorder=z)


            else:
                plt.plot([x1,x2],
                         [y1,y2], '-',linewidth=lw,color=color_ij,zorder=z)
                plt.plot([x1,x2],
                         [y1,y2], '-',linewidth=lw, color=color_ij,zorder=z)

def crear_T(texto0,texto1):
    data = np.loadtxt(texto0)
    topology = np.loadtxt(texto1)
    x = data[:,1]
    y = data[:,2]
    vi = topology[:,0]
    vj = topology[:,1]
    Stateij = topology[:,6]

    #CREACION VERTICES, CELULAR Y TEJIDO.
    pos0 = []
    t = []

    for i in range(len(x)):
        pos0.append((x[i],y[i],0))
        t.append(Vertex(pos0[i]))

    T1 = Tissue(t,vi,vj,Stateij)
    T1.draw_tissue()


plt.close('all')


factor = np.sqrt(2/(3*np.sqrt(3)))
L_x = 46 * np.sqrt(3)*factor
L_y = 60 * np.sqrt(3)*np.sqrt(3)/2*factor
plt.figure(figsize=(2,2))
dt = 0.05
saveResolutionFactor = 5
d_pasos = 25 // saveResolutionFactor

for ii in [99]:##
    fig, ax = plt.subplots(nrows=1, ncols=1)
    i = int(ii)
    crear_T('./results/'+str(i)+'_vertices.txt','./lengths_tensions/'+str(i)+'.txt')
    plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,
    left=False,
    labelleft=False,    # ticks along the top edge are off
    labelbottom=False)
    plt.axis('scaled')
    fig.set_size_inches(5, 5)
    plt.xlim([0,L_x])
    plt.ylim([0,L_y])
    plt.tight_layout(pad=1.05)
    size = fig.get_size_inches()*400
    plt.savefig('./saved_snapshots/'+str(i)+'.png', dpi=400)#,bbox_inches='tight')

    plt.close('all')
