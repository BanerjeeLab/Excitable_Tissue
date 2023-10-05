import numpy as np
from shutil import copy
from shapely.geometry import Polygon

class Vertex:
    def __init__(self,n, position, velocity, cells,adjacents, dt, axis, pair_4fold):
        self.id = n
        self.r = position
        self.v = velocity
        self.cells = cells
        self.adj = adjacents
        self.dt = dt
        self.eje = axis
        self.par = pair_4fold

    @property
    def x(self):
        return self.r[0]

    @property
    def y(self):
        return self.r[1]

    @property
    def vx(self):
        return self.v[0]

    @property
    def vy(self):
        return self.v[1]

class Cell:
    def __init__(self, c_id, id_vertices, A):
        self.n = c_id
        self.ind = id_vertices
        self.area = A
        self.center = (0,0,0)
        self.pol = Polygon([ ])

    def cal_area(self):
        self.area = self.pol.area

    def cal_centro(self):
        xx, yy = self.pol.exterior.coords.xy
        xx = xx.tolist()[:-1]
        yy = yy.tolist()[:-1]
        self.center = (np.mean(xx),np.mean(yy),0)

class Tissue:
    def __init__(self, vertices, cells, KA, area0, friction, L_x, L_y, dt, lT1):
        self.vs = vertices
        self.cs = cells
        self.KA = KA
        self.area0_c = area0
        self.mu = friction
        self.Lx = L_x
        self.Ly = L_y
        self.dt = dt
        self.lT1 = lT1

    def vm_forces(self, v):
        velocity = np.array([0, 0 , 0])

        for j in range(len(self.vs[v].cells)):
            c = self.vs[v].cells[j]
            vi = int(self.vs[v].adj[j][0])
            vd = int(self.vs[v].adj[j][1])

            l_ii = self.modulo(self.resta_vectores(self.vs[vi].r,self.vs[v].r))
            l_id = self.modulo(self.resta_vectores(self.vs[vd].r,self.vs[v].r))
            T_ii = (T00 + Tact[v][vi]* l_ii) * 0.5
            T_id = (T00 + Tact[v][vd]* l_id) * 0.5

            a = self.resta_vectores(self.vs[vi].r,self.vs[v].r)
            b = self.resta_vectores(self.vs[vd].r,self.vs[v].r)

            velocity = velocity + (T_ii/l_ii)*a + (T_id/l_id)*b

            v_out_of_plane = np.array([0.0, 0.0, 1.0])
            rest = self.resta_vectores(self.vs[vd].r, self.vs[vi].r)
            v_tot = np.cross(rest, v_out_of_plane)
            pressure = 0.5* self.KA * (self.cs[c].area - self.area0_c)

            velocity =  velocity +  pressure * v_tot

        return (1/self.mu) * velocity

    def cal_polygon(self, v):
        id_vs = self.cs[v].ind
        vectors_PBC = []
        vfixed = self.vs[id_vs[0]].r
        for i in range(len(id_vs)):
            vi = self.vs[id_vs[i]].r
            vectors_PBC.append(self.resta_vectores(vi,vfixed) + vfixed)
        self.cs[v].pol = Polygon(vectors_PBC)

    def cal_polygons(self):
        for i in range(len(self.cs)):
            self.cal_polygon(i)

    def evol_Tact(self, now):
        for i in range(n_vertices):
            for j in range(i):
                if j in self.vs[i].adj:
                    L = self.modulo(self.resta_vectores(self.vs[i].r,self.vs[j].r))
                    L0 = LargosNat[i][j]
                    strain = L/L0 - 1
                    # Case for activating junction
                    if (strain > strainOn and not junctionActive[i][j] and
                            now > timeDeactivated[i][j] + cooldownTime and JunctionState[i][j]==3):
                        Tact[i][j] = Gamma_act
                        Tact[j][i] = Gamma_act
                        junctionActive[i][j] = True
                        junctionActive[j][i] = True
                        timeActivated[i][j] = now
                        timeActivated[j][i] = now
                        JunctionState[i][j] =1
                        JunctionState[j][i] =1
                    # Case for deactivating junction: From Active to Refractory
                    elif (strain < -strainOff or now - timeActivated[i][j] > turnoverTime
                            and junctionActive[i][j]):
                        Tact[i][j] = 0
                        Tact[j][i] = 0
                        junctionActive[i][j] = False
                        junctionActive[j][i] = False
                        timeDeactivated[i][j] = now
                        timeDeactivated[j][i] = now
                        JunctionState[i][j] =2
                        JunctionState[j][i] =2

                    # From Refractory to Passive
                    elif (not junctionActive[i][j] and now > timeDeactivated[i][j] + cooldownTime
                            and JunctionState[i][j]==2):
                        JunctionState[i][j] =3
                        JunctionState[j][i] =3

    def T1_tissue(self, now):
        flats = []
        for i in range(n_vertices):
            if len(self.vs[i].adj) == 0:
                flat = []
            else:
                flat = [self.vs[i].adj[0][0],self.vs[i].adj[0][1],self.vs[i].adj[1][0]]
            flats.append(flat)

        for i in range(n_vertices):
            for j in flats[i]:
                cells_i = self.vs[i].cells
                cells_j = self.vs[j].cells
                c_int_i = 0
                c_int_j = 0
                for ci in cells_i:
                    if ci not in cells_borde:
                        c_int_i = c_int_i + 1
                for cj in cells_j:
                    if cj not in cells_borde:
                        c_int_j= c_int_j + 1

                if j>i and j in np.array(self.vs[i].adj).flatten() and c_int_i==3 and c_int_j==3:
                    lij = self.modulo(self.resta_vectores(self.vs[i].r, self.vs[j].r))
                    if lij < l_T1:
                        # Write T1 to file
                        T1s.write(f"{i} {j} {now}\n")
                        print(f"{round(Kc/Kl,2)} {i} {j} {round(now,2)}")

                        # Turn off active tension
                        if junctionActive[i][j]:
                            print("off")
                            Tact[i][j] = 0
                            Tact[j][i] = 0
                            junctionActive[i][j] = False
                            junctionActive[j][i] = False
                            JunctionState[i][j] = 3
                            JunctionState[j][i] = 3

                        self.delete_v1_original(i,j)
                        nT12 = T00
                        nl120 = 1.5*l_T1
                        deleted_info = [[],[],[],[]]
                        indices_k = []
                        for k in range(len(i_act)):
                            if i_act[k] ==i and j_act[k] ==j:
                                deleted_info[0].append(i_act[k])
                                deleted_info[1].append(j_act[k])
                                deleted_info[2].append(ti_act[k])
                                deleted_info[3].append(delta_act[k])
                                indices_k.append(k)
                        for v in range(len(indices_k)):
                                i_act.pop(indices_k[0])
                                j_act.pop(indices_k[0])
                                ti_act.pop(indices_k[0])
                                delta_act.pop(indices_k[0])
                        self.evol_Tact(now)
                        self.create_v1_perpendicular(i,nT12,nl120)

    def evol_vertex(self,now):
        self.pos_nuevas_vertex()
        self.cal_polygons()
        self.cal_areas()
        self.T1_tissue(now)
        self.new_restlegth()
        self.cal_polygons()
        self.cal_areas()

    def pos_nuevas_vertex(self):
        for n in range(len(self.vs)):
            if len(self.vs[n].cells)==0:
                self.vs[n].v = np.array([0,0,0])
            else:
                self.vs[n].v  = self.vm_forces(n)

        for n in range(len(self.vs)):
            self.vs[n].r = self.vs[n].r + self.dt * self.vs[n].v
            self.vs[n].r = self.resta_inicial(self.vs[n].r)

    def data(self):
        pos = [] #vertex position
        A = [] #cell areas
        C = [] #cell center

        celda0 = []
        celda1 = []

        new01 = []
        for i in range(len(self.vs)):
            pos.append(self.vs[i].r)
            new01.append(i)

        new02 = []
        for i in range(len(self.cs)):
            A.append(self.cs[i].area)
            C.append(self.cs[i].center)
            new02.append(i)
            for k in range(len(T1.cs[i].ind)-1):
                celda0.append(i)
                celda1.append(T1.cs[i].ind[k])

        np.savetxt('./results/data_vertices.txt',  np.c_[new01, pos], fmt='%1.5f')
        np.savetxt('./results/data_cells.txt',  np.c_[new02,A,C], fmt='%1.10f')
        np.savetxt('./results/data_celda.txt',  np.c_[celda0, celda1], fmt='%1.0f')

    def cal_areas(self):
        for i in range(len(self.cs)):
            self.cs[i].cal_area()

    def new_restlegth(self):
        for i in range(len(self.vs)):
            for j in [T1.vs[i].adj[0][0],T1.vs[i].adj[0][1],T1.vs[i].adj[1][0]]:
                if j>i:
                    lij = self.modulo(self.resta_vectores(T1.vs[i].r, T1.vs[j].r))
                    e_ij = (lij - LargosNat[i][j])/LargosNat[i][j]

                    v_Lij = np.round(-Kl*(LargosNat[i][j]-lij),6)
                    new_Lij = LargosNat[i][j] + self.dt * v_Lij
                    LargosNat[i][j] = new_Lij
                    LargosNat[j][i] = new_Lij

    def resta_inicial(self,pos):
        x,y = pos[0], pos[1]
        if pos[0] < 0.0:
            x = self.Lx + pos[0]
        if pos[0] > self.Lx:
            x = -self.Lx + pos[0]
        if pos[1] < 0.0:
            y = self.Ly + pos[1]
        if pos[1] > self.Ly:
            y = -self.Ly + pos[1]
        return np.array([x,y,0])

    def resta_vectores(self,a,b):
        if abs(a[0]-b[0]) < self.Lx/2.:
            q = a[0]-b[0]
        if a[0]-b[0] >= self.Lx/2.:
            q = a[0] - self.Lx -b[0]
        if a[0]-b[0] <= -self.Lx/2.:
            q = a[0] + self.Lx -b[0]
        if abs(a[1]-b[1]) < self.Ly/2.:
            n = a[1]-b[1]
        if a[1]-b[1] >= self.Ly/2.:
            n = a[1] - self.Ly -b[1]
        if a[1]-b[1] <= -self.Ly/2.:
            n = a[1] + self.Ly -b[1]
        return np.array([q, n, 0])

    def modulo(self,x):
        return np.sqrt(x.dot(x))

    def create_v1_perpendicular(self,vi1,nT12,nl120):
        cells_i1 = self.vs[vi1].cells
        cellA = cells_i1[0]
        cellB = cells_i1[1]
        cellC = cells_i1[2]
        cellD = cells_i1[3]

        vert_A = self.cs[cellA].ind
        vert_B = self.cs[cellB].ind
        vert_C = self.cs[cellC].ind
        vert_D = self.cs[cellD].ind


        for i in vert_A:
            if i in vert_B and i not in vert_C and i not in vert_D:
                vert2 = i
        for i in vert_A:
            if i in vert_D and i not in vert_B and i not in vert_C:
                vert3 = i
        for i in vert_C:
            if i in vert_B and i not in vert_A and i not in vert_D:
                vert6 = i
        for i in vert_C:
            if i in vert_D and i not in vert_A and i not in vert_B:
                vert5 = i
        vert4 = vi1
        vi2 = self.vs[vi1].par
        vert1 = self.vs[vi1].par

        #MODIFY adj
        for i in range(len(self.vs[vert3].adj)):
            for j in range(len(self.vs[vert3].adj[i])):
                if self.vs[vert3].adj[i][j] == vert4:
                    self.vs[vert3].adj[i][j] = vert1
        for i in range(len(self.vs[vert5].adj)):
            for j in range(len(self.vs[vert5].adj[i])):
                if self.vs[vert5].adj[i][j] == vert4:
                    self.vs[vert5].adj[i][j] = vert1

        self.vs[vert4].cells = [cellA, cellB, cellC]
        self.vs[vert4].adj = [[vert2,vert1],[vert6,vert2],[vert1,vert6]]
        self.vs[vert1].cells = [cellC, cellD, cellA]
        self.vs[vert1].adj = [[vert5,vert4],[vert3,vert5],[vert4,vert3]]

        #MODIFY REST LENGTHS
        LargosNat[vert3][vert1]= LargosNat[vert3][vert4]
        LargosNat[vert1][vert3]= LargosNat[vert3][vert4]
        LargosNat[vert5][vert1]= LargosNat[vert5][vert4]
        LargosNat[vert1][vert5]= LargosNat[vert5][vert4]
        LargosNat[vert3][vert4]= 0
        LargosNat[vert4][vert3]= 0
        LargosNat[vert5][vert4]= 0
        LargosNat[vert4][vert5]= 0

        #MODIFY IND cell A
        add_A = []
        q = 0
        for i in range(len(self.cs[cellA].ind)-1):
            if self.cs[cellA].ind[i] == vert4 and q==0:
                add_A.append(vert4)
                add_A.append(vert1)
                q=1
            else:
                add_A.append(self.cs[cellA].ind[i])
        self.cs[cellA].ind = add_A
        if self.cs[cellA].ind[0] != self.cs[cellA].ind[len(T1.cs[cellA].ind)-1]:
            self.cs[cellA].ind.append(self.cs[cellA].ind[0])

        #MODIFY IND cell C
        add_C = []
        q = 0
        for i in range(len(self.cs[cellC].ind)-1):
            if self.cs[cellC].ind[i] == vert5 and q==0:
                add_C.append(vert5)
                add_C.append(vert1)
                q=1
            else:
                add_C.append(T1.cs[cellC].ind[i])
        self.cs[cellC].ind = add_C
        if self.cs[cellC].ind[0] != self.cs[cellC].ind[len(T1.cs[cellC].ind)-1]:
            self.cs[cellC].ind.append(self.cs[cellC].ind[0])

        #MODIFY IND cell D
        add_D = []
        q = 0
        for i in range(len(self.cs[cellD].ind)-1):
            if self.cs[cellD].ind[i] == vert4 and q==0:
                add_D.append(vert1)
                q=1
            else:
                add_D.append(self.cs[cellD].ind[i])
        self.cs[cellD].ind = add_D
        if self.cs[cellD].ind[0] != self.cs[cellD].ind[len(self.cs[cellD].ind)-1]:
            self.cs[cellD].ind.append(self.cs[cellD].ind[0])

        c1y2 = []
        for c in [cellA,cellB,cellC,cellD]:
            if c not in self.vs[vert4].eje:
                c1y2.append(c)
        c1 = c1y2[0]
        c2 = c1y2[1]

        self.vs[vert4].eje = [-1,-1]
        self.vs[vert4].par = -1

        r_T1 = self.resta_vectores(self.cs[c1].center,self.cs[c2].center)
        rnorm_T1 =(1/self.modulo(r_T1)) * r_T1

        pos_vi1p = self.vs[vert4].r + (3/4)*self.lT1 * rnorm_T1
        pos_vi2p = self.vs[vert4].r - (3/4)*self.lT1 * rnorm_T1
        self.vs[vi1].r = pos_vi1p
        self.vs[vi2].r = pos_vi2p

        for c in [cellA,cellB,cellC,cellD]:
            self.cal_polygon(c)
            self.cs[c].calculate_area()
            self.cs[c].calculate_centro()

    def delete_v1_original(self,vi1,vi2):
        not_repeated = []
        cells_i1 = self.vs[vi1].cells #A,B,D    o    B,D,A    o    D,A,B
        cells_i2 = self.vs[vi2].cells #B,C,D    o    C,D,B    o    D,B,C
        for i in range(len(cells_i1)):
            if cells_i1[i] in cells_i2:
                qq=0
            else:
                not_repeated.append(cells_i1[i]) #A
        for i in range(len(cells_i2)):
            if cells_i2[i] in cells_i1:
                qq=0
            else:
                not_repeated.append(cells_i2[i]) #C
        if cells_i1[0]== not_repeated[0]: #A,B,D
            cells_i1_new = [cells_i1[0],cells_i1[1],not_repeated[1]]
        elif cells_i1[1]== not_repeated[0]: #D,A,B
            cells_i1_new = [cells_i1[1],cells_i1[2],not_repeated[1]]
        else:  #B,D,A
            cells_i1_new = [cells_i1[2],cells_i1[0],not_repeated[1]]
        if cells_i2[0]== not_repeated[1]: #C,D,B
            cells_i2_new = [cells_i2[0],cells_i2[1],not_repeated[0]]
        elif cells_i2[1]== not_repeated[1]: #B,C,D
            cells_i2_new = [cells_i2[1],cells_i2[2],not_repeated[0]]
        else: #D,B,C
            cells_i2_new = [cells_i2[2],cells_i2[0],not_repeated[0]]

        cellA = cells_i1_new[0]
        cellB = cells_i1_new[1]
        cellC = cells_i1_new[2]
        cellD = cells_i2_new[1]
        vert_A = self.cs[cellA].ind[:len(self.cs[cellA].ind)-1]
        vert_B = self.cs[cellB].ind[:len(self.cs[cellB].ind)-1]
        vert_C = self.cs[cellC].ind[:len(self.cs[cellC].ind)-1]
        vert_D = self.cs[cellD].ind[:len(self.cs[cellD].ind)-1]

        for i in vert_A:
            if i in vert_B and i not in vert_C and i not in vert_D:
                vert2 = i
        for i in vert_A:
            if i in vert_D and i not in vert_B and i not in vert_C:
                vert3 = i
        for i in vert_C:
            if i in vert_B and i not in vert_A and i not in vert_D:
                vert6 = i
        for i in vert_C:
            if i in vert_D and i not in vert_A and i not in vert_B:
                vert5 = i
        vert4 = vi1
        vert1 = vi2

        #MODIFY adj
        for i in range(len(self.vs[vert6].adj)):
            for j in range(len(self.vs[vert6].adj[i])):
                if self.vs[vert6].adj[i][j] == vert1:
                    self.vs[vert6].adj[i][j] = vert4

        for i in range(len(self.vs[vert5].adj)):
            for j in range(len(self.vs[vert5].adj[i])):
                if self.vs[vert5].adj[i][j] == vert1:
                    self.vs[vert5].adj[i][j] = vert4
        self.vs[vert4].cells = [cellA, cellB, cellC, cellD]
        self.vs[vert4].adj = [[vert2,vert3],[vert6,vert2],[vert5,vert6],[vert3,vert5]]
        self.vs[vert1].cells = []
        self.vs[vert1].adj = []

        #MODIFY TENSIONES AND NATURAL LENGTHS
        Tensiones[vert6][vert4]= Tensiones[vert6][vert1]
        Tensiones[vert4][vert6]= Tensiones[vert6][vert1]
        Tensiones[vert5][vert4]= Tensiones[vert5][vert1]
        Tensiones[vert4][vert5]= Tensiones[vert5][vert1]
        Tensiones[vert6][vert1] = 0
        Tensiones[vert1][vert6] = 0
        Tensiones[vert5][vert1] = 0
        Tensiones[vert1][vert5] = 0

        LargosNat[vert6][vert4]= LargosNat[vert6][vert1]
        LargosNat[vert4][vert6]= LargosNat[vert6][vert1]
        LargosNat[vert5][vert4]= LargosNat[vert5][vert1]
        LargosNat[vert4][vert5]= LargosNat[vert5][vert1]
        LargosNat[vert6][vert1] = 0
        LargosNat[vert1][vert6] = 0
        LargosNat[vert5][vert1] = 0
        LargosNat[vert1][vert5] = 0

        #MODIFY IND cell B
        for i in T1.cs[cellB].ind:
            if i == vi2:
                self.cs[cellB].ind.remove(vi2)
        if self.cs[cellB].ind[0] != self.cs[cellB].ind[len(self.cs[cellB].ind)-1]:
            self.cs[cellB].ind.append(self.cs[cellB].ind[0])


        #MODIFY IND cell C
        add_C = []
        q = 0
        for i in range(len(self.cs[cellC].ind)-1):
            if self.cs[cellC].ind[i] == vi2 and q==0:
                add_C.append(vi1)
                q=1
            else:
                add_C.append(self.cs[cellC].ind[i])
        self.cs[cellC].ind = add_C
        if self.cs[cellC].ind[0] != self.cs[cellC].ind[len(self.cs[cellC].ind)-1]:
            self.cs[cellC].ind.append(self.cs[cellC].ind[0])

        #MODIFY IND cell D
        for i in self.cs[cellD].ind:
            if i == vi2:
                self.cs[cellD].ind.remove(vi2)
        if self.cs[cellD].ind[0] != self.cs[cellD].ind[len(T1.cs[cellD].ind)-1]:
            self.cs[cellD].ind.append(self.cs[cellD].ind[0])

        self.vs[vert4].eje = [cellA,cellC]
        self.vs[vert4].par = vert1

        mid = self.vs[vert4].r + 0.5 * self.resta_vectores(self.vs[vert1].r,self.vs[vert4].r)
        mid = self.resta_inicial(mid)

        self.vs[vi1].r = mid
        self.vs[vi2].r = mid

        for c in [cellA,cellB,cellC,cellD]:
            self.cal_polygon(c)
            self.cs[c].calculate_area()
            self.cs[c].calculate_centro()
