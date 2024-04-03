#Creation of vertices
pos0 = []
vel0 = []
cells0 = []
t = []

for i in range(n_vertices):
    pos0.append(np.array([x[i],y[i],0]))
    vel0.append(np.array([0,0,0]))
    cells_i_horario = cells_tipo1[i]
    ady_i_horario = adjs_tipo1[i]
    t.append(Vertex(i, pos0[i],vel0[i],np.array(cells_i_horario), np.array(ady_i_horario), dt, [-1,-1], -1))


#Creations of cells
celulas = []
for i in range(n_cells):
    b = []
    for j in range(len(celda)):
        if celda[j] == i:
            b.append(int(v_celda[j]))
    celulas.append(Cell(i,b, 0))

    b.append(b[0])


#Creation of tissue
T1 = Tissue(t, celulas, KA, A0, mu, L_x, L_y, dt, l_T1)

LargosNat = np.zeros((n_vertices, n_vertices))
Tact = np.zeros((n_vertices, n_vertices))

for i in range(len(lista_i)):
    v1 = int(lista_i[i])
    v2 = int(lista_j[i])
    LargosNat[v1][v2] = lij0_save[i]
    LargosNat[v2][v1] = lij0_save[i]

# Initialize state of each junction
JunctionState = 3*np.ones((n_vertices, n_vertices))
# Initialize TactTemp to hold active tension during calculation
TactTemp = np.zeros((n_vertices, n_vertices))
# Rule-based active tension propagation
junctionActive = np.zeros((n_vertices, n_vertices)).astype(bool)
# Records time of junction activation, for active tension propagation model
timeActivated = np.zeros((n_vertices, n_vertices))
# Records time of most recent deactivation, initialize so that everything can be
# activated immediately.
timeDeactivated = np.full((n_vertices, n_vertices), -cooldownTime * 1.1)

for k in range(len(i_act)):
    i = int(i_act[k])
    j = int(j_act[k])
    Tact[i][j] = Gamma_act
    Tact[j][i] = Gamma_act
    JunctionState[i][j] =1
    JunctionState[j][i] =1
    junctionActive[i][j] = True
    junctionActive[j][i] = True
    timeActivated[i][j] = 0.0
    timeActivated[j][i] = 0.0

cells_borde = []
