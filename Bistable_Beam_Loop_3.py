

import numpy as np
import math
from ansys.mapdl.core import launch_mapdl
import os



H1_Array = [5]

H2_Array = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]

L2_Array = [5, 50/9, 50/8, 50/7,  25/3, 10, 12.5, 50/3, 25]

t_Array = [0.5, 1, 1.5, 2]


mapdl = launch_mapdl(nproc=16, memory=128)

# JSON_CURLY_BEAM = np.zeros((np.size(H1_Array)*np.size(H2_Array)*np.size(L2_Array)*np.size(t_Array), 7))
# JSON_Counter = -1
for H1i in H1_Array:
    for L2i in L2_Array:
        for H2i in H2_Array:
            for ti in t_Array:

                mapdl.clear()

                Directory_Counter = 'H1_' + str(H1i) + '_H2_' + str(H2i) + '_L2_' + str(L2i) + '_t_' + str(ti)
                path_APDL = 'C:/Users/User_WS/.conda/envs/AnsysApdl/APDL/SET_5/' + str(Directory_Counter)

                if not os.path.exists(path_APDL):

                    os.mkdir(path_APDL)
                    mapdl.cwd(path_APDL)

                    os.chdir(path_APDL)
                    path_Python = os.getcwd()

                    mapdl.prep7()
                    mapdl.units('USER')  # User system (mm, MPa, N,kg, s, K).

                    Pi = np.pi

                    # Mooney Rivlin 3
                    # C10 = 0.0308307
                    # C01 = 0
                    # C11 = 0.0269727


                    # Mooney Rivlin 5
                    # C10 = -3196.35028854309
                    # C01 = 4242.29059644814
                    # C20 = 624.131775843925
                    # C11 = -2632.47362255165
                    # C02 = 4367.82258390202
                    # D1 = 0

                    Tol1 = 0.00001
                    Tol2 = ti/4
                    Tol3 = 2
                    Tol4 = 2.5

                    Delta = 2 * H1i

                    N1 = 5  # Number of divisions for first SIN height
                    N2 = 5  # Number of divisions for second SIN height
                    N3 = 5  # Number of divisions for second SIN length
                    N4 = 5  # Number of divisions for thickness
                    N5 = 101  # Number of Keypoints for geometry

                    L1 = 50
                    L2 = L2i
                    H1 = H1i
                    H2 = H2i
                    t = ti/4

                    Z = 0

                    mapdl.et(1, "PLANE183", kop3=3)
                    mapdl.r(1, 5)  # thickness of 0.001 meters)


                    # mapdl.tb("HYPER", 1, "", npts=5, tbopt="MOONEY")
                    # mapdl.tbdata(1, C10, C01, C20, C11, C02, D1)

                    # mapdl.tb("HYPER", 1, "", npts=3, tbopt="MOONEY")
                    # mapdl.tbdata(1, C10, C01, C11)

                    mapdl.mp('EX', 1, 2300)
                    mapdl.mp('NUXY', 1, 0.3)

                    X_array = np.linspace(0, L1, N5)
                    i = 0

                    for KP_X1 in X_array:
                        i = i + 1
                        X1_1 = KP_X1
                        X1_1 = KP_X1 + (H1 * (-t) * Pi * np.sin((2 * Pi * KP_X1) / L1)) / (L1 * ((H1 ** 2 * Pi ** 2 * np.sin((2 * KP_X1 * Pi) / L1) ** 2) / L1 ** 2 + 1) ** (1 / 2))
                        Y1_1 = - (H1 * (np.cos((2 * Pi * KP_X1) / L1) - 1)) / 2 - (-t) / ((H1 ** 2 * Pi ** 2 * np.sin((2 * KP_X1 * Pi) / L1) ** 2) / L1 ** 2 + 1) ** (1 / 2)

                        X2_1 = KP_X1 + (H2 * (-t) * Pi * np.sin((2 * Pi * KP_X1) / L2)) / (L2 * ((H2 ** 2 * Pi ** 2 * np.sin((2 * KP_X1 * Pi) / L2) ** 2) / L2 ** 2 + 1) ** (1 / 2))
                        Y2_1 = - (H2 * (np.cos((2 * Pi * KP_X1) / L2) - 1)) / 2 - (-t) / ((H2 ** 2 * Pi ** 2 * np.sin((2 * KP_X1 * Pi) / L2) ** 2) / L2 ** 2 + 1) ** (1 / 2)

                        X1 = (X1_1 + X2_1)/2
                        Y1 = Y1_1 + Y2_1

                        mapdl.k(i, X1, Y1, Z)

                    for KP_X2 in X_array:
                        i = i + 1
                        X1_2 = KP_X2 + (H1 * t * Pi * np.sin((2 * Pi * KP_X2) / L1)) / (L1 * ((H1 ** 2 * Pi ** 2 * np.sin((2 * KP_X2 * Pi) / L1) ** 2) / L1 ** 2 + 1) ** (1 / 2))
                        Y1_2 = - (H1 * (np.cos((2 * Pi * KP_X2) / L1) - 1)) / 2 - t / ((H1 ** 2 * Pi ** 2 * np.sin((2 * KP_X2 * Pi) / L1) ** 2) / L1 ** 2 + 1) ** (1 / 2)

                        X2_2 = KP_X2 + (H2 * t * Pi * np.sin((2 * Pi * KP_X2) / L2)) / (L2 * ((H2 ** 2 * Pi ** 2 * np.sin((2 * KP_X2 * Pi) / L2) ** 2) / L2 ** 2 + 1) ** (1 / 2))
                        Y2_2 = - (H2 * (np.cos((2 * Pi * KP_X2) / L2) - 1)) / 2 - t / ((H2 ** 2 * Pi ** 2 * np.sin((2 * KP_X2 * Pi) / L2) ** 2) / L2 ** 2 + 1) ** (1 / 2)

                        X2 = (X1_2 + X2_2)/2
                        Y2 = Y1_2 + Y2_2

                        mapdl.k(i, X2, Y2, Z)

                    mapdl.allsel()

                    SPLN = int((N5-6)/5+1)

                    for SPLi in range(1, SPLN+1):
                        mapdl.spline(5*(SPLi-1)+1, 5*(SPLi-1)+2, 5*(SPLi-1)+3, 5*(SPLi-1)+4, 5*(SPLi-1)+5, 5*(SPLi-1)+6)

                    for SPLi in range(1, SPLN+1):
                        mapdl.spline(N5+5*(SPLi-1)+1, N5+5*(SPLi-1)+2, N5+5*(SPLi-1)+3, N5+5*(SPLi-1)+4, N5+5*(SPLi-1)+5, N5+5*(SPLi-1)+6)

                    mapdl.allsel()

                    mapdl.l(1, N5+1)
                    mapdl.l(N5, 2*N5)

                    mapdl.allsel()
                    mapdl.lsel("S", "LINE", "", 1, N5-1,)
                    mapdl.lcomb("All")

                    mapdl.allsel()
                    mapdl.lsel("S", "LINE", "", N5, 2*(N5-1),)
                    mapdl.lcomb("All")

                    mapdl.allsel()

                    mapdl.numstr("line", 1)
                    mapdl.numcmp("line")

                    mapdl.allsel()
                    mapdl.al(1, 3, 2, 4)

                    mapdl.allsel()
                    #mapdl.smrtsize("OFF", 0.2, 1, 2, 7.5, 15, 1.4, "ON", "ON", 4, "OFF")
                    mapdl.desize(7, 7, 3000)
                    mapdl.amesh("All")

                    mapdl.allsel()

                    #mapdl.eplot()

                    mapdl.finish()

                    mapdl.run('/SOLU')

                    mapdl.antype(0)
                    mapdl.nlgeom("ON")
                    mapdl.nsubst(700, 4000, 500, "OFF")
                    mapdl.autots("OFF")
                    mapdl.outres("All", "All")

                    mapdl.allsel()
                    mapdl.nsel("S", "LOC", "X", -Tol1, Tol1,)
                    mapdl.d("All", "UX", 0)
                    mapdl.d("All", "UY", 0)

                    mapdl.allsel()

                    mapdl.nsel("S", "LOC", "X", L1-Tol1, L1+Tol1)
                    mapdl.d("All", "UX", 0)
                    mapdl.d("All", "UY", 0)

                    mapdl.allsel()

                    mapdl.nsel("S", "LOC", "X", L1/2-Tol3, L1/2+Tol3)

                    if L2==50/3 or L2==50/5 or L2==50/7 or L2==50/9:

                        mapdl.nsel("R", "LOC", "Y", H1 + H2 + t - Tol2, H1 + H2 + t + Tol2)

                    else:
                        mapdl.nsel("R", "LOC", "Y", H1 + t - Tol2, H1 + t + Tol2)


                    mapdl.d("All", "UY", -Delta)

                    mapdl.allsel()
                    output = mapdl.solve()

                    NUMSUBSTEP = int(mapdl.get('NUMSUBSTEP', 'ACTIVE', 0, 'SET', 'NSET'))

                    mapdl.finish()

                    result = mapdl.result

                    #result.plot_nodal_solution((1, 200), 'y', label='Displacement')

                    mapdl.finish()

                    mapdl.save()

                    mapdl.post1()

                    mapdl.allsel()

                    mapdl.nsel("S", "LOC", "X", 5, L1 - 5)
                    mapdl.nsel("U", "LOC", "X", L1 / 2 - Tol4, L1 / 2 + Tol4)

                    NUMSELECTEDNODES = int(mapdl.get('NUMSELECTEDNODES', 'node', '', 'count'))

                    MINSELECTEDNODES = int(mapdl.get('MINSELECTEDNODES', 'node', '', 'num', 'min'))

                    SELECTEDNODES = np.zeros(NUMSELECTEDNODES, dtype=int)
                    SELECTEDNODES[0] = MINSELECTEDNODES

                    SELECTEDNODES_X = np.zeros(NUMSELECTEDNODES)
                    SELECTEDNODES_Y = np.zeros(NUMSELECTEDNODES)

                    for SELECTEDNODEi in range(1, NUMSELECTEDNODES):
                        SELECTEDNODES[SELECTEDNODEi] = int(
                            mapdl.get('Next', 'node', SELECTEDNODES[SELECTEDNODEi - 1], 'NXTH'))

                    counter_node = 0
                    for SELECTEDNODEi in SELECTEDNODES:
                        SELECTEDNODES_X[counter_node] = (mapdl.get('node_x', 'node', SELECTEDNODEi, 'LOC', 'X'))
                        SELECTEDNODES_Y[counter_node] = (mapdl.get('node_y', 'node', SELECTEDNODEi, 'LOC', 'Y'))
                        counter_node += 1

                    SELECTEDNODES_DATA = np.zeros((NUMSELECTEDNODES, 3))

                    SELECTEDNODES_DATA[:, 0] = SELECTEDNODES
                    SELECTEDNODES_DATA[:, 1] = SELECTEDNODES_X
                    SELECTEDNODES_DATA[:, 2] = SELECTEDNODES_Y

                    mapdl.allsel()
                    mapdl.lsel("S", "LINE", "", 3)
                    mapdl.nsll("S", 1)
                    NUMSIDENODES = int(mapdl.get('NUMSIDENODES','node', '', 'count'))

                    MINSIDENODES = int(mapdl.get('MINSIDENODES','node', '', 'num', 'min'))

                    SIDENODES = np.zeros(NUMSIDENODES, dtype=int)
                    SIDENODES[0] = MINSIDENODES
                    for SIDENODEi in range(1, NUMSIDENODES):
                        SIDENODES[SIDENODEi] = int(mapdl.get('Next', 'node', SIDENODES[SIDENODEi-1], 'NXTH'))

                    mapdl.allsel()
                    mapdl.post1()

                    Von_Mises_strain = np.zeros((NUMSELECTEDNODES, NUMSUBSTEP))
                    Side_Nodal_Reaction_Force = np.zeros((NUMSIDENODES, NUMSUBSTEP))

                    for SUBSTEPi in range(1, NUMSUBSTEP):

                        Tot_Node_Numbers, Tot_Nodal_Elastic_Strain = result.nodal_elastic_strain((1, SUBSTEPi))
                        reaction_forces = result.nodal_reaction_forces((1, SUBSTEPi))
                        reaction_forces_matrix = np.zeros((3, np.size(reaction_forces[0])))
                        reaction_forces_matrix[0, :] = reaction_forces[0]
                        reaction_forces_matrix[1, :] = reaction_forces[1]
                        reaction_forces_matrix[2, :] = reaction_forces[2]

                        Counter_i = 0
                        for NODEi in SELECTEDNODES:

                            if not math.isnan(Tot_Nodal_Elastic_Strain[NODEi, 6]):

                                Von_Mises_strain[Counter_i, SUBSTEPi-1] = Tot_Nodal_Elastic_Strain[NODEi, 6]
                                Counter_i = Counter_i + 1

                        Counter_j = 0
                        for NODEj in SIDENODES:

                            for NODEk in range(1, np.size(reaction_forces[0])):

                                if reaction_forces_matrix[1, NODEk] == NODEj and reaction_forces_matrix[2, NODEk] == 2:
                                    Side_Nodal_Reaction_Force[Counter_j, SUBSTEPi - 1] = reaction_forces_matrix[0, NODEk]
                                    Counter_j = Counter_j + 1
                                    break

                    Max_Von_Mises_strain_Array = np.amax(Von_Mises_strain, axis=0)
                    Max_Von_Mises_strain = np.amax(Von_Mises_strain)
                    Sum_Side_Nodal_Reaction_Force_Array = 2*np.sum(Side_Nodal_Reaction_Force, axis=0)
                    Max_Side_Nodal_Reaction_Force = np.amax(Sum_Side_Nodal_Reaction_Force_Array)
                    Min_Side_Nodal_Reaction_Force = np.amin(Sum_Side_Nodal_Reaction_Force_Array)

                    np.savetxt("Von_Mises_strain.csv", Von_Mises_strain, delimiter=",")
                    np.savetxt("Side_Nodal_Reaction_Force.csv", Sum_Side_Nodal_Reaction_Force_Array, delimiter=",")
                    np.savetxt("Selectednodes_coordinate.csv", SELECTEDNODES_DATA, delimiter=",")

                    mapdl.finish()

mapdl.exit()


