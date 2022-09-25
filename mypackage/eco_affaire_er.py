
import numpy as np
from cmath import log10
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# npax=300


class eco_affaire_er:
    def __init__(self, npax, raction, coefremp, mach, allongV, srefe):
        # Initiation
        self.sref = srefe
        self.npax = npax
        self.mach = mach
        self.P = 30088
        self.T = 228.56
        self.Cz = 0.5
        self.Cs = 1.75*(10**(-5))
        self.r_action = raction  # Nm
        self.Nmot = 2
        self.altitude = 30000  # ft
        self.Coef_rempli = coefremp  # % Coef de remplissage avion
        self.prix_billet_eco = 0.25  # $/NM/pax eco
        self.prix_billet_aff = 0.5  # $/NM/pax aff
        self.prix_billet_er = 0.75  # $/NM/pax aff

        # naca4415
        self.m = 0.04  # Max camber
        self.p = 0.4  # Max camber position
        self.t = 0.15  # Thickness
        self.c = 1.0

        # Voilure:
        self.lamda = allongV
        self.phi25 = 30
        self.epsilon1 = 0.38

        # Empenage H:
        self.VH = 1
        self.phi25H = self.phi25+3
        self.lamdaH = 5.5  # -------a entrer
        self.epsilonH = 0.3
        self.e_sur_l_H = 10/100

        # Empenage V

        self.Vv = 0.09
        self.phi25V = self.phi25+10
        self.lamdaV = 1.3  # -------a entrer
        self.epsilonV = 0.5
        self.e_sur_l_V = 11/100

    def initialisation(self):
        MTOW = 1000  # 160000
        MZFW = 100  # 120000

        return MTOW, MZFW

    def Sref(self):
        if 100 <= self.npax < 246:
            Sref = 122.6
        else:
            Sref = 246  # 361.6
        return Sref

    def npf(self):

        if 100 <= self.npax <= 119:
            np_front = [5]

        if 119 < self.npax <= 150:
            np_front = [5, 6]

        if 150 < self.npax <= 210:
            np_front = [6]

        if 210 < self.npax <= 239:
            np_front = [6, 7]

        if 240 <= self.npax <= 320:
            np_front = [7, 8]

        return np_front

    def l_cylindre(self, Npf, nb_range_co, nb_rang_affaire, nb_rang_er, rang_AV, rang_AR):

        A = [42, 72, 36]
        B = [32, 72, 36]
        C = [30, 48, 20]
        I = [24, 48, 20]
        III = [20, 36, 0]
        Lpax_inch = []
        Lcyl_m = []

        if 80 < self.npax <= 109:
            L = [I[2]+2*III[2]]

        if 110 <= self.npax <= 139:
            L = [2*I[2]+III[2]]

        if 140 <= self.npax <= 179:
            L = [2*I[2]+2*III[2]]

        if 180 <= self.npax <= 299:
            if Npf <= 6:
                L = [(self.npax//75)*B[2], (self.npax//55)*C[2],
                     (self.npax//45)*I[2], (self.npax//35)*III[2]]

            else:
                L = [(self.npax//110)*A[2], (self.npax//55)*C[2],
                     (self.npax//45)*I[2], (self.npax//35)*III[2]]

        if self.npax > 299:
            L = [(self.npax//110)*A[2], (self.npax//45)
                 * I[2], (self.npax//55)*C[2]]

        for i in L:
            Lpax_inch = Lpax_inch+[nb_range_co*34 +
                                   nb_rang_affaire*36+nb_rang_er*40+i]

        for i in Lpax_inch:
            Lcyl_m = Lcyl_m+[(i-(rang_AR*34+rang_AV*40))*0.0254]

        return Lcyl_m, Lpax_inch

    def rayon(self, Npf):

        Marge = 1
        L_couloir = 10
        L_simple = 22
        L_double = 42
        L_triple = 62
        L_quadruple = 82
        L_quintuple = 102
        rayon_intern = []
        rayon_externe = []
        L_av = []
        L_ar = []
        L_ar_av = []
        i = 0

        if Npf == 5:
            rayon_intern = [((2*Marge+L_triple+L_couloir+L_double)*0.0254)/2]
            rayon_externe = [rayon_intern[0]*1.1]

        if Npf == 6:
            rayon_intern = [((2*Marge+2*L_triple+L_couloir)*0.0254)/2]
            rayon_externe = [rayon_intern[0]*1.1]

        if Npf == 7:
            rayon_intern = [
                ((2*Marge+2*L_double+2*L_couloir+L_triple)*0.0254)/2]
            rayon_externe = [rayon_intern[0]*1.1]

            if Npf == 8:
                rayon_intern = [((2*Marge+2*L_double+2*L_couloir+L_quadruple)*0.0254)/2,
                                ((2*Marge+2*L_triple+2*L_couloir+L_quadruple)*0.0254)/2]

        for i in rayon_intern:
            rayon_externe = rayon_externe+[i*1.1]

        for i in rayon_externe:

            L_ar = L_ar+[2*i*3.6]
            L_av = L_av+[2*i*1.7]

        for (item1, item2) in zip(L_ar, L_av):
            L_ar_av.append(item1 + item2)

        return L_ar_av, rayon_externe, L_av, L_ar

    def L_sur_D(self, l_cylindre, l_ar_av, rayon_extern):
        l_fuselage = l_cylindre+l_ar_av
        LsurD = l_fuselage/(2*rayon_extern)
        return LsurD, l_fuselage

    def closest(self,lst, K):

        return lst[min(range(len(lst)), key=lambda i: abs(lst[i]-K))]

    def Coordonee(self, s_ref, r_extern, l_cyl_m1, lav, lar, L_fuselage1):

        def symetrie(listes):
            liste_sym = []
            for i in listes:
                liste_sym = liste_sym+[-i]
            return liste_sym

        B = math.sqrt(self.lamda*s_ref)
        Y4 = B/2
        Y2 = r_extern
        Y3 = 0.4*Y4
        BF = 2*Y2

        L1 = (s_ref-(Y3-Y2)*(Y3+Y2)*math.tan(self.phi25*math.pi/180)) / \
            ((1+self.epsilon1)/2*(B-BF)+BF-3 *
             (1-self.epsilon1)*(Y3-Y2)*(Y3+Y2)/(2*(B-BF)))
        L2 = L1+(Y3-Y2)*(math.tan(self.phi25*math.pi/180) -
                         3/2*(1-self.epsilon1)/(B-BF)*L1)
        L4 = self.epsilon1*L1
        L3 = L4+(L1-L4)*(Y4-Y3)/(Y4-Y2)

        # Point 1-->5  et 1'--->2' et X3-X4
        X_1_prime = [-L2/4]

        Y_2_prime = [Y4]
        X_2_prime = [-(Y_2_prime[0]-Y2) *
                     math.tan(self.phi25*math.pi/180)+X_1_prime[0]]

        X = [0, X_2_prime[0]+L4/4, (X_2_prime[0]+L4/4)-L4, -L2, -L2, 0]
        Y = [Y2, Y4, Y4, Y3, Y2, Y2]

        X_sym = symetrie(X)
        Y_sym = symetrie(Y)
        x_a = [0, (X_2_prime[0]+L4/4)-L4]
        y_a = [0, 0]
        Y_1_prime = [Y[0]]

        X3 = [X[3]+L3]
        X4 = [X[1]]

        # Fleche BA et effilement

        Fleche_BA = -math.atan(X[1]/(Y[1]-Y[0]))*180/math.pi
        effilement = L4/L2

        # Calcule CMA

        L0 = (3*Y2*math.pow(L2, 2)+(Y3-Y2)*(math.pow(L2, 2)+math.pow(L3, 2) +
              L2*L3)+(Y4-Y3)*(math.pow(L3, 2)+math.pow(L4, 2)+L4*L3))*2/(3*s_ref)
        X0 = (X3[0]*((Y3-Y2)*(2*L3+L2)+(Y4-Y3)*(2*L3+L4)) +
              X4[0]*(Y4-Y3)*(2*L4+L3))/(3*s_ref)
        Y0 = (3*L2*math.pow(Y2, 2)+(Y3-Y2)*(L3*(Y2+2*Y3)+L2*(Y3+2*Y2)) +
              (Y4-Y3)*(L4*(Y3+2*Y4)+L3*(Y4+2*Y3)))/(3*s_ref)

        X_debut = X0
        Y_debut = Y0
        X_fin = X_debut-L0
        Y_fin = Y_debut

        X_CMA = [X_debut, X_fin]
        Y_CMA = [Y_debut, Y_fin]

        # epaisseur relative

        Aero = (0.89-(self.mach+0.02) *
                math.sqrt(math.cos(self.phi25*math.pi/180)))/L0
        Emplanture = 1.24*Aero
        Cassure = 0.94*Aero
        Extremite = 0.86*Aero

        # Calcule du volume Carburant

        A = math.pow(s_ref, 1.5)*math.pow(self.lamda, -0.4) * \
            (0.6*Emplanture+0.4*Extremite)
        MFW = 224*A+1570

        # Foyer

        X_Foyer = [X_debut-0.25*L0]
        Y_Foyer = [Y_fin]

        # Fuselage:

        X_pt_av = X_Foyer[0]+l_cyl_m1/2+lav
        Y_pt_av = 0
        X_pt_ar = X_Foyer[0]-l_cyl_m1/2-lar
        Y_pt_ar = 0
        X_debut_cyl = X_pt_av-lav
        Y_debut_cyl = r_extern
        X_fin_cyl = X_pt_ar+lar
        Y_fin_cyl = r_extern

        X_fuselage = [X_debut_cyl, X_fin_cyl]
        Y_fuselage = [Y_debut_cyl, Y_fin_cyl]

        X_fuselage_sym = symetrie(X_fuselage)
        Y_fuselage_sym = symetrie(Y_fuselage)

        # Empenage H:

        Lfa = X_pt_av-X0+0.25*L0

        LpH = L_fuselage1*0.91-Lfa

        SH = s_ref*L0*self.VH/LpH

        Y3H = math.sqrt(self.lamdaH*SH)/2
        L1H = SH/Y3H/(1+self.epsilonH)
        L3H = L1H*self.epsilonH
        X3H = L1H/4+Y3H*math.tan(self.phi25H*math.pi/180)-L3H/4

        LOH = (2*Y3H/(3*SH))*(L1H**2+L3H**2+L1H*L3H)
        XOH = X3H*Y3H*(2*L3H+L1H)/(3*SH)
        YOH = (Y3H**2)*(2*L3H+L1H)/(3*SH)

        XOH_prime = X_Foyer[0]-LpH

        # CMA Empenage H

        X_CMA_EH = [XOH_prime+0.25*LOH, XOH_prime-0.75*LOH]
        Y_CMA_EH = [YOH, YOH]

        # Coordonee Empenage H

        X_EH = [X_CMA_EH[0]+XOH, X_CMA_EH[0]+XOH-X3H,
                X_CMA_EH[0]+XOH-X3H-L3H, X_CMA_EH[0]+XOH-L1H]
        Y_EH = [0, Y3H, Y3H, 0]

        X_EH_sym = symetrie(X_EH)
        Y_EH_sym = symetrie(Y_EH)

        # Empenage V

        LpV = L_fuselage1*0.88-Lfa
        SV = (self.Vv*s_ref*B)/LpV

        Z2V = math.sqrt(self.lamdaV*SV)
        L1V = (SV*2)/(Z2V*(1+self.epsilonV))
        L2V = self.epsilonV*L1V

        L0V = Z2V*(L1V**2+L2V**2+L1V*L2V)/(3*SV)
        X2V = L1V/4+Z2V*math.tan(self.phi25V*math.pi/180)-L2V/4
        X0V = X2V*Z2V*(2*L2V+L1V)/(6*SV)
        Z0V = (Z2V**2)*(2*L2V+L1V)/(6*SV)

        XOH_prime = X_Foyer[0]-LpV
        X_CMA_EV = [XOH_prime+0.25*L0V, XOH_prime-0.75*L0V]
        Z_CMA_EV = [Z0V, Z0V]

        X_EV = [X_CMA_EV[0]+X0V, X_CMA_EH[0]+X0V-X2V,
                X_CMA_EV[0]+X0V-X2V-L2V, X_CMA_EV[0]+X0V-L1V]
        Y_EV = [Y_fin_cyl, Z2V+Y_fin_cyl, Z2V+Y_fin_cyl, Y_fin_cyl]

        # moteur

        L_nac = 3.201
        D_nac = 1.602

        X_nac = [X_CMA[0]+L_nac/2, X_CMA[0]+L_nac/2, X_CMA[0] -
                 L_nac/2, X_CMA[0]-L_nac/2, X_CMA[0]+L_nac/2]
        Y_nac = [Y_CMA[0]-D_nac/2, Y_CMA[0]+D_nac/2, Y_CMA[0] +
                 D_nac/2, Y_CMA[0]-D_nac/2, Y_CMA[0]-D_nac/2]
        Y_nac_sym = symetrie(Y_nac)

        SPF = (L2+L3)*(Y3-Y2)+(L3+L4)*(Y2-Y3)

        return Cassure, L0, L2, LOH, SH, effilement, Y4, Y_debut_cyl, X_fuselage, Y_fuselage, Y_fuselage_sym, X, Y_sym, X_EH, Y_EH, Y_EH_sym, X_nac, Y_nac, Y_nac_sym, Y, MFW, X_debut_cyl, X_fin_cyl, X_pt_av, Y_pt_av, X_pt_ar, Y_pt_ar, X_EV, Y_EV, x_a, y_a, SV, B, Y2, Extremite, Emplanture, SPF

    def draw(self, nfront, front_reduit, er_AV, aff_AV, er_cyl, aff_cyl, eco_cyl, rang_nonreduitAR, rang_reduitAR, X_debut_cyl, X_fin_cyl, Y_debut_cyl, X_fuselage, Y_fuselage, Y_fuselage_sym, X, Y, Y_sym, X_EH, Y_EH, Y_EH_sym, X_nac, Y_nac, Y_nac_sym, X_av, Y_av, X_ar, Y_ar, XEV, YEV, xa, ya):

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7, 3))

        j = er_AV
        i = 0
        bas = 0
        haut1 = 0
        haut = 0

        if j != 0:
            if nfront-2 == 3:
                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+1)/39.3700787)),
                                    40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+1)/39.3700787)),
                                    (40)/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+1) /
                                    39.3700787)), (40)/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))

                    i = i+40/39.3700787
                    j = j-1
            if nfront-2 == 4:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-haut1-(
                        ((42/2)+1+(62/6))/39.3700787)), 40/39.3700787, 21/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-haut1-(((42/2)+(
                        42/2)+1+(62/6))/39.3700787)), (40)/39.3700787, (21)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((42/2)+19+19+(
                        42/2)+(42/2)+1+(62/6))/39.3700787)), 40/39.3700787, 21/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((42/2)+19+(42/2)+(
                        42/2)+1+(42/2)+19+(62/6))/39.3700787)), (40)/39.3700787, (21)/39.3700787, fc='none', color='chocolate'))
                    i = i+40/39.3700787
                    j = j-1
                    haut1 = haut1+(62/12)/39.3700787

            if nfront-2 == 5:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)-haut), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((62/3)+1+(62/3)+(
                        62/3)+(62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((62/3)+(62/3)+1+(62/3)+(
                        62/3)+(62/3)+19)/39.3700787)), (40)/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((42/2)+(42/2)+19+(62/3)+(62/3)+(
                        62/3)+19+(42/2)+(42/2)+1)/39.3700787)+bas), (40)/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    i = i+40/39.3700787
                    j = j-1
                    bas = bas+(62/6)/39.3700787
                    haut = haut+(62/6)/39.3700787

            if nfront-2 == 6:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((62/3)+1+(
                        62/3))/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((62/3)+1+(
                        62/3)+(62/3))/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((62/3)+(62/3)+1+(
                        62/3)+(62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+1+(
                        62/3)+(62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl+(aff_AV*36/39.3700787)+i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+(
                        62/3)+1+(62/3)+(62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    j = j-1
                    i = i+40/39.3700787

        j = rang_nonreduitAR
        i = 34/39.3700787

        if j != 0:

            if nfront == 5:
                while j > 0:
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+1)/39.3700787)),
                                    34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+1+(62/3))/39.3700787)),
                                    34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+1+(62/3)+(
                        62/3))/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((42/2)+(62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((42/2)+(42/2)+(62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    i = i+34/39.3700787
                    j = j-1

            if nfront == 6:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+1)/39.3700787)),
                                    34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+1+(62/3))/39.3700787)),
                                    34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+1+(62/3)+(
                        62/3))/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+(62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))

                    j = j-1
                    i = i+34/39.3700787

            if nfront == 7:

                while j > 0:

                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((42/2)+1)/39.3700787)),
                                    (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((42/2)+(42/2)+1)/39.3700787)),
                                    (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+19+(42/2)+(
                        42/2)+1)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+19+(42/2)+(
                        42/2)+1)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+19+(42/2)+(
                        42/2)+1)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((42/2)+19+(62/3)+(62/3)+(62/3)+19+(
                        42/2)+(42/2)+1)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((42/2)+(42/2)+19+(62/3)+(62/3)+(62/3)+19+(
                        42/2)+(42/2)+1)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))

                    i = i+34/39.3700787
                    j = j-1

            if nfront == 8:

                while j > 0:

                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((42/2)+1)/39.3700787)),
                                    (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((42/2)+(42/2)+1)/39.3700787)),
                                    (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-(((82/4)+19+(42/2)+(
                        42/2)+1)/39.3700787)), 34/39.3700787, (82/4)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-((((82/4))+((82/4))+19+(42/2)+(
                        42/2)+1)/39.3700787)), 34/39.3700787, (82/4)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-((((82/4))+((82/4))+((82/4))+19+(
                        42/2)+(42/2)+1)/39.3700787)), 34/39.3700787, (82/4)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-((((82/4))+(82/4)+((82/4))+((82/4))+19+(
                        42/2)+(42/2)+1)/39.3700787)), 34/39.3700787, (82/4)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-((((82/4))+(82/4)+((82/4))+((82/4))+19+(
                        42/2)+19+(42/2)+(42/2)+1)/39.3700787)), 34/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-i, Y_debut_cyl-((((82/4))+(82/4)+((82/4))+((82/4))+19+(42/2)+19+(
                        42/2)+(42/2)+(42/2)+1)/39.3700787)), 34/39.3700787, (42/2)/39.3700787, fc='none', color='black'))

                    i = i+34/39.3700787
                    j = j-1

        j = rang_reduitAR
        i = 34/39.3700787
        z = 0

        if rang_reduitAR != 0:

            if nfront-front_reduit == 1:
                if nfront == 6:
                    while j > 0:
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (62/3+62/6+22+1)/39.3700787)-z), 34/39.3700787, 22/39.3700787, fc='none', color='black'))

                        i = i+34/39.3700787
                        j = j-1
                        z = 4/39.3700787
                if nfront == 7:
                    while j > 0:
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (42/4+42/2+1)/39.3700787)-z), 34/39.3700787, 22/39.3700787, fc='none', color='black'))

                        i = i+34/39.3700787
                        j = j-1
                        z = 4/39.3700787
                else:
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*32/39.3700787)-i, Y_debut_cyl-(
                        ((62/6)+(42/2)+1)/39.3700787)), 32/39.3700787, 21/39.3700787, fc='none', color='black'))
                    i = i+32/39.3700787
                    j = j-1
            if nfront-front_reduit == 2:
                if nfront == 6:
                    while j > 0:
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (62/6+(42/2)+1)/39.3700787)-z), 34/39.3700787, 21/39.3700787, fc='none', color='black'))
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (62/6+2*(42/2)+1)/39.3700787)-z), (34)/39.3700787, (21)/39.3700787, fc='none', color='black'))
                        i = i+34/39.3700787
                        j = j-1
                        z = 4/39.3700787+z
                if nfront == 7:
                    while j > 0:
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (42/4+(42/2)+1)/39.3700787)-z), 34/39.3700787, 21/39.3700787, fc='none', color='black'))
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (42/4+2*(42/2)+1)/39.3700787)-z), (34)/39.3700787, (21)/39.3700787, fc='none', color='black'))
                        i = i+34/39.3700787
                        j = j-1
                        z = 4/39.3700787+z
            if nfront-front_reduit == 3:
                while j > 0:
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+(62/3)+1)/39.3700787)-z), (34)/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((62/3)+(
                        62/3)+1+(62/3))/39.3700787)-z), (34)/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    i = i+34/39.3700787
                    j = j-1
            if nfront-front_reduit == 4:
                if nfront == 6:
                    while j > 0:
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (62/3+(42/2)+1)/39.3700787)-z), 34/39.3700787, 21/39.3700787, fc='none', color='black'))
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (62/3+2*(42/2)+1)/39.3700787)-z), (34)/39.3700787, (21)/39.3700787, fc='none', color='black'))
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (62+19+(42/2)+(62/6))/39.3700787)-z), 34/39.3700787, 21/39.3700787, fc='none', color='black'))
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (62+19+2*(42/2)+(62/6))/39.3700787)-z), (34)/39.3700787, (21)/39.3700787, fc='none', color='black'))

                        i = i+34/39.3700787
                        j = j-1
                if nfront == 7:
                    while j > 0:
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            ((42)+1)/39.3700787)-z), 34/39.3700787, 21/39.3700787, fc='none', color='black'))
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (42+19+62/6+42/2)/39.3700787)-z), (34)/39.3700787, (21)/39.3700787, fc='none', color='black'))
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (42+19+62/6+42/2+(42/2))/39.3700787)-z), 34/39.3700787, 21/39.3700787, fc='none', color='black'))
                        ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                            (62+2*19+(42/2)+42)/39.3700787)-z), (34)/39.3700787, (21)/39.3700787, fc='none', color='black'))

                        i = i+34/39.3700787
                        j = j-1

            if nfront-front_reduit == 5:
                while j > 0:
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1+(62/3))/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1+(62/3)+(62/3))/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((42/2)+(62/3)+1+(
                        62/3)+(62/3)+19)/39.3700787)-z), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((42/2)+(42/2)+(
                        62/3)+1+1+(62/3)+(62/3)+19)/39.3700787)-z), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    i = i+34/39.3700787
                    j = j-1

            if nfront-front_reduit == 6:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1+(62/3))/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1+(62/3)+(62/3))/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((62/3)+(
                        62/3)+1+(62/3)+(62/3)+19)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((62/3)+(62/3)+(
                        62/3)+1+(62/3)+(62/3)+19)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+(
                        62/3)+1+(62/3)+(62/3)+19)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))

                    j = j-1
                    i = i+34/39.3700787
            if nfront-front_reduit == 7:

                while j > 0:

                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((42/2)+1)/39.3700787)-z), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(
                        ((42/2)+(42/2)+1)/39.3700787)-z), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((62/3)+19+(
                        42/2)+(42/2)+1)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((62/3)+(
                        62/3)+19+(42/2)+(42/2)+1)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((62/3)+(62/3)+(
                        62/3)+19+(42/2)+(42/2)+1)/39.3700787)-z), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((42/2)+19+(62/3)+(62/3)+(
                        62/3)+19+(42/2)+(42/2)+1)/39.3700787)-z), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_fin_cyl-(rang_nonreduitAR*34/39.3700787)-i, Y_debut_cyl-(((42/2)+(42/2)+19+(62/3)+(
                        62/3)+(62/3)+19+(42/2)+(42/2)+1)/39.3700787)-z), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))

                    i = i+34/39.3700787
                    j = j-1

        j = er_cyl
        i = 40/39.3700787

        if er_cyl != 0:

            if nfront-2 == 3:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+1)/39.3700787)),
                                    40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+1)/39.3700787)),
                                    (40)/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+1+19+(62/3)+(
                        62/3))/39.3700787)), (40)/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    i = i+40/39.3700787
                    j = j-1

            if nfront-2 == 4:
                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((42/2)+1)/39.3700787)),
                                    40/39.3700787, 21/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((42/2)+(42/2)+1)/39.3700787)),
                                    (40)/39.3700787, (21)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((42/2)+(62/3)+(62/3)+19+(42/2)+(
                        42/2)+1)/39.3700787)), 40/39.3700787, 21/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((42/2)+19+(42/2)+(62/3)+(42/2)+1+(
                        42/2)+19)/39.3700787)), (40)/39.3700787, (21)/39.3700787, fc='none', color='chocolate'))

                    i = i+40/39.3700787
                    j = j-1

            if nfront-2 == 5:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+1)/39.3700787)),
                                    40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+1+(62/3)+(62/3)+(
                        62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+1+(62/3)+(62/3)+(
                        62/3)+19)/39.3700787)), (40)/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((42/2)+(42/2)+19+(62/3)+(62/3)+(62/3)+19+(
                        42/2)+(42/2)+1)/39.3700787)), (40)/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    i = i+40/39.3700787
                    j = j-1

            if nfront-2 == 6:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+1)/39.3700787)),
                                    40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+1+(62/3))/39.3700787)),
                                    40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+1+(62/3)+(
                        62/3))/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+(62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 40/39.3700787, (62/3)/39.3700787, fc='none', color='chocolate'))

                    i = i+40/39.3700787
                    j = j-1

        j = aff_cyl
        i = 36/39.3700787

        if aff_cyl != 0:

            if nfront-1 == 4:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(
                        62/3)+1)/39.3700787)), (36)/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(62/3)+(
                        62/3)+1)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+(42/2)+(62/3)+1+(
                        62/3)+(62/3)+19)/39.3700787)), (36)/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    i = i+36/39.3700787
                    j = j-1

            if nfront-1 == 5:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+1+(
                        62/3))/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+1+(
                        62/3)+(62/3))/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+(62/3)+1+(62/3)+(
                        62/3)+(62/3)+19)/39.3700787)), (36)/39.3700787, (42/2)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+(42/2)+(62/3)+(
                        62/3)+1+(62/3)+(62/3)+19)/39.3700787)), (36)/39.3700787, (42/2)/39.3700787, fc='none', color='firebrick'))
                    i = i+36/39.3700787
                    j = j-1

            if nfront-1 == 6:

                while j > 0:
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+1+(
                        62/3))/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+1+(62/3)+(
                        62/3)+19)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(62/3)+1+(
                        62/3)+(62/3)+19)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+1+(
                        62/3)+(62/3)+19)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+(42/2)+19+(62/3)+(62/3)+(
                        62/3)+19+(42/2)+(42/2)+1)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))

                    i = i+36/39.3700787
                    j = j-1

            if nfront-1 == 7:

                while j > 0:

                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((42/2)+1)/39.3700787)), (36)/39.3700787, (42/2)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+(
                        42/2)+1)/39.3700787)), (36)/39.3700787, (42/2)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+19+(
                        42/2)+(42/2)+1)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(62/3)+19+(
                        42/2)+(42/2)+1)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(62/3)+(62/3)+19+(
                        42/2)+(42/2)+1)/39.3700787)), 36/39.3700787, (62/3)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+19+(62/3)+(62/3)+(
                        62/3)+19+(42/2)+(42/2)+1)/39.3700787)), (36)/39.3700787, (42/2)/39.3700787, fc='none', color='firebrick'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+(42/2)+19+(62/3)+(62/3)+(
                        62/3)+19+(42/2)+(42/2)+1)/39.3700787)), (36)/39.3700787, (42/2)/39.3700787, fc='none', color='firebrick'))

                    i = i+36/39.3700787
                    j = j-1
        j = eco_cyl
        i = 34/39.3700787

        if eco_cyl != 0:

            if nfront == 5:

                while X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i > X_fin_cyl:
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1+(62/3))/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1+(62/3)+(62/3))/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((42/2)+(62/3)+1+(62/3)+(62/3)+19)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+(
                        42/2)+(62/3)+1+(62/3)+(62/3)+19)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    i = i+34/39.3700787
                    j = j-1

            if nfront == 6:

                while X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i > X_fin_cyl:
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1+(62/3))/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+1+(62/3)+(62/3))/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+(62/3)+1+(62/3)+(62/3)+19)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(
                        62/3)+(62/3)+1+(62/3)+(62/3)+19)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(
                        62/3)+(62/3)+(62/3)+1+(62/3)+(62/3)+19)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))

                    i = i+34/39.3700787
                    j = j-1

            if nfront == 7:

                while X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i > X_fin_cyl:

                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((42/2)+1)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((42/2)+(42/2)+1)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+19+(42/2)+(42/2)+1)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(
                        ((62/3)+(62/3)+19+(42/2)+(42/2)+1)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(((62/3)+(
                        62/3)+(62/3)+19+(42/2)+(42/2)+1)/39.3700787)), 34/39.3700787, (62/3)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+19+(
                        62/3)+(62/3)+(62/3)+19+(42/2)+(42/2)+1)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))
                    ax[0].add_patch(Rectangle((X_debut_cyl-(aff_cyl*36/39.3700787+er_cyl*40/39.3700787)-i, Y_debut_cyl-(((42/2)+(42/2)+19+(
                        62/3)+(62/3)+(62/3)+19+(42/2)+(42/2)+1)/39.3700787)), (34)/39.3700787, (42/2)/39.3700787, fc='none', color='black'))

                    i = i+34/39.3700787
                    j = j-1

        c = ((Y_debut_cyl-Y_av)**2)/(4*((X_debut_cyl-X_av)))
        y = np.linspace(-Y_debut_cyl, Y_debut_cyl, 1000)
        x = ((y-Y_av)**2)/(4*c) + X_av

        c2 = ((Y_debut_cyl-Y_ar)**2)/(4*((X_fin_cyl-X_ar)))
        y2 = np.linspace(-Y_debut_cyl, Y_debut_cyl, 1000)
        x2 = ((y-Y_ar)**2)/(4*c2) + X_ar

        ax[0].plot(x, y, x2, y2, lw=1, color='blue')

        Y_ar_nouveau = +(X_fin_cyl-X_ar)*math.tan(15*(math.pi/180))-Y_debut_cyl
        x_1 = np.linspace(XEV[0], XEV[3], 100)
        list1x = x_1.tolist()
        y_1 = np.linspace(YEV[0], YEV[3], 100)
        list1y = y_1.tolist()
        x_2 = np.linspace(XEV[3], X_ar, 100)
        list2x = x_2.tolist()
        y_2 = np.linspace(YEV[3], Y_ar_nouveau, 100)
        list2y = y_2.tolist()
        x_3 = np.linspace(X_ar, X_fin_cyl, 100)
        list3x = x_3.tolist()
        y_3 = np.linspace(Y_ar_nouveau, -Y_debut_cyl, 100)
        list3y = y_3.tolist()
        X1 = [X_fin_cyl]+list1x+list2x+list3x
        Y1 = [Y_debut_cyl]+list1y+list2y+list3y

        # NACA profil:

        def camber_line(x, m, p, c):
            return np.where((x >= 0) & (x <= (c*p)),
                            m * (x / np.power(p, 2)) * (2.0 * p - (x / c)),
                            m * ((c - x) / np.power(1-p, 2)) * (1.0 + (x / c) - 2.0 * p))

        def dyc_over_dx(x, m, p, c):
            return np.where((x >= 0) & (x <= (c*p)),
                            ((2.0 * m) / np.power(p, 2)) * (p - x / c),
                            ((2.0 * m) / np.power(1-p, 2)) * (p - x / c))

        def thickness(x, t, c):
            term1 = 0.2969 * (np.sqrt(x/c))
            term2 = -0.1260 * (x/c)
            term3 = -0.3516 * np.power(x/c, 2)
            term4 = 0.2843 * np.power(x/c, 3)
            term5 = -0.1015 * np.power(x/c, 4)
            return 5 * t * c * (term1 + term2 + term3 + term4 + term5)

        def naca4(x, m, p, t, c=1):
            dyc_dx = dyc_over_dx(x, m, p, c)
            th = np.arctan(dyc_dx)
            yt = thickness(x, t, c)
            yc = camber_line(x, m, p, c)
            return ((x - yt*np.sin(th), yc + yt*np.cos(th)),
                    (x + yt*np.sin(th), yc - yt*np.cos(th)))

        x4 = np.linspace(0, 1, 1000)

        ax[0].plot(X_fuselage, Y_fuselage, X_fuselage, Y_fuselage_sym, X, Y, X, Y_sym,
                   X_EH, Y_EH, X_EH, Y_EH_sym, X_nac, Y_nac, X_nac, Y_nac_sym, lw=1, color='blue')
        ax[0].axis('scaled')
        ax[0].xaxis.set_visible(False)
        ax[0].yaxis.set_visible(False)

        ax[1].plot(X_fuselage, Y_fuselage, X_fuselage,
                   Y_fuselage_sym, color='blue', lw=1)
        ax[1].plot(XEV, YEV)
        ax[1].plot(xa, ya)
        ax[1].plot(X1, Y1, color='blue', lw=1)
        #ax[1].plot(x3, y3, color='blue', lw=1)
        ax[1].plot(x, y, color='blue', lw=1)
        ax[1].axis('scaled')
        ax[1].xaxis.set_visible(False)
        ax[1].yaxis.set_visible(False)

        for spine in ['top', 'right', 'left', 'bottom']:
            ax[0].spines[spine].set_visible(False)
            ax[1].spines[spine].set_visible(False)

        plt.savefig('static/resultats/aircraft.png', dpi=1200)

        fig2, ax = plt.subplots(nrows=1, ncols=1)

        for item in naca4(x4, self.m, self.p, self.t, self.c):
            ax.plot(item[0], item[1], 'b', lw=1)

        ax.axis('scaled')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        for spine in ['top', 'right', 'left', 'bottom']:
            ax.spines[spine].set_visible(False)

        #ax[1,0].set_xlim((-0.05, 1.05))
        ax.set_ylim((-0.2, 0.2))

        plt.savefig('static/resultats/airfoil.png', dpi=1200)
        # plt.show()

    def Aero_GV(self, cassure, L0, L2, rextern, sref, LOH, SH, Lfus, Lav, Lar, lcyl, effilement, Y4):

        Re_sur_m = 47899*((self.P*self.mach*((1+0.126*math.pow(self.mach, 2))*self.T+110.4))
                          )/(math.pow(self.T, 2)*math.pow((1+0.126*self.mach), 5/2))

        # Cx0 Voilure

        # Facteur d'epaisseure relative
        Kev = 4.688*(cassure**2)+3.146*(cassure)

        # Facteur de cambrure:
        Cz_sur_cosphi_carre = self.Cz/(math.cos(self.phi25*(math.pi)/180))
        Kc = 2.859*math.pow(Cz_sur_cosphi_carre, 3)-1.849 * \
            math.pow(Cz_sur_cosphi_carre, 2)+0.382*(Cz_sur_cosphi_carre)+0.06

        # Facteur de fleche:
        Kphi = 1-0.000178 * \
            math.pow((self.phi25*math.pi/180), 2)-0.0065*self.phi25*math.pi/180

        # Facteur d'iteration voilure/fuselage:
        Kiv = 0.04

        # coef de frottement:
        cfv = 0.455/(1+0.126*(self.mach**2)) * \
            (1/((math.log10(Re_sur_m*L0))**2.58))

        # surface mouillee de la voilure
        SMA = 2*(sref-rextern*L2)

        # CxOA
        CxOA = ((Kev+Kc)*Kphi+Kiv+1)*cfv*SMA/sref

        # ---------------------------------------------
        # Cx0 Empennage

        # Facteur d'epaisseure relative
        Kee = 4.688*((10/100)**2)+3.146*(10/100)

        # Facteur de fleche:
        Kphie = 1-0.000178 * \
            math.pow((self.phi25H*math.pi/180), 2) - \
            0.0065*self.phi25H*math.pi/180

        # Facteur d'iteration voilure/fuselage:
        Kie = 0.04

        # coef de frottement:
        cfe = 0.455/(1+0.126*self.mach**2)/(math.log10(Re_sur_m*LOH))**2.58

        SMH = 2*SH

        # CxOH
        CxOH = ((Kee*Kphie+Kie/4+1)*cfe*SMH/sref)

        # ---------------------------------------

        # fuselage

        # coef de frottement:
        cff = 0.455/(1+0.126*self.mach**2)/(log10(Re_sur_m*Lfus))**2.58

        SMF = 2.45*2*rextern*Lav+math.pi*2*rextern*lcyl+2.3*2*rextern*Lar

        # CxOF
        CxOF = (0.98+0.754*2*rextern/sref)*cff*SMF/sref
        # -------------------------------------------------

        # Parasite

        # KT
        Kt = 0.8

        # Smtot
        Smtot = SMA+SMH+SMF

        # Kp
        Kp = -0.00000000000239*(Smtot**3)+0.0000000258 * \
            (Smtot**2)-0.000089*Smtot+0.163

        # Cx Parasite
        Cx_para = Kt*Kp*(CxOA+CxOH+CxOF)

        # Derive

        # estimation
        estimation = CxOH

        # Nacelle+mats
        forfait = 0.002

        Cx0 = CxOA+CxOH+CxOF+Cx_para+estimation+forfait

        # trainee induite
        delta0 = 0.233*(effilement**2)-0.068*effilement+0.012
        fv = (0.7/self.Cz-1)**2
        Cxi = (1.03+delta0*fv+(2/rextern/(2*Y4))**2) * \
            (self.Cz**2)/(math.pi*self.lamda)

        # trainee d'equilibrage
        delta_cx_eq = 0.000589*self.Cz

        # trainee de compressibilite
        delta_cx_c = 0.002

        Cx_tot = Cx0+Cxi+delta_cx_eq+delta_cx_c
        finesse_croisiere = self.Cz/Cx_tot

        return Cx_tot, finesse_croisiere, Cx0, Cxi, SMF, Smtot

    def Aero_BV(self, Cx0):

        Cx0BV = Cx0
        delta_Cx0hyper = 0.006
        KBV = 1.05
        Czmax = 2.523

        Cz = Czmax/(1.2**2)
        Cx = Cx0BV+delta_Cx0hyper+KBV*(Cz**2)/(math.pi*self.lamda)
        finesse_D_sur_L = Cz/Cx

        return finesse_D_sur_L, Cz

    def Masse(self, MTOW, MZFW, MFW, nfront, lcyl, SMF, rext, SV, SH, B, Lfus, Y4, Y2, Smtot, L2, extremite, Cassure, emplanture, sref, SPF, finesse_D_sur_L, finesse_croisiere, Lpax):

        MLW = 1.06*MZFW

        # masse voilure:
        nextr = 1.5*2.5
        m1 = 1.05*MZFW
        nm1 = nextr*m1
        m2 = MTOW
        MCV = min(0.8*MFW, (MTOW-MZFW))
        nm2 = nextr*(m2-0.55*MCV)

        max_nm = max(nm1, nm2)

        # Masse voilure
        Kvoil = 1.05
        Kvmot = 1.39
        Lemp = L2
        e_sur_l = 0.5*emplanture+(1/3)*Cassure+(1/6)*extremite
        flexion = 0.00005922*Kvoil * \
            ((max_nm/(Lemp*e_sur_l))*((B/math.cos(self.phi25*math.pi/180))**2))**0.9
        cisaillement = 0.0005184*Kvoil * \
            (max_nm*B/math.cos(self.phi25*math.pi/180))**0.9
        nervure = Kvoil*(1.7009*sref+0.001*max_nm)
        renfort = 0.0044*Kvoil*((MLW)**1.0169)
        secondaire = 0.3285*Kvoil*((MZFW)**0.35)*SPF*Kvmot
        A1 = flexion+cisaillement+nervure+renfort+secondaire

        # masse fuselage
        ktr = 1.05
        kfus = 1
        A2 = SMF*(10+1.2*rext+0.00019*nm1/(rext**1.7))*kfus*ktr

        # Masse empennage
        A3H = SH*(14.4+0.155*SH)
        A3V = SV*(15.45+0.202*SV)
        A3 = A3H+A3V

        # Masse Commande de vol

        A4 = 0.000085*max_nm * \
            (Lfus**0.66+(B/math.cos(self.phi25*math.pi/180))**0.66)

        # Masse Atterisseure
        HTR = 0.9*0.18*Y4
        A5 = 0.028*MTOW+0.12*(HTR**1.5)*(MLW**0.66)

        envergure = B

        # Contrainte decollage
        pente_min = 2.4
        finesse_decollage = finesse_D_sur_L
        FN0_deco = (MTOW*9.81*pente_min/100+MTOW*9.81/finesse_decollage)

        # Contrainte Croisiere
        m_mission = (MTOW+MZFW)/2

        rho = self.P/(287*self.T)
        Vcr = self.mach*math.sqrt(1.4*self.T*287)
        Vzmin = 500*0.3048/60
        pente_CR = math.asin(Vzmin/Vcr)
        finesse_CR = finesse_croisiere
        FN = m_mission*9.81*pente_CR+m_mission*9.81/finesse_CR
        FN0_croisiere = FN*1.225/rho/2

        FN0 = max(FN0_croisiere.real, FN0_deco.real)
        FN0_lbf = FN0*0.224809

        B1 = 0.22*self.Nmot*(FN0_lbf**0.939)

        # Masse Mats
        SMNAC = 0.4*(10**(-3))*FN0_lbf+11
        SMmat = 0.35*SMNAC
        A6 = 1.2*((SMmat*0.6)**0.5) * \
            (23+0.588*((B1/self.Nmot)**0.708))*self.Nmot

        # Masse peinture
        A7 = 0.18*(Smtot+(SMNAC+SMmat)*self.Nmot)

        # Masse tot
        A = A1+A2+A3+A4+A5+A6+A7
        #print(A , A1,A2,A3,A4,A5,A6,A7)
        B2 = 0.02*B1+2*envergure / \
            math.cos(self.phi25*math.pi/180)+0.35*(MZFW**0.66)
        B3 = 25*self.Nmot+0.0035*MFW

        B = B1+B2+B3

        K = 1
        C11 = 11.3*(self.npax**0.64)
        C12 = (0.444*(MTOW**0.66)+2.54*self.npax+0.254*A4)*K
        C13 = (0.256*(MTOW**0.66)+1.46*self.npax+0.254*A4)*K

        C1 = C11+C12+C13

        LCP = Lfus*0.8
        DF = rext*2
        if self.r_action < 3000:
            a = 21
            b = 40
        if 3000 <= self.r_action < 4500:
            a = 17
            b = 35
        else:
            a = 16
            b = 30
        PNC = (self.npax+a)//b
        PNT = 2
        Lsoute = lcyl+0.864*(nfront-5)-0.8*L2

        if nfront > 6:
            Vsoute = 0.541*(self.npax-38)
            amenagement_soute = 4.66*(Vsoute**0.67)+8*LCP
            D1 = 0  # a voir

        if nfront == 6:
            Vsoute = 0.351*(self.npax-38)
            amenagement_soute = 15.4*(Vsoute**0.67)+16*LCP
            D1 = Vsoute

        if nfront <= 5:
            Vsoute = 0.351*(self.npax-38)
            amenagement_soute = 11.3*(Vsoute**0.67)+5*LCP
            D1 = Vsoute

        C21 = (3.6+5.7)*DF*LCP

        if self.r_action <= 4500:
            C22 = 200+(27*(self.npax**0.46))+(7.2*(self.Nmot**0.7)
                                              * (self.npax**0.64))+self.npax+0.0029*(self.npax**0.64)
        else:
            C22 = 450+(51*(self.npax**0.46))+(7.2*(self.Nmot**0.7)
                                              * (self.npax**0.64))+self.npax+0.0029*(self.npax**0.64)

        C23 = 53+2*1.9*9.5*(Y4-Y2)/math.cos(1.2*self.phi25*math.pi/180)
        C24 = 1.4*LCP*DF
        C25 = 27*PNT+18*PNC
        C26 = 80+1.3*self.npax
        C27 = 0.01*B1+2.3*self.npax
        if self.r_action < 3000:
            C28 = 14.5*Lpax*0.0254+40 * \
                (self.npax//50)+0.2*self.npax+1.6*Lpax * \
                0.0254*LCP+1*self.npax+amenagement_soute

        if 3000 <= self.r_action < 4500:
            C28 = 14.5*Lpax*0.0254+65 * \
                (self.npax//50)+0.25*self.npax+1.6*Lpax * \
                0.0254*LCP+3.2*self.npax+amenagement_soute
        else:
            C28 = 14.5*Lpax*0.0254+65 * \
                (self.npax//50)+0.25*self.npax+1.6*Lpax * \
                0.0254*LCP+5.9*self.npax+amenagement_soute

        C2 = C21+C22+C23+C24+C25+C26+C27+C28

        if self.r_action <= 1500:
            a = 150
            k_D5 = 0.1
            C4 = 100
            C6 = 10
            D2 = 9*self.npax
            nourriture = 4*self.npax
            eau = 0.5*(self.npax+PNC)

        if 1500 < self.r_action <= 3000:
            a = 450
            k_D5 = 0.5
            C4 = 200
            C6 = 45
            D2 = 10*self.npax
            nourriture = 4*self.npax
            eau = 0.5*(self.npax+PNC)
        if 3000 < self.r_action <= 4500:
            a = 700
            k_D5 = 1
            C4 = 250
            C6 = 45
            D2 = 10*self.npax
            nourriture = 6*self.npax
            eau = 1.12*(self.npax+PNC)
        else:
            a = 800
            k_D5 = 1.5
            C4 = 350
            C6 = 45
            D2 = 11*self.npax
            nourriture = 9*self.npax
            eau = 1.75*(self.npax+PNC)

        C3 = a+0.033*Lfus*(Y4*2/math.cos(self.phi25*math.pi/180))

        C5 = 100+23.4*Lsoute*2

        C = C1+C2+C3+C4+C5+C6
        D3 = nourriture+eau
        D4 = 1.5*self.npax
        D5 = k_D5*self.npax

        D = D1+D2+D3+D4+D5

        E = 85*PNT+75*PNC

        G = 120*self.npax

        F = MTOW-(A+B+C+D+E+G)
        #print(A1, A2, A3, A4, A5, A6,A7)
        return F, A, B, C, D, E, G, FN0_lbf, B1, nourriture

    def perfo(self, MTOW, Cx0, Cxi):

        Mi = MTOW
        ki = Cxi/((self.Cz)**2)

        Czopt = math.sqrt(Cx0.real/(3*ki.real))
        Cxopt = Cx0+ki*(Czopt**2)

        V_croisiere = self.mach*math.sqrt(1.4*287*self.T)
        s_voil = (Mi*9.81)/(0.7*self.P*(self.mach**2)*Czopt)

        Finesse = Czopt/Cxopt
        Dfran = self.r_action*1852

        Mcarb = Mi*(1-math.exp(-Dfran*self.Cs*9.81/(V_croisiere*Finesse.real)))

        return V_croisiere, Finesse, Mcarb, Czopt, Cxopt, s_voil

    def boucle_conv(self, A, B, C, D, E, G, B1, nourriture):

        MV = A+B+C  # Masse a vide
        MVE = MV+D  # masse a vide equipe
        MVOE = MVE+E  # masse a vide en ordre d'exploitation
        ZFW = MVOE+G  # masse a vide sans carburant
        MC = MVOE-(E+B1+nourriture)
        return MC, ZFW

    def TOC(self, FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW):

        Tbloc = (10+1+2+(self.r_action/v_croi)*60+5)/60
        Gbloc = self.Nmot*FN0*0.002+self.Nmot*FN0*0.007+4.274*(10**(-7))*MTOW*self.altitude+(self.r_action/v_croi)*2.15*(10**(-5))*math.sqrt(
            (288.15-6.5*(10**(-3))*self.altitude/3.281)/288.15)*math.pow(self.mach, 3.51*(10**(-2))*8.4-1.27*(10**(-5))*self.altitude+31)+self.Nmot*FN0*0.001
        Kdoc = 1.5
        LW = 1.06*ZFW
        MLW = 1.06*MZFW
        prix_billet = 0.25  # $/NM/pax

        if self.r_action < 4500:
            U = (3750*Tbloc)/(Tbloc+0.5)
        else:
            U = 4600

        # Prix cellule:
        PC = 2541*Kdoc*((MC)**0.887)

        # Prix des moteur:
        PM = self.Nmot*(10**6)*Kdoc*(0.885+6.75*(10**(-5))*FN0)

        # Prix de l'avion:
        PA = PM+PC

        # Investissment tot:
        IT = 0.1*PC+PA+0.3*PM

        # Calcule du DOC:
        # ---------------

        # Amortissement et Interet:
        AM_plus_I = (IT*Tbloc)/(14*U) + 0.0296*(IT*Tbloc)/(U)

        # Assurance:
        AS = (0.005*PA*Tbloc)/U

        # Carburant:
        C = 0.3344*((Gbloc/0.79)/3.7854)*0.8

        # Taxes:
        # ------

        if self.r_action < 4500:
            TA = 0.0045*Kdoc*MTOW
            TN = 0.003*Kdoc*self.r_action*math.sqrt(MTOW)
        else:
            TA = 0.005*Kdoc*MTOW
            TN = 0.0061*Kdoc*self.r_action*math.sqrt(MTOW)

        TX = TA+TN

        # Equipage:
        if self.r_action < 4500:
            ET = 310*Kdoc*Tbloc
            EC = 52*(math.ceil(self.npax/35))*Kdoc*Tbloc
        else:
            ET = 380*Kdoc*Tbloc
            EC = 60*(math.ceil(self.npax/35))*Kdoc*Tbloc
        E = ET+EC
        # Maintenance:
        # ------------

        MCMO = ((9*(10**(-5))*MC+6.7-630/(0.0018*MC+135))
                * (1+0.59*(Tbloc-0.25)))*40*Kdoc
        MCP = (6.24+3.08*(Tbloc-0.25)*PC*(10**(-6)))*Kdoc

        # Moteur:
        # -------

        MMMO = self.Nmot*40*Kdoc * \
            (1.44*FN0/49000 + 0.42+(0.8*FN0/49000 + 0.24)*Tbloc)
        MMP = self.Nmot*40*Kdoc * \
            (73.4*FN0/49000 + 5.31+(41.3*FN0/49000 + 3)*Tbloc)

        M = MCMO+MCP+MMMO+MMP

        # DOC

        DOC = AM_plus_I+AS+C+TX+E+M  # $/vol
        DOC_par_Nm_pax = (DOC/(self.r_action*self.npax))*100

        # I.O.C:
        # -----

        # Instalation au sol
        IS = 0.00081*MLW

        # Commisatiat:
        CO = 2.58*Tbloc*self.npax*self.Coef_rempli

        # service pax:
        if self.r_action < 3000:
            SP = 7.3*self.npax*self.Coef_rempli
        else:
            SP = 8*self.npax*self.Coef_rempli

        # Publicite-reservation:
        PU = 0.0083*(self.npax+self.npax*40/200)*self.r_action*self.Coef_rempli

        # Frais generaux:
        FG = 0.053*((DOC-AM_plus_I-AS)+0.6*IS+CO+SP+PU)

        IOC = (IS+CO+SP+PU+FG)*Kdoc

        # Revenu Par vol
        RP = self.Coef_rempli*self.r_action * ((self.npax-math.ceil(0.2*self.npax)-math.ceil(0.3*self.npax))*(1+40/200)*self.prix_billet_eco + (
            (math.ceil(0.3*self.npax))*(1+40/200))*self.prix_billet_aff+((math.ceil(0.2*self.npax))*(1+40/200))*self.prix_billet_er)

        Prix_Cellule_d = PC/(10**6)
        Prix_Moteur_d = PM/(10**6)
        Prix_Avion_d = PA/(10**6)
        Investissement_totale_d = IT/(10**6)
        DOC_d = DOC
        DOC_Nm_pax_d = DOC_par_Nm_pax
        IOC_d = IOC
        Revenu_par_vol_d = RP
        # Coef_de_remplissage_au_pt_d_equilibre=coe_pt_eq*100
        Prix_Cellule_e = PC*0.94/(10**6)
        Prix_Moteur_e = PM*0.94/(10**6)
        Prix_Avion_e = PA*0.94/(10**6)
        Investissement_totale_e = IT*0.94/(10**6)
        DOC_e = DOC*0.94
        IOC_e = IOC*0.94
        Revenu_par_vol_e = RP*0.94

        return Prix_Cellule_d, Prix_Moteur_d, Prix_Avion_d, Investissement_totale_d, DOC_d, DOC_Nm_pax_d, IOC_d, Revenu_par_vol_d, Prix_Cellule_e, Prix_Moteur_e, Prix_Avion_e, Investissement_totale_e, DOC_e, IOC_e, Revenu_par_vol_e

    def diagrame_tapis(self):

        def Lp(S, FN):
            m = 2*(FN*self.Nmot*self.r_action)/(self.Cz/math.pow(1.2, 2))
            P = math.pow(101325*(1-1000/44330.78), 5.25588)
            T = 288.1-0.0065*1000
            rho_0 = 1.225
            rho = 0.0034837*(P/T)
            sigma = rho/rho_0
            alpha = (m**2)/(self.Cz*self.Nmot*FN*S*(sigma**0.8))
            return 41.43*((alpha/100)**2)+12.713*alpha+287

        n = 4
        s = np.linspace(180, 240, n)
        fn = np.linspace(40000, 55000, n)
        S, FN = np.meshgrid(s, fn)

        Lp = Lp(S, FN)

        fig1 = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(S, FN, Lp)
        plt.title('Lp(S,FN)')
        plt.savefig('static/resultats/lp.png', dpi=1200)

        def V_app(S, FN):
            m = 2*(FN*self.Nmot*self.r_action)/(self.Cz/math.pow(1.2, 2))
            rho_0 = 1.225
            g = 9.81
            alpha = 2*m*g/(rho_0*S*self.Cz)
            return 1.3*(alpha**0.5)

        V_app = V_app(S, FN)

        fig2 = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(S, FN, V_app)
        plt.title('V app(S,FN)')
        plt.savefig('static/resultats/v_app.png', dpi=1200)

        def L_att(S, FN):
            m = 2*(FN*self.Nmot*self.r_action)/(self.Cz/math.pow(1.2, 2))
            rho_0 = 1.225
            g = 9.81
            alpha = 2*m*g/(rho_0*S*self.Cz)
            Vapp = 1.3*(alpha**0.5)
            if self.r_action < 3000:
                K = 0.9
            if 3000 <= self.r_action < 4500:
                K = 1
            if self.r_action >= 4500:
                K = 1.1

            return 260*K*np.expm1(0.0265*Vapp)

        Latt = L_att(S, FN)

        fig3 = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(S, FN, Latt)
        plt.title('L atterissage(S,FN)')
        plt.savefig('static/resultats/atterissage.png', dpi=1200)

        def Consommation(S, FN_phi):
            FN = 0.97*FN_phi*(1-self.mach+(self.mach**2)/2)
            P = math.pow(101325*(1-1000/44330.78), 5.25588)
            T = 288.1-0.0065*1000
            P0 = 101325
            T0 = 288.1
            teta = T/T0
            sigma = P/P0
            z = (FN/(FN_phi*sigma))
            Cs = ((0.331+0.3*(self.mach-0.8))*(z**2)-(0.542+0.75 *
                  (self.mach-0.8))*z+0.908+0.85*(self.mach-0.8))*(teta**0.5)
            return Cs

        Cs = Consommation(S, FN)

        fig3 = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(S, FN, Cs)
        plt.title('Consommation(S,FN)')
        

        plt.savefig('static/resultats/consommation.png', dpi=1200)

    def terminer(self):

        best_combinaison = []
        l = []
        lx = []
        mini1 = 20
        maxi = 0
        maxl = []
        nb_places_tot_er = math.ceil(0.1*self.npax)
        nb_places_tot_affaire = math.ceil(0.17*self.npax)
        state = False
        MTOW = self.initialisation()[0]
        MZFW = self.initialisation()[1]
        sref = self.sref
        i = 0
        while i < 20:

            mini1 = 20
            maxi = 0
            nb_places_tot_er = math.ceil(0.2*self.npax)
            nb_places_tot_affaire = math.ceil(0.3*self.npax)

            for nfront in self.npf():
                x = 0
                while nb_places_tot_affaire-x > nfront-1:
                    y = 0
                    while nb_places_tot_er-y > nfront-2:

                        nb_rang_AV = nfront-4
                        nb_rang_AR = nfront

                        Rang_reduit_AV = nfront-1
                        Rang_redui_AR = nfront-5

                        # nfront_reduit_AV=(nfront-1)
                        nfront_reduit_AV = (nfront-2)

                        if Rang_redui_AR < 0 or Rang_reduit_AV < 0:

                            Rang_redui_AR = 0
                            Rang_reduit_AV = 0

                        place_AV = Rang_reduit_AV*nfront_reduit_AV + \
                            ((nb_rang_AV-Rang_reduit_AV)*(nfront-2))  # (nfront-1)

                        places_er_AV = place_AV
                        places_affaire_AV = 0
                        places_eco_AV = 0

                        i_ar = 1

                        while nfront-i_ar >= 0:

                            place_AR = (nb_rang_AR)*nfront + \
                                Rang_redui_AR*(nfront-i_ar)

                            places_cyl_er = nb_places_tot_er-y-places_er_AV
                            places_cyl_affaire = nb_places_tot_affaire-x-places_affaire_AV
                            places_cyl_eco = self.npax - \
                                (place_AR+places_cyl_affaire+places_affaire_AV +
                                 places_cyl_er+places_er_AV+places_eco_AV)

                            place_cylindre = places_cyl_eco+places_cyl_affaire+places_cyl_er
                            nb_rang_cylindre = (
                                places_cyl_eco/nfront)+(places_cyl_affaire/(nfront-1))+(places_cyl_er/(nfront-2))

                            if (places_cyl_eco % nfront)+(places_cyl_affaire % (nfront-1))+(places_cyl_er % (nfront-2)) == 0 and places_cyl_affaire > 0 and places_cyl_er > 0 :

                                rang_eco = places_cyl_eco/nfront + \
                                    nb_rang_AR+places_eco_AV/(nfront-1)
                                places_eco_tot = places_cyl_eco+place_AR+places_eco_AV

                                rang_affaire = places_affaire_AV / \
                                    (nfront-1)+places_cyl_affaire/(nfront-1)
                                places_affaire_tot = places_affaire_AV+places_cyl_affaire

                                rang_er = places_er_AV / \
                                    (nfront-2)+(places_cyl_er/(nfront-2))
                                places_er_tot = places_er_AV+places_cyl_er

                                place_totale = place_AV+place_AR+place_cylindre
                                rang_totale = rang_eco+rang_affaire+rang_er

                                l_cylindre_m = self.l_cylindre(
                                    nfront, rang_eco, rang_affaire, rang_er, nb_rang_AV, nb_rang_AR)[0]
                                lpax_inch = self.l_cylindre(
                                    nfront, rang_eco, rang_affaire, rang_er, nb_rang_AV, nb_rang_AR)[1]

                                for longeur in l_cylindre_m:
                                    for Lpax in lpax_inch:
                                        for larav in self.rayon(nfront)[0]:
                                            for rext in self.rayon(nfront)[1]:
                                                rapport = self.L_sur_D(
                                                    longeur, larav, rext)[0]
                                                l = l+[rapport]

                                                if rapport > 10 and rapport < mini1:

                                                    mini1 = rapport
                                                    nfront1 = nfront

                                                    rang_eco1 = rang_eco
                                                    rang_affaire1 = rang_affaire
                                                    rang_er1 = rang_er

                                                    rang_er_AV1 = places_er_AV / \
                                                        (nfront1-2)
                                                    rang_affaire_AV1 = places_affaire_AV / \
                                                        (nfront1-1)
                                                    rang_eco_AV1 = places_eco_AV / \
                                                        (nfront1-1)

                                                    places_eco_tot1 = places_eco_tot
                                                    places_affaire_tot1 = places_affaire_tot
                                                    places_er_tot1 = places_er_tot

                                                    rang_er_cycl1 = places_cyl_er / \
                                                        (nfront1-2)
                                                    rang_affaire_cycl1 = places_cyl_affaire / \
                                                        (nfront1-1)
                                                    rang_eco_cycl1 = places_cyl_eco/nfront1

                                                    nfront_rang_reduit1 = nfront-i_ar

                                                    rang_eco_nonreduit1 = (
                                                        nb_rang_AR)
                                                    rang_eco_reduit1 = Rang_redui_AR
                                                    rang_totale1 = rang_totale
                                                    nb_rang_AV1 = nb_rang_AV
                                                    nb_rang_cylindre1 = nb_rang_cylindre
                                                    nb_rang_AR1 = nb_rang_AR
                                                    rayonextern1 = rext
                                                    lcyl1 = longeur
                                                    lav1 = self.rayon(
                                                        nfront1)[2][0]
                                                    lar1 = self.rayon(
                                                        nfront1)[3][0]
                                                    lfus1 = self.L_sur_D(
                                                        longeur, larav, rext)[1]
                                                    Cassure1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[0]
                                                    L01 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[1]
                                                    L21 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[2]
                                                    LOH1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[3]
                                                    SH1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[4]
                                                    effilement1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[5]
                                                    Y41 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[6]
                                                    Y_debut_cyl1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[7]
                                                    X_fuselage1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[8]
                                                    Y_fuselage1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[9]
                                                    Y_fuselage_sym1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[10]
                                                    X1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[11]
                                                    Y_sym1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[12]
                                                    X_EH1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[13]
                                                    Y_EH1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[14]
                                                    Y_EH_sym1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[15]
                                                    X_nac1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[16]
                                                    Y_nac1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[17]
                                                    Y_nac_sym1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[18]
                                                    Y1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[19]
                                                    MFW1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[20]
                                                    X_debut_cyl1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[21]
                                                    X_fin_cyl1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[22]
                                                    X_av1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[23]
                                                    Y_av1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[24]
                                                    X_ar1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[25]
                                                    Y_ar1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[26]
                                                    X_ev1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[27]
                                                    Y_ev1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[28]
                                                    xa1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[29]
                                                    ya1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[30]
                                                    SV1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[31]
                                                    B1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[32]
                                                    Y21 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[33]
                                                    extremite1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[34]
                                                    emplanture1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[35]
                                                    SPF1 = self.Coordonee(
                                                        sref, rayonextern1, lcyl1, lav1, lar1, lfus1)[36]

                                                    Cxtot1 = self.Aero_GV(
                                                        Cassure1, L01, L21, rayonextern1, sref, LOH1, SH1, lfus1, lav1, lar1, lcyl1, effilement1, Y41)[0]
                                                    finesse_croisiere1 = self.Aero_GV(
                                                        Cassure1, L01, L21, rayonextern1, sref, LOH1, SH1, lfus1, lav1, lar1, lcyl1, effilement1, Y41)[1]
                                                    CX01 = self.Aero_GV(
                                                        Cassure1, L01, L21, rayonextern1, sref, LOH1, SH1, lfus1, lav1, lar1, lcyl1, effilement1, Y41)[2]
                                                    CXi1 = self.Aero_GV(
                                                        Cassure1, L01, L21, rayonextern1, sref, LOH1, SH1, lfus1, lav1, lar1, lcyl1, effilement1, Y41)[3]
                                                    SMF1 = self.Aero_GV(
                                                        Cassure1, L01, L21, rayonextern1, sref, LOH1, SH1, lfus1, lav1, lar1, lcyl1, effilement1, Y41)[4]
                                                    Smtot1 = self.Aero_GV(
                                                        Cassure1, L01, L21, rayonextern1, sref, LOH1, SH1, lfus1, lav1, lar1, lcyl1, effilement1, Y41)[5]
                                                    Cz1 = self.Aero_BV(CX01)[1]
                                                    finesse_D_sur_L1 = self.Aero_BV(CX01)[
                                                        0]

                                                if rapport < 10 and rapport > maxi:

                                                    maxi = rapport
                                                    rang_eco2 = rang_eco
                                                    rang_affaire2 = rang_affaire
                                                    rang_er2 = rang_er
                                                    nfront2 = nfront
                                                    nfront_rang_reduit2 = nfront-i_ar
                                                    rang_er_AV2 = places_er_AV / \
                                                        (nfront2-2)
                                                    rang_affaire_AV2 = places_affaire_AV / \
                                                        (nfront2-1)
                                                    rang_eco_AV2 = places_eco_AV / \
                                                        (nfront2-1)

                                                    places_eco_tot2 = places_eco_tot
                                                    places_affaire_tot2 = places_affaire_tot
                                                    places_er_tot2 = places_er_tot

                                                    rang_er_cycl2 = places_cyl_er / \
                                                        (nfront2-2)
                                                    rang_affaire_cycl2 = places_cyl_affaire / \
                                                        (nfront2-1)
                                                    rang_eco_cycl2 = places_cyl_eco/nfront2

                                                    rang_eco_nonreduit2 = (
                                                        nb_rang_AR)
                                                    rang_eco_reduit2 = Rang_redui_AR

                                                    rang_totale2 = rang_totale
                                                    nb_rang_AV2 = nb_rang_AV
                                                    nb_rang_cylindre2 = nb_rang_cylindre
                                                    nb_rang_AR2 = nb_rang_AR
                                                    rayonextern2 = rext
                                                    lcyl2 = longeur
                                                    lav2 = self.rayon(
                                                        nfront2)[2][0]
                                                    lar2 = self.rayon(
                                                        nfront2)[3][0]
                                                    lfus2 = self.L_sur_D(
                                                        longeur, larav, rext)[1]
                                                    Cassure2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[0]
                                                    L02 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[1]
                                                    L22 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[2]
                                                    LOH2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[3]
                                                    SH2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[4]
                                                    effilement2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[5]
                                                    Y42 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[6]
                                                    Y_debut_cyl2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[7]
                                                    X_fuselage2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[8]
                                                    Y_fuselage2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[9]
                                                    Y_fuselage_sym2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[10]
                                                    X2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[11]
                                                    Y_sym2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[12]
                                                    X_EH2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[13]
                                                    Y_EH2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[14]
                                                    Y_EH_sym2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[15]
                                                    X_nac2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[16]
                                                    Y_nac2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[17]
                                                    Y_nac_sym2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[18]
                                                    Y2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[19]
                                                    MFW2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[20]
                                                    X_debut_cyl2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[21]
                                                    X_fin_cyl2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[22]
                                                    X_av2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[23]
                                                    Y_av2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[24]
                                                    X_ar2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[25]
                                                    Y_ar2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[26]

                                                    X_ev2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[27]
                                                    Y_ev2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[28]
                                                    xa2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[29]
                                                    ya2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[30]
                                                    SV2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[31]
                                                    B2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[32]
                                                    Y22 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[33]
                                                    extremite2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[34]
                                                    emplanture2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[35]
                                                    SPF2 = self.Coordonee(
                                                        sref, rayonextern2, lcyl2, lav2, lar2, lfus2)[36]

                                                    Cxtot2 = self.Aero_GV(
                                                        Cassure2, L02, L22, rayonextern2, sref, LOH2, SH2, lfus2, lav2, lar2, lcyl2, effilement2, Y42)[0]
                                                    finesse_croisiere2 = self.Aero_GV(
                                                        Cassure2, L02, L22, rayonextern2, sref, LOH2, SH2, lfus2, lav2, lar2, lcyl2, effilement2, Y42)[1]
                                                    CX02 = self.Aero_GV(
                                                        Cassure2, L02, L22, rayonextern2, sref, LOH2, SH2, lfus2, lav2, lar2, lcyl2, effilement2, Y42)[2]
                                                    CXi2 = self.Aero_GV(
                                                        Cassure2, L02, L22, rayonextern2, sref, LOH2, SH2, lfus2, lav2, lar2, lcyl2, effilement2, Y42)[3]
                                                    SMF2 = self.Aero_GV(
                                                        Cassure2, L02, L22, rayonextern2, sref, LOH2, SH2, lfus2, lav2, lar2, lcyl2, effilement2, Y42)[4]
                                                    Smtot2 = self.Aero_GV(
                                                        Cassure2, L02, L22, rayonextern2, sref, LOH2, SH2, lfus2, lav2, lar2, lcyl2, effilement2, Y42)[5]
                                                    Cz2 = self.Aero_BV(CX02)[1]
                                                    finesse_D_sur_L2 = self.Aero_BV(CX02)[
                                                        0]

                            i_ar += i_ar

                        y = y+1
                    x = x+1
                    lx = lx+[nb_places_tot_affaire-x]

            while state == False:
                if maxi!=0:
                    V_croisiere2 = self.perfo(MTOW, CX02, CXi2)[0]
                    finesse2 = self.perfo(MTOW, CX02, CXi2)[1]
                    Mcarb2 = self.perfo(MTOW, CX02, CXi2)[2]
                    Czopt2 = self.perfo(MTOW, CX02, CXi2)[3]

                    F2 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                    extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[0]
                    A2 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                    extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[1]
                    B_masse2 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                          extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[2]
                    C2 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                    extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[3]
                    D2 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                    extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[4]
                    E2 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                    extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[5]
                    G2 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                    extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[6]
                    FN02 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                      extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[7]
                    B12 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                     extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[8]
                    nourriture2 = self.Masse(MTOW, MZFW, MFW2, nfront2, lcyl2, SMF2, rayonextern2, SV2, SH2, B2, lfus2, Y42, Y22, Smtot2, L22,
                                             extremite2, Cassure2, emplanture2, sref, SPF2, finesse_D_sur_L2, finesse_croisiere2, Lpax)[9]

                    if (0 <= (Mcarb2-F2)/(F2)*100 <= 0.5):
                        state = True

                    else:
                        # print((Mcarb+F))
                        MTOW = MTOW+100
                        MZFW = MTOW-F2
                        sref2 = ((MTOW-F2/2)*9.81)/((1.4)/2 *
                                                    self.P*math.pow(self.mach, 2)*Czopt2)
                    sref = self.sref
                else:
                    
                    V_croisiere1 = self.perfo(MTOW, CX01, CXi1)[0]
                    finesse1 = self.perfo(MTOW, CX01, CXi1)[1]
                    Mcarb1 = self.perfo(MTOW, CX01, CXi1)[2]
                    Czopt1 = self.perfo(MTOW, CX01, CXi1)[3]

                    F1 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                    extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[0]
                    A1 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                    extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[1]
                    B_masse1 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                          extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[2]
                    C1 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                    extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[3]
                    D1 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                    extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[4]
                    E1 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                    extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[5]
                    G1 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                    extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[6]
                    FN01 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                      extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[7]
                    B11 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                     extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[8]
                    nourriture1 = self.Masse(MTOW, MZFW, MFW1, nfront1, lcyl1, SMF1, rayonextern1, SV1, SH1, B1, lfus1, Y41, Y21, Smtot1, L21,
                                             extremite1, Cassure1, emplanture1, sref, SPF1, finesse_D_sur_L1, finesse_croisiere1, Lpax)[9]

                    if (0 <= (Mcarb1-F1)/(F1)*100 <= 0.5):
                        state = True

                    else:
                        # print((Mcarb+F))
                        MTOW = MTOW+100
                        MZFW = MTOW-F1
                        sref1 = ((MTOW-F1/2)*9.81)/((1.4)/2 *
                                                    self.P*math.pow(self.mach, 2)*Czopt1)
                    sref = self.Sref()
            i = i+1

        if maxi != 0:

            best_combinaison = ["L_surD: ", maxi, "npfront: ", nfront2, "Rang tot: ", rang_totale2, "rang AV: ", nb_rang_AV2, "rang cyl: ", nb_rang_cylindre2, "rang Ar:", nb_rang_AR2, 'rang er av', rang_er_AV2, 'rang affaire av', rang_affaire_AV2, 'rang eco av', rang_eco_AV2, 'rang cyl er',
                                rang_er_cycl2, 'rang cyl aff', rang_affaire_cycl2, 'rang cyl eco', rang_eco_cycl2, 'rang eco non reduit', rang_eco_nonreduit2, 'rang eco reduit', rang_eco_reduit2, '----', places_eco_tot2, places_affaire_tot2, places_er_tot2, '-----', rang_eco2, rang_affaire2, rang_er2]
            self.draw(nfront2, nfront_rang_reduit2, rang_er_AV2, rang_affaire_AV2, rang_er_cycl2, rang_affaire_cycl2, rang_eco_cycl2, rang_eco_nonreduit2,
                      rang_eco_reduit2, X_debut_cyl2, X_fin_cyl2, Y_debut_cyl2, X_fuselage2, Y_fuselage2, Y_fuselage_sym2, X2, Y22, Y_sym2, X_EH2, Y_EH2, Y_EH_sym2, X_nac2, Y_nac2, Y_nac_sym2, X_av2, Y_av2, X_ar2, Y_ar2, X_ev2, Y_ev2, xa2, ya2)

        else:

            best_combinaison = ["L_surD: ", mini1, "npfront: ", nfront1, "Rang tot: ", rang_totale1, "rang AV: ", nb_rang_AV1, "rang cyl: ", nb_rang_cylindre1, "rang Ar:", nb_rang_AR1, 'rang er av', rang_er_AV1, 'rang affaire av', rang_affaire_AV1, 'rang eco av', rang_eco_AV1, 'rang cyl er',
                                rang_er_cycl1, 'rang cyl aff', rang_affaire_cycl1, 'rang cyl eco', rang_eco_cycl1, 'rang eco non reduit', rang_eco_nonreduit1, 'rang eco reduit', rang_eco_reduit1, '----', places_eco_tot1, places_affaire_tot1, places_er_tot1, '-----', rang_eco1, rang_affaire1, rang_er1]
            self.draw(nfront1, nfront_rang_reduit1, rang_er_AV1, rang_affaire_AV1, rang_er_cycl1, rang_affaire_cycl1, rang_eco_cycl1, rang_eco_nonreduit1,
                      rang_eco_reduit1, X_debut_cyl1, X_fin_cyl1, Y_debut_cyl1, X_fuselage1, Y_fuselage1, Y_fuselage_sym1, X1, Y1, Y_sym1, X_EH1, Y_EH1, Y_EH_sym1, X_nac1, Y_nac1, Y_nac_sym1, X_av1, Y_av1, X_ar1, Y_ar1, X_ev1, Y_ev1, xa1, ya1)

        self.diagrame_tapis()
