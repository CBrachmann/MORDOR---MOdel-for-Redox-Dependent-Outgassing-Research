# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 09:52:10 2025

@author: cabra
"""
'''This is my Version of the Shi and Saxena fugacity coeffcieint calculations. I had to leave many parameters out becuase the results did not make sense. 
So these fugacity coefficients can only be used as an approximation. The original parameters are stored in FugacityC_ALL_XX'''

import numpy as np
from matplotlib import pyplot as plt

def FNF(V,TK,A,B,P):
    R = 83.14321
    return R * TK / (V - B) - A / ((V * V + B * V) * TK**0.5) - P

def Fugacity_C(P,TK, X_1, species): 
    if P > 5000:
        print('ERROR: not recommended for P>5000') 
    
    if species == 'H2O':
        R = 83.14321
        FNA = (166800000 - 193080 * (TK - 273.15) + 186.4 * (TK - 273.15)**2 - 0.071288 * ((TK - 273.15)**3)) * 1.01325
        FNB = 1.01325 * (73030000 - 71400 * (TK - 273.15) + 21.57 * (TK - 273.15)**2)
        FNC = 1.01325 * (np.exp(-11.071 + 5953 / TK - 2746000 / TK**2 + 464600000 / TK**3) * 0.5 * R * R * TK**2.5 / 1.02668 + 40123800)
        #FNC wahrscheinlich K von H2O + CO2 = H2CO3
        B_1 = 14.6
        B_2 = 29.7
        
        B = X_1 * B_1 + (1 - X_1) * B_2
        A = X_1**2 * FNA + 2 * X_1 * (1 - X_1) * FNC + (1 - X_1)**2 * FNB
        Temp2 = B + 5
        Q = 1
        Temp1 = 0
        while abs(Temp2 - Temp1) >= 0.00001:
            Temp1 = Temp2
            F_1 = (FNF(Temp1 + 0.01, TK, A, B, P) - FNF(Temp1, TK, A, B, P)) / 0.01
            Temp2 = Temp1 - Q * FNF(Temp1, TK, A, B, P) / F_1
            F_2 = (FNF(Temp2 + 0.01, TK, A, B, P) - FNF(Temp2, TK, A, B, P)) / 0.01
            if F_2 * F_1 <= 0:
                Q = Q / 2.
            if abs(Temp2 - Temp1) > 0.00001:
                F_1 = F_2
        V = Temp2
        
        FC_H2O = np.log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * FNA + (1 - X_1) * FNC) * np.log((V + B) / V) / (R * TK**1.5 * B)
        FC_H2O = FC_H2O + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
        FC_H2O = np.exp(FC_H2O)
        FC_CO2 = np.log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * FNC + (1 - X_1) * FNB) * np.log((V + B) / V) / (R * TK**1.5 * B)
        FC_CO2 = FC_CO2 + (np.log((V + B) / V) - B / (V + B)) * A * B_2 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
        FC_CO2 = np.exp(FC_CO2)
        if species == 'H2O':
            FC =  FC_H2O
        # if species == 'CO2':
        #     FC = FC_CO2
        # if species == 'H2O':
        #     Tcr = 647.25
        #     Pcr = 221.1925
    if species in ['CO2', 'CH4', 'CO', 'O2', 'H2', 'S2', 'SO2', 'COS', 'H2S']:
        if species == 'CO2':
            Tcr = 304.15
            Pcr = 73.8659
        if species == 'CH4':
            Tcr = 191.05
            Pcr = 46.4069
        elif species == 'CO':
            Tcr = 133.15
            Pcr = 34.9571
        elif species == 'O2':
            Tcr = 154.75
            Pcr = 50.7638
        elif species == 'H2':
            Tcr = 33.25
            Pcr = 12.9696
        elif species == 'S2':
            Tcr = 208.15
            Pcr = 72.954
        elif species == 'SO2':
            Tcr = 430.95
            Pcr = 78.7295
        elif species == 'COS':
            Tcr = 377.55
            Pcr = 65.8612
        elif species == 'H2S':
            Tcr = 373.55
            Pcr = 90.0779      
        
        if species in ['O2','CO2', 'CO', 'CH4', 'S2', 'COS' ]:
    
            #if P <= 5000: 
                Q1_A = 1
                Q2_A = 0
                Q3_A = 0
                Q4_A = 0
                Q5_A = -0.5917
                Q6_A = 0
                Q7_A = 0
                Q8_A = 0
                Q1_B = 0
                Q2_B = 0
                Q3_B = 0.09122
                Q4_B = 0
                Q5_B = 0
                Q6_B = 0
                Q7_B = 0
                Q8_B = 0
                Q1_C = 0
                Q2_C = 0
                Q3_C = 0
                Q4_C = 0
                Q5_C = -0.0001416
                Q6_C = 0
                Q7_C = 0
                Q8_C = -2.835e-06
                Q1_D = 0
                Q2_D = 0
                Q3_D = 0
                Q4_D = 0
                Q5_D = 0
                Q6_D = 0
                Q7_D = 0
                Q8_D = 0
            # else:
            #     Q1_A = 2.0614
            #     Q2_A = 0
            #     Q3_A = 0
            #     Q4_A = 0
            #     Q5_A = -2.235
            #     Q6_A = 0
            #     Q7_A = 0
            #     Q8_A = -0.3941
            #     Q1_B = 0
            #     Q2_B = 0
            #     Q3_B = 0.05513
            #     Q4_B = 0
            #     Q5_B = 0.03934
            #     Q6_B = 0
            #     Q7_B = 0
            #     Q8_B = 0
            #     Q1_C = 0
            #     Q2_C = 0
            #     Q3_C = -1.8935e-06
            #     Q4_C = 0
            #     Q5_C = -1.1092e-05
            #     Q6_C = 0
            #     Q7_C = -2.189e-05
            #     Q8_C = 0
            #     Q1_D = 0
            #     Q2_D = 0
            #     Q3_D = 5.053E-11
            #     Q4_D = 0
            #     Q5_D = 0
            #     Q6_D = -6.303e-21
            #     Q7_D = 0
            #     Q8_D = 0
                
        elif species == 'H2':
            #if P < 1000:
                Q1_A = 1
                Q2_A = 0
                Q3_A = 0
                Q4_A = 0
                Q5_A = 0
                Q1_B = 0
                Q2_B = 0.09827
                Q3_B = 0
                Q4_B = -0.2709
                Q5_B = 0
                Q1_C = 0
                Q2_C = 0
                Q3_C = -0.00103
                Q4_C = 0
                Q5_C = 0.01427
                Q1_D = 0
                Q2_D = 0
                Q3_D = 0
                Q4_D = 0
                Q5_D = 0
            # else: 
            #     Q1_A = 2.2615
            #     Q2_A = 0
            #     Q3_A = -68.712
            #     Q4_A = 0
            #     Q5_A = -1.0573e+04
            #     Q6_A = 0
            #     Q7_A = 0
            #     Q8_A = -0.16936
            #     Q1_B = -2.6707e-04
            #     Q2_B = 0
            #     Q3_B = 0.20173
            #     Q4_B = 0
            #     Q5_B = 4.5759
            #     Q6_B = 0
            #     Q7_B = 0
            #     Q8_B = 3.1452e-05
            #     Q1_C = -2.3376E-09
            #     Q2_C = 0
            #     Q3_C = 3.4091e-07
            #     Q4_C = 0
            #     Q5_C = 1.4188E-03
            #     Q6_C = 0
            #     Q7_C = 0
            #     Q8_C = 3.0117e-10
            #     Q1_D = 3.2606E-15
            #     Q2_D = 0
            #     Q3_D = 2.4402E-12
            #     Q4_D = 0
            #     Q5_D = -2.4027E-09
            #     Q6_D = 0
            #     Q7_D = 0
            #     Q8_D = 0
        elif species == 'SO2':
                Q1_A = 0.92854
                Q2_A = 0.043269
                Q3_A = -0.24671
                Q4_A = 0
                Q5_A = 0.24999
                Q6_A = 0
                Q7_A = -0.53182
                Q8_A = -0.016461
                Q1_B = 0.84866e-03
                Q2_B = -0.0018379
                Q3_B = 0.066787
                Q4_B = 0
                Q5_B = -0.029427
                Q6_B = 0
                Q7_B = 0.029003
                Q8_B = 0.0054808
                Q1_C = -0.35456e-03
                Q2_C = 0.23316e-04
                Q3_C = 0.94159e-03
                Q4_C = 0
                Q5_C = -0.81653e-03
                Q6_C = 0
                Q7_C = 0.23154e-03
                Q8_C = 0.55542e-04
                Q1_D = 0
                Q2_D = 0
                Q3_D = 0
                Q4_D = 0
                Q5_D = 0
                Q6_D = 0
                Q7_D = 0
                Q8_D = 0
        elif species == 'H2S':
            if P < 500:
                Q1_A = 1.4721
                Q2_A = 1.1177
                Q3_A = 3.9657
                Q4_A = 0
                Q5_A = -10.028
                Q6_A = 0
                Q7_A = 4.5484
                Q8_A = -3.82
                Q1_B = 0.16066
                Q2_B = 0.10887
                Q3_B = 0.29014
                Q4_B = 0
                Q5_B = -0.99593
                Q6_B = 0
                Q7_B = -0.18627
                Q8_B = -0.45515
                Q1_C = -0.28933
                Q2_C = -0.070522
                Q3_C = 0.39828
                Q4_C = 0
                Q5_C = -0.050533
                Q6_C = 0
                Q7_C = 0.1176
                Q8_C = 0.33972
                Q1_D = 0
                Q2_D = 0
                Q3_D = 0
                Q4_D = 0
                Q5_D = 0
                Q6_D = 0
                Q7_D = 0
                Q8_D = 0    
            else:
                Q1_A = 0.59941
                Q2_A = -0.001557
                Q3_A = 0.04525
                Q4_A = 0
                Q5_A = 0.36687
                Q6_A = 0
                Q7_A = -0.79248
                Q8_A = 0.26058
                Q1_B = 0.022545
                Q2_B = 0.0017473
                Q3_B = 0.048253
                Q4_B = 0
                Q5_B = -0.01989
                Q6_B = 0
                Q7_B = 0.032794
                Q8_B = -0.010985
                Q1_C = 0.57375e-03
                Q2_C = -0.20944e-05
                Q3_C = -0.0011894
                Q4_C = 0
                Q5_C = 0.0014661
                Q6_C = 0
                Q7_C = -0.75605e-03
                Q8_C = -0.27985e-03
                Q1_D = 0
                Q2_D = 0
                Q3_D = 0
                Q4_D = 0
                Q5_D = 0
                Q6_D = 0
                Q7_D = 0
                Q8_D = 0  
        Tr = TK/Tcr
        Pr = P/Pcr
        P0r = 1/Pcr
        
        if species in ['O2','CO2', 'CO', 'CH4', 'S2', 'COS', 'H2' , 'H2O']:    
            #if P >= 1000:
                # A = Q1_A + Q2_A*Tr + Q3_A*Tr**-1 + Q4_A*Tr**2 + Q5_A*Tr**-2 + Q6_A*Tr**3 + Q7_A*Tr**-3 + Q8_A*np.log(Tr)
                # B = Q1_B + Q2_B*Tr + Q3_B*Tr**-1 + Q4_B*Tr**2 + Q5_B*Tr**-2 + Q6_B*Tr**3 + Q7_B*Tr**-3 + Q8_B*np.log(Tr)
                # C = Q1_C + Q2_C*Tr + Q3_C*Tr**-1 + Q4_C*Tr**2 + Q5_C*Tr**-2 + Q6_C*Tr**3 + Q7_C*Tr**-3 + Q8_C*np.log(Tr)
                # D = Q1_D + Q2_D*Tr + Q3_D*Tr**-1 + Q4_D*Tr**2 + Q5_D*Tr**-2 + Q6_D*Tr**3 + Q7_D*Tr**-3 + Q8_D*np.log(Tr)
            # elif P < 1000:
                A = Q1_A + Q2_A*Tr**-1 + Q3_A*Tr**-1.5 + Q4_A*Tr**-3 + Q5_A*Tr**-4
                B = Q1_B + Q2_B*Tr**-1 + Q3_B*Tr**-1.5 + Q4_B*Tr**-3 + Q5_B*Tr**-4
                C = Q1_C + Q2_C*Tr**-1 + Q3_C*Tr**-1.5 + Q4_C*Tr**-3 + Q5_C*Tr**-4
                D = Q1_D + Q2_D*Tr**-1 + Q3_D*Tr**-1.5 + Q4_D*Tr**-3 + Q5_D*Tr**-4
        if species in ['SO2', 'H2S']:
                A = Q1_A + Q2_A*Tr + Q3_A*Tr**-1 + Q4_A*Tr**2 + Q5_A*Tr**-2 + Q6_A*Tr**3 + Q7_A*Tr**-3 + Q8_A*np.log(Tr)
                B = Q1_B + Q2_B*Tr + Q3_B*Tr**-1 + Q4_B*Tr**2 + Q5_B*Tr**-2 + Q6_B*Tr**3 + Q7_B*Tr**-3 + Q8_B*np.log(Tr)
                C = Q1_C + Q2_C*Tr + Q3_C*Tr**-1 + Q4_C*Tr**2 + Q5_C*Tr**-2 + Q6_C*Tr**3 + Q7_C*Tr**-3 + Q8_C*np.log(Tr)
                D = Q1_D + Q2_D*Tr + Q3_D*Tr**-1 + Q4_D*Tr**2 + Q5_D*Tr**-2 + Q6_D*Tr**3 + Q7_D*Tr**-3 + Q8_D*np.log(Tr)
        
        #Z = A + B*Pr +C*Pr**2 +D*Pr**3
        INT_Z = A*np.log(Pr/P0r) + B*(Pr- P0r) + (C/2)*(Pr**2- P0r**2) + (D/3)*(Pr**3 - P0r**3)
        #R = 83.14
        #V = Z*R*TK/P
        #FC = np.exp(((V*P-V) - R*TK*np.log(P) + R*TK*np.log(1))/(R*TK))
        FC = np.exp(INT_Z - np.log(P))        
    return(FC)
    
