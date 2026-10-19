# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 08:51:29 2026

@author: cabra
"""


import numpy as np
from solve_gases_CHONS import solve_gases_CHONS
from solve_gases_CHONS_extrusive import solve_gases_CHONS_ext
import Molecular_weights as MW


'__________General parameters__________'


'Gas Behavior'
#set to 'ideal' or 'real'
gas_behavior= 'ideal'

'Partitioning'
#Set Partitioning to True if you want MORDOR to calculate volatile partitioning from rock to melt, set to false if you want to define volatiles in the melt directly
Partitioning = False

'the following parameters are only needed if partitioning is set to True'
#set melting style to bacth melting or fractional melting
melting_style = 'batch_melting'
# F_melt - fraction of the rock that is molten
F_melt = 0.25


'__________Extrusive vs Intrusive____________'

# Degassing style: Intrusive or Extrusive, 
#Intrusive: - solubility is defined by vapor phase in the melt which is calculated in the code, 
#Extrusive: solubility is defined by partial pressure of gases in the atmosphere
Degassing_style = 'Intrusive'

'If Degassing_style is set to extrusive, define partial pressures of volatiles in the atmopshere (in bar)'
PH2O = 0
PH2 = 0
PCO2 = 0
PCO = 0
PCH4 = 0
PH2S = 0
PSO2 = 0
PS2 = 0
PNH3 = 0
PN2 = 0



'__________Melt Composition____________'

'Define the melt composition in weight fraction. Some examples are given below'
#basalt
mSiO2 = 0.483
mTiO2 = 0.013
mAl2O3 = 0.179
mFeO = 0.088
mMgO = 0.081
mCaO = 0.103
mNa2O = 0.024
mK2O = 0.002
mMnO = 0
        
#alkali basalt
#mSiO2 = 0.479
# mTiO2 = 0.028
# mAl2O3 = 0.18
# mFeO = 0.096
# mMgO = 0.033
# mCaO = 0.077
# mNa2O = 0.06
# mK2O = 0.027
        
#rhyolite
# mSiO2 = 0.714
# mTiO2 = 0.0025
# mAl2O3 = 0.1238
# mFeO = 0.0459
# mMgO = 0.0016
# mCaO = 0.0077
# mNa2O = 0.0446
# mK2O = 0.0428



'Define the volatile content in weight fraction '
'if Partitioning is set to True, define volatile content of rock'
'if Partitioning is set to False, define volatile content of the melt'

'H2O'
mH2O = 1000E-6
'CO2 (if Partitioning is set to TRUE set here the maximum Co2 content of the melt. The actual content will be calculated depending on fO2'
mCO2 = 500E-6
'S' 
mS= 1000E-6
'N'
mN = 100E-6

x = mSiO2/MW.SiO2 + mTiO2/MW.TiO2 + mAl2O3/MW.Al2O3 + mFeO/MW.FeO + mMgO/MW.MgO + mCaO/MW.CaO \
        + mNa2O/MW.Na2O + mK2O/MW.K2O + mH2O/MW.H2O + mCO2/MW.CO2 +mS/MW.S + mN/MW.N


'__________Melt Properties____________'

'Temperature of the melt in Kelvin'
T = 1675

'Pressure in the melt in bar'
P = 100

#If Degassing sytle is extrusive, the pressure is defined by the pressure in the atmopshere
if Degassing_style == 'Extrusive':
    P = PH2 + PH2O + PCO + PCO2 + PCH4 + PS2 + PSO2 + PH2S + PN2 + PNH3
    
# Overall P cannot be 0, for the system to work, hence for a non existing atmopshere, P is set to 1E-10 bar
if P == 0:
    P = 1E-10


'''
f_O2 = oxygen fugacity of the melt
    possible inputs:
        'IW'    - Iron Wustite
        'QFM'   - Quartz Fayalite Magnetite
        'NiNiO' - Nickel Nickel Oxide
buffer = if you want to change the oxygen fugacity 
    e.g. for IW - 2 put in: -2'''
    
f_O2 = 'IW'
buffer = 3

#buffer = 8




'''Sol = solubility law that you want to use, 
    possible inputs:
        Sol_H2O:
            IM - Iacono Marziano, 2012
            PW_basalt - Parfit and Wilson, 2008 (basaltic)
            PW_rhyolite - Parfit and Wilson, 2008 (rhyolitic)
            LI - Lichtenberg, 2021
            BO - Bower, 2019
            DI - Dixon 1995, 1997
            BRA_rhyolite - fit from Brachmann 2026 (for rhyolite)
            BRA_basalt - fit from Brachmann 2026 (for basalt)
            
        Sol_CO2:
            IM - Iacono Marziano, 2012
            PW_basalt - Parfit and Wilson, 2008 (basaltic)
            PW_rhyolite - Parfit and Wilson, 2008 (rhyolitic)
            LI - Lichtenberg, 2021
            BO - Bower, 2019
            DI - Dixon 1995, 1997
            BRA_rhyolite -fit from Brachmann 2026 (for rhyolite)
            BRA_basalt - fit from Brachmann 2026 (for basalt)
            
        Sol_S2:
            GAI - Gaillard, 2026
            BOU - Boulling and Wood, 2023
            ONE - O'neill 2021
        
        Sol_S6:
            SCA - Scaillet and Pichavant 2009
            BOU - Boulling and Wood 2022
            ONE - O'neill and Mavrogenes 2022
            
        Sol_N
            BER - Bernadou 2021
            LIB_1 - Libourel 2003
            LIB_2 - Libourel 2025 (no NBO)
            LIB_3 - Libourel 2025 (with NBO)
            BRA   - fit from Brachmann 2026 '''

Sol_H2O = 'BRA_basalt'
Sol_CO2 = 'BRA_basalt'
Sol_S2 = 'GAI'
Sol_S6 = 'BOU'
Sol_N = 'BRA'


'_______________IMPORTANT NOTE_________________'
'''In the intrusive version of the code, the gases are normalized to 1, which is necessary to calcalulate the solubility, here alphaG gas denotes the fraction of total volatile moles residing in the gas phase,
 in the extrusive version, solubility is defined by the partial pressures of the gases in the atmopshere, hence no alphaG is needed. Here the gases produced are given in molar fraction of the entire gases'''

if Degassing_style == 'Intrusive':

    XH2, XH2O, XCH4, XCO, XCO2, XS2, XSO2, XH2S, XN2, XNH3,  XH2O_melt, XCO2_melt, XS_melt, XN_melt,  alphaG = \
            solve_gases_CHONS(P, T, f_O2, buffer, mH2O, mCO2, mS, mN, Sol_H2O ,Sol_CO2, Sol_S2, Sol_S6, Sol_N, \
                              mSiO2, mTiO2, mAl2O3, mFeO, mMgO, mCaO, mNa2O, mK2O, mMnO, Partitioning, melting_style, F_melt, gas_behavior)

elif Degassing_style == 'Extrusive':
    
    XH2, XH2O, XCH4, XCO, XCO2, XS2, XSO2, XH2S, XN2, XNH3,  XH2O_melt, XCO2_melt, XS_melt, XN_melt = \
            solve_gases_CHONS_ext(P, T, f_O2, buffer, mH2O, mCO2, mS, mN, Sol_H2O ,Sol_CO2, Sol_S2, Sol_S6, Sol_N, \
                                  mSiO2, mTiO2, mAl2O3, mFeO, mMgO, mCaO, mNa2O, mK2O, mMnO, Partitioning, melting_style, \
                                  F_melt, gas_behavior, PH2, PH2O, PCO, PCO2, PCH4, PS2, PSO2, PH2S, PN2)