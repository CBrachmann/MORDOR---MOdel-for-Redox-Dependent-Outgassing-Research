# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 12:15:07 2025

@author: cabra
"""


import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import sys
sys.path.insert(1, 'C:/Users/cabra/Nextcloud/PHD/my_codes/Gibbs_energies')
from Gibbs_E import Gibbs_E


def Fe_FeO_Fe2O3(T,P, fugacity, buffer, XFeO_init):

    
    R = 8.31447215
    
    Gf_FeO = - 244118 + 115.56*T - 8.474*T*np.log(T) 
    #Gf_FeO = Gibbs_E(T, 'FeO')*10**3
    Gf_Fe2O3 = Gibbs_E(T, 'Fe2O3')*10**3
    
    dG_10 = Gf_Fe2O3 - 2*Gf_FeO
    

    
    K9 = np.exp(-Gf_FeO/(R*T))
    K10 = np.exp(-dG_10/(R*T))
                           #molar mass of CO2

    
    if fugacity == 'QIF':
        f_O2 = 10**(7.679 - 29637/T + 0.05*(P-1)/T +buffer)     #QIF
    elif fugacity == 'IW':
        f_O2 = 10**(6.899 - 27714/T + 0.05*(P-1)/T +buffer)     #IW
    elif fugacity == 'QFM':
        f_O2 = 10**(5.5976 - 24505/T + 0.8099*np.log10(T) + 0.0937*(P-1)/T +buffer)    #QFM
    elif fugacity == 'NiNiO':
        f_O2 = 10**(8.951 - 24556/T + 0.46*(P-1)/T +buffer)     #NiNiO
   
    
    XFe_init = 1    
    
    def sys1(y):
        XFe, XFeO, XFe2O3 = y
        return(
               K9 * XFe*(f_O2)**0.5 - XFeO,
               K10 * XFeO**2 *(f_O2)**0.5 - XFe2O3,
               XFe + XFeO + 2 * XFe2O3 - XFe_init)
       

    y_init = (XFe_init*0.5, XFe_init*0.4, XFe_init *0.1)
    #sol = optimize.least_squares(sys1, y_init, bounds=(0, np.inf))
    sol = optimize.root(sys1, y_init,options={'maxiter': 100000}, method='lm') 
    XFe, XFeO, XFe2O3 = sol['x']                       

    xFe = XFe*XFeO_init
    xFeO = XFeO*XFeO_init
    xFe2O3 = XFe2O3*XFeO_init
    aFe2O3 = XFe2O3*XFeO_init**2

    return(XFe, XFeO, XFe2O3, xFe, xFeO, xFe2O3, aFe2O3)

'''XFe, XFeO, XFe2O3 - values normalized so that xFe +xFeO +2*xFe2O3 = 1 is true,
XFe, XXFeo, XFe2O3 - values so that XFe +XFeO +2*XFe2O3 = XFeO_init is true,
aFe2O3 activity of Fe2O3 so that fO2 = aFe2O3/(XFeO**2*K10))**2 is correct
'''

