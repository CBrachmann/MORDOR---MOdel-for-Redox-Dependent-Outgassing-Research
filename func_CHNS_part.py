#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 11:16:46 2023

@author: caroline
"""

import numpy as np



def C_H_N_S_partitioning(F_melt, m_H2O_mantle, m_C_mantle, m_S_mantle, m_N_mantle, T, P, fugacity, buffer, x, melting_style):
    fwm = 36.594
    


    
    
    'oxygen fugacity calculation'
    if fugacity == 'QIF':
        f_O2 = 10**(7.679 - 29637/T + 0.05*(P-1)/T +buffer)     #QIF
    elif fugacity == 'IW':
        f_O2 = 10**(6.899 - 27714/T + 0.05*(P-1)/T +buffer)     #IW
    elif fugacity == 'QFM':
        f_O2 = 10**(5.5976 - 24505/T + 0.8099*np.log10(T) + 0.0937*(P-1)/T +buffer)  #QFM
    elif fugacity == 'NiNiO':
        f_O2 = 10**(12.9850 - 25026/T -1.1786*np.log10(T) + 0.0458*(P-1)/T + buffer)    
    elif fugacity == 'MH':
        f_O2 = 10**(-25.7139 - 19375/T +11.4191*np.log10(T) + 0.0185*(P-1)/T +buffer)
    elif fugacity == 'CUO':
        f_O2 = 10**(-19.6938 - 15266/T -2.7902*np.log10(T) + 0.0105*(P-1)/T +buffer)   
        
    
    'C partitioning after Holloway et al 1992'    
    K1 = 10**(40.07639-((2.53932*10.0**(-2))*T)+(5.27096*10.0**(-6))*(T**2)+0.0267*(P-1.0)/T)
    K2 = 10**(-6.24763-(282.56/T)-(0.119242*(P-1000.0)/T))
    
    m_carbonate_melt = (K1*K2*f_O2)/(1.0+(K1*K2*f_O2))
    
    m_C_melt_Hol = (44.01/fwm)*m_carbonate_melt / ( 1.0-(1.0-(44.01/fwm))*m_carbonate_melt )
    C_ppm_1 = m_C_melt_Hol*10**6
    


    D_C = 5.5E-4
    if melting_style == 'batch_melting':   
        m_C_melt_Ros = m_C_mantle/(D_C + F_melt* (1- D_C))
    else:
        m_C_melt_Ros = m_C_mantle/F_melt*(1-(1-F_melt)**(1/D_C))
  
    #m_C_melt = min(C_ppm_2*10**-6, m_C_melt_Ros)
    
    m_C_melt = min(C_ppm_1*10**-6, m_C_melt_Ros)
    
    #calculate Co2 from C
    XC_melt = m_C_melt/12.011/x
    XCO2_melt = XC_melt
    m_CO2_melt = XCO2_melt*44.01*x
    

        
    'N partitioning interpolated after Pal and Dasgupta 2024'
    N_part_data = [[10**(6.899 - 27714/T + 0.05*(P-1)/T -2),0.136], [10**(6.899 - 27714/T + 0.05*(P-1)/T),0.0065], [10**(6.899 - 27714/T + 0.05*(P-1)/T + 4),0.00056]]
    x, y = np.asarray(N_part_data).T
    fit = np.polyfit(x, y, 2)
    D_N = np.polyval(fit, f_O2)
    
    if D_N < 0:
        D_N = 0
        
    if melting_style == 'batch_melting':    
        m_N_melt =  m_N_mantle/(D_N + F_melt* (1- D_N)) 
    else:
        m_N_melt = m_N_mantle/F_melt*(1-(1-F_melt)**(1/D_N))
    
    
    'H partiotioning'
    D_H2O = 0.01  
    if melting_style == 'batch_melting':   
        m_H2O_melt =  m_H2O_mantle/(D_H2O + F_melt* (1- D_H2O)) 
    else:
        m_H2O_melt = m_H2O_mantle/F_melt*(1-(1-F_melt)**(1/D_H2O))
    
    'S partitioning'
    D_S = 0.01 
    if melting_style == 'batch_melting':   
        m_S_melt =  m_S_mantle/(D_S + F_melt* (1- D_S)) 
    else:
        m_S_melt = m_S_mantle/F_melt*(1-(1-F_melt)**(1/D_S))
    
    
    
    return(m_H2O_melt, m_CO2_melt, m_S_melt, m_N_melt,)


#print(C_H_N_S_partitioning(0.1, 400e-6, 10E-6, 50E-6, 100E-6, 2000, 3000, 'IW', 5))