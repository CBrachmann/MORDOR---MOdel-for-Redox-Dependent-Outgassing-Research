# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 13:04:33 2026

@author: cabra
"""




#TODO , improve initial guess, add extrusive



import numpy as np
from scipy import optimize
from scipy.optimize import least_squares
from math import erf
import sys
sys.path.insert(1, 'C:/Users/cabra/Nextcloud/PHD/my_codes/Gibbs_energies')
from Gibbs_E import Gibbs_E
from FugacityC_ALL import Fugacity_C
from initial_guess_ext import initial_guess
import Molecular_weights as MW
from func_CHNS_part import C_H_N_S_partitioning
from func_Fe_FeO_Fe2O3 import Fe_FeO_Fe2O3


def solve_gases_CHONS_ext(P_melt, T_melt, fugacity, buffer, mH2O_init, mCO2_init, mS_init, mN_init, H_law, C_law, S2_law, S6_law, N_law,\
                      mSiO2, mTiO2, mAl2O3, mFeO, mMgO, mCaO, mNa2O, mK2O, mMnO, Partitioning, melting_style, F_melt, gas_behavior,\
                      PH2, PH2O, PCO, PCO2, PCH4, PS2, PSO2, PH2S, PN2):
    
    if gas_behavior == 'ideal':
        phi_H2  = 1
        phi_H2O = 1
        phi_CO  = 1
        phi_CO2 = 1
        phi_CH4 = 1
        phi_S2  = 1
        phi_SO2 = 1
        phi_H2S = 1
        
    else: 
        phi_H2  = Fugacity_C(P_melt, T_melt, 1, 'H2')
        phi_H2O = Fugacity_C(P_melt, T_melt, 1, 'H2O')
        phi_CO  = Fugacity_C(P_melt, T_melt, 1, 'CO')
        phi_CO2 = Fugacity_C(P_melt, T_melt, 1, 'CO2')
        phi_CH4 = Fugacity_C(P_melt, T_melt, 1, 'CH4')
        phi_S2  = Fugacity_C(P_melt, T_melt, 1, 'S2')
        phi_SO2 = Fugacity_C(P_melt, T_melt, 1, 'SO2')
        phi_H2S = Fugacity_C(P_melt, T_melt, 1, 'H2S')
    
    if fugacity == 'IW':
        f_O2 = 10**(6.899 - 27714/T_melt + 0.05*(P_melt-1)/T_melt +buffer)     #IW
    elif fugacity == 'QFM':
        f_O2 = 10**(5.5976 - 24505/T_melt + 0.8099*np.log10(T_melt) + 0.0937*(P_melt-1)/T_melt +buffer)  #QFM
    elif fugacity == 'NiNiO':
        f_O2 = 10**(8.951 - 24556/T_melt + 0.46*(P_melt-1)/T_melt +buffer)     #NiNiO
        
    NNO = 10**(8.951 - 24556/T_melt + 0.46*(P_melt-1)/T_melt)    
    delta_NNO = abs(np.log10(f_O2) - np.log10(NNO))
    FFS = 10**((2/(1987*T_melt)*(-0.2655*(P_melt-1)-35910 + 12.56*T_melt))/2.302581)
    

    'inverse molar weight of melt'
    x = mSiO2/MW.SiO2 + mTiO2/MW.TiO2 + mAl2O3/MW.Al2O3 + mFeO/MW.FeO + mMgO/MW.MgO + mCaO/MW.CaO \
        + mNa2O/MW.Na2O + mK2O/MW.K2O + mMnO/MW.MnO
        
    'melt composition mole fraction'
    XSiO2  = mSiO2/MW.SiO2/x
    XTiO2  = mTiO2/MW.TiO2/x
    XAl2O3 = mAl2O3/MW.Al2O3/x
    XFeO_init   = mFeO/MW.FeO/x 
    XMgO   = mMgO/MW.MgO/x
    XCaO   = mCaO/MW.CaO/x
    XNa2O  = mNa2O/MW.Na2O/x
    XK2O   = mK2O/MW.K2O/x
    XMnO   = mMnO/MW.MnO/x
    
    
    xFe, xFeO, xFe2O3, XFe, XFeO, XFe2O3, aFe2O3 = Fe_FeO_Fe2O3(T_melt,P_melt, fugacity, buffer, XFeO_init)
    
    mFe = XFe*MW.Fe*x
    mFeO = XFeO*MW.FeO*x
    mFe2O3 = XFe2O3*MW.Fe2O3*x
    
    
    'Normalize'
    sum_melt = mSiO2 + mTiO2 + mAl2O3 + mFe + mFeO + mFe2O3 + mMgO + mCaO + mNa2O + mK2O + mMnO
    mSiO2  = mSiO2/sum_melt
    mTiO2  = mTiO2/sum_melt
    mAl2O3 = mAl2O3/sum_melt
    mFe    = mFe/sum_melt
    mFeO   = mFeO/sum_melt
    mFe2O3 = mFe2O3/sum_melt
    mMgO   = mMgO/sum_melt
    mCaO   = mCaO/sum_melt
    mNa2O  = mNa2O/sum_melt
    mK2O   = mK2O/sum_melt
    mMnO   = mMnO/sum_melt

    
    if Partitioning == True:
        mH2O, mCO2, mS, mN, = C_H_N_S_partitioning(F_melt, mH2O_init, mCO2_init, mS_init, mN_init, T_melt, P_melt, fugacity, buffer, x, melting_style)
        
        
    
    else: 
        mH2O = mH2O_init
        mCO2 = mCO2_init
        mS = mS_init
        mN = mN_init

    
    'Recalculate inverse molar weight of melt'
    x = mSiO2/MW.SiO2 + mTiO2/MW.TiO2 + mAl2O3/MW.Al2O3 +mFe/MW.Fe + mFeO/MW.FeO +mFe2O3/MW.Fe2O3 + mMgO/MW.MgO + mCaO/MW.CaO \
        + mNa2O/MW.Na2O + mK2O/MW.K2O + mH2O/MW.H2O + mCO2/MW.CO2 +mS/MW.S + mN/MW.N
        
    
    'Recalculate molar fractions'
    
    XSiO2  = mSiO2/MW.SiO2/x
    XTiO2  = mTiO2/MW.TiO2/x
    XAl2O3 = mAl2O3/MW.Al2O3/x
    XFeO_init   = mFeO/MW.FeO/x 
    XMgO   = mMgO/MW.MgO/x
    XCaO   = mCaO/MW.CaO/x
    XNa2O  = mNa2O/MW.Na2O/x
    XK2O   = mK2O/MW.K2O/x
    XMnO   = mMnO/MW.MnO/x
        
    XH2O   = mH2O/MW.H2O/x
    XCO2   = mCO2/MW.CO2/x
    XS     = mS/MW.S/x
    XN     = mN/MW.N/x
    
    
    NBO = 2*(XNa2O + XK2O + XCaO + XMgO + XFeO - XAl2O3)
    Tetrah = XSiO2 + XTiO2 + XAl2O3
    O = 2* XSiO2 + 2* XTiO2 + 3* XAl2O3 + XMgO + XFeO + XCaO + XNa2O + XK2O
    NBO_T = NBO/Tetrah
    NBO_O = NBO/O  

    
    R =  	8.31446261815324   

    
    #needed bacuse sol law of Boulling and Wood becomes weird at pressures higher than 400 bar
    P_max = min(400, P_melt)
    
    threshhold = 0

    
    'values equal to 0 need to be treated as 1 for the logarithmic system'
    XH2O_init = 1 if XH2O == 0 else XH2O
    XCO2_init = 1 if XCO2 == 0 else XCO2
    XS_init = 1 if XS == 0 else XS
    XN_init = 1 if XN == 0 else XN


    'Sulfide Capacity'
    if S2_law == 'BOU':
        CS2 = 10**(0.225 + (25237*XFeO_init + 5214*XCaO + 12705*XMnO + 19829*XK2O*2 - 1109*XSiO2 - 8879)/T_melt)

    
    if S2_law == 'ONE':
        CS2 = np.exp(-23590/T_melt + 8.77 + (1673/T_melt)*(6.7*(XNa2O*2 + XK2O*2) + 1.8*XAl2O3*2 + 4.9*XMgO + 8.1*XCaO + 5.07*XTiO2 + 8.9*XFeO - 22.2*XFeO*XTiO2 + 7.2*XFeO*XSiO2) - 2.06* erf(-7.2*XFeO) )
    
        
    'Sulfate Capacity'
    if S6_law == 'BOU':
        #CS6 = 10**(-12.806+(14502*XCaO +30375)/T_melt)
        CS6 = 10**(-12.948 + (15602*XCaO + 28649*XNa2O*2 - 9496*XMgO + 4194*XAl2O3*2 + 16016*XMnO + 29244)/T_melt)
    
    if S6_law == 'ONE':
        CS6 = np.exp(-8.02 + (21100 + 44000*XNa2O*2 + 18700*XMgO +4300*XAl2O3*2 + 44200*XK2O*2 +35600*XCaO + 12600*XMnO + 16500*XFeO)/T_melt)
    
    
    'variables for Iacono Marziano Solubility Law'                              
   
    a_H2O   = 0.54    
    b_H2O   = 1.24
    B_H2O   = -2.95
    C_H2O   = 0.02   
                                
    a_CO2   = 1                                 
    b_CO2   = 15.8
    B_CO2   = -5.3
    C_CO2   = 0.14                                    

    d_H2O   = 2.3  
    d_Al2O3 = 3.8
    d_FeO   = -16.3
    d_Na2O  = 20.1
    
    
    'Variables Scaillet and Pichavant'

    a_SP  = 7.28E-6
    b_SP  = 0.00084107
    c_SP  = 0.00488128
    d_SP  = 0.0818873
    e_SP  = -0.0224068
    f_SP  = 0.22801636
    g_Si  = -0.012467
    g_Al  = -0.0015766
    g_Fe3 = 0.37362348
    g_Fe2 = 0.0674383
    g_Mg  = 0.01121929
    g_Ca  = 0.02000831
    g_Na  = 0.05644745
    g_K   = -0.0248037
    g_Ti  = 0.00672403
    g_H   = 0.06868295
    
    giWi = g_Si*mSiO2*100 + g_Al*mAl2O3*100 + g_Fe3*mFe2O3*100 + g_Fe2*mFeO*100 + g_Mg*mMgO*100 + g_Ca*mCaO*100 + g_Na*mNa2O*100 + g_K*mK2O*100 + g_Ti*mTiO2*100
    xidi = (XAl2O3/(XCaO + XK2O +XNa2O))*d_Al2O3 + (XFeO+ XMgO)* d_FeO + (XNa2O + XK2O)*d_Na2O
    'Non-bridging-oxygen /oxygen'

    F1 = np.log(1/(MW.CO2*x*10**6)) + C_CO2 * P_melt /T_melt + xidi + B_CO2 + b_CO2*NBO_O  #np.log(1/(MW_CO2*x*10**6)) +
    F2 = np.log(1/(MW.H2O*x*100)) +C_H2O * P_melt / T_melt + B_H2O + b_H2O*NBO_O
    
    'Gibbs energies and equilibrium constants'
    Gf_CH4   = Gibbs_E(T_melt, 'CH4')*10**3#-77437   + 22.5098 * T_melt *np.log10(T_melt) + 29.5967 * T_melt            #2C (graphite) + 2H2 = CH4;  
    Gf_CO    = Gibbs_E(T_melt, 'CO') *10**3#(-214104 + 25.2183 * T_melt *np.log10(T_melt) -262.1545 * T_melt)/2         #2C (graphite) + O2 = 2CO
    Gf_CO2   = Gibbs_E(T_melt, 'CO2')*10**3#-392647  +  4.5855 * T_melt *np.log10(T_melt) - 16.9762 * T_melt            #2C (graphite) + O2 = CO2
    Gf_H2O   = Gibbs_E(T_melt, 'H2O')*10**3#(-483095 + 25.3687 * T_melt *np.log10(T_melt) + 21.9563 * T_melt)/2         #2H2 + O2 = 2H2O
    Gf_H2S = Gibbs_E(T_melt, 'H2S')*10**3
    Gf_SO2 = Gibbs_E(T_melt, 'SO2')*10**3
    Gf_NH3 = Gibbs_E(T_melt, 'NH3')*10**3

    dG_2 = Gf_CO - Gf_CO2                        #CO2 = CO + 1/2O2
    dG_3 = Gf_CH4 -(Gf_CO2 + 2*Gf_H2O)           #CO2 + 3H2O = CH4 + 2O2
    dG_8 = (Gf_H2O-Gf_H2S)
    dG_10 = 2*Gf_NH3 - 3*Gf_H2O
    
    'eqilibrium constants'
    K1 = np.exp(((Gf_H2O))/(R*T_melt))
    K2 = np.exp(-dG_2 /(R*T_melt))
    K3 = np.exp(-dG_3 /(R*T_melt))
    K4 = np.exp(-Gf_SO2/(R*T_melt))       
    K8 = np.exp(-dG_8/(R*T_melt))
    K10 = np.exp((-dG_10)/(R*T_melt))
    


    
    def _solve_XOH(XH2O_mol, guess=1e-6):

        def f(XOH):
            return -np.log(XOH**2 / (XH2O_mol * (1.0 - XOH - XH2O_mol))) \
                   - (0.403 + 15.33 * XOH + 10.894 * XH2O_mol)
        XOH, = optimize.fsolve(f, guess, xtol=1e-12)

        if XOH < 0.0:
            XOH = 0.0
        if XOH > 1.0 - XH2O_mol:
            XOH = 1.0 - XH2O_mol
        return(XOH)
    
    
    def safe_log(z, tiny=1e-300):
        return np.log(np.clip(z, tiny, None))
    
    
    def ln_sulfur_trioxide_from_s2(XS2, T, P, f_O2):
        K_SO3 = 55921.0/T - 25.07 + 0.6465*np.log(T)
        return safe_log(K_SO3) + 0.5*safe_log(XS2) + 1.5*safe_log(f_O2) - 0.5*safe_log(P)
    
    
    
    'S and N solubility - are calculated before speciation beacuse of their fO2 dependence C and H solubilities are calculated later'

    'Sulfur'
    if fugacity == 'IW' and buffer < 5 or fugacity == 'QFM' and buffer < 1.3 or fugacity == 'NiNiO' and buffer < 0.3:
        if S2_law == 'GAI':
            lnX_S = 0.5 * np.log(PS2) - 0.5*np.log(f_O2) + 0.12*mFeO \
            + (13.84 - 26480.0/T_melt) - np.log(MW.S * x) - 6.0*np.log(10.0)
            
            
        if S2_law == 'BOU':
            lnX_S = np.log(CS2) + 0.5 * np.log(PS2) - 0.5*np.log(f_O2) - np.log(100) - np.log(MW.S) - np.log(x)
    
            
        if S2_law == 'ONE':
            lnX_S = (np.log(CS2) - 0.5 * np.log(PS2) - 0.5*np.log(f_O2)) - 6*np.log(10) - np.log(MW.S) - np.log(x)
    
    
    elif fugacity == 'IW' and buffer >= 5 or fugacity == 'QFM' and buffer >= 1.3 or fugacity == 'NiNiO' and buffer >= 0.3:
        
        
        if S6_law == 'SCA':
            deltaFFS2 = np.log10(PS2) - np.log10(FFS)
            XS_melt = 10**(7.28E-6*P_melt + 0.00084107*(T_melt-273.15) + 0.00488128*delta_NNO**3 
                           + 0.0818873*delta_NNO**2 + -0.0224068*delta_NNO*deltaFFS2 + 0.22801636*deltaFFS2 + giWi) *10**-6/MW.S/x
            lnX_S = np.log(XS_melt)
    
        
        if S6_law == 'BOU':
            lnX_S = np.log(CS6) +np.log(PS2) + 0.5*np.log(f_O2) - (18887*np.log(10))/T_melt + 3.8064*np.log(10) \
            - (1.5237*np.log(10))/T_melt*(P_max - 0.1) - np.log(MW.S * x) - 2.0*np.log(10.0)
    
            
        if S6_law == 'ONE':
            lnXSO3_aux = ln_sulfur_trioxide_from_s2(PS2, T_melt, P_melt, f_O2)
            lnX_S = (np.log(CS6) + lnXSO3_aux + np.log(P_melt)- 6 * np.log(10) - np.log(MW.S)- np.log(x))
    
    XS_melt = np.exp(lnX_S)
    XS_gas = max(0, XS_init - XS_melt)
    
    
    if N_law == 'BER':
        dG_13 = 29344 - T_melt * -121 + P_melt * 4
        dG_14 = 183733 - T_melt * -172 + P_melt * -5
        K13 = np.exp(-dG_13/(R*T_melt))
        K14 = np.exp(-dG_14/(R*T_melt))

        
        ln_term1 = np.log(K13) + np.log(PN2)
        ln_term2 = np.log(K14) - 0.75*np.log(f_O2) + 0.5*np.log(PN2)
        lnX_N = np.logaddexp(ln_term1, ln_term2)
            
    if N_law == 'LIB_1':
        lnX_N = -6 *np.log(10) - np.log(MW.N) - np.log(x) + np.log(0.0611*( PN2) 
                                      + 5.97 *10**-10 *f_O2**(-3/4)*np.sqrt(PN2) )

    
    
    if N_law == 'LIB_2':
        lnX_N = -6 *np.log(10) - np.log(MW.N) - np.log(x) + np.log(0.0611*( PN2) 
                                      + 0.8 *10**-7 *f_O2**(-0.55)*np.sqrt(PN2) )

    
    
    if N_law == 'LIB_3':
        lnX_N = -6 *np.log(10) - np.log(MW.N) - np.log(x) + np.log((0.23 - 0.16*NBO_T)*( PN2) 
                                     + f_O2**(-3/4)*(10**-10 + 2.2*10**-10*NBO_T)*np.sqrt(PN2) )


    
    
    if N_law == 'DAS':
        IW = 10**(6.899 - 27714/T_melt + 0.05*(P_melt-1)/T_melt)
        delta_IW = np.log10(f_O2) - np.log10(IW)
        mN_melt = ((PN2*10**-4)**0.5*np.exp(5908 *(P_melt*10**-4)**0.5/T_melt - 1.6 *delta_IW )
                   + (PN2*10**-4)*np.exp(4.67 + 7.11*XSiO2 - 13.06* XAl2O3 -120*XTiO2))*10**-6
        XN_melt = mN_melt/MW.N/x
        lnX_N = np.log(XN_melt)
 
    
    if N_law == 'BRA':

        lnX_N = -2 *np.log(10) - np.log(MW.N) - np.log(x) + np.log(3.102140e-06*( PN2) + 6.685873e-13
                *f_O2**(-6.519758e-01)*np.sqrt(PN2) *np.exp(2.342337*np.sqrt(P_melt*1E-4)))
        
        
    XN_melt = np.exp(lnX_N)
    XN_gas = max(0, XN_melt)
    
    def sys_log_raw(y):
        (
         yXH2, yXH2O, yXCH4, yXCO, yXCO2,
         yXS2, yXSO2, yXH2S, yXN2, yXNH3) = y
    
    
        XH2   = inv_logit(yXH2)
        XH2O  = inv_logit(yXH2O)
        XCH4  = inv_logit(yXCH4)
        XCO   = inv_logit(yXCO)
        XCO2  = inv_logit(yXCO2)
        XS2   = inv_logit(yXS2)
        XSO2  = inv_logit(yXSO2)
        XH2S  = inv_logit(yXH2S)
        XN2   = inv_logit(yXN2)
        XNH3  = inv_logit(yXNH3)

    
    
        lnXH2   = safe_log(XH2)
        lnXH2O  = safe_log(XH2O)
        lnXCH4  = safe_log(XCH4)
        lnXCO   = safe_log(XCO)
        lnXCO2  = safe_log(XCO2)
        lnXS2   = safe_log(XS2)
        lnXSO2  = safe_log(XSO2)
        lnXH2S  = safe_log(XH2S)
        lnXN2   = safe_log(XN2)
        lnXNH3  = safe_log(XNH3)



        'Mass balance'
        
        if mH2O > threshhold and mCO2 > threshhold and XS_gas > threshhold and XN_gas > threshhold:
            MB_HCSN = (XH2O_init- (2*XCH4 + XH2O + XH2 + XH2S + 1.5* XNH3) )
            
        elif mH2O > threshhold and mCO2 > threshhold and XS_gas > threshhold and XN_gas <= threshhold:
            MB_HCS  = (XH2O_init- (2*XCH4 + XH2O + XH2 + XH2S) )
            
        elif mH2O > threshhold and mCO2 > threshhold and XS_gas <= threshhold and XN_gas <= threshhold:
            MB_HC  = (XH2O_init- (2*XCH4 + XH2O + XH2) )
        
        elif mH2O > threshhold and mCO2 > threshhold and XS_gas <= threshhold and XN_gas > threshhold:
            MB_HCN = (XH2O_init- (2*XCH4 + XH2O + XH2 + 1.5* XNH3) )
        
        elif mH2O > threshhold and mCO2 <= threshhold and XS_gas > threshhold and XN_gas <= threshhold:
            MB_HS  = (XH2O_init- ( XH2O + XH2 + XH2S) )
        
        elif mH2O > threshhold and mCO2 <= threshhold and XS_gas > threshhold and XN_gas > threshhold:
            MB_HSN = (XH2O_init- ( XH2O + XH2 + XH2S + 1.5* XNH3 ))
            
        elif mH2O > threshhold and mCO2 <= threshhold and XS_gas <= threshhold and XN_gas > threshhold:
            MB_HN = (XH2O_init- (XH2O + XH2 + 1.5* XNH3) )    
            
        elif mH2O > threshhold and mCO2 <= threshhold and XS_gas <= threshhold and XN_gas <= threshhold:
            MB_H  = (XH2O_init- - (XH2O + XH2) )
         
            
        if mH2O > threshhold and mCO2 > threshhold:
            MB_CH = (XCO2_init- (XCH4 + XCO + XCO2) )
        elif mH2O <= threshhold and mCO2 > threshhold:
            MB_C = (XCO2_init- (XCO + XCO2) )
           
        if mH2O > threshhold and XS_gas > threshhold:    
            MB_SH = (XS_gas- (2*XS2 + XSO2 + XH2S) )
        if mH2O <= threshhold and XS_gas > threshhold:  
            MB_S = (XS_gas - (2*XS2 + XSO2) )
            
        if mH2O > threshhold and XN_gas > threshhold:
            MB_NH = XN_gas  - (XNH3 +2*XN2) 
        if mH2O <= threshhold and XN_gas > threshhold:
            MB_N = XN_gas - (2*XN2) 
            
        # -------------------
        # 7–9: Gas equilibria
        # -------------------
        
        lnP = np.log(P_melt)
        lnfH2 = np.log(phi_H2) + lnXH2 + lnP
        lnfH2O = np.log(phi_H2O) + lnXH2O + lnP
        lnfCO = np.log(phi_CO) + lnXCO + lnP
        lnfCO2 = np.log(phi_CO2) + lnXCO2 + lnP
        lnfCH4 = np.log(phi_CH4) + lnXCH4 + lnP
        lnfS2 = np.log(phi_S2) + lnXS2 + lnP
        lnfSO2 = np.log(phi_SO2) + lnXSO2 + lnP
        lnfH2S = np.log(phi_H2S) + lnXH2S + lnP
        lnfO2 = np.log(f_O2)
        
        EQ_H2 = lnfH2 - (np.log(K1) - 0.5*lnfO2 + lnfH2O)
        
        EQ_CO = lnfCO - (np.log(K2) - 0.5*lnfO2 + lnfCO2)
        
        EQ_CH4 = np.log(K3) - 2*lnfO2 + lnfCO2 + 2*lnfH2O - lnfCH4
        
        EQ_S2 = np.log(K4) + lnfO2 + 0.5*lnfS2 - lnfSO2
        
        EQ_H2S = np.log(K8) + lnfH2S + 0.5*lnfO2 - 0.5*lnfS2 - lnfH2O
    
        EQ_N2 = ( np.log(K10) + lnXN2 + 3 * lnXH2O + 2 * np.log(P_melt) - 2 * lnXNH3 - 1.5 * np.log(f_O2))
    

        
        
        
        if mH2O > threshhold and mCO2 > threshhold and XS_gas > threshhold and XN_gas > threshhold:
            return np.array([ MB_HCSN, MB_CH, MB_SH, MB_NH, EQ_H2, EQ_CO, EQ_CH4, EQ_S2, EQ_H2S, EQ_N2])
        
        if mH2O > threshhold and mCO2 > threshhold and XS_gas > threshhold and XN_gas <= threshhold:
            return np.array([ MB_HCS, MB_CH, MB_SH, 0, EQ_H2, EQ_CO, EQ_CH4, EQ_S2, EQ_H2S, 0])
        
        if mH2O > threshhold and mCO2 > threshhold and XS_gas <= threshhold and XN_gas > threshhold:
            return np.array([ MB_HCN, MB_CH, 0, MB_NH, EQ_H2, EQ_CO, EQ_CH4, 0, 0, EQ_N2])
        
        if mH2O > threshhold and mCO2 > threshhold and XS_gas <= threshhold and XN_gas <= threshhold:
            return np.array([ MB_HC, MB_CH, 0, 0, EQ_H2, EQ_CO, EQ_CH4, 0, 0, 0])
        
        if mH2O > threshhold and mCO2 <= threshhold and XS_gas <= threshhold and XN_gas > threshhold:
            return np.array([ MB_HN, 0, 0, MB_NH, EQ_H2, 0, 0, 0, 0, EQ_N2])
        
        if mH2O > threshhold and mCO2 <= threshhold and XS_gas <= threshhold and XN_gas <= threshhold:
            return np.array([ MB_H, 0, 0, 0, EQ_H2, 0, 0, 0, 0, 0])
            
        if mH2O > threshhold and mCO2 <= threshhold and XS_gas > threshhold and XN_gas > threshhold:
            return np.array([ MB_HSN, 0, MB_SH, MB_NH, EQ_H2, 0, 0, EQ_S2, EQ_H2S, EQ_N2])
        
        if mH2O > threshhold and mCO2 <= threshhold and XS_gas > threshhold and XN_gas <= threshhold:
            return np.array([ MB_HS, 0, MB_SH, 0, EQ_H2, 0, 0, EQ_S2, EQ_H2S, 0])
            
        if mH2O <= threshhold and mCO2 > threshhold and XS_gas > threshhold and XN_gas > threshhold:
            return np.array([ 0, MB_C, MB_S, MB_N, 0, EQ_CO, 0, EQ_S2, 0, 0])
        
        if mH2O <= threshhold and mCO2 > threshhold and XS_gas > threshhold and XN_gas <= threshhold:
            return np.array([  0, MB_C, MB_S, 0, 0, EQ_CO, 0, EQ_S2, 0, 0])
            
        if mH2O <= threshhold and mCO2 > threshhold and XS_gas <= threshhold and XN_gas > threshhold:
            return np.array([  0, MB_C, 0, MB_N, 0, EQ_CO, 0, 0, 0, 0])
        
        if mH2O <= threshhold and mCO2 > threshhold and XS_gas <= threshhold and XN_gas <= threshhold:
            return np.array([ 0, MB_C, 0, 0, 0, EQ_CO, 0, 0, 0, 0])
            
        if mH2O <= threshhold and mCO2 <= threshhold and XS_gas > threshhold and XN_gas > threshhold:
            return np.array([  0, 0, MB_S, MB_N,  0, 0, 0, EQ_S2, 0, 0])
        
        if mH2O <= threshhold and mCO2 <= threshhold and XS_gas > threshhold and XN_gas <= threshhold:
            return np.array([  0, 0, MB_S, 0,  0, 0, 0, EQ_S2, 0, 0])
        
        if mH2O <= threshhold and mCO2 <= threshhold and XS_gas <= threshhold and XN_gas > threshhold:
            return np.array([  0, 0, 0, MB_N,  0, 0, 0, 0, 0, 0])
        
        if mH2O <= threshhold and mCO2 <= threshhold and XS_gas <= threshhold and XN_gas <= threshhold:
            return np.array([  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0])
        
        
    def logit(p):
        eps = np.finfo(float).eps
        p = np.clip(p, eps, 1 - eps)
        return np.log(p/(1-p))

    def inv_logit(y):
        return 1.0/(1.0 + np.exp(-y))
      
    XH2_g, XH2O_g, XCH4_g, XCO_g, XCO2_g, XS2_g, XSO2_g, XH2S_g, XN2_g, XNH3_g   = initial_guess(P_melt, XH2O_init, XCO2_init, XS_gas, XN_gas, fugacity, buffer, threshhold, Partitioning)
    y_init = np.log([ XH2_g, XH2O_g, XCH4_g, XCO_g, XCO2_g, XS2_g, XSO2_g, XH2S_g, XN2_g, XNH3_g])
    y_init = np.append(y_init, 0.1)
     

    guess = np.array([
                      XH2_g, XH2O_g, XCH4_g, XCO_g, XCO2_g,
                      XS2_g, XSO2_g, XH2S_g, XN2_g, XNH3_g])
    guess = np.clip(guess, 1e-12, 1-1e-12)           
    y0 = logit(guess)                                              

    # solve the syestem
    sol = least_squares(sys_log_raw, y0, method='trf',
                      max_nfev=1_000_000, ftol=1e-12, xtol=1e-12)
    
    #sol = optimize.root(sys_log_raw, y0, options={'maxiter': 100000}, method='lm')

    if not sol.success:
        raise RuntimeError(f"System does not converge: {sol.message}")

    # transform back
    (
     yXH2, yXH2O, yXCH4, yXCO, yXCO2,
     yXS2, yXSO2, yXH2S, yXN2, yXNH3) = sol.x



    XH2   = inv_logit(yXH2)
    XH2O  = inv_logit(yXH2O)
    XCH4  = inv_logit(yXCH4)
    XCO   = inv_logit(yXCO)
    XCO2  = inv_logit(yXCO2)
    XS2   = inv_logit(yXS2)
    XSO2  = inv_logit(yXSO2)
    XH2S  = inv_logit(yXH2S)
    XN2   = inv_logit(yXN2)
    XNH3  = inv_logit(yXNH3)

    
    'Solubility'
    'H2O'
    if H_law == 'IM':
        lnX_H2O = (a_H2O * np.log(PH2O) + F2)

    if H_law == 'PW_basalt':
        lnX_H2O = (np.log(0.1078)+ 0.7 * (np.log(PH2O) + np.log(0.1))- np.log(100) - np.log(MW.H2O)- np.log(x) )


    if H_law == 'PW_rhyolite':
        lnX_H2O = (np.log(0.4111)+ 0.5 * (np.log(PH2O)+ np.log(0.1)) - np.log(100) - np.log(MW.H2O)- np.log(x))


    if H_law == 'LI':
        lnX_H2O = (np.log(1.033)+ (1 / 1.747) * (np.log(PH2O) + 5 * np.log(10))- 6 * np.log(10)- np.log(MW.H2O) - np.log(x))


    if H_law == 'BO':
        lnX_H2O = (np.log(0.068) + 0.7 * (np.log(PH2O) + 5 * np.log(10)) - 6 * np.log(10) - np.log(MW.H2O) - np.log(x))


    if H_law == 'DI':
        XH2O_coef = -3.04E-5 + 1.29E-6 * mSiO2 * 100
        FC_H2O = Fugacity_C(P_melt, T_melt, 1, 'H2O')
        lnXH2O_mol = (np.log(XH2O_coef) + np.log(PH2O) + np.log(FC_H2O) + (-12 * (P_melt - 1)) / (83.15 * T_melt))
        XOH_mol = _solve_XOH(np.exp(lnXH2O_mol))
        lnXOH_mol = safe_log(XOH_mol)
        lnX_H2O =  np.logaddexp(lnXH2O_mol, lnXOH_mol)

    if H_law == 'BRA_basalt':
        lnX_H2O = (-3.191550 - np.log(100 * MW.H2O * x)+ 6.343862e-01 * np.log(PH2O)  -2.454957e-03 * np.log(PH2O)**2 ) 
        #lnX_H2O = (-2.850134 - np.log(100 * MW.H2O * x)+ 5.959924e-01 * lnPH2O  -3.601632e-03 * lnPH2O**2 ) 


    if H_law == 'BRA_rhyolite':
        lnX_H2O = ( -2.464232 - np.log(100 * MW.H2O * x) + 5.880529e-01 * np.log(PH2O)-4.336921e-03 * np.log(PH2O)**2)


    
    
    'CO2'
    if C_law == 'IM':
        lnX_CO2 = (d_H2O * np.exp(lnX_H2O) + a_CO2 * np.log(PCO2) + F1)

    if C_law == 'PW_basalt' or C_law == 'PW_rhyolite':
        lnX_CO2 = np.log(2.3e-6) + np.log(PCO2) - np.log(MW.CO2) - np.log(x)


    if C_law == 'LI':
        lnX_CO2 = np.log(1.937e-9) + (1/0.714) * (np.log(PCO2) + 5*np.log(10)) - 6*np.log(10) - np.log(MW.CO2) - np.log(x)


    if C_law == 'BO':
        lnX_CO2 = np.log(4.4e-7) + np.log(PCO2) - np.log(MW.CO2) - np.log(x)


    if C_law == 'DI':
        XCO2_coef = 8.7E-6 - 1.7E-7 * mSiO2 * 100
        FC_CO2 = Fugacity_C(P_melt, T_melt, 1, 'CO2')
        lnXCO32_melt = np.log(XCO2_coef) + np.log(PCO2) + np.log(FC_CO2) + (-23 * (P_melt - 1)) / (83.15 * T_melt)
        XCO32_melt = np.exp(lnXCO32_melt)
        lnX_CO2 = np.log(4400) + lnXCO32_melt - np.log(36.6 - 44 * XCO32_melt) - 2*np.log(10) - np.log(MW.CO2) - np.log(x)

        
    if C_law == 'BRA_rhyolite' or C_law == 'BRA_basalt':

        lnX_CO2= ( -1.061231e+01 - np.log(100 * MW.CO2 * x) + 1.033317 * np.log(PCO2) + 5.087957e-03 * np.log(PCO2)**2)

    
    
    XH2O_melt = np.exp(lnX_H2O)
    XCO2_melt = np.exp(lnX_CO2)
    XH2_melt  = 2.572*10**-6 *(PH2 *10**5)**(1/1)     *10**-6/MW.H2/x
    XCH4_melt = 9.937*10**-8 *(PCH4*10**5)**(1/1)     *10**-6/MW.CH4/x
    XCO_melt  = 1.6  *10**-7 *(PCO *10**5)**(1/1)     *10**-6/MW.CO/x
    
    
    XH2  = max(XH2  - XH2_melt,  0)
    XH2O = max(XH2O - XH2O_melt, 0)
    XCH4 = max(XCH4 - XCH4_melt, 0)
    XCO  = max(XCO  - XCO_melt,  0)
    XCO2 = max(XCO2 - XCO2_melt, 0)
 
    

        
    
    return(XH2, XH2O, XCH4, XCO, XCO2, XS2, XSO2, XH2S, XN2, XNH3,  XH2O_melt, XCO2_melt, XS_melt, XN_melt)
    

