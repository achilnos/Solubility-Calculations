#Solubility\ Model.py
#solubility model for Ar in rhyolite based on composition and thermodynamic state
#Using model from Iacono_Marziano 2010

import csv
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

Temperature = 1000. #Celcius
Pressure = 200. #MPa
SiO2 = 48.86
TiO2 = 1.73
Al2O3 = 16.77
FeO = 9.71
Fe2O3 = 0.00
MgO = 6.65
CaO = 9.86
Na2O = 3.62
K2O = 1.93

def sol_calc(time, Temperature, Pressure):
    comp_vec = (SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, CaO, Na2O, K2O)
    nor_maj_ele_comp = ( (SiO2*100./np.sum(comp_vec)), (TiO2*100./np.sum(comp_vec)), (Al2O3*100./np.sum(comp_vec)), (FeO*100./np.sum(comp_vec)), (Fe2O3*100./np.sum(comp_vec)), (MgO*100./np.sum(comp_vec)), (CaO*100./np.sum(comp_vec)), (Na2O*100./np.sum(comp_vec)), (K2O*100./np.sum(comp_vec)) )
    mol_frac = (    (  ( nor_maj_ele_comp[0]/60.9/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) ), ( nor_maj_ele_comp[1]/79.9/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) ), ( nor_maj_ele_comp[2]/101.96/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) ), ( nor_maj_ele_comp[3]/71.85/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) ), ( nor_maj_ele_comp[4]/159.7/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) ), ( nor_maj_ele_comp[5]/40.32/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) ), ( nor_maj_ele_comp[6]/56.08/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) ), ( nor_maj_ele_comp[7]/61.98/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) ), ( nor_maj_ele_comp[8]/94.2/ (nor_maj_ele_comp[0]/60.9 + nor_maj_ele_comp[1]/79.9 + nor_maj_ele_comp[2]/101.96 + nor_maj_ele_comp[3]/71.85 + nor_maj_ele_comp[4]/159.7 + nor_maj_ele_comp[5]/40.32 +  nor_maj_ele_comp[6]/56.08 + nor_maj_ele_comp[7]/61.98 + nor_maj_ele_comp[8]/94.2) )  )   )
    mol_num_ox = (   (0.25*( mol_frac[0] - 0.5 * (mol_frac[3] + mol_frac[5] + mol_frac[6])  - mol_frac[7] - mol_frac[8] )), (0.25*mol_frac[1]), (0.375*mol_frac[2]), (0.25*mol_frac[3]), (0.375*mol_frac[4]), (0.25*mol_frac[5]), (0.25*mol_frac[6]), (0.375*mol_frac[7]), (0.375*mol_frac[8]) )
    mol_frac_ox = ( mol_num_ox[0]/np.sum(mol_num_ox), mol_num_ox[1]/np.sum(mol_num_ox), mol_num_ox[2]/np.sum(mol_num_ox), mol_num_ox[3]/np.sum(mol_num_ox), mol_num_ox[4]/np.sum(mol_num_ox), mol_num_ox[5]/np.sum(mol_num_ox), mol_num_ox[6]/np.sum(mol_num_ox), mol_num_ox[7]/np.sum(mol_num_ox), mol_num_ox[8]/np.sum(mol_num_ox) )
    Vca = (13.2983670477867, 13.8124494135477, 20.180231396784, 7.04655952564533, 20.4744046748107, 7.56821655072, 9.14865277864533, 14.4990257250773, 23.9907159149707)
    Wmol = (60.09, 79.9, 101.96, 71.85, 159.7, 40.32, 56.08, 61.98, 94.2)
    Wmol_ox = (240.36, 319.6, 271.9, 407.58, 425.866666666667, 281.46, 344.5, 325.52, 411.44)
    Vmol = (26.92, 22.43, 36.8, 13.35, 41.44, 11.24, 16.27, 28.02, 44.61)
    dPMVdT = (0, 0.00724, 0.00262, 0.00909, 0.00292, 0.00262, 0.00292, 0.00741, 0.01191)
    dPMVdP = (-0.000173, -0.000222, -0.000182, -0.000054, -0.000054, 0.00002, 0.000008, -0.000259, -0.000879)
    dPMVdPdT = (0.00000004, -0.00000003, 0.00000034, -0.00000013, -0.00000013, -0.00000033, -0.00000012, -0.00000041, -0.00000107)
    dPMVdP = ( (dPMVdT[0]+dPMVdPdT[0]*(Temperature+273.-1673.)), (dPMVdT[1]+dPMVdPdT[1]*(Temperature+273.-1673.)), (dPMVdT[2]+dPMVdPdT[2]*(Temperature+273.-1673.)), (dPMVdT[3]+dPMVdPdT[3]*(Temperature+273.-1673.)), (dPMVdT[4]+dPMVdPdT[4]*(Temperature+273.-1673.)), (dPMVdT[5]+dPMVdPdT[5]*(Temperature+273.-1673.)), (dPMVdT[6]+dPMVdPdT[6]*(Temperature+273.-1673.)), (dPMVdT[7]+dPMVdPdT[7]*(Temperature+273.-1673.)), 
    (dPMVdT[8]+dPMVdPdT[8]*(Temperature+273.-1673.)) )
    ref_Vmol = ( (Vmol[0] + dPMVdT[0] * (Temperature + 273. - 1573.) + dPMVdP[0]), (Vmol[1] + dPMVdT[1] * (Temperature + 273. - 1573.) + dPMVdP[1]), (Vmol[2] + dPMVdT[2] * (Temperature + 273. - 1573.) + dPMVdP[2]), (Vmol[3] + dPMVdT[3] * (Temperature + 273. - 1573.) + dPMVdP[3]), (Vmol[4] + dPMVdT[4] * (Temperature + 273. - 1573.) + dPMVdP[4]), (Vmol[5] + dPMVdT[5] * (Temperature + 273. - 1573.) + dPMVdP[5]), (Vmol[6] + dPMVdT[6] * (Temperature + 273. - 1573.) + dPMVdP[6]), (Vmol[7] + dPMVdT[7] * (Temperature + 273. - 1573.) + dPMVdP[7]), (Vmol[8] + dPMVdT[8] * (Temperature + 273. - 1573.) + dPMVdP[8]) )
    melt_weight = np.dot(mol_frac, Wmol)
    melt_weight_ox = np.dot(mol_frac, Wmol_ox)
    melt_vol = np.dot(mol_frac, ref_Vmol)
    melt_vol_ox = melt_vol * melt_weight_ox / melt_weight
    ind_var = (mol_frac[0], mol_frac[1], mol_frac[2], mol_frac[3], mol_frac[4], mol_frac[5], mol_frac[6], mol_frac[7], mol_frac[8], np.power(mol_frac[7], 2))
    model_para = (-7.18411363109547, -15.531671574456, -10.1237252230073, -1.96914424991609, -9.03006905177394, -2.20447137229592, -2.09887726655598, -8.30617337714253, -13.3745366316107, 6.36552265543969)
    alpha = 0.511931205324756
    beta = -53.1042466088579
    fug_coef = 1.37841592470347
    comp_temp_vec = ( mol_frac[0]*(1./Temperature-1./1300.), mol_frac[2]*(1./Temperature-1./1300.), (mol_frac[5]+mol_frac[6])*(1./Temperature-1./1300.), (mol_frac[7]+mol_frac[8])*(1/Temperature-1./1300.), mol_frac[4]*(1./Temperature-1./1300.), np.power(mol_frac[7],2)*(1./Temperature-1./1300.) )
    comp_pressure_vec = ( mol_frac[0]*(Pressure-1./10.)*10., mol_frac[2]*(Pressure-0.1)*10., (mol_frac[5]+mol_frac[6])*(Pressure-0.1)*10., mol_frac[7]*(Pressure-0.1)*10., (mol_frac[3]+mol_frac[4])*(Pressure-0.1)*10., mol_frac[8]*(Pressure-0.1)*10.)
    gamma = (-770.126004476107, -2811.09773370873, 2118.90407498594, 6100.42831453808, -4139.59812962638, -7459.97546972548)
    kappa = (-0.00000206509719744227, 0.0000389611190309771, 0.0000466832446720345, 0.000161765036283955, 0.000102858625737339, 0.000071126014158019)
    lnK = alpha * (100. - 100. / (melt_vol_ox*melt_weight/melt_weight_ox) * (np.dot(Vca, mol_frac) + np.dot(model_para, ind_var) + np.dot(gamma, comp_temp_vec) + np.dot(comp_pressure_vec, kappa)) ) + beta
    wt_per = np.exp(lnK) * 100. * 39.948 * Pressure * 10. * fug_coef / melt_weight
    ccSTPg = np.exp(lnK) * 22414 * Pressure * 10. * fug_coef / melt_weight
    return (lnK, wt_per, ccSTPg, time, Pressure, Temperature)

def data_read():
    time=[]
    pressure=[]
    temperature=[]
    with open('TZMW4_data.csv', 'rU') as csvfile:
        reader=csv.DictReader(csvfile)
        for row in reader:
            #print(row['time'], row['Pressure'], row['Temperature'])
            time.append(row['time'])
            pressure.append(row['Pressure']) 
            temperature.append(row['Temperature'])
        return (time, pressure, temperature)

def solubility_iter():
    t_P_T = []
    time_pressure_temperature = data_read()
    time_pressure_temperature = np.array(time_pressure_temperature)
    for i in time_pressure_temperature.T:
        time, Pressure, Temperature = i
        time = float(time)
        Pressure = float(Pressure)
        Temperature = float(Temperature)
        t_P_T.append(np.array(sol_calc(time, Temperature, Pressure)))
    return t_P_T
    
def plotter():
    t_P_T = solubility_iter()
    count = 0
    lnK_wt_cc_t_P_T = np.empty([1,6])
    while count < len(t_P_T):
        row = t_P_T[count]
        lnK_wt_cc_t_P_T = np.vstack((lnK_wt_cc_t_P_T,row))
        count = count + 1
    lnK_wt_cc_t_P_T = np.delete(lnK_wt_cc_t_P_T,(0),axis=0)
    x = lnK_wt_cc_t_P_T[:,3]
    y = lnK_wt_cc_t_P_T[:,0]
    figure()
    plot(x, y, 'r')
    xlabel('runtime')
    ylabel('solubility (wt%)')
    title('solubility over TZM4')
    show()#problem with the pressure dependance in the model calculation. Need to go through it to find the bug!

    
plotter()