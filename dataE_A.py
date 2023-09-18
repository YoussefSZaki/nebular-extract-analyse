#!/usr/bin/env python3

import sys
import csv
import re
import pyneb as pn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename=sys.argv[1]
newfile=sys.argv[2]

def elements_extraction(filename):
    elements=[]
    with open(filename, 'r') as observation_data:
        obs_reader=csv.reader(observation_data)
        for row in observation_data.readlines() :
            if obs_reader.line_num == 1:
                continue
            pattern=r'\[\w+ \w+\]|(He|H|C) \w+'
            match=re.search(pattern,row)
            if match == None:
                continue
            elif match.group() in elements:
                continue
            elif match.group() == '[Kr IV]':
                continue
            elif match.group() == '[Fe III]':
                continue
            elements.append(match.group())
    elements.sort()
    return elements 


def wavelengths_extract(filename):
    wavelengths_data={}
    with open(filename, 'r') as observation_data:
        obs_reader=csv.reader(observation_data)
        for row in obs_reader:
            for element in elements_extraction(filename):  
                    if element in row:
                        if  element in wavelengths_data:
                            wavelengths_data[element].append(row[0]+'A')

                        else:  
                            wavelengths_data[element] = [row[0]+'A']
    return wavelengths_data

def recandcoll_mod(data):
    for ion in list(data.keys()):
        pattern=r'^(\w+) (\w+)$'        
        match=re.search(pattern,ion)
        if match == None:
            continue
        elif match[2] == 'I':
            data['{}1r'.format(match[1])]=data.pop(ion)
        elif match[2] == 'II':
            data['{}2r'.format(match[1])]=data.pop(ion)
    
    for ion in list(data.keys()):
        pattern=r'\[(\w+) (\w+)\]'        
        match=re.search(pattern,ion)
        if match == None:
            continue
        elif match[2] == 'I':
            data['{}1'.format(match[1])]=data.pop(ion)
        elif match[2] == 'II':
            data['{}2'.format(match[1])]=data.pop(ion)
        elif match[2] == 'III':
            data['{}3'.format(match[1])]=data.pop(ion)
        elif match[2] == 'IV':
            data['{}4'.format(match[1])]=data.pop(ion)
        elif match[2] == 'V':
            data['{}5'.format(match[1])]=data.pop(ion)
    return data


def combine_name(filename):
    pyneb_labels=[]
    for ion, wavelength in recandcoll_mod(wavelengths_extract(filename)).items():
            for label,values in pn.LINE_LABEL_LIST.items():
                    for wave in wavelength:
                        for x in values:
                            if ion == label and wave == x:
                                pyneb_labels.append('{}_{}'.format(ion,wave))
    return pyneb_labels

def flux_err_extract(filename):
    flux_err_data={}
    with open(filename, 'r') as observation_data:
        obs_reader=csv.reader(observation_data)
        for row in obs_reader:
            for ion in combine_name(filename):
                pattern=r'^([A-Za-z]+|[A-Z]+)'
                match=re.search(pattern,ion)
                if ion[-5:-1] in row[0] and match[1] in row[2]:
                    if  ion in flux_err_data:
                        continue
                    else:  
                        flux_err_data[ion] = [row[3],row[4]]
        return flux_err_data

def file_output(filename,newfile):
    with open(newfile,'w+') as file:
        file.write('LINE    Flux    err_Flux\n')
        for ion, data in flux_err_extract(filename).items():
            if ion == 'O2_7319A' or ion == 'O2_7330A':
                file.write('{}+   {}  {} \n'.format(ion,data[0],data[1]))
            else:
                file.write('{}   {}  {} \n'.format(ion,data[0],data[1]))
    return newfile


def Compute_obs(filename,newfile):
    
    #First we read the data file we created 
    
    obs_data = file_output(filename,newfile)
    obs = pn.Observation(obs_data, fileFormat='lines_in_rows_err_cols',corrected=False, errIsRelative=False)
    
    #Then we obtain the cHbeta & E(B-V)
    
    six=float((obs.getIntens(returnObs=True)['H1r_6563A']))
    four=float((obs.getIntens(returnObs=True)['H1r_4861A']))
    result=six/four
    rc = pn.RedCorr() 
    rc.law = 'S79 H83 CCM89'
    rc.setCorr(obs_over_theo=result / 2.86, wave1=6563., wave2=4861.)
    obs.extinction.cHbeta = rc.cHbeta
    obs.extinction.E_BV= rc.E_BV
    
    #Afterwards we correct the data and generate a report containing all the corrected intensities of observed ions
    
    obs.correctData(normWave=4861.)
    
    with open('Corrected Data.csv','w+') as file:
        file.write('Line,Corrected Flux\n')
        for line in obs.getSortedLines():
            file.write('{},{:.3E}\n'.format(line.label, line.corrIntens[0]))
    
    #Then we compute Te and Ne using the specified diagnostics
    
    diags = pn.Diagnostics()
    diags.addDiag(['[OIII] 4363/5007+',
                   '[NII] 5755/6548',
                   '[OII] 3726/3729',
                   '[SII] 6731/6716', 
                   '[ClIII] 5538/5518',
                   '[ArIV] 4740/4711',
                   '[ArIII] 5192/7300+',
                   '[SII] 4072+/6720+',
                   '[OI] 5577/6300+'])
    
    low_density=['[SII] 6731/6716','[OII] 3726/3729']
    low_temprature=['[NII] 5755/6548','[SII] 4072+/6720+','[OI] 5577/6300+']
    medium_density=['[ClIII] 5538/5518','[ArIV] 4740/4711']
    medium_temprature=['[ArIII] 5192/7300+','[OIII] 4363/5007+']


    low_Tmean, low_Nmean, med_Tmean, med_Nmean, nS_II, tS_II, nO_II, tO_II, nCl_III, tCl_III, nAr_IV, tAr_IV = ([] for i in range(12))

    data1 = {'Temp/Den':['N II','S II', 'O I']}
    data2 = {'Temp/Den':['Ar III','O III']}

    #To create denisty & Temprature values table for low ionizing regions

    for x in low_density:
        for y in low_temprature:
            Te, Ne = diags.getCrossTemDen(y,x,obs=obs)
            if x == '[SII] 6731/6716':
                nS_II.append('{:.2f} cm-1'.format(Ne))
                low_Nmean.append(Ne)
                tS_II.append('{:.0f} K'.format(Te))
                low_Tmean.append(Te)
            elif x == '[OII] 3726/3729':
                low_Nmean.append(Ne)
                nO_II.append('{:.2f} cm-1'.format(Ne))
                low_Tmean.append(Te)
                tO_II.append('{:.0f} K'.format(Te))

    #Density Table
    df1 = pd.DataFrame(data1, columns = ['Temp/Den'])
    df1['S II'] = nS_II
    df1['O II'] = nO_II
    fig, ax = plt.subplots(1, 1)
    plt.title("Density Values for low ionizing regions")
    ax.axis('off')
    ax.axis('tight')
    ax.table(cellText=df1.values, colLabels=df1.columns, loc='center')
    plt.savefig('Low_Density_table.png')

    #Temprature Table
    df2 = pd.DataFrame(data1, columns = ['Temp/Den'])
    df2['S II'] = tS_II
    df2['O II'] = tO_II
    fig, ax = plt.subplots(1, 1)
    plt.title("Temprature Values for low ionizing regions")
    ax.axis('off')
    ax.axis('tight')
    ax.table(cellText=df2.values, colLabels=df2.columns, loc='center')
    plt.savefig('Low_Temprature_table.png')
    #To create denisty & Temprature values table for medium ionizing regions

    for x in medium_density:
        for y in medium_temprature:
            Te, Ne = diags.getCrossTemDen(y,x,obs=obs)
            if x == '[ClIII] 5538/5518':
                nCl_III.append('{:.2f} cm-1'.format(Ne))
                med_Nmean.append(Ne)
                tCl_III.append('{:.0f} K'.format(Te))
                med_Tmean.append(Te)
            elif x == '[ArIV] 4740/4711':
                med_Nmean.append(Ne)
                nAr_IV.append('{:.2f} cm-1'.format(Ne))
                med_Tmean.append(Te)
                tAr_IV.append('{:.0f} K'.format(Te))

    #Density Table
    df3 = pd.DataFrame(data2, columns = ['Temp/Den'])
    df3['Cl III'] = nCl_III
    df3['Ar IV'] = nAr_IV
    fig, ax = plt.subplots(1, 1)
    plt.title("Density Values for medium ionizing regions")
    ax.axis('off')
    ax.axis('tight')
    ax.table(cellText=df3.values, colLabels=df3.columns, loc='center')
    plt.savefig('Medium_density_table.png')

    #Temprature Table
    df4 = pd.DataFrame(data2, columns = ['Temp/Den'])
    df4['Cl III'] = tCl_III
    df4['Ar IV'] = tAr_IV
    fig, ax = plt.subplots(1, 1)
    plt.title("Temprature Values for medium ionizing regions")
    ax.axis('off')
    ax.axis('tight')
    ax.table(cellText=df4.values, colLabels=df4.columns, loc='center')
    plt.savefig('Medium_temprature_table.png')

    #To compute the mean values for Te and Ne
    
    Low_Te_mean = np.mean(low_Tmean)
    Low_Ne_mean = np.mean(low_Nmean)
    Medium_Te_mean = np.mean(med_Tmean)
    Medium_Ne_mean = np.mean(med_Tmean)
    
    
    # We will now compute all the ionic abundances for our observations
        
    all_atoms = pn.getAtomDict(atom_list=obs.getUniqueAtoms())
    line_ab = {}
    ion_ab = {}
    temp = Low_Te_mean
    dens = Low_Ne_mean
    for line in obs.getSortedLines():
        line_ab[line.label] = all_atoms[line.atom].getIonAbundance(line.corrIntens, temp, dens, to_eval=line.to_eval)
        if line.atom not in ion_ab:
            ion_ab[line.atom] = []
        ion_ab[line.atom].append(line_ab[line.label][0])
            
    # Now creating a report containing the ions and the mean value from all the lines observed for each ion
    
    with open('Ion_abun_report.csv','w+') as file:
        atom_abun={}
        file.write('Ion,Ionic Abundance\n')
        for atom in ion_ab:
            mean = np.mean(np.asarray(ion_ab[atom]))
            atom_abun[atom] = mean
            file.write('{},{:.3E}\n'.format(atom, mean))

    
    # performs all the possible element abundance computations given the input ionic abundance set and available icf  
    # and then creates a report containing detailed info about each icf used and all the elemental abundances computed
    
    icf = pn.ICF()
    elem_abun = icf.getElemAbundance(atom_abun,icf_list=['direct_O.23',
                                                        'direct_S.23',
                                                        'direct_Ne.345',
                                                        'direct_Cl.23',
                                                        'direct_Cl.34',
                                                        'direct_Cl.234',
                                                        'PHCD07_12',
                                                        'PHCD07_13',
                                                        'PHCD07_16',
                                                        'PHCD07_17',
                                                        'direct_Ar.23',
                                                        'direct_Ar.234',
                                                        'direct_Ar.345',
                                                        'KB94_A30.0',
                                                        'KB94_A30.10',
                                                        'KB94_A30.10b', 
                                                        'KB94_A32', 
                                                        'KH01_4g', 
                                                        'KH01_4txt', 
                                                        'DIMS14_35', 
                                                        'DIMS14_36',
                                                        'KH01_4f',
                                                        'DIMS14_29', 
                                                        'DIMS14_29b', 
                                                        'DIMS14_32',
                                                        'direct_He.23', 
                                                        'KH01_4a', 
                                                        'DIMS14_10',
                                                        'direct_N.23', 
                                                        'direct_N.234', 
                                                        'direct_N.2345',
                                                        'KB94_A0', 
                                                        'KB94_A1.6', 
                                                        'KB94_A1.8', 
                                                        'KB94_A1.10', 
                                                        'KB94_A1.10b',
                                                        'KH01_4c', 
                                                        'DIMS14_14', 
                                                        'DIMS14_14b',
                                                        'direct_Ne.23', 
                                                        'direct_Ne.235', 
                                                        'direct_Ne.2345', 
                                                        'direct_Ne.2356',
                                                        'KB94_A27',
                                                        'KB94_A28.6', 
                                                        'KB94_A28.8', 
                                                        'KB94_A28.10', 
                                                        'KB94_A28.10b', 
                                                        'KH01_4d',
                                                        'DIMS14_17a',
                                                        'DIMS14_17b', 
                                                        'DIMS14_17c', 
                                                        'DIMS14_20', 
                                                        'direct_O.234', 
                                                        'direct_O.2345',
                                                        'KB94_A6',
                                                        'KB94_A8', 
                                                        'KB94_A10', 
                                                        'KB94_A10b', 
                                                        'KH01_4b',
                                                        'DIMS14_12',
                                                        'direct_S.234', 
                                                        'direct_S.2345', 
                                                        'KB94_A36.6',
                                                        'KB94_A36.8', 
                                                        'KB94_A36.10', 
                                                        'KB94_A36.10b', 
                                                        'KB94_A38.6', 
                                                        'KB94_A38.8', 
                                                        'KB94_A38.10', 
                                                        'KH01_4e', 
                                                        'DIMS14_23', 
                                                        'DIMS14_26'])

    with open('Elem_abun_report.csv','w+') as file:
        file.write('Ion,Abundance,Name,Expression,Refrence\n')
        for atom,abun in elem_abun.items():
            if np.isnan(abun):
                continue
            else:
                file.write('{},{:.3E},{},{},{}\n'.format(icf.all_icfs[atom]['elem'], abun, atom, icf.getExpression(atom), icf.getReference(atom)))

if __name__ == '__main__':
    print(Compute_obs(filename,newfile))
