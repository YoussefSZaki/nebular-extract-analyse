#!/usr/bin/env python3

import sys
import csv
import re
import pyneb as pn
import numpy as np

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
    
    with open('Corrected Data.txt','w+') as file:
        for line in obs.getSortedLines():
            file.write('{}   {}\n'.format(line.label, line.corrIntens[0]))
    
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
    density_sensitive=['[ClIII] 5538/5518','[ArIV] 4740/4711','[SII] 6731/6716','[OII] 3726/3729']
    temprature_sensitive=['[NII] 5755/6548','[ArIII] 5192/7300+','[OIII] 4363/5007+','[SII] 4072+/6720+','[OI] 5577/6300+']
    
    # We will generate a report containing the obtained results from all the diagnostics 
    # and the values of cHbeta and E(B-V) as we compute Te & Ne
    # We will also compute the mean values for Te and Ne 
    
    Tmean= []
    Nmean= []
    with open('diag_report.txt','w+') as file:
        for x in density_sensitive:
            file.write('Using {} as the density sensitive line ratio:\n'.format(x))
            for y in temprature_sensitive:
                Te, Ne = diags.getCrossTemDen(y,x,obs=obs)
                Tmean.append(Te)
                Nmean.append(Ne)
                file.write(' Using {} as the temperature sensitive line ratio: Te = {} K, Ne = {} cm-1\n'.format(y,Te,Ne))
        file.write('The value of cHbeta = {}, The value of E(B-V) = {}'.format(rc.cHbeta,rc.E_BV))
        
    #To compute the mean values for Te and Ne
    
    Te_mean = np.mean(Tmean)
    Ne_mean = np.mean(Nmean)
    
    # We will now compute all the ionic abundances for our observations
        
    all_atoms = pn.getAtomDict(atom_list=obs.getUniqueAtoms())
    line_ab = {}
    ion_ab = {}
    temp = Te_mean
    dens = Ne_mean
    for line in obs.getSortedLines():
        if line.atom != 'H1' and line.atom != 'He1' and line.atom != 'He2':
            line_ab[line.label] = all_atoms[line.atom].getIonAbundance(line.corrIntens, temp, dens, to_eval=line.to_eval)
            if line.atom not in ion_ab:
                ion_ab[line.atom] = []
            ion_ab[line.atom].append(line_ab[line.label][0])
            
    # Now creating a report containing the ions and the mean value from all the lines observed for each ion
    
    with open('Ion_abun_report.txt','w+') as file:
        atom_abun={}
        for atom in ion_ab:
            mean = np.mean(np.asarray(ion_ab[atom]))
            atom_abun[atom] = mean
            file.write('{} : {}\n'.format(atom, mean))

    
    # performs all the possible element abundance computations given the input ionic abundance set and available icf  
    # and then creates a report containing detailed info about each icf used and all the elemental abundances computed
    
    icf = pn.ICF()
    elem_abun = icf.getElemAbundance(atom_abun) 
    
    with open('Elem_abun_report.txt','w+') as file:
        for atom,abun in elem_abun.items():
            if np.isnan(abun):
                continue
            if atom == 'He':
                file.write('{} {} = {} using {}\n'.format(atom, icf.all_icfs[atom]['elem'],abun,icf.getExpression(atom)))

            else:
                file.write('{} {} = {} using {}\n'.format(atom, icf.all_icfs[atom]['elem'],abun,icf.getExpression(atom)))

if __name__ == '__main__':
    print(Compute_obs(filename,newfile))
