#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 14:20:02 2022

@author: pappel
"""
def minform_droop(f,w,x,t):  
  total_cations = sum(f[k] for k in f)
  if 'Fe' in f:
      f['Fe2'] = f.pop('Fe')
  if (total_cations > t):    
      f3 = 2 * x * (1 - t / total_cations)  # Droop Formula
      new_form = {k: f[k] * t / total_cations for k in f}
      new_form['Fe2'] = new_form['Fe2'] - f3
      new_form['Fe3'] = f3
      new_wts = w.copy()
      new_wts['FeO'] = w['FeO'] * new_form['Fe2']/(new_form['Fe2'] + new_form['Fe3'])
      new_wts['Fe2O3'] = 1.1113 * w['FeO'] * new_form['Fe3']/(new_form['Fe2'] + new_form['Fe3'])
      formula = new_form
  else:
      new_wts = w
      formula = f
      formula['Fe3'] = float('NaN')
  return(formula, new_wts)
 
    
def minform(wts, ox):
  # only make sure that keys in mmasses, catmult and oxmult are the same
  mmasses = {'SiO2':60.083, 'TiO2': 79.8658, 'Al2O3':101.959, 'FeO':71.844, 'Fe2O3': 159.69, 'MnO':70.9374, 'MgO':40.304, 'CaO':56.077, 'Na2O':61.979, 'K2O':94.196}
  catmult = {'SiO2':1, 'Al2O3':2, 'FeO':1, 'Fe2O3':2, 'MnO':1, 'MgO':1, 'TiO2': 1, 'CaO':1, 'Na2O':2, 'K2O':2}
  oxmult = {'SiO2':2, 'Al2O3':3, 'FeO':1, 'Fe2O3':3, 'MnO':1, 'MgO':1, 'CaO':1, 'TiO2': 2, 'Na2O':1, 'K2O':1}
  key_map_dict = {'Al2O3':'Al', 'SiO2':'Si', 'TiO2':'Ti', 'CaO':'Ca', 'Na2O':'Na', 'FeO':'Fe2','Fe2O3':'Fe3', 'MgO':'Mg', 'MnO':'Mn', 'K2O':'K'}
  moleratios = {k:wts[k]/mmasses[k] for k in wts if k in mmasses}
  molepercat = {k:moleratios[k]*catmult[k] for k in moleratios}
  moleperox = {k:moleratios[k]*oxmult[k] for k in moleratios} 
  normalize_fac = ox/sum(list(moleperox.values()))
  # result keeps a dict with key/values of the stoichiometric coefficients
  result = {k:molepercat[k]*normalize_fac for k in molepercat}
  # replace the oxide form of the keys by the atomic
  formula = {(key_map_dict[k] if k in key_map_dict else k):v for (k, v) in result.items() }
  return(formula)


def glcferric(f, w, t, s):
        

        
    
        

    return(formula, new_wts)


def tcba_assigmnet(formula_data, ferric):
    tcba = {'TSi':0, 'TAlIV':0, 'CAlVi':0, 'CFe3':0, 'CCr':0, 'CMg':0, 'CFe2':0, 'CMn':0, 'BMg':0, 'BFe2':0, 'BMn':0, 'BCa':0, 'BNa':0, 'ANa':0, 'AK':0}
    cb = [['CMg', 'Mg', 'BMg'], ['CFe2', 'Fe2', 'BFe2'],['CMn', 'Mn', 'BMn']]
    
    # T-site assignment 
    if formula_data['Si'] <= 8:
        tcba['TSi'] = formula_data['Si']
        tcba['TAlIV'] = 8-tcba['TSi']
    else:
        print('Si > 8 - Caution!')
        tcba['TSi'] = formula_data['Si']
        tcba['TAlIV'] = 0
        
    # # assigning C
    tcba['CAlVi'] = formula_data['Al']-tcba['TAlIV']
    tcba['CFe3'] = 0 if (ferric == 0) else ferric
        

    """
    if ferric == 0:
        tcba['CFe3'] = 0
    else:
        tcba['CFe3'] = ferric
    """    

    tcba['CCr'] = formula_data['Cr'] if "Cr" in formula_data else 0
    C = tcba['CAlVi']+tcba['CFe3']+tcba['CCr']
    
    # C-Site assignment 
    for c,n,b in cb:
        if (C +formula_data[n]) <= 5:
            tcba[c]= formula_data[n]
            tcba[b]= 0
        elif C<5: 
            tcba[b] = (C+ formula_data[n])-5
            tcba[c] = formula_data[n]-tcba[b]
        else:
            tcba[b] = formula_data[n]
        C = C +tcba[c]

    # B- and A-site assignment 
    tcba['BCa'] = formula_data['Ca'] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    B = tcba['BMn'] + tcba['BFe2']+tcba['BMg'] + tcba['BCa']

    if (B+ formula_data['Na'])<= 2:
        tcba['BNa'] = formula_data['Na']
    elif B < 2:
        tcba['BNa'] = 2-B
        tcba['ANa'] = formula_data['Na']-tcba['BNa']
    else:
        tcba['BNa'] = 0 
        tcba['ANa'] = formula_data['Na']
        
    tcba['AK'] = formula_data['K'] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    return tcba

def amph_ferric_minform(wts):
    ox = 23
    formula = minform(wts, ox)

    tcba = tcba_assigmnet(formula, 0)
    slii = {'Al':1.5, 'Si':2, 'Ti':2, 'Ca':1, 'Na':0.5, 'Fe2':1, 'Mg':1, 'K':0.5, 'Mn':1, 'Cr':1.5}

    # and the formula for the all -ferric case, needed for the maximum estimate
    allferric_wts = wts
    allferric_wts['Fe2O3'] = wts['FeO']*1.1113
    allferric_wts['FeO'] = 0
    allferric_formel = minform(allferric_wts, ox) 

    # formula for min ferric 
    _8Si =8/tcba['TSi']
    _16CAT = 16/sum(tcba.values())
    all_ferrous = 1
    _15eNK = 15/(sum(tcba.values())-tcba['BNa']-tcba['ANa']-tcba['AK'])

    if ((_8Si > 1) and (_16CAT > 1) and (_15eNK > 1)):
        print("no Fe3+ for minimum estimate")
        # ... do something more ...
        # return the all-ferrous formula with Fe3 = 0
        formula['Fe3'] = 0
        return (formula, tcba_assigmnet(formula, 0))
    else:
        minFerric = min(_8Si,_16CAT,all_ferrous,_15eNK)
        minFerricFormula = {k: formula[k]*minFerric for k in formula}
        oxygensPerCations = {k: minFerricFormula[k]*slii[k] for k in formula}
        # calculate the ferric iron and the new ferrous iron of the formula
        minFerricFormula['Fe3'] = (ox - sum(list(oxygensPerCations.values()))) * 2
        minFerricFormula['Fe2'] = minFerricFormula['Fe2'] - minFerricFormula['Fe3']

    # minFerricResult = tcba_assigmnet(minFerricFormula, minFerricFormula['Fe3'])
    # ferricResult = pd.DataFrame([minFerricResult])

    # formula for max and average ferric
    _13eCNK = 13/(sum(tcba.values())-tcba['BNa']-tcba['ANa']-tcba['AK']- tcba['BMn'])
    _15eK = 15/(sum(tcba.values())- tcba['AK'])
    allFerric = ox/(ox+(0.5* formula['Fe2'])) 
    _8SiAl = 8/(tcba['CAlVi']+tcba['TAlIV']+tcba['TSi'])
    maxFerric = max(_13eCNK, _15eK, _8SiAl, allFerric)
    averageFerric = (minFerric + maxFerric)/2   
  
    for x in [maxFerric, averageFerric]:
        ferricFormula = {k: formula[k]*x for k in formula}  
        oxygensPerCations = {k: ferricFormula[k]*slii[k] for k in formula}
        ferricFormula['Fe3'] = (ox - sum(list(oxygensPerCations.values()))) * 2
        ferricFormula['Fe2'] = ferricFormula['Fe2'] - ferricFormula['Fe3']
        
    # ferricResult contains in the following order (minFerric, maxFerric, averageFerric) 
#    return ferricResult
    return (ferricFormula, tcba_assigmnet(ferricFormula, ferricFormula['Fe3']))
