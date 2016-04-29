import os, sys, numpy as np, matplotlib.pyplot as plt
import re

import sympy as sp
from math import log10, sqrt
from scipy.optimize import fsolve,root,brute,fmin,minimize,newton_krylov, anderson, newton,fmin_tnc,broyden1,root
from sympy import Symbol as sym
from decimal import Decimal as dec
if not os.path.exists('Data') and os.path.exists('../Data'):
    sys.path.insert(1,os.path.abspath('..'))
from collections import Counter
import time
import condensation
import get_data
#hack to allow scripts to be placed in subdirectories next to condensation:
if not os.path.exists('condensation') and os.path.exists('../condensation'):
    sys.path.insert(1,os.path.abspath('..'))

fugacities = ['H','Cl','F','O','N']



def Mass_balance_fun(guess_list, *args):
    R = 8.3144621e-2
    #units of L bar mol^-1 K^-1
    Element_dict = args[0]
    K = args[1]
    Par_Pres = args[2]
    gasses = args[3]
    Name = args[4]
    RT = R*args[5]
    condensing_solids = args[6]
    guess_dict = {}

    test = {}
    counter = 0
    for i in Name:
        guess_dict[i]=guess_list[counter]
        counter+=1



    final = [] # Sum of all activity terms for individual element
    outs = 1
    out = []
    final_dict = {}
    finals = []
    coefs ={}
    out_dict = {}
    #print guess_dict
    for i,j in dict.iteritems(Element_dict):
        out = 0.
        for k in j:
            #print k
            coefs = gasses.get(k) #open an individual gas species
            entry=0.
            for x,y in dict.iteritems(coefs):
                if guess_dict[x] <= 0.:
                    entry = 1.e999
                    break
                else:
                    if x == i:
                        if x in fugacities:
                            entry = entry+((y/2.)*log10(guess_dict.get(x)*RT)+log10(y))
                        else:
                            entry = entry+(y*log10(guess_dict.get(x)*RT)+log10(y))
                    else:
                        if x in fugacities:
                            entry = entry+((y/2.)*log10(guess_dict.get(x)*RT))
                        else:
                            entry = entry+(y*log10(guess_dict.get(x)*RT))
            #print k
            entry = entry+(K.get(k))-log10(RT)
            entry = pow(10.,entry)
            out+=entry

        out_dict.update({i:out})
    for i,j in dict.iteritems(Par_Pres):

            final_dict.update({i:log10(j)-log10(out_dict.get(i))})

    if len(condensing_solids)>0:

        solid_molecule_dict,Y,solid_element_dict,Z = get_data.get_equations(condensing_solids)
        for i,j in dict.iteritems(solid_element_dict): #open an element
            Sum=0.
            for k in j: #open molecule that contains that element
                if guess_dict[k]<=0.:
                    Sum = 1.e999
                    break
                else:
                    for x in solid_molecule_dict.get(k):
                        if x[0]==i:
                            if x[0] in fugacities:
                                Sum+=(float(x[1])/2.)*guess_dict[k]
                            else:
                                Sum+=float(x[1])*guess_dict[k]

            final_dict.update({i:1.-((log10(out_dict[i]+Sum)/log10(Par_Pres.get(i))))})




    counter = 0
    for i in Name:
        if i not in condensing_solids:
            finals.append(final_dict.get(i))
    #condensing_solids = ['Al2O3']
    #solid_molecule_dict,Y,solid_element_dict,Z = get_data.get_equations(condensing_solids)
    if len(condensing_solids)>0:

        for i,j in dict.iteritems(solid_molecule_dict):

            pressures = 0.
            for k in j:
                if guess_dict.get(k[0]) < 0.:
                    pressures = 1.e999
                    break
                else:
                    if k[0] in fugacities:
                        pressures += (float(k[1])/2.)*(log10(guess_dict.get(k[0])*RT))
                    else:
                        pressures += (float(k[1]))*(log10(guess_dict.get(k[0])*RT))

            if K.get(i) ==0.:
                finals.append(-K.get(i) - pressures)
            else:
                finals.append((1.-(pressures/-K.get(i))))

    #Check whether solid equilibria even matters...

    return finals

def get_denominator(gasses,K,T,guess_dict):
    R = 8.3144621e-2
    total = 0.

    for i,j in dict.iteritems(gasses): #opens a gas
        K_gas = K.get(i)
        entry = 0.
        for k,l in dict.iteritems(j):
            entry += (l*log10(guess_dict.get(k)*R*T))
        if i == 'H1' or i == 'H2':
            print i, entry
        total += (pow(10.,entry)*K_gas)

    return total




