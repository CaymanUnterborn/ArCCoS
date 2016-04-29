
import numpy as np
import os
import matplotlib.pyplot as plt
from math import log10, sqrt
import math
import time
from math import log as log
from math import exp as exp
import bisect
from collections import Counter
#hack to allow scripts to be placed in subdirectories next to condensation:
if not os.path.exists('condensation') and os.path.exists('../condensation'):
    sys.path.insert(1,os.path.abspath('..'))
import operator

from condensation import fun

import itertools
import re
from scipy.optimize import root,brentq



fugacities = ['H','Cl','F','O','N']
#fugacities = ['x','y']
def get_equations(condensing_solids):
    Names = condensing_solids
    X_dict = {}
    Total = []
    Element_dict = {}
    Molecule_dict={}

    for i in Names:
        Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', i)
        X_dict.update(dict(Mol_name))
        No_numbers = []
        for j in Mol_name:
            No_numbers.append(j[0])
            Total.append(j[0])
        Molecule_dict.update({i : Mol_name})
    X = []

    for i,j  in dict.iteritems(X_dict):
        X.append(i)
        Element_dict.update({i:[]})

    for i in Names:
        No_numbers=[]
        Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', i)
        for j in Mol_name:
            No_numbers.append(j[0])
        for k in No_numbers:
            Z = Element_dict.get(k)
            Z.append(i)
            Element_dict[k] = Z
    Z = Counter(Total)
    for i,j in dict.iteritems(Molecule_dict):
        entry=[]
        for k in j:
            entry.append((k[0],float(k[1])))
        Molecule_dict[i]=entry


    return Molecule_dict,X,Element_dict,Z

def get_abundance(solar_filename,abundance_filename,type):
    Abundance_path = open("Data/Stellar/abundance_"+abundance_filename+".dat","rU")

    Abundance_rows = Abundance_path.readlines()
    element_names = []
    element_multipliers = []

    for i in Abundance_rows:
        split = i.split('\t')
        element_names.append(split[0])
        element_multipliers.append(float(split[1].rstrip('\n')))
    Abundance_path.close()
    abundance_dict  = dict(zip(element_names,element_multipliers))
    Solar_path = open("Data/Solar/abundance_"+solar_filename+".dat","rU")
    Solar_rows = Solar_path.readlines()
    name = []
    abundance =[]
    if type == 'absolute':
        for aRow in Solar_rows:
            values = aRow.split('\t')

            name.append(values[0])
            solar_abun = float(values[1].rstrip('\n'))

            abun_multiplier = abundance_dict.get(values[0])
            abundance.append(abun_multiplier*solar_abun)
        Solar_path.close

        return name, abundance
    if type == 'relative':
        for aRow in Solar_rows:
            values = aRow.split('\t')

            name.append(values[0])
            solar_abun = float(values[1].rstrip('\n'))
            abun_multiplier = pow(10.,abundance_dict.get(values[0]))
            abundance.append(abun_multiplier*solar_abun*1.e10)
        return name, abundance
    return 0
def get_ratio(Abundance_dict,T,P_tot,Name):
    Sum = 0.
    R = 8.3144621e-2
    ratio = {}
    H = Abundance_dict['H']
    He = Abundance_dict['He']
    Ne = Abundance_dict['Ne']
    Ar = Abundance_dict['Ar']
    A =  2.9836
    B = -11452.
    K = 10.**(get_K_gasses('H1',T))
    X = (pow(K,2.))
    H2_coef_a = 4.+X
    H2_coef_b = -((4.*H)+(X*H))
    H2_coef_c = pow(H,2.)

    H2 = ((-H2_coef_b)-np.sqrt(pow(H2_coef_b,2.)-(4.*H2_coef_a*H2_coef_c)))/(2.*H2_coef_a)
    monoH= H - 2.*H2

    Sum = monoH+H2+He+Ne+Ar

    counter = 0
    for i,j in dict.iteritems(Abundance_dict):
        ratio.update({i:j/Sum})
        counter+=1

    return ratio

def get_ParPres(Abun_norm_dict,P_tot,Temp,Name):
    R = 8.3144621e-2
    Par_Pres = {}
    RT = R*Temp
    counter=0
    #print "Abun",Abun_norm_dict
    for i,j in dict.iteritems(Abun_norm_dict):
        Par_Pres.update({i:j*(P_tot/(RT))})
        counter+=1
    #print "par",Par_Pres
    return Par_Pres

def make_list_solids(Name,filename):

    Master_List=[]
    counter = 0
    Master_dict={}
    entry_dict = {}

    path = filename

    for j in path:
        entry = []
        counter = 0
        real_temp = []
        if not j.startswith("#"):
            temp=j.split('\t')
            temp[-1] = temp[-1].strip()

            Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', temp[0])
            entry = []
            for i in temp:
                if temp.index(i) == 0:
                    entry.append(i)
                else:
                    entry.append(float(i))
            counter = 0
            bool = False
            for k in Name:
                for i in Mol_name:
                    if i[0] == k:
                        entry.append(float(i[1]))
                        bool = True
                if bool == False:
                    entry.append(0.)
                counter +=1
                bool = False

            Master_List.append(entry)
    for i in Master_List:
        entry_dict = {}
        counter = 0
        for j in i[1:]:
            if j > 0:
                entry_dict.update({Name[counter]:j})
            counter+=1
        Master_dict.update({i[0]:entry_dict})
    return Master_dict
def make_list(list):
    return 0

def make_list_gas(Name,filename):

    Master_List=[]
    Master_dict = {}
    entry_dict = {}
    counter = 0

    path = filename

    for j in path:

        entry = []
        counter = 0
        real_temp = []
        if not j.startswith("#"):
            temp=j.split('\t')
            temp[-1] = temp[-1].strip()

            Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', temp[0])
            entry = []
            for i in temp:
                if temp.index(i) == 0:
                    entry.append(i)
                else:
                    entry.append(float(i))
            counter = 0
            bool = False
            for k in Name:
                for i in Mol_name:
                    if i[0] == k:
                        entry.append(float(i[1]))
                        bool = True
                if bool == False:
                    entry.append(0.)
                counter +=1
                bool = False
            Master_List.append(entry)


    for i in Master_List:
        entry_dict = {}
        counter = 0
        for j in i[1:]:
            if j > 0:
                entry_dict.update({Name[counter]:j})
            counter+=1
        Master_dict.update({i[0]:entry_dict})

    string = "Name	K   "+'	'.join(Name)
    """
    for i in Name:
        output = []
        for j in Master_List:
            if j[Name.index(i)+1] !=0:
                output.append(j)

        output_filename = "Data/"+i+".dat"
        list = []
        for x in output:
            #print x
            line = tuple(x)
            #print line
            list.append(line)
        #for x in output:
            #output_filename.write(x)
        #print out

        np.savetxt(output_filename,list,"%s","\t",newline='\n', header=string, footer='', comments='# ')
        #f = open(output_filename, 'wb')
        #data = zip(output)
        #np.savetxt(f, data, fmt='%.10e', delimiter='\t')
    """
    return Master_dict

def FindK(molecule,T):
    K = pow(10.,(molecule[0] + (molecule[1]/T)))

    return K



def get_K_data(solid,T):
    molecule=[]
    path = open("Data/Solids/Solids_Cp.dat","rU")
    for temp in path:
        if not temp.startswith("#"):
            values = temp.split('\t')
            if values[0]+"_s"==solid :
                molecule = values
    Name = molecule[0]
    Method = molecule[1]

    Temp_ref = []
    del_H_div_T = []
    S_ref = []
    H_ref = []
    K_ref = []
    del_H_ref = []
    del_G_ref = []

    if Method != 'J' and Method != 'ref':

        del_H = float(molecule[2])
        S0 = float(molecule[3])
        k0 = float(molecule[4])
        k1 = float(molecule[5])
        k2 = float(molecule[6])
        k3 = float(molecule[7])


    if Method == 'ref':
        return 0.

    if Method == 'B8' or Method == 'B5' or Method=='JB':
        C0 = k0*((T-298.15)-(T*(log(T)-log(298.15))))
        C1 = 2.*k1*(((T**0.5)-(298.15**0.5))+(T*(((T**-0.5))-((298.15**-0.5)))))
        C2 = -k2*(((T**-1.)-(298.15**-1.))-((T/2.)*(((T**-2.)-(298.15**-2.)))))
        C3 = -k3*((((T**-2.)-(298.15**-2.))/2.)-((T/3.)*(((T**-3.)-(298.15**-3.)))))


    if Method == 'B3':
        C0 = k0*((T-298.15)-(T*(log(T)-log(298.15))))
        C1 = -k1*(((T**-1.)-(298.15**-1.))-((T/2.)*(((T**-2.)-(298.15**-2.)))))
        C2 = 2.*k2*(((T**0.5)-(298.15**0.5))+(T*(((T**-0.5))-((298.15**-0.5)))))
        C3 = k3*((log(T)-log(298.15))+(T*((T**-1.)-(298.15**-1.))))


    if Method == 'R':
        C0 = k0*((T-298.15)-(T*(log(T)-log(298.15))))
        C1 = k1*((0.5*((T**2)-(298.15**2)))-(T*(T-298.15)))
        C2 = 2.*k2*(((T**0.5)-(298.15**0.5))+(T*(((T**-0.5))-((298.15**-0.5)))))
        C3 = -k3*(((T**-1.)-(298.15**-1.))-((T/2.)*(((T**-2.)-(298.15**-2.)))))

    if Method == 'J':
        path = open("Data/Janaf/"+Name+".dat","rU")

        for aRow in path:
            values = aRow.split('\t')
            if not aRow.startswith('#'):
                Temp_ref.append(float(values[0]))
                del_H_ref.append(float(values[4]))
                H_ref.append(float(values[5]))
                S_ref.append(float(values[2]))
                del_G_ref.append(float(values[6]))
                K_ref.append(float(values[7]))

        if T > Temp_ref[-1]:
            del_H = float(molecule[2])
            S0 = float(molecule[3])
            k0= float(molecule[4])
            k1 = float(molecule[5])
            k2 = float(molecule[6])
            k3 = float(molecule[7])
            C0 = k0*((T-298.15)-(T*(log(T)-log(298.15))))
            C1 = 2.*k1*(((T**0.5)-(298.15**0.5))+(T*(((T**-0.5))-((298.15**-0.5)))))
            C2 = -k2*(((T**-1.)-(298.15**-1.))-((T/2.)*(((T**-2.)-(298.15**-2.)))))
            C3 = -k3*((((T**-2.)-(298.15**-2.))/2.)-((T/3.)*(((T**-3.)-(298.15**-3.)))))
        else:
            K = lookup_and_interpolate(Temp_ref,K_ref,T)
            return K


    Mol_name_split = re.findall(r'([A-Z][a-z]*)(\d*)', Name)
    reference_correction = 0.

    for i in Mol_name_split:
        Temp_ref = []
        del_H_ref = []
        S_ref = []
        G_ref =[]
        H_ref = []


        path = open("Data/Reference/"+i[0]+"_ref.dat","rU")
        for aRow in path:
            if not aRow.startswith("#"):
                values = aRow.split('\t')
                Temp_ref.append(float(values[0]))
                del_H_ref.append(float(values[4]))
                H_ref.append(float(values[5]))
                S_ref.append(float(values[2]))

        H0 = lookup_and_interpolate(Temp_ref,H_ref,298.15)*1000.
        change_H = lookup_and_interpolate(Temp_ref,del_H_ref,T)*1000.
        change_S = lookup_and_interpolate(Temp_ref,S_ref,T)
        #print "change_S",i[0],i[1],change_S
        if i[0] in fugacities:
            reference_correction += (float(i[1])/2.)*(H0+change_H-(T*(change_S)))
        else:
            reference_correction += (float(i[1]))*(H0+change_H-(T*(change_S)))


    del_G_formation_solid = del_H -(T*S0)+C0+C1+C2+C3



    #print "ref",reference_correction/1e3
    #print "solid",del_G_formation_solid/1e3
    #print del_G_formation_solid/1e3
    #print reference_correction/1e3


    del_G_formation_reaction= del_G_formation_solid-reference_correction

    K = log10(exp(1))*((-del_G_formation_reaction/(8.3144621*T)))


    return K


def get_K_gasses(mol_name,T):
    path = open("Data/Gasses.dat","\rU")
    for temp in path:
        if not temp.startswith("#"):
            values = temp.split('\t')
            if values[0]==mol_name :
                molecule = values
    Name = molecule[0]
    Method = molecule[1]
    Temp_ref = []
    del_G_ref = []
    del_H_ref = []
    S_ref = []
    H_ref = []

    if Method == 'J':
        path_molecule = open("Data/Gasses/"+Name+".dat","\rU")
        for aRow in path_molecule:
            values = aRow.split('\t')
            if not aRow.startswith('#'):
                Temp_ref.append(float(values[0]))
                del_H_ref.append(float(values[4]))
                H_ref.append(float(values[5]))
                S_ref.append(float(values[2]))
                del_G_ref.append(float(values[6]))
        H0 = lookup_and_interpolate(Temp_ref,H_ref,298.15)*1000.
        change_H = lookup_and_interpolate(Temp_ref,del_H_ref,T)*1000.
        change_S = lookup_and_interpolate(Temp_ref,S_ref,T)

        #del_G_gas = H0+change_H-(T*(change_S))
        del_G_gas = lookup_and_interpolate(Temp_ref,del_G_ref,T)*1000

        if T > Temp_ref[-1]:
            return 0
        else:
            K=((-del_G_gas/(8.3144621*T))*log10(exp(1.)))
        #if mol_name =='H2O1':

            #print K

        return K

    if Method == 'K':
        H0 = float(molecule[2])
        S0 = float(molecule[3])
        A = float(molecule[4])
        B = float(molecule[5])*1e-3
        C = float(molecule[6])*1e6

        delH = H0 + A*(T-298.15) + ((B/2.)*((T**2)-(298.15)**2)) - C*((1./T)-(1./298.15))
        S = S0 + A*(log(T)-log(298.15)) + B*(T-298.15) - (C/2.)*(((T**-2)-(298.15)**-2))

        delG = delH - T*S
        Mol_name_split = re.findall(r'([A-Z][a-z]*)(\d*)', Name)
        reference_correction = 0.

        for i in Mol_name_split:
            Temp_ref = []
            del_H_ref = []
            S_ref = []
            G_ref =[]
            H_ref = []


            path = open("Data/Reference/"+i[0]+"_ref.dat","rU")
            for aRow in path:
                if not aRow.startswith("#"):
                    values = aRow.split('\t')
                    Temp_ref.append(float(values[0]))
                    del_H_ref.append(float(values[4]))
                    H_ref.append(float(values[5]))
                    S_ref.append(float(values[2]))

            H0 = lookup_and_interpolate(Temp_ref,H_ref,298.15)*1000.
            change_H = lookup_and_interpolate(Temp_ref,del_H_ref,T)*1000.
            change_S = lookup_and_interpolate(Temp_ref,S_ref,T)
            if i[0] in fugacities:
                reference_correction += (float(i[1])/2.)*(H0+change_H-(T*(change_S)))
            else:
                reference_correction += (float(i[1]))*(H0+change_H-(T*(change_S)))

        del_G_corrected = delG - reference_correction

        K = ((-del_G_corrected/(8.3144621*T))*log10(exp(1.)))

        return K

    if Method == 'P':
        path_molecule = open("Data/Gasses/"+Name+".dat","\rU")
        for aRow in path_molecule:
            values = aRow.split('\t')
            if not aRow.startswith('#'):
                Temp_ref.append(float(values[0]))
                del_G_ref.append(float(values[1]))
        G0 = lookup_and_interpolate(Temp_ref,del_G_ref,T)

        del_G_gas = lookup_and_interpolate(Temp_ref,del_G_ref,T)


        Mol_name_split = re.findall(r'([A-Z][a-z]*)(\d*)', Name)
        reference_correction = 0.

        for i in Mol_name_split:
            Temp_ref = []
            del_H_ref = []
            S_ref = []
            G_ref =[]
            H_ref = []


            path = open("Data/Reference/"+i[0]+"_ref.dat","rU")
            for aRow in path:
                if not aRow.startswith("#"):
                    values = aRow.split('\t')
                    Temp_ref.append(float(values[0]))
                    del_H_ref.append(float(values[4]))
                    H_ref.append(float(values[5]))
                    S_ref.append(float(values[2]))

            H0 = lookup_and_interpolate(Temp_ref,H_ref,298.15)*1000.
            change_H = lookup_and_interpolate(Temp_ref,del_H_ref,T)*1000.
            change_S = lookup_and_interpolate(Temp_ref,S_ref,T)
            if i[0] in fugacities:
                reference_correction += (float(i[1])/2.)*(H0+change_H-(T*(change_S)))
            else:
                reference_correction += (float(i[1]))*(H0+change_H-(T*(change_S)))
        del_G_corrected = del_G_gas - reference_correction
        K = ((-del_G_corrected/(8.3144621*T))*log10(exp(1.)))


        return K


def lookup_and_interpolate(table_x, table_y, x_value):

    idx = bisect.bisect_left(table_x, x_value) - 1
    if (idx < 0):
        return table_y[0]
    elif (idx < len(table_x)-1):
        return linear_interpol(x_value, table_x[idx], table_x[idx+1],table_y[idx], table_y[idx+1])
    else:
        return table_y[idx]



def get_K(gasses,solids,T):
    K_dict = {}
    K_gas = {}
    for i,j in dict.iteritems(gasses):
        gas_K=(get_K_gasses(i,T))
        K_gas.update({i:gas_K})
    K_solids = {}
    for i,j in dict.iteritems(solids):
        solid_K = get_K_data(i,T)
        K_solids.update({i:solid_K})
    K_dict = K_gas.copy()
    K_dict.update(K_solids)
    return K_dict

def get_gas_fractions(n_x,Element_dict,K_dict,gasses,RT,Par_Pres):
    gas_number_dict_by_element = {}
    gas_number_dict_by_molecule = {}


    for i,j in dict.iteritems(Element_dict):
        molecules = {}
        for k in j:
            number = 0.
            Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', k)
            for l in Mol_name:
                if l[0] == i:
                    if l[0] in fugacities:
                       number += ((float(l[1])/2.)*(log10(n_x.get(l[0]))+log10(RT)))
                    else:
                        number += (float(l[1]))*(log10(n_x.get(l[0]))+log10(RT))+(log10(float(l[1])))
                else:
                    if l[0] in fugacities:
                       number += (float(l[1])/2.)*(log10(n_x.get(l[0]))+log10(RT))
                    else:
                        number += (float(l[1]))*(log10(n_x.get(l[0]))+log10(RT))
            number += K_dict.get(k) - log10(RT)
            molecules.update({k:(10.**number)})
        gas_number_dict_by_element.update({i:molecules})

    for i in gasses:
        Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', i)
        number = 0.
        for l in Mol_name:
                if l[0] in fugacities:
                   number += (float(l[1])/2.)*(log10(n_x.get(l[0]))+log10(RT))
                else:
                    number += (float(l[1]))*(log10(n_x.get(l[0]))+log10(RT))
        number += K_dict.get(k) - log10(RT)
        gas_number_dict_by_molecule.update({i:10.**number})



    return gas_number_dict_by_element,gas_number_dict_by_molecule



def check_out(condensing_solids,num_solids,num_solids_old,T,T_old,condensed_solids_old,condensed_solids):
    Temp_out = {}
    dT = abs(T_old-T)
    #print "out loop", condensed_solids_old
    for i in condensing_solids:
        #if i == 'Ca1Ti1O3_s':
        #    print "perov ", T, "%.2f" % (num_solids.get(i)/1.e-10)
        j = num_solids.get(i)
        old = num_solids_old.get(i)
        if old == None:
            old = -1.

        if j <= 1.e-10 and old >= 1.e-10:


            Temp_out.update({i:get_out_temp(T,T_old,num_solids,num_solids_old,dT,i)})
            """
            T = Temp_out
            print
            condensing_solids.remove(i)
            del num_solids[i]
            Name.remove(i)
            del guess_dict[i]
            guess = []
            for z in Name:
                for x,y in dict.iteritems(guess_dict):
                    if x == z:
                        if x in condensing_solids:
                            guess.append(1.e-7*sum(Par_Pres))
                        else:
                            guess.append(y)
            Abun_norm = get_data.get_ratio(Abundance_dict,T,P_tot,Name)

            for i,j in dict.iteritems(Abundance_dict):
                if i not in condensing_solids:
                    Abun_norm_dict.update({i:Abun_norm[Abundance.index(j)]})
            Par_Pres = get_data.get_ParPres(Abun_norm_dict,P_tot,T,Name)


            for i,j in dict.iteritems(Abundance_dict):
                if i not in condensing_solids:
                    Par_Pres_dict.update({i:Par_Pres[Abundance.index(j)]})

            K_dict = get_data.get_K(gasses,solids,T)
            guess_dict = dict(zip(Name,guess))
            args = (Element_dict,K_dict,Par_Pres_dict,gasses,Name,T,condensing_solids)
            gas_activity = root(fun.Mass_balance_fun,guess,args=args,method='lm')
            n_x = dict(zip(Name, gas_activity.x))
            errors = fun.Mass_balance_fun(gas_activity.x,*args)
            error_dict = dict(zip(Name,errors))

            num_solids = {}
            for q in condensing_solids:
                num_solids.update({q:n_x.get(q)/sum(Par_Pres)})
            """
    out_temp = 0.
    for i,j in dict.iteritems(Temp_out):
        if j > out_temp:
            out_temp = j
            name = i

    if out_temp>0:
        return name,out_temp
    else:
        return False,False
def get_out_temp(T,T_old,num_solids,num_solids_old,dT,i):
    y_then = log10(num_solids_old.get(i))
    y_now = log10(num_solids.get(i))
    x_now = T
    x_then = T_old

    test_point = T
    diff_int = y_now
    divisor = 1.
    sign = 1.
    diff_int_old = y_then
    #print "out ",i, T, T_old

    #print "y_down",y_now,"at",T
    #print "y_up",y_then,"at",T_old
    test_point = brentq(linear_interpol_out,x_now,x_then,args=(x_now,x_then,y_now,y_then),xtol=0.0000000001)
    
    return test_point
def check_in(solids,n_x,T,K_dict,condensing_solids,T_old,K_dict_old,n_x_old,removed_solids_old,removed_solids,InOut,wrong_solid):

    R = 8.3144621e-2

    new_solids = []

    for i,j in dict.iteritems(solids):
        cool = True
        Pressure_sum=0.
        diff_sum = 0.
        for x,y in dict.iteritems(j):
            if x in fugacities:
                Pressure_sum += (y/2.)*(log10(n_x.get(x)*R*T))
            else:
                Pressure_sum += y*(log10(n_x.get(x)*R*T))

        diff = (-K_dict.get(i)-Pressure_sum)



        if i not in condensing_solids and T != T_old:

            if diff < 0.:
                y_now = diff
                x_now = T
                x_then = T_old

                Pressure_sum_old=0.

                for x,y in dict.iteritems(j):
                    if x in fugacities:
                        Pressure_sum_old += (y/2.)*(log10(n_x_old.get(x)*R*(T_old)))
                        #print "oxy",x,y
                    else:
                        Pressure_sum_old += (y)*(log10(n_x_old.get(x)*R*(T_old)))
                        #print "not oxy",x,y
                y_then = (-K_dict_old.get(i)-Pressure_sum_old)
                #if i == 'Mg2Si1O4_s':
                #    print "OLIVINE then", y_then
                #    print "Olivine Now",y_now
                #-1.e-5
                if y_then <= 0.:
                    #print "in this weird loop"
                    #print removed_solids_old,T
                    if i not in removed_solids_old and i not in removed_solids:
                        #print "Okay to come in:",removed_solids_old,removed_solids,i
                        new_solids.append([i,T_old])
                    else:
                        new_solids.append([i,T_old])
                    """
                    if y_now > y_then:
                        print i,"weird"
                        print x_now,x_then,y_now,y_then
                        test_point = linear_extrapolate(x_now,x_then,y_now,y_then)
                        if test_point >= T:
                            new_solids.append([i,test_point])
                            cool = False
                    else:
                        print "error"

                        print x_now,x_then,y_now,y_then
                        #print i,x_now,x_then,y_now,y_then
                        #print test_point
                        test_point = linear_extrapolate(x_now,x_then,y_now,y_then)
                        new_solids.append([i,test_point])
                    """
                    cool = False

                if cool == True:
                    test_point = brentq(linear_interpol_in,x_now,x_then,args=(x_now,x_then,y_now,y_then),xtol=0.0000000001)
                    """

                    test_point = T
                    diff_int = y_now
                    divisor = 1.
                    sign = 1.
                    diff_int_old = y_then
                    #print "in ",i, T, T_old
                    dT = abs(T-T_old)

                    while diff_int < -1.e-4 or diff_int >= 0.:

                        #print i,x_now-test_point
                        #print i,x_now,test_point,x_then,diff_int
                        #print " ",y_now,"   ",y_then
                        diff_int = linear_interpol(test_point,x_now,x_then,y_now,y_then)

                        if diff_int_old > 0. and diff_int < 0.: #flipped from positive to negative
                            sign = 1.
                            divisor = divisor*10.

                        if diff_int_old < 0. and diff_int > 0.: #flipped from negative to positive
                            sign = -1.
                            divisor = divisor*10.

                        diff_int_old = diff_int
                        test_point += (sign)*dT/divisor
                        if test_point >= x_then:
                            test_point -= (sign)*dT/divisor
                            divisor = divisor*10.
                            test_point += (sign)*dT/divisor

                    """
                    new_solids.append([i,test_point])

    RT = R*T
    #os.system('cls' if os.name == 'nt' else 'clear')
    #print new_solids
    #print "T_old ",T_old
    if len(new_solids)>0:
        Test_temp = 0.

        for a in new_solids:
            #print a,"starting loop and Test temp",Test_temp
            if a[1] > Test_temp:
                if a[1]>T_old:
                    #new_solids.remove(a)
                    #print "new solids in loop", new_solids
                    print a[0], "prevented from coming in"
                else:
                    #print "i in else",a
                    Test_temp = a[1]
                    solid = a[0]
        if len(new_solids)>0 and Test_temp != 0.:
            #print new_solids
            return solid,Test_temp
        else:
            #print "well that's good"
            return False,False
    else:
        return False,False

def linear_extrapolate(x_now,x_then,y_now,y_then):

    m = (y_now-y_then)/(x_now-x_then)

    T = (-y_then/m)+x_then

    return T


def get_num(n_x,coefs,K,RT):
    total=0.

    for x,y in dict.iteritems(coefs):
        if x in fugacities:
            total += (y/2.)*(log10(n_x.get(x)*RT))
        else:
            total += y*(log10(n_x.get(x)*RT))

    total += K - log10(RT)

    number = pow(10.,total)
    return number

def get_num_solid(n_x,coefs,K,RT):
    total=0.

    for x,y in dict.iteritems(coefs):
        if x in fugacities:
            total += (y/2.)*(log10(n_x.get(x)*RT))
        else:
            total += y*(log10(n_x.get(x)*RT))

    total+= -K - log10(RT)

    number = pow(10.,total)
    return number

def linear_interpol_out(x,*args):
    x1 = args[0]
    x2 = args[1]
    y1 = args[2]
    y2 = args[3]
    """
    Linearly interpolate to point x, between
    the points (x1,y1), (x2,y2)
    """
    #print "x1",x1
    #print "x2",x2
    #print "x",x
    assert(x1<=x)
    assert(x2>=x)
    assert(x1<=x2)

    alpha = (x - x1) / (x2-x1)
    return -10.-((1.-alpha)*y1 + alpha*y2)

def linear_interpol_in(x,*args):
    x1 = args[0]
    x2 = args[1]
    y1 = args[2]
    y2 = args[3]
    """
    Linearly interpolate to point x, between
    the points (x1,y1), (x2,y2)
    """
    #print "x1",x1
    #print "x2",x2
    #print "x",x
    assert(x1<=x)
    assert(x2>=x)
    assert(x1<=x2)

    alpha = (x - x1) / (x2-x1)
    return ((1.-alpha)*y1 + alpha*y2)
def linear_interpol(x,x1,x2,y1,y2):
    """
    Linearly interpolate to point x, between
    the points (x1,y1), (x2,y2)
    """
    #print "x1",x1
    #print "x2",x2
    #print "x",x
    assert(x1<=x)
    assert(x2>=x)
    assert(x1<=x2)

    alpha = (x - x1) / (x2-x1)
    return ((1.-alpha)*y1 + alpha*y2)

def get_total_atoms(Element_dict,n_x,gasses,condensing_solids,K,RT,Element_solid_dict,solids):
    n_el = {}
    total = 0.
    for i,j in dict.iteritems(Element_dict):
        interim = 0.
        for k in j:
            coefs = gasses.get(k)
            K_coef = K.get(k)
            interim += (coefs.get(i))*get_num(n_x,gasses.get(k),K_coef,RT)

        n_el.update({i:interim})
    for i,j in dict.iteritems(n_el):
        total += j

    if len(condensing_solids)>0:
        total = 0.
        for i,j in dict.iteritems(Element_solid_dict):
            interim = 0.
            for k in j:
                if k in condensing_solids:
                    coefs = solids.get(k)
                    interim += (coefs.get(i))*n_x.get(k)
            n_el.update({i:n_el.get(i)+interim})
        for i,j in dict.iteritems(n_el):
            total += j

    return total

def get_total_atoms_final(Element_dict,n_x,gasses,condensing_solids,K,RT,Element_solid_dict,solids):
    n_el = {}
    total = 0.
    for i,j in dict.iteritems(Element_dict):
        interim = 0.
        for k in j:
            coefs = gasses.get(k)
            K_coef = K.get(k)
            interim += (coefs.get(i))*get_num(n_x,gasses.get(k),K_coef,RT)

        n_el.update({i:interim})
    for i,j in dict.iteritems(n_el):
        total += j

    if len(condensing_solids)>0:
        total = 0.
        for i,j in dict.iteritems(Element_solid_dict):
            interim = 0.
            for k in j:
                if k in condensing_solids:
                    coefs = solids.get(k)
                    interim += (coefs.get(i))*n_x.get(k)
            n_el.update({i:n_el.get(i)+interim})
        for i,j in dict.iteritems(n_el):
            total += j

    return n_el

def get_total_atoms_solids(Element_dict,n_x,gasses,condensing_solids,K,RT,Element_solid_dict,solids):
    n_elements = get_total_atoms_final(Element_dict,n_x,gasses,condensing_solids,K,RT,Element_solid_dict,solids)

    n_el = {}
    for i,j in dict.iteritems(Element_solid_dict): #open an element
        if len(condensing_solids)>0:

            solids_in_element = {}
            for k in j:
                interim = 0.
                if k in condensing_solids: #if it's a solid
                    coefs = solids.get(k)
                    interim = 100.*(coefs.get(i))*n_x.get(k)/n_elements.get(i)
                solids_in_element.update({k:interim})
            n_el.update({i:solids_in_element})
        else:
            n_el.update({i:0.})



    return n_el