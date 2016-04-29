import os, sys, numpy as np, matplotlib.pyplot as plt
import sympy as sp
from math import log10, sqrt
from scipy.optimize import fsolve, root, brute, fmin, minimize, newton_krylov, anderson, newton, fmin_tnc, broyden1, \
    root
from sympy import Symbol as sym
from decimal import Decimal as dec

if not os.path.exists('Data') and os.path.exists('../Data'):
    sys.path.insert(1, os.path.abspath('..'))
import time
import condensation
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to condensation:
if not os.path.exists('condensation') and os.path.exists('../condensation'):
    sys.path.insert(1, os.path.abspath('..'))
from condensation import get_data, fun, write
import re

fugacities = ['H', 'Cl', 'F', 'O', 'N']
R = 8.3144621e-2

#if __name__ == "__main__":
def arccos(Abundance_filename,Abundance_type,Solar_abun_filename,P_tot,sequence_output_filename,abundance_output_filename):
    # starting temperature
    T_up = 2500.

    # stopping temperature
    T_down = 500.

    # Temperature step, best to leave this about 2 K
    dT = 2.

    # starting total pressure
    #P_tot = 1.e-4
    #Abundance_type = 'relative'
    #Abundance_filename = 'relative'
    #Solar_abun_filename = 'grossman'

    # Matrix of partial pressures for each element. P_X = P_tot * X_i
    Par_Pres = []

    # calls function that opens each individual element's data file. Returns an array of element names and unnormalized abundance
    Name, Abundance = get_data.get_abundance(Solar_abun_filename, Abundance_filename, Abundance_type)
    Abundance_dict = {}
    Abun_norm_dict = {}

    for i in range(len(Name)):
        Abundance_dict.update({Name[i]: Abundance[i]})

    # Calculate the normalized abundance of elements
    Abun_norm_dict = get_data.get_ratio(Abundance_dict, T_up, P_tot, Name)

    # Python likes to switch the order of things and unfortunately order matters in some of these calculations so we
    # rearrange Name to reflect this
    Name = []
    for i, j in dict.iteritems(Abun_norm_dict):
        Name.append(i)

    # list of all molecules considered in
    solids = []
    gasses = []
    Solids_list = []
    Gasses_list = []

    # Function to create master_list of all molecules considered
    path = open("Data/Gasses.dat", "\rU")
    for i in path:
        if not i.startswith("#"):
            values = i.split('\t')
            Gasses_list.append(values[0])

    # list of gasses for each element
    gasses = get_data.make_list_gas(Name, Gasses_list)

    path = open("Data/Solids/Solids_Cp.dat", "\rU")
    for i in path:
        if not i.startswith("#"):
            values = i.split('\t')
            Solids_list.append(values[0] + "_s")

    # list of solids for each element
    solids = get_data.make_list_solids(Name, Solids_list)

    Element_dict = {}
    Element_solid_dict = {}
    for i in Name:
        Element_dict.update({i: []})
        Element_solid_dict.update({i: []})

    for i, j in dict.iteritems(gasses):
        No_numbers = []
        Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', i)
        for j in Mol_name:
            No_numbers.append(j[0])

        for k in No_numbers:
            Z = Element_dict.get(k)
            Z.append(i)
            Element_dict[k] = Z

    for i, j in dict.iteritems(solids):
        No_numbers = []
        Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', i)
        for j in Mol_name:
            No_numbers.append(j[0])
        for k in No_numbers:
            Z = Element_solid_dict.get(k)

            Z.append(i)
            Element_solid_dict[k] = Z

    Par_Pres_dict = {}
    guess_dict = {}

    Par_Pres_dict = get_data.get_ParPres(Abun_norm_dict, P_tot, T_up, Name)

    #initial gueses
    guess_dict = Par_Pres_dict

    # Equation solver must have a guess to begin solving. Same size as Name array.
    guess = []

    # The solver requires a list for the guess, not a dictionary
    for i, j in dict.iteritems(guess_dict):
        guess.append(j)

    output = []
    condensing_solids = []
    K_dict_old = {}
    n_x_old = {}
    output_dict = {}
    got_names = False
    T = T_up
    num_solids_old = {}
    any_in = True
    any_out = True
    out_temp = 0.
    T_old = T_up
    Percent_element_condensed = dict(zip(Abun_norm_dict, [0. for i in Name]))
    solids_names =[]
    numbers = []
    pressures =[]
    temp =[]

    #Percentage condensed for each element
    per_O = []
    per_Al = []
    per_Ca = []
    per_Mg = []
    per_Si = []
    per_Fe = []
    per_Ti = []
    per_Na = []
    per_Ni = []
    per_C = []
    per_N = []
    per_Co = []
    per_Cr = []
    per_K = []
    per_Mn = []

    Output = []
    counter_out = 0
    removed_solids = []
    removed_solids_old = []
    condensed_solids = []
    condensed_solids_old = []
    in_solid = ''
    out_solid = ''
    T_old = T_up
    InOut = False
    Out_loop = False
    num_solid_atom = {}
    num_solid_atom_old = {}

    Redo_solids = False
    elements = ['Ca','Mg','Al','Si','O','Fe','Ti']

    output_dict = {}
    for i in elements:
        output_dict.update({i: []})
    counter=0.
    trigger = False

    #################################################################################
    while T >= T_down:
        Total_elements_condensed = dict(zip(Abun_norm_dict, [0. for i in Name]))

        RT = R * T
        K_gas = {}

        Abun_norm_dict = get_data.get_ratio(Abundance_dict, T, P_tot, Name)

        Par_Pres = []
        Par_Pres_dict = get_data.get_ParPres(Abun_norm_dict, P_tot, T, Name)


        K_dict = get_data.get_K(gasses, solids, T)
        args = (Element_dict, K_dict, Par_Pres_dict, gasses, Name, T, condensing_solids)
        gas_activity = root(fun.Mass_balance_fun, guess, args=args, method='lm',options={'maxiter': 100000000,'ftol':1.e-15})

        n_x = dict(zip(Name, gas_activity.x))
        num_solids = {}
        guess = gas_activity.x
        guess_dict = dict(zip(Name, guess))
        Sum_Par_pres2 = 0.
        for i, j in dict.iteritems(Par_Pres_dict):
            Sum_Par_pres2 += j

        T_lower_bound = T
        Sum_Par_pres = get_data.get_total_atoms(Element_dict, n_x, gasses, condensing_solids, K_dict, RT,
                                                Element_solid_dict, solids)


        for i in condensing_solids:
            num_solids.update({i: n_x.get(i) / Sum_Par_pres})
        num_solid_atom = get_data.get_total_atoms_solids(Element_dict, n_x, gasses, condensing_solids, K_dict, RT,
                                                         Element_solid_dict, solids)



        lst_sqr = 0.
        errors = fun.Mass_balance_fun(gas_activity.x, *args)
        error_dict = dict(zip(Name, errors))
        for i,j in dict.iteritems(error_dict):
            if i not in condensing_solids:
                lst_sqr += pow(j, 2.)

        lst_sqrs = pow(lst_sqr, 0.5)

        while any_in == True and any_out == True:
            in_solid, in_temp = get_data.check_in(solids, n_x, T_lower_bound, K_dict, condensing_solids, T_old,
                                                  K_dict_old, n_x_old, removed_solids_old, removed_solids, InOut,
                                                  in_solid)
            if lst_sqrs > 1.e-13:
                in_solid = False
                in_temp = 0.

            if in_solid == False:
                any_in = False
                in_temp = 0.

            new_solid = [in_solid]
            out_solid, out_temp = get_data.check_out(condensing_solids, num_solids, num_solids_old, T_lower_bound,
                                                     T_old, condensed_solids_old, condensed_solids)
            if out_solid == False:
                any_out = False
                out_temp = 0.

            if in_temp > out_temp:

                print "solid in", in_solid, in_temp
                n_el = get_data.get_total_atoms_final(Element_dict, n_x, gasses, condensing_solids, K_dict, RT,
                                              Element_solid_dict, solids)
                condensed_solids.append(in_solid)
                in_solid_old = in_solid
                any_out = False
                any_in = True

                condensing_solids.append(in_solid)
                elements_in_solid =[]
                #guess = np.append(guess, 1.e-20)

                for i,j in dict.iteritems(solids.get(in_solid)):
                    if i not in elements_in_solid:
                        elements_in_solid.append(i)

                tester = False

                counter =0


                for i in elements_in_solid:
                    if Percent_element_condensed[i] < 0.95 and i != 'O':
                        tester = True
                        break

                if tester == True:
                    if in_solid == 'Fe1_s' or in_solid == 'Ni1_s' or in_solid == 'Co1_s':
                        name_solid = in_solid.strip('1_s')
                        
                        if in_solid == 'Fe1_s':
                            new_guess =  1.e-3*Par_Pres_dict.get(name_solid)
                            guess = np.append(guess, new_guess)
                            

                        if in_solid == 'Ni1_s':
                            new_guess = 0.5e-1*Par_Pres_dict.get(name_solid)
                            guess = np.append(guess,new_guess)
                            
                        if in_solid == 'Co1_s':
                            new_guess = 3.e-1*Par_Pres_dict.get(name_solid)
                            guess = np.append(guess,8.e-14)
                            
                    else:
                        #if in_solid == 'Ca2Mg1Si2O7_s' and T < 1300.:
                            #print "in Ca2Mg1Si2O7_s"
                            #new_guess = 5   .e-7*Sum_Par_pres
                            #guess = np.append(guess,new_guess)
                        #else:
                        if in_solid != 'Mg1Si1O3_ortho_s':
                            if in_solid == 'Ca1Al4O7_s':
                                #new_guess = 3.e-8*Sum_Par_pres
                                guess = np.append(guess,3.e-14)
                            elif in_solid == 'Si1O2_s':
                                new_guess = 1.e-6*Sum_Par_pres
                                guess = np.append(guess,new_guess)
                            elif in_solid == 'Ca2Mg1Si2O7_s':
                                new_guess = 1.e-8*Sum_Par_pres
                                print "here here"
                                guess = np.append(guess,1.e-13)
                            elif in_solid == 'Ca2Al2Si1O7_s':
                            	guess = np.append(guess,8.e-14)
                            else:
                                new_guess = 3.e-7*Sum_Par_pres
                                guess = np.append(guess,5.e-14)

                        else:
                            new_guess = 3.e-6*Sum_Par_pres
                            guess = np.append(guess,2.e-14)
                else:
                    if T < 1300. and 'Ti' in elements_in_solid:

                        if 'Fe' in elements_in_solid:
                            guess = np.append(guess,1.e-5*n_el.get('Fe'))
                        else:
                            
                            guess = np.append(guess,1.e-2*n_el.get('Ti'))

                    elif T > 1300. and 'Ti' in elements_in_solid:
                        guess = np.append(guess, 3.e-1*n_el.get('Ti'))
                    else:
                        guess = np.append(guess,9.e-8*Sum_Par_pres)

                Name.append(in_solid)

                guess_dict = dict(zip(Name, guess))



                elements_in_solid = []
                for i,j in dict.iteritems(solids.get(in_solid)):
                    elements_in_solid.append(i)
                elements_done = []

                guess = []
                for z in Name:
                    for x, y in dict.iteritems(guess_dict):
                        if x == z:
                            if x in condensing_solids and x != in_solid and in_solid != 'Fe1_s'and in_solid != 'Co1_s' and in_solid != 'Ni1_s' and x != 'Fe1_s' and x != 'Ni1_s':
                                guess.append(.88*guess_dict.get(x))
                            else:
                                guess.append(y)


                Output.append(str(in_temp) + " In " + in_solid)
                T_old = T
                K_dict_old = K_dict
                n_x_old = n_x
                num_solids_old = num_solids
                T = in_temp

                Abun_norm_dict = get_data.get_ratio(Abundance_dict, T, P_tot, Name)

                Par_Pres_dict = get_data.get_ParPres(Abun_norm_dict, P_tot, T, Name)

                num_elements = len(guess) - len(condensing_solids)


                guess_dict = dict(zip(Name, guess))

            if out_temp > in_temp:
                print "solid out", out_solid, out_temp
                Output.append(str(out_temp) + " Out " + out_solid)

                any_out = True
                any_in = False
                T = out_temp
                out_solid_old = out_solid
                condensing_solids.remove(out_solid)
                del num_solids[out_solid]
                Name.remove(out_solid)
                del guess_dict[out_solid]
                removed_solids.append(out_solid)
                guess = []
                Abun_norm_dict = get_data.get_ratio(Abundance_dict, T, P_tot, Name)
                Par_Pres_dict = get_data.get_ParPres(Abun_norm_dict, P_tot, T, Name)
 
                Sum_Par_pres = get_data.get_total_atoms(Element_dict, n_x, gasses, condensing_solids, K_dict, RT,
                                                        Element_solid_dict, solids)

                elements_done=[]
                elements_in_solid=[]
                for i,j in dict.iteritems(solids.get(out_solid)):
                    elements_in_solid.append(i)

                for z in Name:
                    for x, y in dict.iteritems(guess_dict):
                        if x == z:
                            if x in condensing_solids:
                                guess.append(y)
                            else:
                                guess.append(y)

            if any_out == True or any_in == True:
                counter_out = 0
                K_dict = get_data.get_K(gasses, solids, T)
                Par_Pres_dict = get_data.get_ParPres(Abun_norm_dict, P_tot, T, Name)
                args = (Element_dict, K_dict, Par_Pres_dict, gasses, Name, T, condensing_solids)

                gas_activity = root(fun.Mass_balance_fun, guess, args=args, method='lm',options={'maxiter': 100000000,'ftol':1.e-15})
                guess = gas_activity.x
                guess_dict = dict(zip(Name, guess))

                errors = []
                error_dict = {}
                lst_sqr = 0.
                errors = fun.Mass_balance_fun(gas_activity.x, *args)
                error_dict = dict(zip(Name, errors))

                for i,j in dict.iteritems(error_dict):
                    if i not in condensing_solids:
                        lst_sqr += pow(j, 2.)

                lst_sqrs = pow(lst_sqr, 0.5)
                print
                #print "after in/out",lst_sqrs
                #print error_dict
                print
                # root(fun.Mass_balance_fun,guess,args=args,method='lm')
                n_x = dict(zip(Name, gas_activity.x))
                Sum_Par_pres = get_data.get_total_atoms(Element_dict, n_x, gasses, condensing_solids, K_dict, RT,
                                                        Element_solid_dict, solids)

                guess = gas_activity.x
                guess_dict = dict(zip(Name, guess))
            num_solids = {}

            for i in condensing_solids:
                num_solids.update({i: n_x.get(i) / Sum_Par_pres})
            num_solid_atom = get_data.get_total_atoms_solids(Element_dict, n_x, gasses, condensing_solids, K_dict, RT,
                                                             Element_solid_dict, solids)

            K_dict_old = K_dict
            n_x_old = n_x
            num_solids_old = num_solids

            guess = gas_activity.x
            guess_dict = dict(zip(Name, guess))
            T_old = T

            entry = [T]
            for i, j in dict.iteritems(solids):
                if got_names == False:
                    solids_names.append(i)
                if n_x.get(i) == None:
                    entry.append(0.)
                else:
                    entry.append(n_x.get(i))
                got_names = True

        numbers.append(entry)
        name_entry = []
        entry = []
        for i, j in dict.iteritems(Par_Pres_dict):
            entry.append(j)
            # name_entry.append(i+'_P_'+str((abs(int(log10(P_tot))))))
            name_entry.append(i + '_P_Lod')

        temp.append(T)
        guess = gas_activity.x
        for i in condensing_solids:
            elements_in_solid = solids.get(i)
            for j, k in dict.iteritems(elements_in_solid):
                Total_elements_condensed[j] = Total_elements_condensed[j] + (k * n_x.get(i))

        for i, j in dict.iteritems(Total_elements_condensed):
            #if j / (Par_Pres_dict[i]) > 1.:
            #    Percent_element_condensed[i] = 1.
            #else:
            Percent_element_condensed[i] = j / (Par_Pres_dict.get(i))

        per_O.append(Percent_element_condensed['O'] * 100.)
        per_Cr.append(Percent_element_condensed['N'] * 100.)
        per_Co.append(Percent_element_condensed['Co'] * 100.)
        per_N.append(Percent_element_condensed['Cr'] * 100.)
        per_Al.append(Percent_element_condensed['Al'] * 100.)
        per_Ca.append(Percent_element_condensed['Ca'] * 100.)
        per_Mg.append(Percent_element_condensed['Mg'] * 100.)
        per_Si.append(Percent_element_condensed['Si'] * 100.)
        per_Fe.append(Percent_element_condensed['Fe'] * 100.)
        per_Ti.append(Percent_element_condensed['Ti'] * 100.)
        per_Na.append(Percent_element_condensed['Na'] * 100.)
        per_Ni.append(Percent_element_condensed['Ni'] * 100.)
        per_C.append(Percent_element_condensed['C'] * 100.)
        per_K.append(Percent_element_condensed['K'] * 100.)
        per_Mn.append(Percent_element_condensed['Mn'] * 100.)



# (counter_out < 40) and
        if T != T_up:
            if (lst_sqrs >1.e-12) and counter_out < 1000:
                if (counter_out%50) == 0:
                	print
                	print "tried: ", counter_out+1, "times"
                	print  "mass balance error:", lst_sqrs, "mol"
                elements_in_solid = []

                #################
                del per_O[-1]
                del per_Al[-1]
                del per_Ca[-1]
                del per_Mg[-1]
                del per_Si[-1]
                del per_Fe[-1]
                del per_Ti[-1]
                del per_Na[-1]
                del per_C[-1]
                del per_Ni[-1]
                del per_N[-1]
                del per_Co[-1]
                del per_Cr[-1]
                del per_K[-1]
                del per_Mn[-1]
                del temp[-1]

                K_dict_old = K_dict
                n_x_old = n_x
                num_solids_old = num_solids
                guess = gas_activity.x
                guess_dict = dict(zip(Name, guess))
                T -= 0.0001

                any_in = False
                any_out = False
                counter_out += 1


                condensed_solids = []

            else:
                #########################


                if len(condensing_solids)>0:
                    print "Temperature:" , "%.2f" % T
                    print "Element", "% Condensed"
                    print "O", "%.2f" % per_O[-1]
                    print "Si", "%.2f" % per_Si[-1]
                    print "Mg", "%.2f" % per_Mg[-1]
                    print "Fe", "%.2f" % per_Fe[-1]
                    print "Ni", "%.2f" % per_Ni[-1]
                    print "Ca", "%.2f" % per_Ca[-1]
                    print "Al", "%.2f" % per_Al[-1]
                    print "Ti", "%.2f" % per_Ti[-1]
                    print "Na", "%.2f" % per_Na[-1]
                    print "C", "%.2f" % per_C[-1]
                    print


                text_file = open("output/sequence/Sequence_" + sequence_output_filename + ".dat", "w")

                for i in Output:
                    text_file.write(i + "\n")
                text_file.close()
                if counter_out > 0:
                    T-= .0001
                    K_dict_old = K_dict
                    n_x_old = n_x
                    num_solids_old = num_solids
                    condensed_solids_old = condensed_solids
                    removed_solids_old = removed_solids
                else:
                    if lst_sqrs > 1.e-14:
                        T-= .0001
                        K_dict_old = K_dict
                        n_x_old = n_x
                        num_solids_old = num_solids
                        condensed_solids_old = condensed_solids
                        removed_solids_old = removed_solids
                    else:
                        output_dict = write.write(T, n_x, gasses, solids, condensing_solids, Element_solid_dict, K_dict,
                                          Abundance_dict, abundance_output_filename, Element_dict, output_dict,Par_Pres_dict)

                        T -= dT

                counter_out = 0.

                any_in = True
                any_out = True

                K_dict_old = K_dict
                n_x_old = n_x
                num_solids_old = num_solids
                condensed_solids_old = condensed_solids
                removed_solids_old = removed_solids

                condensed_solids = []
                removed_solids = []
                Out_loop = False
                Redo_solids = False
                num_solid_atom_old = num_solid_atom
                counter +=1.

                if per_Si[-1] >= 99.99 and per_Mg[-1] >= 99.99 and per_Fe[-1] >= 99.99 and per_Ca[-1] >= 99.99 and per_Al[-1] >=99.99:
                    print "DONE w/ Refractories at",  "%.2f" % temp[-1], " K"
                    T_down = T
                    break
        else:
            T -= dT
            any_in = True
            any_out = True


    print "O", "%.2f" % per_O[-1]
    print "Si", "%.2f" % per_Si[-1]
    print "Mg", "%.2f" % per_Mg[-1]
    print "Fe", "%.2f" % per_Fe[-1]
    print "Ni", "%.2f" % per_Ni[-1]
    print "Ca", "%.2f" % per_Ca[-1]
    print "Al", "%.2f" % per_Al[-1]
    print "Ti", "%.2f" % per_Ti[-1]
    print "Na", "%.2f" % per_Na[-1]
    print "C", "%.2f" % per_C[-1]

    print
    print "C/O", "%.2f" % (Abundance_dict.get('C') / Abundance_dict.get('O'))

    print
    print "Mg/Si = ", "%.2f" % (Abundance_dict.get('Mg') / Abundance_dict.get('Si'))
    print
    print "Mg/Fe = ", "%.2f" % (Abundance_dict.get('Mg') / Abundance_dict.get('Fe'))
    print
    print "Percent Oxygen Abundance Condensed = ", "%.2f" % (per_O[-1])

    Output.append("")
    Output.append("DONE w/ Refractories at " + str(temp[-1]) + " K")
    Output.append("")
    Output.append("Percent O Condensed = " + str(per_O[-1]))



    text_file = open("output/sequence/Sequence_"+sequence_output_filename+".dat", "w")
    for i in Output:
        text_file.write(i + "\n")
    text_file.close()
    print "Condensation Sequence list output to: Condensation_" + sequence_output_filename + ".dat"
	
    plt.plot(temp, per_O, color='b', linestyle='-', marker='o', \
             markerfacecolor='b', markersize=4, label='O')
    plt.plot(temp, per_Al, color='g', linestyle='-', marker='o', \
             markerfacecolor='g', markersize=4, label='Al')
    plt.plot(temp, per_Mg, color='r', linestyle='-', marker='o', \
             markerfacecolor='r', markersize=4, label='Mg')
    plt.plot(temp, per_Ca, color='k', linestyle='-', marker='o', \
             markerfacecolor='k', markersize=4, label='Ca')
    plt.plot(temp, per_Si, color='y', linestyle='-', marker='o', \
             markerfacecolor='y', markersize=4, label='Si')
    plt.plot(temp, per_Fe, color='k', linestyle='-', marker='x', \
             markerfacecolor='k', markersize=4, label='Fe')
    plt.plot(temp, per_Ti, color='b', linestyle='-', marker='x', \
             markerfacecolor='b', markersize=4, label='Ti')
    plt.plot(temp, per_Na, color='r', linestyle='-', marker='x', \
             markerfacecolor='r', markersize=4, label='Na')
    plt.legend(loc='upper right')
    plt.xlabel("Temperature (K)")
    plt.ylabel("% Element Condensed")
    plt.xlim(temp[-1], 1850)
    plt.ylim(0, 100)
    plt.savefig("Ratios_" + Abundance_filename + ".jpg")
    #plt.show()
