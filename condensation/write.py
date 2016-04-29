import os
#hack to allow scripts to be placed in subdirectories next to condensation:
if not os.path.exists('condensation') and os.path.exists('../condensation'):
    sys.path.insert(1,os.path.abspath('..'))
from condensation import get_data
import numpy as np
R = 8.3144621e-2
def write(T,n_x,gasses,solids,condensing_solids,Element_solid_dict,K_dict,Abundance_dict,Abundance_filename,Element_dict,output_dict,Par_Pres_dict):


    elements = ['Ca','Mg','Al','Si','O','Fe','Ti']
    RT = R * T
    for name_for_file in elements:
        final_output = output_dict.get(name_for_file)
        line_name = []
        line_name.append('T')

        for i, j in dict.iteritems(Element_solid_dict):
            if i == name_for_file:
                for k in j:
                    line_name.append(k)
        for i, j in dict.iteritems(Element_dict):
            if i == name_for_file:
                # get_data.get_number(i,n_x)
                for k in j:
                    line_name.append(k)

        line_item = [T]

        ##Correct for coefficients
        baseline = [0]
        for i in line_name:
            if i is not 'T':
                if gasses.get(i) is not None:
                    for j, k in dict.iteritems(gasses.get(i)):
                        if j == name_for_file:
                            line_item.append(k * get_data.get_num(n_x, gasses.get(i), K_dict.get(i), RT))
                else:
                    if i in condensing_solids:
                        elements_in_solid = solids.get(i)
                        for j, k in dict.iteritems(elements_in_solid):
                            if j == name_for_file:
                                line_item.append(k * n_x.get(i))
                    else:
                        line_item.append(0.)

        n_el = get_data.get_total_atoms_final(Element_dict, n_x, gasses, condensing_solids, K_dict, RT,
                                              Element_solid_dict, solids)

        for i in range(len(line_item)):
            if line_item[i] is not T:
                if (line_item[i] / n_el.get(name_for_file)) > 1.e-10:
                    line_item[i] = line_item[i] / Par_Pres_dict.get(name_for_file)
                    #line_item[i] = line_item[i]

                else:
                    line_item[i] = 0.

        line_item_strip = [0] + line_item[1::]

        for i in range(len(line_item_strip)):
            inter = []

            for j in range(len(line_item_strip)):
                if i == 0:
                    inter = []
                else:
                    if j < i - 1.:
                        inter.append(line_item_strip[j])

            baseline.append(sum(inter))

        for i in range(len(line_item)):
            if line_item[i] is not T:
                if line_item[i] > 0.:
                    line_item[i] = line_item[i] + baseline[i + 1]


        final_output.append(line_item)

        string_element = '	'.join(line_name)
        np.savetxt("output/abundance/"+Abundance_filename+"_"+ name_for_file + ".dat", final_output, '%10.10e', "\t", newline='\n',
                   header=string_element, footer='', comments='# ')

        output_dict.update({i:final_output})

    return output_dict

