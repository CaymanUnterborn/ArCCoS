import os
if not os.path.exists('condensation') and os.path.exists('../condensation'):
    sys.path.insert(1, os.path.abspath('..'))
from condensation import arccos
from math import log10

if __name__ == "__main__":

    #input pressure in bars
    P_tot = 1.e-3

    #Input the filename of the abundances, can be relative to Solar ([X/H] notation) or in pure moles
    Abundance_filename = 'Solar'

    #Is this filename relative file name ([X/H] format) or absolute (fraction change relative to solar in moles)?
    Abundance_type = 'relative'

    #Which Solar model would you like to use? asplund09 (Asplund et al., 2009), lodders (Lodders, 2003), 
    #or andersgrev (Anders & Grevesse, 1989)
    Solar_abun_filename = 'asplund09'

    #output filename of appearance/disappearance temperatures text file. No need to include .dat
    sequence_output_filename = Abundance_filename+"_"+Solar_abun_filename+"_"+str(log10(P_tot)).strip('.0')

    #output filename of individual element phase diagram filenames
    abundance_output_filename = Abundance_filename+"_"+Solar_abun_filename+"_"+str(log10(P_tot)).strip('.0')

    arccos.arccos(Abundance_filename,Abundance_type,Solar_abun_filename,P_tot,sequence_output_filename,abundance_output_filename)