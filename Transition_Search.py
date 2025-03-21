import numpy
import matplotlib.pyplot as pyplot
from arc import *
pyplot.rcParams.update(pyplot.rcParamsDefault) #reset matplotlib defaults as ARC overrides them!

def get_all_trans(freq_range, atom_type = 'Cs', ground = 'default',\
min_dme = 5, nmax = 80, ryd_laser_range = (700, 1000), save = True):
    '''
    Creates an array of all possible transitions within a certain frequency range 
    for a specified atom type. This function only considers dipole allowed transitions 
    (\Delta l = +-1) and the value of the dipole matrix element given is for pi transitions.

    Args:
        freq_range (tuple of floats): Frequency range (in GHz) in which to 
            search for possible Rydberg transitions. 

        atom_type (str 'Cs' or 'Rb', optional): The atomic species to consider. 
            Default is Cs (caesium), Rb is rubidium 87.
        ground (tuple of floats, optional): The state from which we excite to 
            Rydberg in the form (n,l,j). Default is `'default'` which corresponds 
            to 7S_{1/2} in Caesium, and 6S_{1/2} for Rubidium.
        min_dme (float, optional): the minimum value of the dipole matrix element 
            (DME) to be stored (in units a_0e). Transitions with DME less than 
            this will not appear in the output array. Default is 5 a_0e.
        nmax (float, optional): the maximum value of principal quantum number n
            that will be considered. Default is 80.
        ryd_laser_range (tuple of floats, optional): wavelength range of the laser 
            used to excite from the specified 'ground' state to the R_1 Rydberg 
            state (in nm). Default is (700,1000). 
        save (bool, optional): whether the output array should be saved. Filename 
            will be "Transition_Search_{frequency_range}GHz_{atom_type}\_from_{ground}.csv". 
            Default is True.
 
    Returns:
        data_arr (2darray): details of all transitions within the specified parameter ranges. 
            In order [n_1, l_1, j_1, n_2, l_2, j_2, Transition frequency (GHz), 
            DME (a_0e), Wavelength from chosen 'ground' state to Rydberg state R_1 (nm), 
            DME of transition from 'ground' state to Rydberg state R_1 ($a_0e$)]'''
    from_freq, to_freq = freq_range[0], freq_range[1]
    atom_dict = {'Cs':Caesium(), 'Rb':Rubidium87()}
    atom = atom_dict[atom_type]
    notation = ['S','P','D','F']

    ## This sets the 'ground state', the state from which we excite to Rydberg
    ground_dict = {'Cs':(7,0,0.5), 'Rb':(6,0,0.5)}
    if ground == 'default':
        ground_state = ground_dict[atom_type]
    else:
        ground_state = ground
    print('Exciting to Rydberg from {}{}{:.0f}/2 (in {})'.format(ground_state[0],\
    notation[int(ground_state[1])], ground_state[2]*2, atom_type))
    
    Rydberg_ns = numpy.arange(ground_state[0],nmax)
    max_delta_n = 10 
    ## This is the maximum change in n that will be considered. 
    ## Only needs to be changed if you are wanting to find very weak transitions.
    delta_ns = numpy.arange(-max_delta_n, max_delta_n +1)
    freqs = []
    dmes = []
    ryd_wvls = []
    ryd_dmes = []
    froms = []
    tos = []

    ## Identify all the possible values of l we can reach from the specified ground state.
    ## Note: only considering dipole allowed transitions (delta l = +- 1)
    ls = [ground_state[1] - 1, ground_state[1] + 1]
    for l in ls:
        if l<0 or l>ground_state[0] - 1:
            ls.remove(l)

    for n in Rydberg_ns:
        for l_from in ls:
            ## Create an array of j values, check that j = l +- 1/2 and delete if not
            js = numpy.arange(ground_state[2] - 1, ground_state[2] + 2, 1)
            for j,j_from in enumerate(js):
                if j_from < l_from - 0.5 or j_from > l_from + 0.5:
                    js = numpy.delete(js, j)
                    
                else:
                    s_from = (n, l_from, j_from)
                    ryd_trans = ground_state + s_from
                    ryd_wvl = abs(atom.getTransitionWavelength(*ryd_trans))
                    if ryd_wvl*1e9>ryd_laser_range[0] and ryd_wvl*1e9<ryd_laser_range[1]:
                        ## If the transition is out of reach of the final step laser, stop
                        for del_n in delta_ns:
                            n_to = n+del_n
                            if n_to<ground_state[0]:
                                n_to = ground_state[0]
                            for delta_l in (-1,1):
                                l_to = abs(l_from+delta_l)
                                if l_to<n_to:
                                    for delta_j in (-1,0,1):
                                        j_to = abs(j_from + delta_j)
                                        if j_to>=l_to-0.5 and j_to<=l_to+0.5:
                                            s_to = (n_to, l_to, j_to)
                                            trans = s_from+s_to
                                            dme = s_from + (-0.5,) +s_to + (-0.5,0)
                                            thz_trans = atom.getTransitionFrequency(*trans)
                                            if from_freq<abs(thz_trans)/1e9<to_freq:
                                                if abs(j_from-j_to)<1.5:
                                                    thz_dme = atom.getDipoleMatrixElement(*dme)
                                                    ryd_trans_dme = ground_state + (-0.5,) + s_from + (-0.5, 0)
                                                    ryd_dme = atom.getDipoleMatrixElement(*ryd_trans_dme)
                                                    if abs(thz_dme)>min_dme:
                                                        if len(dmes)<1:
                                                            freqs.append(abs(thz_trans*1e-9)) #take absolute frequency as direction of transition is known
                                                            dmes.append(thz_dme)
                                                            froms.append(s_from)
                                                            tos.append(s_to)
                                                            ryd_wvls.append(ryd_wvl)
                                                            ryd_dmes.append(ryd_dme)
                                                        elif thz_dme != dmes[-1]:
                                                            freqs.append(abs(thz_trans*1e-9))
                                                            dmes.append(thz_dme)
                                                            froms.append(s_from)
                                                            tos.append(s_to)
                                                            ryd_wvls.append(ryd_wvl)
                                                            ryd_dmes.append(ryd_dme)
    froms_array = numpy.array(froms, dtype='float')
    tos_array = numpy.array(tos, dtype='float')
    freqs_array = numpy.array(freqs)
    dmes_array = numpy.array(numpy.abs(dmes))
    wvls_array = numpy.array(ryd_wvls)*1e9 #convert to nm
    ryd_dmes_array = numpy.asarray(numpy.abs(ryd_dmes), dtype = 'float')
    fromto = numpy.hstack((froms_array, tos_array))
    data_arr = numpy.c_[fromto, freqs_array, dmes_array, wvls_array, ryd_dmes_array]
    if save:
        numpy.savetxt('Transition_Search_{}-{}GHz_{}_from_{}{}{}_2.csv'.format(from_freq, to_freq, atom_type, ground_state[0], notation[ground_state[1]], int(2*ground_state[2])), data_arr, delimiter = ',')
    return data_arr