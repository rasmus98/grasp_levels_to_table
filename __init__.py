import numpy as np
from .library import * 


def read_output(grasp_out_path,num_orbitals_shown=np.inf,display_inner_terms=True):
    csf_array,orb_labels,num_nrcsf,num_orbs,mode = create_csf_array_from_output(grasp_out_path)
    sorted_orbital_array,csf_sorted_by_orbital = sort_orbitals(orb_labels,csf_array)
    if mode == 0:
        print('-------------')
        print('your GRASP printing is in mode 0 - pressing on but it is not likely your grasp.out contains the jj2ls data i want.')
        print('-------------')

    csf_strings_prepared = make_csf_strings(num_orbitals_shown,csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf)
    map,num_rcsf = find_relativistic_csfs(grasp_out_path,num_nrcsf)

    if display_inner_terms:
        inner_terms = find_inner_occupation_terms(grasp_out_path,mode,num_rcsf)
    else:
        inner_terms = []
    #print(inner_terms)
    states = find_levels(grasp_out_path,inner_terms,csf_strings_prepared,map)
    
    if display_inner_terms:
        states = add_inner_occupation_strings_to_eigenclass(grasp_out_path,mode,states,display_inner_terms,csf_strings_prepared,map)

    return states
