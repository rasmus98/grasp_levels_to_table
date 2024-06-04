from library import * 
from library import * 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file',  help='Specify path of GRASP.OUT')
args = parser.parse_args()

if args.file:

#less tested - not perfect.
    graspout = args.file
    #graspout = '/Users/leomulholland/YIII/GRASP.OUT'
    #graspout = 'legacy_testing/GRASP.OUT'


    csf_array,orb_labels,num_nrcsf,num_orbs,mode = create_csf_array_from_output(graspout)
    sorted_orbital_array,csf_sorted_by_orbital = sort_orbitals(orb_labels,orbitals_order,csf_array)
    csf_strings_prepared = make_csf_strings(2,csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf)
    map = find_relativistic_csfs(graspout,num_nrcsf)[0]
    states = find_levels(graspout,False,csf_strings_prepared,map)



    shifted_energies = find_shifted_energies(graspout) 
    if len(shifted_energies) > 0:
        states = add_shifted_energies_to_many_eigenstates(states,shifted_energies)

    data = find_transition_probabilities(graspout)

    data = convert_raw_data_to_transition_class(data,states)

    print_out_a_values(data)
else:
    print('no grasp out')