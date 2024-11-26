from as_lib import * 

path = '/Users/leomulholland/he_opt/olg'

das_file_numpy,orbital_strings,num_csfs,lambda_array = read_das('/Users/leomulholland/SeIII/as/das5/das_5')
csf_strings,pseudo_array = make_csf_strings(das_file_numpy,orbital_strings,num_csfs,5,lambda_array)

states = read_oic_into_list_of_eigenstates('/Users/leomulholland/SeIII/as/das5/oic',csf_strings,num_levels=2**63)


transitions = get_transitions(path,states)

print_out_a_values(transitions)
