from as_lib import * 

path = '/Users/leomulholland/he_opt/olg'

das_file_numpy,orbital_strings,num_csfs = read_das('/Users/leomulholland/he_opt/das')
csf_strings = make_csf_strings(das_file_numpy,orbital_strings,num_csfs,4)

states = read_oic_into_list_of_eigenstates('/Users/leomulholland/he_opt/oic',csf_strings,num_levels=2**63)


transitions = get_transitions(path,states)

print_out_a_values(transitions)
