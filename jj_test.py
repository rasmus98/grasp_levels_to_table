from jj_library import * 
import library
grasp_out_path = 'legacy_testing/GRASP.OUT'

num_orbitals_shown = 4

csf_array,orb_labels,num_nrcsf,num_orbs,mode = library.create_csf_array_from_output(grasp_out_path)
sorted_orbital_array,csf_sorted_by_orbital = library.sort_orbitals(orb_labels,library.orbitals_order,csf_array)

csf_strings_prepared = library.make_csf_strings(num_orbitals_shown,csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf)

print(csf_strings_prepared)

map,num_rcsf = library.find_relativistic_csfs(grasp_out_path,num_nrcsf)

#print(csf_sorted_by_orbital)

index = locate_eigenvectors(grasp_out_path)

#print(index)
num_components_to_be_shown = 3
eigenvectors = read_in_eigenvectors(index,num_rcsf,grasp_out_path)

x = find_levels_jj(grasp_out_path)

levels = read_in_eigenlevels(grasp_out_path,x,num_rcsf)

pos = locate_cfout_jj(grasp_out_path)
jj_terms_raw = read_cfout_jj(grasp_out_path,pos,num_rcsf)
jj_terms = decode_jj_terms(jj_terms_raw)
#print(jj_terms)
eigenstate_array = produce_eigenstate_array(levels,eigenvectors)

eigenstate_array = make_many_eigenstate_strings(eigenstate_array,num_components_to_be_shown,jj_terms,map,csf_strings_prepared)
output_table(eigenstate_array,20)

