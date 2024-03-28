from jj_library import * 

grasp_out_path = 'legacy_testing/GRASP.OUT'

index = locate_eigenvectors(grasp_out_path)

#print(index)
num_to_be_shown = 10
num_rcsf = 31 

eigenvectors = read_in_eigenvectors(index,31,grasp_out_path)

x = find_levels_jj(grasp_out_path)

levels = read_in_eigenlevels(grasp_out_path,x,num_rcsf)

eigenstate_array = produce_eigenstate_array(levels,eigenvectors)

eigenstate_array = make_many_eigenstate_strings(eigenstate_array,num_to_be_shown)
output_table(eigenstate_array,10)