from library import * 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file',  help='Specify path of GRASP.OUT')
parser.add_argument('-n', '--num',  help='Requested number of states (all states by default)',type=int)
parser.add_argument('-o', '--num_orbitals',  help='number of orbitals to show in NRCSF (2 default)',type=int)
parser.add_argument('-i', '--inner_terms',  help='show inner terms',action='store_true')

args = parser.parse_args()

def main(grasp_out_path,num_orbitals_shown,user_num_levels,display_inner_terms):
    csf_array,orb_labels,num_nrcsf,num_orbs,mode = create_csf_array_from_output(grasp_out_path)
    sorted_orbital_array,csf_sorted_by_orbital = sort_orbitals(orb_labels,orbitals_order,csf_array)

    csf_strings_prepared = make_csf_strings(num_orbitals_shown,csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf)
    map,num_rcsf = find_relativistic_csfs(grasp_out_path,num_nrcsf)

    if args.inner_terms:
        inner_terms = find_inner_occupation_terms(grasp_out_path,mode,num_rcsf)
    else:
        inner_terms = []
    print(inner_terms)
    states = find_levels(grasp_out_path,inner_terms,csf_strings_prepared,map)
    
    if display_inner_terms:
        states = add_inner_occupation_strings_to_eigenclass(grasp_out_path,mode,states,display_inner_terms,csf_strings_prepared,map)

    output_table(states,user_num_levels)

    return 0

if not args.file:
    print('no grasp out given. stopping')
else:
    if args.num:
        user_num_levels = args.num 
    else:
        user_num_levels = 0

    num_orbitals_shown = 2
    if args.num_orbitals:
        num_orbitals_shown = max(num_orbitals_shown,args.num_orbitals)
    display_inner_terms = False
    if args.inner_terms:
        display_inner_terms = True
        
    main(args.file,num_orbitals_shown,user_num_levels,display_inner_terms)



