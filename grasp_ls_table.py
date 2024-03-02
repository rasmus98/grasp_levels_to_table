from library import * 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file',  help='Specify path of GRASP.OUT')
parser.add_argument('-n', '--num',  help='Requested number of states (all states by default)',type=int)
parser.add_argument('-o', '--num_orbitals',  help='number of orbitals to show in NRCSF (2 default)',type=int)
args = parser.parse_args()

def main(grasp_out_path,num_orbitals_shown,user_num_levels):
    csf_array,orb_labels,num_nrcsf,num_orbs = create_csf_array_from_output(grasp_out_path)
    sorted_orbital_array,csf_sorted_by_orbital = sort_orbitals(orb_labels,orbitals_order,csf_array)
    csf_strings_prepared = make_csf_strings(num_orbitals_shown,csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf)
    map = find_relativistic_csfs(grasp_out_path,num_nrcsf)
    states = find_levels(grasp_out_path)
    output_table(csf_strings_prepared,map,states,user_num_levels)

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
        
    main(args.file,num_orbitals_shown,user_num_levels)



