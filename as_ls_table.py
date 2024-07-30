from as_lib import * 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-d', '--das',  help='Specify path of das')
parser.add_argument('-t', '--terms',  help='Specify path of TERMS')
parser.add_argument('-l', '--oic',  help='Specify path of oic')

parser.add_argument('-n', '--num',  help='Requested number of states (all states by default)',type=int)
parser.add_argument('-o', '--num_orbitals',  help='number of orbitals to show in NRCSF (2 default)',type=int)

args = parser.parse_args()


if not args.das:
    print('no das given. stopping')
else:

    if (not args.terms) and (not args.oic):
        print('no terms or oic given. stopping')
    else:

        if args.num:
            user_num_levels = args.num 
        else:
            user_num_levels = 2**63

        num_orbitals_shown = 2
        if args.num_orbitals:
            num_orbitals_shown = max(num_orbitals_shown,args.num_orbitals)


        das_file_numpy,orbital_strings,num_csfs,lambda_array = read_das(args.das)
        csf_strings,pseudo_array = make_csf_strings(das_file_numpy,orbital_strings,num_csfs,num_orbitals_shown,lambda_array)

        if args.terms:
            read_terms_and_output(args.terms,csf_strings,user_num_levels,pseudo_array)
        if args.oic:
            states = read_oic_into_list_of_eigenstates(args.oic,csf_strings,user_num_levels)
            header = 'Index       Energy(Ry)     CSF(TERM)        J     LV'
            print(header)
            for state in states:
                state.display_state()


