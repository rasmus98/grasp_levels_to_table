from library import * 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file',  help='Specify path of GRASP.OUT')
parser.add_argument('-n', '--num',  help='Requested number of states (all states by default)',type=int)

#parser.add_argument('-o', '--num_orbitals',  help='number of orbitals to show in NRCSF (2 default)',type=int)
#parser.add_argument('-i', '--inner_terms',  help='show inner terms',action='store_true')

args = parser.parse_args()
num_disp = np.inf

def main(grasp_out_path):
    num_orbitals_shown = 4
    display_inner_terms = False
    csf_array,orb_labels,num_nrcsf,num_orbs,mode = create_csf_array_from_output(grasp_out_path)
    sorted_orbital_array,csf_sorted_by_orbital = sort_orbitals(orb_labels,orbitals_order,csf_array)

    if mode == 0:
        print('-------------')
        print('your GRASP printing is in mode 0 - pressing on but it is not likely your grasp.out contains the jj2ls data i want.')
        print('-------------')

    csf_strings_prepared = make_csf_strings(num_orbitals_shown,csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf)
    map,num_rcsf = find_relativistic_csfs(grasp_out_path,num_nrcsf)

    csf_strings_adas = make_csfs_strings_for_adas(csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf)

    inner_terms = []
    states = find_levels(grasp_out_path,inner_terms,csf_strings_adas,map)
    string = '{:>1}{:11} {:4} {:>12.3f} {:>4}'
    nummmm = min(len(states),num_disp)
    shifted_energies = find_shifted_energies(grasp_out_path) 
    if len(shifted_energies) > 0:
        states = add_shifted_energies_to_many_eigenstates(states,shifted_energies)
        nummmm = min(nummmm,len(shifted_energies))
        
        new_order = np.argsort(shifted_energies)
        #print(new_order)

    #print(shifted_energies)
    
    print("&ADASEX NLEVS= {} NUMTMP=19 IRDTMP=1 ITCC=1 IBORN=0 IELAS=0 IEL='SYMBOL_HERE' FIPOT=ION_POT_CM_HERE IONTRM='TERM_HERE'/".format(nummmm))
    print('1.00+03 1.50+03 1.80+03 2.00+03 2.50+03 5.00+03 7.50+03 1.00+04 1.50+04 1.80+04 2.00+04 3.00+04 4.00+04 5.00+04 6.00+04 7.00+04 8.00+04 9.00+04 1.00+05')
    
    for index in range(0,nummmm):
        index_of_interest = index 
        if len(shifted_energies) > 0:
            index_of_interest = new_order[index]

        state = states[index_of_interest]
        adas_label = state.leading_csf+state.leading_term_only_adas
        if len(adas_label) > 13:
            adas_label = adas_label[len(adas_label)-13:len(adas_label)]
        else:
            adas_label = ' ' * (13-len(adas_label)) + adas_label
        energy = state.wavenumber
        if len(shifted_energies) > 0:
            energy = state.shifted_energy_cm
        x = string.format(state.parity_character,adas_label,state.angular_momentum_float,energy,index+1)
        print(x)
        #print(string.format(state.parity_character),state.leading_csf,state.multiplicity,state.angular_momentum_orbital_int,state.angular_momentum_float,state.wavenumber)

    #if display_inner_terms:
    #    states = add_inner_occupation_strings_to_eigenclass(grasp_out_path,mode,states,display_inner_terms,csf_strings_prepared,map)
    print('NAME:')
    print('DATE:')
    print('')
    print('')
    print('')
    print('                      ENTER DETAILS OF CALCULATION')                  
    print('')
    print('')
    print('')
    print('')
    print('.')
    return 0
if args.num:
    num_disp = int(args.num)
if not args.file:
    print('no grasp out given. stopping')
else:
    main(args.file)

    


