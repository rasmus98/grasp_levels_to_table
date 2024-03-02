import numpy as np 
import sys

LINE_READING_LIMIT = 10**7

#preferred orbital order. you might need to change this for your preferred application.
orbitals_order = ['1S', '2S', '2P', '3S', '3P', '4S', '3D', '4P', '4D', '5S', '5P', '6S', '4F', '5D', '6P', '7S', '5F', '6D', '7P', '8S','5G','6G','6H','8P','7D']

class energy_eigenstate:
    def __init__(self,level,terms_strings,angular_momentum,parity,mixing_indices,mixing_amounts,eigenenergy):
        self.mixing_indices = mixing_indices
        self.mixing_amounts = mixing_amounts
        self.terms_strings = terms_strings
        self.angular_momentum = angular_momentum
        self.parity = parity
        self.eigenenergy = eigenenergy
        self.level = level

def check_for_condensed_notation_in_row(occupation_string_array):

    for occupation_string in occupation_string_array:
        x = occupation_string.find('*')
        if x != -1:
            return True 
        
    return False

def decode_condensed_notation(occupation_string_array,num_csfs):
    array = np.zeros(num_csfs)

    offset = 0

    for string in occupation_string_array:
        asterisk_index = string.find('*')

        if asterisk_index == -1:
            array[offset] = int(string)
            offset += 1 
        else:
            num_repeat = int(string[0:asterisk_index])
            repeated_occupation_number = int(string[asterisk_index+1:])
            array[offset:offset+num_repeat] = repeated_occupation_number
            offset = offset + num_repeat


    return array 

def create_csf_array_from_output(grasp_out_path):
    grasp_out = open(grasp_out_path,'r')
    found = False
    while found == False:

        current_line = grasp_out.readline()
        #print(current_line.split())
        split = current_line.split()
        length = len(split)
        #print(length)
        if length > 0: 
            
            #print(length,split)
            if split[0] == 'Input': 
                found = True
                #print('yippee')
        if not current_line:
            print('ran off end of file your GRASP.OUT probably doesnt contain the string (Input) which is what im looking for')
            print('stopping')
            sys.exit() 
    #it's 11pm and im lazy
    grasp_out.readline()
    grasp_out.readline()

    pertinent_details = grasp_out.readline().split()
    num_nrcsf = int(pertinent_details[0])
    num_orbitals = int(pertinent_details[1])
    mode = int(pertinent_details[2])
    if mode == 0:
        print('-------------')
        print('your GRASP printing is in mode 0 - pressing on but it is not likely your grasp.out contains the jj2ls data i want.')
        print('-------------')

    orb_labels = []    

    csf_array = np.zeros([num_nrcsf,num_orbitals])
    for jj in range(0,num_orbitals):
        
        line = grasp_out.readline()
        x = line.find('!')
        if x != -1:
            #print(line[0:x])
            line = line[0:x]

        split_string = line.split()
   
        
        #print(split_string)
        current_orb_string = split_string[0]

        orb_labels.append(current_orb_string)
        occupation_string_array = split_string[1:]
        #print(occupation_string_array)
        length_string_array = len(occupation_string_array)

        condensed = check_for_condensed_notation_in_row(occupation_string_array)

        #print(occupation_string_array)

        #removing comments:
        #print(line)


        if condensed == False: 
            if length_string_array == 1: 
                csf_array[:,jj] = int(occupation_string_array[0])
            elif length_string_array == num_nrcsf:
                csf_array[:,jj] = np.array(occupation_string_array)
            else:
                print("unable to decode this orbital. check orbital no.",jj,current_orb_string)
                sys.exit()
        else:
            csf_array[:,jj] = decode_condensed_notation(occupation_string_array,num_nrcsf)
    #print('found orbitals')
    #print(orb_labels)
    grasp_out.close()
    return csf_array,orb_labels,num_nrcsf,num_orbitals 


def create_csf_array_from_input(grasp_inp_path):

    grasp_inp = open(grasp_inp_path,'r')
    grasp_inp_read = grasp_inp.readlines()
    
    pertinent_details = grasp_inp_read[1].split()
    num_csf = int(pertinent_details[0])
    num_orbs = int(pertinent_details[1])

    orb_labels = []    
    csf_array = np.zeros([num_csf,num_orbs])
    for jj in range(2,num_orbs+2):

        split_string = grasp_inp_read[jj].split()
        
        current_orb_string = split_string[0]

        orb_labels.append(current_orb_string)
        occupation_string_array = split_string[1:]

        length_string_array = len(occupation_string_array)

        condensed = check_for_condensed_notation_in_row(occupation_string_array)


        if condensed == False: 
            if length_string_array == 1: 
                csf_array[:,jj-2] = int(occupation_string_array[0])
            elif length_string_array == num_csf:
                csf_array[:,jj-2] = np.array(occupation_string_array)
            else:
                print("unable to decode this orbital. check orbital no.",jj-2,current_orb_string)
        else:
            csf_array[:,jj-2] = decode_condensed_notation(occupation_string_array,num_csf)
        


    return csf_array,orb_labels,num_csf,num_orbs



def sort_orbitals(orb_labels,orbitals_order,csf_array):
    new_order = np.zeros(len(orb_labels))
    #this is a naive search but the number of orbitals should be small enough for it to not matter

    scores = []

    for jj in range(0,len(orb_labels)):
        found = False 
        for ii in range(0,len(orbitals_order)):
            if orbitals_order[ii] == orb_labels[jj]:
                found = True 
                #print(ii,jj)
                scores.append(ii)
        if found == False:
            print('orbital not found',orb_labels[jj])
            assert(1==0),'orbital not found, not implemented. stopping. add your orbital to the orbitals order global variable'
        
    new_order = np.argsort(scores)
    csf_sorted_by_orbital = np.zeros_like(csf_array)
    sorted_orbital_array = []
    #print(new_order)
    for jj in range(0,len(new_order)):
        new_order[jj]
        sorted_orbital_array.append(orb_labels[int(new_order[jj])])
        csf_sorted_by_orbital[:,jj] = csf_array[:,int(new_order[jj])]


    return sorted_orbital_array,csf_sorted_by_orbital

def make_csf_strings(num_orbitals_to_be_shown,csf_array_resorted_orbitals,sorted_orbital_strings,num_csfs):
    
    num_orbitals_to_be_shown = min(num_orbitals_to_be_shown,len(sorted_orbital_strings))

    csf_strings_in_components = []
    
    csf_strings = []
    
    lengths = []

    #print(sorted_orbital_strings)

    for ii in range(0,num_csfs):
        kept_orbital_strings = []
        current_csf = csf_array_resorted_orbitals[ii]
        occupations_of_interest_indices = np.concatenate(np.argwhere(current_csf>0))
        occupations_of_interest = current_csf[occupations_of_interest_indices]
        #print(occupations_of_interest_indices)

        num_occupations_of_interest = len(occupations_of_interest_indices)

        for kk in range(0,num_occupations_of_interest):
            kept_orbital_strings.append(sorted_orbital_strings[occupations_of_interest_indices[kk]])

        #print(kept_orbital_strings)

        csf_string_broken_down = []
        #the max business is in case all orbitals are requested -- this could go negative, and cause me a headache.
        for kk in range(max(num_occupations_of_interest - num_orbitals_to_be_shown,0),num_occupations_of_interest):
            addition = kept_orbital_strings[kk] + str(int(occupations_of_interest[kk]))
            lengths.append(len(addition))
            csf_string_broken_down.append(addition)

        #csf_string = csf_string.lower()
        
        #print(csf_string)
        csf_strings_in_components.append(csf_string_broken_down)
    
    max_length = max(lengths)

    for jj in range(0,len(csf_strings_in_components)):
        current_string = ''
        for kk in csf_strings_in_components[jj]:
            current_string += kk + (max_length-len(kk))*' ' + ' '
        csf_strings.append(current_string)
        #length_of_this = len(csf_strings[jj])
        #csf_strings[jj] = csf_strings[jj] + (max_length-length_of_this)*' '


    return csf_strings

def find_relativistic_csfs(grasp_out_path,num_csf):
    
    graspout = open(grasp_out_path,'r')

    num_rcsfs_per_csf = np.zeros(num_csf)

    found = False


    while (found == False):
        x = graspout.readline()
        split = x.split()

        if len(split) > 0:
            if split[0] =='NR':
                #print('yes')
                #print(split)

                found = True
                num_rcsfs_per_csf[int(split[2])-1] += int(split[4])
        if not x:
            print('ran off end of file. your GRASP.OUT probably doesnt contain the string (NR CSF) which is what im looking for')
            print('stopping')
            sys.exit()
    while found == True:
        x = graspout.readline()
        split = x.split()
        if len(split) == 0:
            found = False 
        else:
            num_rcsfs_per_csf[int(split[2])-1] += int(split[4])

    #print(num_rcsfs_per_csf)

    num_rcsf_total = int(np.sum(num_rcsfs_per_csf))

    rcsfs_map_to_nrcsfs = np.zeros(num_rcsf_total)

    offset = 0
    
    for kk in range(0,num_csf):
        rcsfs_map_to_nrcsfs[offset:offset+int(num_rcsfs_per_csf[kk])] = kk
        offset = offset+int(num_rcsfs_per_csf[kk])

    graspout.close()
    return rcsfs_map_to_nrcsfs

def find_levels(grasp_out_path):
    graspout = open(grasp_out_path,'r')
    found = False
    

    while found == False:
        lineunplit = graspout.readline()
        line = lineunplit.split()
        if len(line)>0:
            if line[-1] == 'others)':
                #print(line)
                found = True

        if not lineunplit:
            print('ran off end of file. your GRASP.OUT probably doesnt contain the string [ others) ] which is what im looking for')
            print('stopping')
            sys.exit() 

    for jj in range(0,4):
        graspout.readline()

    end = False

    current_length = 0

    states = []
    in_a_state = False
    while end == False:

        line = graspout.readline().split() 
        current_length = len(line)


        if current_length == 0:
            end = True
        else:
            #print(current_length)
            if current_length == 8: #new eigenstate 

                if in_a_state == True:
                    #we are now in a new eigenstate, so save the previous one:
                    #print(current_length)
                    state = energy_eigenstate(level,term_strings,angularmomentum,parity,csf_index,mixing,eigenergy)
                    states.append(state)
                
                in_a_state = True

                #print(line)
                eigenergy = float(line[7])
                level = int(line[0])
                term_strings = [line[1]+line[2]]
                angularmomentum = line[3]
                parity = line[4]
                csf_index = [int(line[5])-1]
                mixing = [float(line[6])]
                
            else: #same eigenstate
                term_strings.append(line[0]+line[1])
                csf_index.append(int(line[4])-1)
                mixing.append(float(line[5]))
                #print(line)
    state = energy_eigenstate(level,term_strings,angularmomentum,parity,csf_index,mixing,eigenergy)
    states.append(state)

    return states

def output_table(csf_strings_prepared,rcsfs_map_to_nrcsfs,eiegenstates_array,user_chosen_num_levels=0):
    num = len(eiegenstates_array)

    if (user_chosen_num_levels != 0) and (user_chosen_num_levels < num):
        num = user_chosen_num_levels 

    header = 'Index,     Energy (Ry),   Parity,        J,      LS Level Composition'
    print(header)
    for jj in range(0,num):

        current_state = eiegenstates_array[jj]
        energy = current_state.eigenenergy
        csf_mixing_coefficients = current_state.mixing_amounts
        csf_mixing_states = current_state.mixing_indices
        current_terms = current_state.terms_strings
        #print(current_state.eigenenergy)
        level = current_state.level
        parity = current_state.parity
        angular_momentum = current_state.angular_momentum

        csf_string = ''

        for kk in range(0,len(csf_mixing_coefficients)):
            current_csf_component_index = csf_mixing_states[kk]

            current_nrcsf_index = rcsfs_map_to_nrcsfs[current_csf_component_index]
            current_nrcsf_string = csf_strings_prepared[int(current_nrcsf_index)].lower()
            
            this_csf_contribution = '{:6.2f}% ({}) {}'.format(round(100*csf_mixing_coefficients[kk]**2,2),current_nrcsf_string,current_terms[kk])

            #print(test)

            csf_string += this_csf_contribution+' '
            if kk < len(csf_mixing_coefficients)-1:
                csf_string +='+ ' 
            
        output_string = '{:5},  {:14E},     {:8},  {:4},    {}'.format(level,energy,parity,angular_momentum,csf_string)
        print(output_string)