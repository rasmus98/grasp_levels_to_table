import numpy as np 
import sys

LINE_READING_LIMIT = 10**7
RYDBERG_CM = 109737.316 
#preferred orbital order. you might need to change this for your preferred application.
orbitals_order = ['1S', '2S', '2P', '3S', '3P', '4S', '3D', '4P', '4D', '5S', '4F', '5P', '6S', '5D', '6P', '7S', '5F', '6D', '7P', '8S','5G','6G','6H','8P','7D']

class energy_eigenstate:

    def make_total_string(self,display_inner_terms,csf_strings_prepared,rcsfs_map_to_nrcsfs):

        csf_mixing_coefficients = self.mixing_amounts
        csf_mixing_states = self.mixing_indices
        current_terms = self.terms_strings
        #print(current_state.eigenenergy)
        csf_string =''
        inner_terms = self.inner_terms_strings

        for kk in range(0,len(csf_mixing_coefficients)):
            current_csf_component_index = csf_mixing_states[kk]

            current_nrcsf_index = rcsfs_map_to_nrcsfs[0][current_csf_component_index]
            #print(current_nrcsf_index)

            current_nrcsf_string = csf_strings_prepared[int(current_nrcsf_index)].lower()
            
            if display_inner_terms and len(inner_terms) >0:
                index = find_place_for_inner_term(current_nrcsf_string)

                current_nrcsf_string = current_nrcsf_string[0:index] + '('+inner_terms[kk]+ ') ' + current_nrcsf_string[index:]

            this_csf_contribution = '{:6.2f}% [{} ({})]'.format(round(100*csf_mixing_coefficients[kk]**2,2),current_nrcsf_string,current_terms[kk])

            csf_string += this_csf_contribution+' '
            if kk < len(csf_mixing_coefficients)-1:
                csf_string +='+ ' 
            self.csf_string = csf_string 
    
    def make_leading_term_string(self,csf_strings_prepared,rcsfs_map_to_nrcsfs):
        

        current_csf_component_index = self.mixing_indices[0]
        current_nrcsf_index = rcsfs_map_to_nrcsfs[0][current_csf_component_index]


        leading_term_string = csf_strings_prepared[int(current_nrcsf_index)].lower()
        leading_term_string += '('+self.terms_strings[0]

        if self.parity == 'odd':
            leading_term_string += '*)' 
        else:
            leading_term_string += ' )' 
        leading_term_string += '_{'+self.angular_momentum+'}'
        self.leading_term_string = leading_term_string



    def __init__(self,level,terms_strings,angular_momentum,parity,mixing_indices,mixing_amounts,eigenenergy,csf_strings_prepared,rcsfs_map_to_nrcsfs,display_inner_terms,inner_terms_strings=[]):
        self.mixing_indices = mixing_indices
        self.mixing_amounts = mixing_amounts
        self.terms_strings = terms_strings
        self.inner_terms_strings = inner_terms_strings
        self.angular_momentum = angular_momentum

        if len(angular_momentum) > 2:

            if angular_momentum[-2] == "/":
                self.angular_momentum_float = float(angular_momentum[0]) / float(angular_momentum[-1])
            else:
                self.angular_momentum_float = float(angular_momentum)


        self.parity = parity
        self.eigenenergy = eigenenergy
        self.level = level
        self.wavenumber = eigenenergy * RYDBERG_CM
        if eigenenergy<0:
            self.wavenumber = 0.0

        self.shifted_energy_ryd = -1.0
        self.shifted_energy_cm= -1.0

        self.make_total_string(display_inner_terms,csf_strings_prepared,rcsfs_map_to_nrcsfs)
        self.make_leading_term_string(csf_strings_prepared,rcsfs_map_to_nrcsfs)

    def set_shifted_energy(self,shifted_energy_ryd):
        self.shifted_energy_ryd = shifted_energy_ryd
        self.shifted_energy_cm = shifted_energy_ryd * RYDBERG_CM


class raw_grasp_data:
    def __init__(self,grasp_line,electric):
        
        self.lower_index_unshifted = min(int(grasp_line[0]),int(grasp_line[1]))
        self.upper_index_unshifted = max(int(grasp_line[0]),int(grasp_line[1]))
        self.wavelength_nm_vac = float(grasp_line[2]) / 10.0
        self.a_ji = float(grasp_line[3])
        self.f_ij = float(grasp_line[4])
        self.lstrength = float(grasp_line[5])

        if electric:
            self.vel_len = float(grasp_line[6])
        else:
            self.vel_len = -np.inf

class transition:
    def __init__(self,eigenstates:list[energy_eigenstate],raw_grasp_data_class: raw_grasp_data,lower_energy_ryd_shifted = -1.0,upper_energy_ryd_shifted = -1.0):
        
        #i hate - myself

        self.lower_index_unshifted = raw_grasp_data_class.lower_index_unshifted
        self.upper_index_unshifted = raw_grasp_data_class.upper_index_unshifted
        
        self.lower_eigenstate = eigenstates[raw_grasp_data_class.lower_index_unshifted - 1]
        self.upper_eigenstate = eigenstates[raw_grasp_data_class.upper_index_unshifted - 1]

        if (self.lower_eigenstate.shifted_energy_cm == -1.0) or (self.upper_eigenstate.shifted_energy_cm == -1.0):
            self.wavelength_nm_vac = raw_grasp_data_class.wavelength_nm_vac
        else:
            self.wavelength_nm_vac = 1e7 / np.abs (self.lower_eigenstate.shifted_energy_cm - self.upper_eigenstate.shifted_energy_cm )
        
        self.wavelength_ang_vac = self.wavelength_nm_vac * 10.0
        self.wavelength_ang_air = vac_to_air(self.wavelength_ang_vac)
        self.wavelength_nm_air = 0.1 * self.wavelength_ang_air



        self.upper_weight = 2.0 * self.upper_eigenstate.angular_momentum_float + 1.0
        self.lower_weight = 2.0 * self.lower_eigenstate.angular_momentum_float + 1.0

        self.upper_ang = self.upper_eigenstate.angular_momentum_float
        self.lower_ang = self.lower_eigenstate.angular_momentum_float

        self.avalue = raw_grasp_data_class.a_ji
        self.f = raw_grasp_data_class.f_ij

        self.gavalue = self.avalue * self.upper_weight
        self.gf = self.f * self.lower_weight
        self.loggf = np.log10(self.gf)
        self.s = raw_grasp_data_class.lstrength
        self.vel_len = raw_grasp_data_class.vel_len
        self.energy_lower_shifted_ryd = lower_energy_ryd_shifted
        self.energy_upper_shifted_ryd = upper_energy_ryd_shifted

        self.energy_lower_shifted_wn = lower_energy_ryd_shifted * RYDBERG_CM
        self.energy_upper_shifted_wn = upper_energy_ryd_shifted * RYDBERG_CM

    def display_transition(self):
        string = '{:5} {:5} {:14.4F}  {:3.1F}   {:14.4F}  {:3.1F}    {:14.4F}    {:14.4F}  {:10.2E} {:10.2E}    {:10.2E}  {:10.2E}  {:10.4f} {:10.4f}    {}    {}'
        
        if (self.lower_eigenstate.shifted_energy_cm != -1.0) and (self.upper_eigenstate.shifted_energy_cm != -1.0):
            lower_energy = self.lower_eigenstate.shifted_energy_cm
            upper_energy = self.upper_eigenstate.shifted_energy_cm
        else:
            lower_energy = self.lower_eigenstate.wavenumber
            upper_energy = self.upper_eigenstate.wavenumber           
        
        print(string.format(self.lower_index_unshifted,\
                            self.upper_index_unshifted,\
                            lower_energy,\
                            self.lower_ang,\
                            upper_energy,\
                            self.upper_ang,\
                            self.wavelength_nm_vac,\
                            self.wavelength_nm_air,\
                            self.avalue,\
                            self.gavalue,\
                            self.f,\
                            self.gf,\
                            self.loggf,\
                            self.vel_len,\
                            self.lower_eigenstate.leading_term_string,\
                            self.upper_eigenstate.leading_term_string))

def vac_to_air(wl_ang_vac):

    #converts a vaccumm wl (ang) to air. there are supposedly many conversion.
    #the commented one is from vald, the other is from https://www.sdss3.org/dr8/spectro/spectra.php
    
    #they give identical results as far as I can see. The other one might be cheaper tbf.

    if (wl_ang_vac > 2000.0 and wl_ang_vac < 100000.0):
        ssq = 1e8 / (wl_ang_vac * wl_ang_vac)
        #n =  1.0 + 0.0000834254 + 0.02406147 / (130.0 - ssq) + 0.00015998  / (38.9 - ssq)
        n = 1.0 + 2.735182E-4 + 131.4182 / wl_ang_vac**2 + 2.76249E8 / wl_ang_vac**4
        return wl_ang_vac / n
    else:
        return wl_ang_vac


def check_for_condensed_notation_in_row(occupation_string_array):   
    #checks for condensed notatation in grasp input
    for occupation_string in occupation_string_array:
        x = occupation_string.find('*')
        if x != -1:
            return True 
        
    return False

def decode_condensed_notation(occupation_string_array,num_csfs):
    #decodes grasp inp condensed notation
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
    
    #makes an np array of input csf occupations from grasp out

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
    #if mode == 0:
    #    print('-------------')
    #    print('your GRASP printing is in mode 0 - pressing on but it is not likely your grasp.out contains the jj2ls data i want.')
    #    print('-------------')

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
    return csf_array,orb_labels,num_nrcsf,num_orbitals,mode 


def create_csf_array_from_input(grasp_inp_path):
    #    #makes an np array of input csf occupations from grasp inp - retired

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
    #makes strings out of the input csfs.
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
    #reads graspout and finds relativistic csfs indicies.
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
    return rcsfs_map_to_nrcsfs,num_rcsf_total

def find_levels(grasp_out_path,inner,csf_strings_prepared,rcsfs_map_to_nrcsfs):
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
                    state = energy_eigenstate(level,term_strings,angularmomentum,parity,csf_index,mixing,eigenergy,csf_strings_prepared,rcsfs_map_to_nrcsfs,inner)
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
    state = energy_eigenstate(level,term_strings,angularmomentum,parity,csf_index,mixing,eigenergy,csf_strings_prepared,rcsfs_map_to_nrcsfs,inner)
    states.append(state)

    return states


def find_place_for_inner_term(csf_string):
    for kk in range(len(csf_string)-1,-1,-1):
        if csf_string[kk] != ' ':
            break 
    
    for ii in range(kk,-1,-1):
        if csf_string[ii] == ' ':
            break 

    return ii    


def output_table(eiegenstates_array,user_chosen_num_levels=0):
    num = len(eiegenstates_array)

    if (user_chosen_num_levels != 0) and (user_chosen_num_levels < num):
        num = user_chosen_num_levels 

    header = 'Index,     Energy (Ry),   Parity,        J,      LS Level Composition'
    print(header)
    for jj in range(0,num):

        current_state = eiegenstates_array[jj]
        energy = current_state.eigenenergy

        level = current_state.level
        parity = current_state.parity
        angular_momentum = current_state.angular_momentum
        csf_string = current_state.csf_string
            
        output_string = '{:5},  {:14E},     {:8},  {:4},    {}'.format(level,energy,parity,angular_momentum,csf_string)
        print(output_string)




def find_inner_occupation_terms(graspoutpath,outputmode,num_rcsf):
    if outputmode < 2:
        print("output mode is less than 2 - inner term strings cannot be found. continuing.")
        return []
    
    grasp_out = open(graspoutpath,'r')

    found = False 
    counter = 0
    
    while found == False: 
        counter +=1
        current_line_read = grasp_out.readline()
        current_line_split = current_line_read.split()
        if len(current_line_split) > 2:
            if current_line_split[1] =='NROUT:':
                #print(current_line_read)
                break
        if not current_line_read:
            print('failure in locating inner terms')
            print('ran off end of file. your GRASP.OUT probably doesnt contain the string (CSFs are defined using format) which is what im looking for')
            print('stopping')
            sys.exit()
    found = False
    while found == False:
        current_line_read = grasp_out.readline().split()
        if len(current_line_read) > 0:
            if current_line_read[0] == 'CSFs':
                found = True
    for jj in range(0,3):
        grasp_out.readline()
    
    inner_term_string_array = []
    

    for kk in range(0,num_rcsf):
        #print(kk)
        inner_term_string = ''
        for ii in range(0,6):
            current_line_read = grasp_out.readline()
            #print(current_line_read)
            #print(current_line_read)
            split = current_line_read.split()
            length = len(split)
            #print(length)
            if ii == 2:
                if length < 4: #only one thing outside a closed shell. therefore the closed shell is 1S
                    inner_term_string = inner_term_string + '1'
                else:
                    inner_term_string = inner_term_string + split[-4]
            if ii ==3:
                if length < 5: #only one thing outside a closed shell. therefore the closed shell is 1S
                    inner_term_string = inner_term_string + 'S'
                else:
                    inner_term_string = inner_term_string + split[-5]
        inner_term_string_array.append(inner_term_string)

    return inner_term_string_array

def add_inner_occupation_strings_to_eigenclass(graspoutpath,outputmode,energy_eigenstate_list,display_inner_terms,csf_strings_prepared,map):

    num_rcsf = len(energy_eigenstate_list)

    inner_term_string_array = find_inner_occupation_terms(graspoutpath,outputmode,num_rcsf)

    for kk in range(0,len(energy_eigenstate_list)):
        inner_term_list = []
        for index in energy_eigenstate_list[kk].mixing_indices:
            inner_term_list.append(inner_term_string_array[index])

        energy_eigenstate_list[kk].inner_terms_strings = inner_term_list 
        energy_eigenstate_list[kk].make_total_string(display_inner_terms,csf_strings_prepared,map)

    #state = energy_eigenstate_list[0]
    #print(state.inner_terms_strings)


    return energy_eigenstate_list


def find_transition_probabilities(grasp_out_path):
    
    graspout = open(grasp_out_path,'r')
    found = False
    while (found == False):
        x = graspout.readline()
        split = x.split()

        if len(split) > 0:
            if split[0] =='Electric':
                #print('yes')
                #print(split)
                found = True
        if not x:
            print('ran off end of file. your GRASP.OUT probably doesnt contain the string (Electric ... emission transition probability (sec-1) ... length) which is what im looking for')
            print('stopping')
            sys.exit()

    for jj in range(0,11):
        graspout.readline()

    electric_data = []
    x = graspout.readline()
    split = x.split()

    while len(split) != 0:
        raw_grasp_transition = raw_grasp_data(split,True)
        electric_data.append(raw_grasp_transition)
        x = graspout.readline()
        split = x.split()




    found = False
    while (found == False):
        x = graspout.readline()
        split = x.split()

        if len(split) > 0:
            if split[0] =='Magnetic':
                #print('yes')
                #print(split)
                found = True
        if not x:
            print('ran off end of file. your GRASP.OUT probably doesnt contain the string (Magnetic ) which is what im looking for')
            print('stopping')
            sys.exit()

    for jj in range(0,11):
        graspout.readline()

    magnetic_data = []
    x = graspout.readline()
    split = x.split()

    while len(split) != 0:
        raw_grasp_transition = raw_grasp_data(split,False)
        magnetic_data.append(raw_grasp_transition)
        x = graspout.readline()
        split = x.split()

    #print(magnetic_data[-1].lstrength)
    graspout.close()


    total_data = electric_data.copy()

    total_data.extend(magnetic_data)
    
    #print(electric_data)

    return total_data

def find_shifted_energies(grasp_out_path):
    
    graspout = open(grasp_out_path,'r')
    found = False
    while (found == False):
        x = graspout.readline()
        split = x.split()

        if len(split) > 0:
            if split[0] =='OSCL':
                #print('yes')
                #print(split)
                found = True
        if not x:
            print('ran off end of file. your GRASP.OUT probably doesnt contain the string (OSCL ... emission transition probability (sec-1) ... length) which is what im looking for')
            print('stopping')
            sys.exit()

    if '17' not in split: #then there is no shifting
        return []
    else:
        shifted_energies = []
        x = graspout.readline()
        x = x.replace(':',' ')
        split = x.split()
        num = int(split[-1])
        for jj in range(0,num):
            x = graspout.readline()
            x = x.replace('c',' ')
            shifted_energies.append(float(x))

        shifted_energies_ryd = np.array(shifted_energies)

        return shifted_energies_ryd

def add_shifted_energies_to_many_eigenstates(eigenstates:list[energy_eigenstate],shifted_energies_ryd):
    counter = 0
    for state in eigenstates:
            state.set_shifted_energy(shifted_energies_ryd[counter])
            counter +=1 
            if counter > len(shifted_energies_ryd)-1:
                break

    return eigenstates

def convert_raw_data_to_transition_class(total_data:list[raw_grasp_data],eigenstates:list[energy_eigenstate]):

    total_data_class = []

    for transition_raw in total_data:
        total_data_class.append( transition(eigenstates,transition_raw) )
    
    return total_data_class 

def print_out_a_values(total_data_class_array:list[transition]):
    
    
    #string = '{:5} {:5} {:14.4F}  {:3.1F}   {:14.4F}  {:3.1F}    {:14.4F}  {:10.2E} {:10.2E}    {:10.2E}  {:10.2E}  {:10.2E}    {}    {}'

    header = '#{:>4s} {:>5s}   {:14} {:3}    {:14} {:3}   {:>14s}    {:>14s}  {:>10s} {:>10s}     {:>10s} {:>10s}  {:>10s} {:>10s}   {:>20s}  {:>20s}'
    print(header.format(
            'Low','Upp'
            ,'Lower (cm-1)','JL'
            ,'Upper (cm-1)','JU'
            ,'λ (vac,nm)'
            ,'λ (air,nm)'
            ,'A (s-1)','gA (s-1)'
            ,'f ','gf','log(gf)'
            ,'vel\len'
            ,'LowDesignation','UppDesignation'
    ))
    for trans in total_data_class_array:
        trans.display_transition()

    return 0