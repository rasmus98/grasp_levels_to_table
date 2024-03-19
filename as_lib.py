import numpy as np 

ANGULAR_SYMBOLS = ['s','p','d','f','g','h','i']

def read_das(path_to_das):

    das_file_opened = open(path_to_das,'r')
    das_file_read = das_file_opened.readlines()

    orbital_data = das_file_read[2].split()

    total_num_orbs = int( len(orbital_data) / 2 ) 
    
    orbital_strings = []

    for kk in range(0,len(orbital_data),2):
        orbital_strings.append(orbital_data[kk]+translate_angular(int(orbital_data[kk+1])))


    #print(orbital_strings)

    for jj in range(0,len(das_file_read)): 
        current_line = das_file_read[jj]

        if current_line.split()[0] == '&SMINIM':
            print(jj-3,'csfs found')
            break    
    
    num_csfs = jj - 3

    das_file_opened.close() 

    das_file_numpy = np.loadtxt(path_to_das,skiprows=3,max_rows=num_csfs,dtype=int)

    #print(das_file_numpy[0])

    return das_file_numpy,orbital_strings,num_csfs

def make_csf_strings(das_file_numpy,orbital_strings,num_csfs,num_requested_orbitals):
    csf_strings = []

    

    for jj in range(0, num_csfs):
        current_csf_ints = das_file_numpy[jj]

        concerned_occupations_locations = np.argwhere(current_csf_ints > 0 )

        #print(concerned_occupations_locations)
        concerned_occupations_numbers = current_csf_ints[concerned_occupations_locations]

        concerned_occupations_orbitals = []

        for kk in range(0,len(concerned_occupations_locations)):
            concerned_occupations_orbitals.append(orbital_strings[concerned_occupations_locations[kk][0]])

        
        num_orbitals = len(concerned_occupations_orbitals)
        current_csf_string = ''

        min = max(0,num_orbitals - num_requested_orbitals)

        for kk in range(min,num_orbitals):
            current_csf_string += concerned_occupations_orbitals[kk] + str(concerned_occupations_numbers[kk][0]) + " "
        #print(current_csf_string)

        csf_strings.append(current_csf_string)
    return csf_strings

def read_terms_and_output(terms,csf_strings,num_levels):

    terms_data = np.loadtxt(terms,skiprows=1) 
    #print(np.shape(terms_data))
    if num_levels > np.shape(terms_data)[0]:
        num_levels = np.shape(terms_data)[0]


            
    output_string = '{:5},   {:10.6f},     {}'


    #print(num_levels)

    for kk in range(0,num_levels-1):
        line = terms_data[kk]
        multiplicity = line[0]
        angular = line[1]
        parity = line[2]
        cf_number = int(line[3]-1)
        energy_ryd = line[5]

        term_string = '('+str(int(multiplicity)) + translate_angular(int(angular))

        if int(parity) == 0: 
            term_string += ' '
        elif int(parity) == 1:
            term_string += '*'
        else:
            print('failure in parity idenficiation at level',kk+1,'parity found',parity)

        term_string += ')'
        print(output_string.format(kk+1,energy_ryd,csf_strings[cf_number] + term_string.upper()))



    return 0

def read_levels_and_output(levels,csf_strings,num_levels):

    terms_data = np.loadtxt(levels,skiprows=1) 
    #print(np.shape(terms_data))
    if num_levels > np.shape(terms_data)[0]:
        num_levels = np.shape(terms_data)[0]


            
    output_string = '{:5},   {:12.8f},     {}  {}'


    #print(num_levels)

    header = 'Index       Energy(Ry)     CSF(TERM)        J'
    print(header)
    for kk in range(0,num_levels-1):
        line = terms_data[kk]

        j = int(line[0]) / 2
        parity = line[1]

        multiplicity = abs(line[2])
        angular = line[3]
        cf_number = int(line[4]-1)
        energy_ryd = line[-1]

        term_string = '('+str(int(multiplicity)) + translate_angular(int(angular))

        if int(parity) == 0: 
            term_string += ' '
        elif int(parity) == 1:
            term_string += '*'
        else:
            print('failure in parity idenficiation at level',kk+1,'parity found',parity)

        term_string += ')'
        print(output_string.format(kk+1,energy_ryd,csf_strings[cf_number] + term_string.upper(),j))



    return 0

def translate_angular(angular_number):
    if angular_number < 7:
        return ANGULAR_SYMBOLS[angular_number]
    else:
        return str(angular_number)
