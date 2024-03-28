import sys
import numpy as np 
import warnings
RYDBERG_CM = 109737.316 
EIGENVECTOR_CUT_OFF = 0.001
class eigenstate:
    def __init__(self,energy_ryd_unshifted,j_value,even_boolean,eigenvector,level_number):

        self.eigenvector = eigenvector
        self.energy_ryd_unshifted = energy_ryd_unshifted
        self.energy_wn_unshifted = energy_ryd_unshifted * RYDBERG_CM
        self.jvalue = j_value
        self.even = even_boolean
        self.level_number =  level_number
        self.label_string = ''
        self.mixing_coefficients = []
        self.mixing_indices = []
            

    def make_csf_string(self,num_states_to_be_shown):
        ind,coef = isolate_most_contributing_csfs_from_eigenvector(self.eigenvector,num_states_to_be_shown)
        self.mixing_coefficients = coef 
        self.mixing_indices = ind
        string = ''
        for jj in range(0,len(ind)):
            if abs(coef[jj]) > EIGENVECTOR_CUT_OFF:
                string+= str(round(coef[jj],3))+'*' + str(ind[jj]) +' '
        self.label_string = string

        


def locate_eigenvectors(graspoutpath):
    index = 0
    grasp_out_read = open(graspoutpath,'r').readlines()

    for jj in range(0,len(grasp_out_read)):
        current_line = grasp_out_read[jj].split()
        #print(current_line)
        if len(current_line) == 1:
            if current_line[0] == 'eigenvectors':
                index = jj 
                return index + 3

    print('ran off end of file')
    sys.exit()
    return 0 

def read_in_eigenvectors(starting_position,num_rcsf,graspoutpath):

    num_blocks = int(np.ceil(num_rcsf / 5 ))
    #print(num_blocks)

    blocks = []

    for jj in range(0,num_blocks):

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            current_block = np.loadtxt(graspoutpath,skiprows=starting_position,max_rows=num_rcsf)
        starting_position +=  3 + num_rcsf
        shape = np.shape(current_block)
        #print(shape)
        if len(shape) > 1:
            num_cols = shape[1]
            for kk in range(0,num_cols):
                blocks.append(current_block[:,kk])
        else:
            blocks.append(current_block)
        
        #print(current_block[:,0])

    blocks = np.transpose(np.array(blocks))

    

    return blocks

def isolate_most_contributing_csfs_from_each_eigenvector(eigenvector_matrix):
    
    #print(np.shape(eigenvector_matrix))

    blocks_temp = np.abs(eigenvector_matrix)
    
    print(eigenvector_matrix[:,4])

    for ii in range(0,12):
        current_eigenvector = blocks_temp[:,ii]
        ind = np.argpartition(current_eigenvector, -4)[-4:][::-1]

        ind = ind[np.argsort(np.abs(current_eigenvector[ind]))]
        ind = ind[::-1]
        string = ''
        for jj in range(0,4):
            current_eigenvector[ind][jj]
            string += str(round(current_eigenvector[ind][jj],2)) +'*'+ str(ind[jj]) +' '
        print(string)
        #print(ind) 
        #print()
    return 0 

def isolate_most_contributing_csfs_from_eigenvector(eigenvector,num_to_be_shown):
        
    ind = np.argpartition(np.abs(eigenvector), -num_to_be_shown)[-num_to_be_shown:][::-1]

    ind = ind[np.argsort(np.abs(eigenvector[ind]))]
    ind = ind[::-1]
    string = ''
    indices = ind 
    coefficients = eigenvector[ind]

    return indices,coefficients



def find_levels_jj(graspout):

    file  = open(graspout,'r')
    grasp_out_read = file.readlines()

    

    for jj in range(0,len(grasp_out_read)):
        current_line = grasp_out_read[jj].split() 
        if len(current_line) > 0 :
            if current_line[-1] == 'lowest':
                file.close()
                return jj + 5 
                break 
    file.close()
    print('ran off end of file looking for eigenlevels')
    sys.exit()

    return 0 

def read_in_eigenlevels(grasp_out_path,position_of_eigen_levels,num_rcsf):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        level_data = np.loadtxt(grasp_out_path,skiprows=position_of_eigen_levels,max_rows=num_rcsf,dtype=str)

    return level_data

def produce_eigenstate_array(level_data,eigenvectors) -> list[eigenstate]:

    j_values = level_data[:,1]
    parities = level_data[:,2]
    energies = (level_data[:,-1]).astype(float)

    num_rcsfs = np.shape(eigenvectors)[0]

    eigenstate_array = []

    for jj in range(0,num_rcsfs):
        energy = energies[jj]
        jval = j_values[jj]
        eigenvector = eigenvectors[:,jj] 
        if parities[jj] =='even':
            even_boolean = True
        else:
            even_boolean = False 

        

        current_eigenstate = eigenstate(energy,jval,even_boolean,eigenvector,jj+1)
        eigenstate_array.append(current_eigenstate)
    
    return eigenstate_array 

def make_many_eigenstate_strings(eigenstate_array:list[eigenstate],num_to_be_shown) -> list[eigenstate]:
    for state in eigenstate_array:
        state.make_csf_string(num_to_be_shown)
    return eigenstate_array

def output_table(eiegenstates_array:list[eigenstate],user_chosen_num_levels=0):
    num = len(eiegenstates_array)

    if (user_chosen_num_levels != 0) and (user_chosen_num_levels < num):
        num = user_chosen_num_levels 

    header = 'Index,     Energy (Ry),   Parity,        J,     jj Level Composition'
    print(header)
    for jj in range(0,num):

        current_state = eiegenstates_array[jj]
        energy = current_state.energy_ryd_unshifted

        level = current_state.level_number
        parity = current_state.even
        if current_state.even:
            parity = 'even'
        else:
            parity = ' odd'
        angular_momentum = current_state.jvalue
        csf_string = current_state.label_string
            
        output_string = '{:5},  {:14E},     {:8},  {:4},    {}'.format(level,energy,parity,angular_momentum,csf_string)
        print(output_string)