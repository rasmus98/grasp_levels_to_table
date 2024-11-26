#import numpy as np  
import f90nml
file = 'DSTG2.OUT'

def read_stg2_input_return_csfs(file):
    
    file_opened = open(file,'r')
    


    found = False
    num_orbitals = 0
    num_nrcsfs = 0

    while found == False:
        lineread = file_opened.readline()
        linereadsplit = lineread.split()


        if len(linereadsplit) > 0:
            if linereadsplit[0] == '&DSTG2':
                found = True
                tmp_namelist = open('tmp.nml','w')
                tmp_namelist.write(lineread)
                tmp_namelist.close()
                namelist = f90nml.read('tmp.nml')
    print(namelist['DSTG2']['NWM'])
    lineread = file_opened.readline().split()
    print(lineread)
                    
    orb = f90nml.read('ORGB.NML')
    print(orb['ORB'])

    

    file_opened.close()
    return 0

read_stg2_input_return_csfs(file)