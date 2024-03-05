# grasp_levels_to_table
usage: python3 (your path to)/grasp_ls_table.py 

-f (your grasp.out, required - in at least output mode 1)

 -o (desired num orbitals in csfs, default 2, cannot print less than 2. not sure why yet) 
-n (num levels out, default all) 

-i (requires output mode 2 or higher in grasp. just a flag, prints term of all the orbitals except the outer one if used. e.g for 3s3p 3d it will print 3s3p(term of all this) 3d(term of total state) without: 3s3p3d (term) )

