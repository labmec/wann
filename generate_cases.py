# to change strings in files:
#https://www.geeksforgeeks.org/how-to-search-and-replace-text-in-a-file-in-python/

import random

# 1) generate all .geo files
# 2) generate input files (parameters)
# 3) call gmsh to generate the meshes

# 1) generate all .geo files
# in a loop:
#   read basemesh
#   copy and save to another file (case i)
#   replace the targets
#   close the file

#//wellbore dimensions as in the basemesh
#Lw = 1;   //[m] wellbore length //400 m
#Dw = 0.1; //[m] wellbore diameter //0.1 m
#
#//reservoir dimensions as in the basemesh
#Lr = 4; //[m] length //2000 m
#Wr = 3; //[m] width  //1000 m
#Hr = 1; //[m] height //100 m

def get_geometric_parameters(case): #{{{
    # range of geometric parameter
    # wellbore:
    #   Lw [400 1200] m
    #   Dw [0.127 0.1778] m // 5" to 7"
    # 
    # reservoir:
    #   Lr [1.5 3.0] x Lw 
    #   Wr [1.0 3.0] x Lw
    #   Hr [15 60] m
    #
    # constraints:
    #
    
    # random.seed(a=0) # set seed to reproduce the same cases

    #Lw = random.uniform(400., 600.)
    #Lw = random.randint(400, 600)
    Lw = random.randrange(100, 250, 50) 
    Dw = random.uniform(0.127, 0.1778)

    Lr = round(random.uniform(1.2, 2.0) * Lw)
    Wr = round(random.uniform(1.0, 2.0) * Lw)
    #Hr = random.uniform(50, 100) 
    #Hr = random.randint(50, 100)
    Hr = random.randrange(50, 110, 10)

    if 0:
        Lw = 2
        Dw = 0.2
        Lr = 5
        Wr = 4
        Hr = 2

    return Lw, Dw, Lr, Wr, Hr
#}}}

number_of_cases = 4

for case in range(1,number_of_cases+1):
    print("generating files for case "+str(case))
    basemesh = "./geo/mesh3D_rev08.geo" 
    casemesh = "./cases/case_"+str(case)+".geo"
    caseparam= "./cases/case_"+str(case)+".txt"
    
    print("   get geometric and physical parameters")
    Lw, Dw, Lr, Wr, Hr = get_geometric_parameters(case)

    print("   open base mesh file and replace geometric parameters")
    with open(basemesh, 'r') as file:
        data = file.read()
        data = data.replace("Lw = 1",   "Lw = "+str(Lw))
        data = data.replace("Dw = 0.1", "Dw = "+str(Dw))
        data = data.replace("Lr = 4",   "Lr = "+str(Lr))
        data = data.replace("Wr = 3",   "Wr = "+str(Wr))
        data = data.replace("Hr = 1",   "Hr = "+str(Hr))
        file.close()

    print("   save .geo")
    with open(casemesh, 'w') as file: 
        file.write(data) 
        file.close()
    
    print("   save .txt")
    with open(caseparam, 'w') as file: 
        dataparam=str(Lw) + "\t" + str(Dw) + "\t" + str(Lr) + "\t" + str(Wr) + "\t" + str(Hr) + "\n"
        file.write(dataparam) 
        file.close()




#@import fileinput
#with fileinput.FileInput(basemesh, inplace=True, backup='.bak') as file:
#    for line in file:
#        print(line.replace("Lw = 1", "Lw = 2"), end='')
