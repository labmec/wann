# to change strings in files:
#https://www.geeksforgeeks.org/how-to-search-and-replace-text-in-a-file-in-python/

import random
import json
import subprocess
import os

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

class SimulationParameters:
    def __init__(self):
        self.MeshData = {
            "file": "mesh3D_rev08.msh",
            "NumUniformRef": 0,
            "NumDirRef": 0,
            "ToCylindrical": 1
        }
        self.WellboreData = {
            "name": "curve_wellbore",
            "perm": 10.0,
            "pOrder": 2,
            "radius": 0.1,
            "length": 1.0,
            "eccentricity": [0.0, 0.0, 0.0],
            "BCs": [
                {"name": "point_heel", "type": 0, "value": 2.0},
                {"name": "point_toe", "type": 1, "value": 0.0}
            ]
        }
        self.ReservoirData = {
            "name": "volume_reservoir",
            "perm": 1.0,
            "porosity": 1.0,
            "pOrder": 1,
            "height": 1,
            "width": 1,
            "length": 1,
            "BCs": [
                {"name": "surface_farfield", "type": 0, "value": 1.0},
                {"name": "surface_wellbore_heel", "type": 1, "value": 0.0},
                {"name": "surface_wellbore_toe", "type": 1, "value": 0.0}
            ]
        }
        self.FluidData = {
            "name": "oil",
            "viscosity": 0.005,
            "density": 800.0
        }
        self.PostProcData = {
            "wellbore_vtk": "postproc_wellbore",
            "reservoir_vtk": "postproc_reservoir",
            "training_data": "training_data.txt",
            "vtk_resolution": 0,
            "training_resolution": 400,
            "nthreads": 0
        }
        self.VerbosityLevel = 0

    def to_dict(self):
        return {
            "MeshData": self.MeshData,
            "WellboreData": self.WellboreData,
            "ReservoirData": self.ReservoirData,
            "FluidData": self.FluidData,
            "PostProcData": self.PostProcData,
            "VerbosityLevel": self.VerbosityLevel
        }

    def write_to_json(self, json_file: str):
        with open(json_file, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)
        print(f"Input written to {json_file}")

def get_geometric_parameters(case): #{{{
    # range of geometric parameter
    # wellbore:
    #   Lw [400 1200] m
    #   Dw [0.127 0.1778] m // 5" to 7"
    #   ecc [-90 90] % wellbore eccentricity 0: no eccentricity, -90%: very close to the reservoir bottom, 90% very close to the reservoir top   
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

def search_root_directory():
    current_dir = os.getcwd()
    for _ in range(2):  # Check current and parent directory
        if "pyscripts" in os.listdir(current_dir):
            print(f"'root' directory found in: {current_dir}")
            break
        current_dir = os.path.dirname(current_dir)
    return current_dir

root_dir = search_root_directory()
input_dir = os.path.join(root_dir, "input")
geo_dir = os.path.join(root_dir, "geo")

number_of_cases = 1

for case in range(1,number_of_cases+1):
    print("generating files for case "+str(case))
    basemesh = geo_dir + "/mesh3D_rev08.geo"
    casegeo = "case_"+str(case)+".geo"
    casemesh = "case_"+str(case)+".msh"
    caseparam = "case_"+str(case)+".json"
    
    print("   get geometric and physical parameters")
    Lw, Dw, Lr, Wr, Hr = get_geometric_parameters(case)
    params = SimulationParameters()
    params.WellboreData["length"] = Lw
    params.WellboreData["radius"] = Dw / 2
    params.ReservoirData["length"] = Lr
    params.ReservoirData["width"] = Wr
    params.ReservoirData["height"] = Hr
    params.MeshData["file"] = casemesh

    print("   save .json")
    params.write_to_json(f"{input_dir}/{caseparam}")

    print("   open base mesh file and replace geometric parameters")
    with open(basemesh, 'r') as file:
        data = file.read()
        # data = data.replace("Lw = 1",   "Lw = "+str(Lw))
        # data = data.replace("Dw = 0.1", "Dw = "+str(Dw))
        # data = data.replace("Lr = 4",   "Lr = "+str(Lr))
        # data = data.replace("Wr = 3",   "Wr = "+str(Wr))
        # data = data.replace("Hr = 1",   "Hr = "+str(Hr))
        file.close()

    print("   save .geo")
    with open(geo_dir+"/"+casegeo, 'w') as file: 
        file.write(data) 
        file.close()
    
    print("   call gmsh to generate the mesh")
    subprocess.run(["gmsh", "-3", "-o", input_dir+"/"+casemesh, geo_dir+"/"+casegeo, "Mesh.Coherence"], check=True)
    
    # print("   save .txt")
    # with open(caseparam, 'w') as file: 
    #     dataparam=str(Lw) + "\t" + str(Dw) + "\t" + str(Lr) + "\t" + str(Wr) + "\t" + str(Hr) + "\n"
    #     file.write(dataparam) 
    #     file.close()




#@import fileinput
#with fileinput.FileInput(basemesh, inplace=True, backup='.bak') as file:
#    for line in file:
#        print(line.replace("Lw = 1", "Lw = 2"), end='')
