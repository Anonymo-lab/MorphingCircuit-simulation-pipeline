# -*- coding: mbcs -*-
# Do not delete the following import lines
import os
import copy
import math as m

from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import __main__


def readTxTToLine(read_path):
    '''
    Read the file with txt format into lines. This function is for 
    material definition extraction and conversion.
    
    Parameters:
        read_path: String. The path of the specified .txt file.
    Returns:
        lines: List of strings. The list of lines of the file's content. 
    '''
    
    with open(read_path, 'rt') as f:
        lines = f.read().splitlines()
    return lines


def genMatDefTuple(mat_file_path):
    '''
    Read the material definition file (in comma delimited file), and convert the data into the tuple structure 
    for definition in Abaqus. 
    
    Parameters:
        mat_file_path: String. The path of the specified material definition data file.
    Returns:
        mat_tuple: Tuple. Containing sub-tuples, which is one single line of the material data. 
    '''
    
    mat_tuple = []
    lines_temp = readTxTToLine(mat_file_path)
    for line in lines_temp:
        one_line_list = line.split(',')
        tuple_temp = [float(item) for item in one_line_list]
        mat_tuple.append(tuple(tuple_temp))
    
    return tuple(mat_tuple)


def importPart(directory, sat_file_name):
    '''
    Import the geometry from the .sat files
    Parameters:
        directory: String. The main directory. 
        sat_file_name: String. The name of the .sat file (one of the CAD_Model_file_list). 
    Return: part_name: String. The name of the part. 
    '''
    real_path = os.path.join(directory, sat_file_name)
    acis = mdb.openAcis(real_path, scaleFromFile=OFF)
    part_name = sat_file_name.split('.')[0]
    mdb.models['Model-1'].PartFromGeometryFile(name=part_name, geometryFile=acis, 
        combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts[part_name]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    
    return part_name

def computeDistance(coord_1, coord_2):
    '''
    Compute the euclidean distance between two points. 
    Parameters: 
        coord_1: Tuple of three float(s). The coordinates (X, Y, Z) of the first point. 
        coord_2: Tuple of three float(s). The coordinates (X, Y, Z) of the second point. 
    Return: 
        dist: the distance result. 
    '''
    
    dist = m.sqrt(abs(coord_1[0] - coord_2[0])**2 + abs(coord_1[1] - coord_2[1])**2 + abs(coord_1[2] - coord_2[2])**2)
    return dist


def findClosestFace(target_point, Assembly, instance_name_list):
    '''
    Used in boundary condition settings. Find the closest face from a specified point. 
    Parameters: 
        point_coord: Tuple of three float(s). The 3d coordinates of the target point (default: (0., 0., 0.))
        Assembly: abaqus object. The target assembly (containing several instances). 
        instance_name_list: List of strings. The names of all the instances in the input assembly. 
    Return:
        closest_face: GeomSequence. The closest face from the assembly to the target point. 
    '''
    
    all_faces, dist_temp, closest_face = None, None, None # GeomSequence, Float, GeomSequence
    for index, name in enumerate(instance_name_list): 
        if index == 0: all_faces = Assembly.instances[name].faces
        else: all_faces += Assembly.instances[name].faces
    
    for instance_name in instance_name_list:
        faces_temp = Assembly.instances[instance_name].faces
        r = faces_temp.getClosest(coordinates=(target_point,))
        if dist_temp == None: 
            if r == {}: continue
            else:
                dist_temp = computeDistance(target_point, r[0][1])
                closest_face = all_faces[all_faces.index(r[0][0]):all_faces.index(r[0][0])+1]
        else: 
            if r == {} or dist_temp < computeDistance(target_point, r[0][1]): continue
            else: 
                dist_temp = computeDistance(target_point, r[0][1])
                closest_face = all_faces[all_faces.index(r[0][0]):all_faces.index(r[0][0])+1]
    
    return closest_face

def computeStress(alpha, residual_stress):
    '''
    Used in the initial condition definition part. 
    Parameters:
        alpha: Float. The direction of the printing (given in the part_name). Unit: radians
        residual_stress: Float. The embedded pre-stress of the material. Unit: MPa
    Return:
        sigma_11: the first item of the stress definition.
        sigma_22: the second item of the stress definition. 
        tau_12: the third item of the stress definition. 
    '''
    
    cos_alpha = m.cos(alpha)
    sin_alpha = m.sin(alpha)
    sigma_11 = residual_stress * (cos_alpha**2)
    sigma_22 = residual_stress - sigma_11
    tau_12 = -residual_stress * sin_alpha * cos_alpha
    
    return sigma_11, sigma_22, tau_12


def modelling(directory, CAD_Model_file_list, mat_file_list, density, Poissons_ratio, residual_stress, gravity, task_name,
              fixing_point=(0.1, 0.1, 0.), mesh_size=1.0):
    '''
    Build the model for the FEA simulation.
    Parameters:
        directory: String. The main folder of the files.
        CAD_Model_file_list: List of strings. The list of file names of all .sat files. 
        mat_file_list: List of strings. The list of file names of all .txt files.
        density: Float. The density of the material. From the experiments. 
        Poissons_ratio: Float. The Poisson's Ratio of the material. From the experiments. 
        residual_stress: Float. The embedded stress of the printing process. 
        gravity: Float. The gravity of the triggering environment. 
        task_name: String. The name of the task. User can personalize this name. 
        fixing_point: Tuple of three float(s). The point of fixation. Based on the operation of the user, this poinbt should be (0., 0., 0.) (as the default value).
        mesh_size: Float. The size of the elemnt. it will control the speed of the simulation. 
    '''
    
    # Import parts in the FEA model. 
    
    part_name_list = []
    for sat_file in CAD_Model_file_list:
        part_name_temp = importPart(directory, sat_file)
        part_name_list.append(part_name_temp)
    
    # Define Materials.
    
    for file in mat_file_list:
        if file.split('_')[0] == "hyper": hyper_tuple = genMatDefTuple(os.path.join(directory, file))
        elif file.split('_')[0] == "vis": visco_tuple = genMatDefTuple(os.path.join(directory, file))
        else: continue
    
    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Density(table=((density, ), ))
    mdb.models['Model-1'].materials['Material-1'].Hyperelastic(materialType=ISOTROPIC, 
        type=MARLOW, volumetricResponse=POISSON_RATIO, poissonRatio=Poissons_ratio, 
        table=())
    mdb.models['Model-1'].materials['Material-1'].hyperelastic.UniaxialTestData(table=hyper_tuple)
    mdb.models['Model-1'].materials['Material-1'].Viscoelastic(domain=FREQUENCY, 
        frequency=TABULAR, type=ISOTROPIC, preload=UNIAXIAL, table=visco_tuple)
    mdb.models['Model-1'].HomogeneousSolidSection(name='Material-1', material='Material-1', 
        thickness=None)
    
    for name in part_name_list:
        p = mdb.models['Model-1'].parts[name]
        session.viewports['Viewport: 1'].setValues(displayedObject=p)
        c = p.cells
        region = p.Set(cells=c, name='Set-1')
        p = mdb.models['Model-1'].parts[name]
        p.SectionAssignment(region=region, sectionName='Material-1', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)
    
    # Assemble the parts
    
    a = mdb.models['Model-1'].rootAssembly
    
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a1 = mdb.models['Model-1'].rootAssembly
    a1.DatumCsysByDefault(CARTESIAN)
    instance_name_list = []
    for name in part_name_list:
        instance_name_temp = name + "-1"
        p = mdb.models['Model-1'].parts[name]
        a1.Instance(name=instance_name_temp, part=p, dependent=ON)
        instance_name_list.append(instance_name_temp)

    # Define step

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        adaptiveMeshConstraints=ON)
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', 
        initialInc=0.01, minInc=1e-20, nlgeom=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON, 
        adaptiveMeshConstraints=OFF)
    
    # Define interaction    
    
    a = mdb.models['Model-1'].rootAssembly
    
    all_faces, finished_faces, non_contacting_faces = None, None, None
    for index, name in enumerate(instance_name_list): 
        if index == 0: all_faces = a.instances[name].faces
        else: all_faces += a.instances[name].faces

    count_temp = 1
    for face_temp in all_faces:
        if finished_faces != None and (face_temp in finished_faces): continue
        other_face_list = [item for item in all_faces if item != face_temp]
        other_centroid_list = [tuple(item.getCentroid()) for item in other_face_list]
        
        centroid_temp = tuple(face_temp.getCentroid())
        if centroid_temp in other_centroid_list: 
            pairing_face_temp = other_face_list[other_centroid_list.index(centroid_temp)]
            face_temp_name, pairing_face_name = "Pair_{}_1".format(str(count_temp)), "Pair_{}_2".format(str(count_temp))
            
            a.Surface(side1Faces=all_faces[all_faces.index(face_temp):all_faces.index(face_temp)+1], name=face_temp_name)
            a.Surface(side1Faces=all_faces[all_faces.index(pairing_face_temp):all_faces.index(pairing_face_temp)+1], name=pairing_face_name)
            
            region_face_temp = a.surfaces[face_temp_name]
            region_pairing_temp = a.surfaces[pairing_face_name]
            
            tie_name_temp = "Tie_{}".format(str(count_temp))
            mdb.models['Model-1'].Tie(name=tie_name_temp, master=region_face_temp, 
                slave=region_pairing_temp, positionToleranceMethod=COMPUTED, adjust=ON, 
                constraintEnforcement=SURFACE_TO_SURFACE)
            
            if count_temp == 1: 
                finished_faces = (all_faces[all_faces.index(face_temp):all_faces.index(face_temp)+1] + 
                                  all_faces[all_faces.index(pairing_face_temp):all_faces.index(pairing_face_temp)+1])
            else:
                finished_faces += all_faces[all_faces.index(face_temp):all_faces.index(face_temp)+1]
                finished_faces += all_faces[all_faces.index(pairing_face_temp):all_faces.index(pairing_face_temp)+1]
            
            count_temp += 1
            
        else: 
            if non_contacting_faces == None: non_contacting_faces = all_faces[all_faces.index(face_temp):all_faces.index(face_temp)+1]
            else: non_contacting_faces += all_faces[all_faces.index(face_temp):all_faces.index(face_temp)+1]
            continue

    # Define boundary condition

    a = mdb.models['Model-1'].rootAssembly
    
    closest_face = findClosestFace(fixing_point, a, instance_name_list) # GeomSequence
    a.Set(faces=closest_face, name='BC')

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, interactions=OFF, constraints=OFF, 
        engineeringFeatures=OFF)
    mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1', 
        comp3=gravity, distributionType=UNIFORM, field='')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
    a = mdb.models['Model-1'].rootAssembly
    region = a.sets['BC']
    mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
        region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    
    # Define initial condition
    
    a = mdb.models['Model-1'].rootAssembly
    
    for index, instance_name in enumerate(instance_name_list):
        IC_name_temp = "IC_{}".format(str(index))
        
        alpha = float(instance_name.split('_')[-1].split('-')[0])
        alpha_radians = alpha * m.pi / 180.0
        
        sigma_11_temp, sigma_22_temp, tau_12_temp = computeStress(alpha_radians, residual_stress)
        
        region = a.instances[instance_name].sets["Set-1"]
        mdb.models['Model-1'].Stress(name=IC_name_temp, region=region, 
            distributionType=UNIFORM, sigma11=sigma_22_temp, sigma22=0.0, sigma33=sigma_11_temp, 
            sigma12=0, sigma13=tau_12_temp, sigma23=0.0)
    
    # Meshing
    
    for part_name in part_name_list:
        p = mdb.models['Model-1'].parts[part_name]
        p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
        p.generateMesh()

    a1 = mdb.models['Model-1'].rootAssembly
    a1.regenerate()
    
    # Generate model file
    
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    mdb.Job(name=task_name, model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numDomains=1, numGPUs=0)
    mdb.jobs[task_name].writeInput(consistencyChecking=OFF)


def simulationFEA(task_name, cpu_processors=8):
    '''
    Call the FEA solver in Abaqus, and conduct the simulation. 
    Parameters:
        task_name: String. The name of the input file generated from the "Modelling" algorithm. Also the user-defined task name.    
        cpu_processors: Int. The number of processors of the local CPU. This argument should be modified if different number of local processors are used. Default: 8.   
    '''
    mdb.JobFromInputFile(name=task_name, inputFileName=task_name, type=ANALYSIS, atTime=None, 
        waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, 
        userSubroutine='', scratch='', multiprocessingMode=DEFAULT, numCpus=cpu_processors, numDomains=cpu_processors, numGPUs=0)
    mdb.jobs[task_name].submit()
    mdb.jobs[task_name].waitForCompletion()


def sepFiles(directory):
    '''
    Categorize files in the directory and return files in different types. 
    Parameters:
        directory: String. The path of main directory
    Return:
        file_list_sat: List of String. Containing all the files with the .sat format. 
        file_list_py: List of String. Containing all the files with the .py format. 
        file_list_txt: List of String. Containing all the files with the .txt format. 
    '''
    
    file_list = os.listdir(directory)
    file_list_sat, file_list_py, file_list_txt = [], [], []
    for file in file_list:
        real_path = os.path.join(directory, file)
        if os.path.isdir(real_path): continue
        
        if file.split('.')[-1] == "sat": file_list_sat.append(file)
        elif file.split('.')[-1] == "py": file_list_py.append(file)
        elif file.split('.')[-1] == "txt": file_list_txt.append(file)
        else: continue
    
    return file_list_sat, file_list_py, file_list_txt


def main():
    # Note: the material definition files, the python code files (this file), and the CAD geometry files (.sat files) need to be in this folder. 
    directory = '/'.join(os.path.abspath("simulation_pipeline.py").split('\\')[:-1]) # Users need to change the directory of the simulation folder. 
    # directory = os.path.dirname(os.path.abspath(__file__)) # Users need to change the directory of the simulation folder.
    print(directory)
    os.chdir(directory)
    CAD_Model_file_list, code_file_list, mat_file_list = sepFiles(directory)
    
    density = 1.25e-9 # The density of PLA. Unit: t*mm^-3
    Poissons_ratio = 0.406 # The Poisson's ratio from the DMA experiments. 
    residual_stress = 0.158 # The residual stress embedded during printing. Unit: MPa. Obtained from the DMA experiment. 
    
    gravity = -9800 # Unit: mm*s^-2
    meshSize = 1.0
    task_name = "Job-1"
    
    modelling(directory, CAD_Model_file_list, mat_file_list, density, Poissons_ratio, 
              residual_stress, gravity, task_name, mesh_size=meshSize)
    simulationFEA(task_name)
    
    file_list = os.listdir(directory)
    for file in file_list:
        if file.split('.')[-1] == "odb": return os.path.join(directory, file)
        else: continue
    
    return "No result is generated. "

if __name__ == "__main__":
    result_file_path = main()
    print(result_file_path)

