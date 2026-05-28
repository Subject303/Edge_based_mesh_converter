from vtk import *
import pyvista
import os
import struct
import gc

from vtkmodules.vtkIOEnSight import *
from vtkmodules.util.misc import *
from vtkmodules.vtkFiltersPoints import *

from vtkmodules.all import *

VTK_DATA_ROOT = vtkGetDataRoot()

# set angle to split boundary edges
feature_angle = 15



print('reading raw case file')
# read the raw data
algo = vtkGenericEnSightReader()
cdp = vtkCompositeDataPipeline()
# Make sure all algorithms use the composite data pipeline
# might not be needed
algo.SetDefaultExecutivePrototype(cdp)
del cdp
algo.SetCaseFileName("../case/star.case") 
# will want to change to somthing more generic
algo.Update()
raw_data = algo.GetOutput()
del algo


# block 0 is always the primary datablock here
polyblock = raw_data.GetBlock(0)  # vtkUnstructuredGrid format


del raw_data
# freeing up the multiblock view


print('cleaning case file')
# general cleanup
algo = vtkStaticCleanUnstructuredGrid()
algo.SetInputData(polyblock)
algo.RemoveUnusedPointsOn()
algo.Update()
polyblock = algo.GetOutput()
del algo


print('converting element types to polyhedra')
# convert all cells to polyhedra
algo = vtkConvertToPolyhedra()
algo.SetInputData(polyblock)
algo.OutputAllCellsOn()
algo.Update()
polyblock = algo.GetOutput()
del algo


# # this splits the points on feature edges
# # we do this before assigning global IDs so the IDs are seperate
# algo = vtkSplitSharpEdgesPolyData()
# algo.SetInputData(polyblock)
# algo.SetFeatureAngle(feature_angle)
# algo.Update()
# polyblock = algo.GetOutput()
# del algo


# print('generating global IDs')
# # create global point IDs that I can reference down the line
# algo = vtkGenerateGlobalIds()
# algo.SetInputData(polyblock)
# algo.Update()
# polyblock = algo.GetOutput()
# del algo


print('extracting point and element counts')
# Pull general mesh info
npoin = polyblock.GetNumberOfPoints()
nelem = polyblock.GetNumberOfCells()
print(npoin,' nodes ', nelem,' cells')

gc.collect()

print('initialising output arrays')
# initialise the output arrays
x_coords = []
y_coords = []
z_coords = []

c_p_sum = 0
f_p_sum = 0
e_p_sum = 0

c_p_index = []
f_p_index = []
e_p_index = [] 

c_p_obj_relation_array = []
f_p_obj_relation_array = []
e_p_obj_relation_array = []

c_f_sum = 0
f_e_sum = 0

c_f_index = []
f_e_index = []

c_f_obj_relation_array = []
f_e_obj_relation_array = []


print('extracting coordinate data')
# extract coordinate data

for p in range(npoin):
    x_coords.append(polyblock.GetPoint(p)[0])
    y_coords.append(polyblock.GetPoint(p)[1])
    z_coords.append(polyblock.GetPoint(p)[2])
    # node: these will include the duplicated points from the boundary data



# extract cell, face and edge datas from the primary data block
faceid = -1
edgeid = -1
print('extracting connectivity and mapping data')
for c in range(nelem):
    
    print(nelem-c)
    
    cell = polyblock.GetCell(c)
    
    temp = []
    for p in range(cell.GetNumberOfPoints()):
        temp.append(cell.GetPointId(p))
    c_p_obj_relation_array.append(sorted(temp))
    temp2=[]
    for f in range(cell.GetNumberOfFaces()): 
        faceid = faceid + 1
        temp2.append(faceid)
        
        face = cell.GetFace(f)
        
        temp = []
        for p in range(face.GetNumberOfPoints()):
            temp.append(cell.GetPointId(p))
        f_p_obj_relation_array.append(sorted(temp))
        
        temp3 = []
        for e in range(face.GetNumberOfEdges()):
            edgeid = edgeid + 1
            temp3.append(edgeid)
            edge = face.GetEdge(e)
            temp = []
            temp.append(edge.GetPointId(0))
            temp.append(edge.GetPointId(1))
            f_p_obj_relation_array.append(sorted(temp))
        f_e_obj_relation_array.append(temp3)
    c_f_obj_relation_array.append(temp2)

del polyblock
del temp2
del temp
del cell
del face
del edge

nface = faceid
nedge = edgeid

gc.collect()

print('sorting relation arrays')
print('   c_p')
c_p_obj_relation_array.sort()
print('   f_p')
f_p_obj_relation_array.sort()
print('   e_p')
e_p_obj_relation_array.sort()

print('   c_f')
c_f_obj_relation_array.sort()
print('   f_e')
f_e_obj_relation_array.sort()

print('generating index arrays')
for obj in c_p_obj_relation_array:
    c_p_index.append(len(obj))
    
for obj in f_p_obj_relation_array:
    f_p_index.append(len(obj))
    
for obj in e_p_obj_relation_array:
    e_p_index.append(len(obj))
    
for obj in c_f_obj_relation_array:
    c_f_index.append(len(obj))
    
for obj in f_e_obj_relation_array:
    f_e_index.append(len(obj))


# ok I should Have all data I need to output now.
# and indexify the index arrays
# and sort the relation arrays and delete duplicated stuff
# I just need to filter out boundary duplicated nodes


print('removing duplicate connectivites')

temp = []
temp2 = []
print('   c_p')
for i in range(len(c_p_index)-1):
    if c_p_index[i] == c_p_index[i+1]:    
        obj1 = c_p_obj_relation_array[i]
        obj2 = c_p_obj_relation_array[i+1]
        for p in range(c_p_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(c_p_index[i])
                break
            
c_p_obj_relation_array = temp
c_p_index = temp2
del temp
del temp2

temp = []
temp2 = []
print('   f_p')
for i in range(len(f_p_index)-1):
    if f_p_index[i] == f_p_index[i+1]:    
        obj1 = f_p_obj_relation_array[i]
        obj2 = f_p_obj_relation_array[i+1]
        for p in range(f_p_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(f_p_index[i])
                break
            
f_p_obj_relation_array = temp 
f_p_index = temp2
del temp
del temp2

temp = []
temp2 = []
print('   e_p')
for i in range(len(e_p_index)-1):
    if e_p_index[i] == e_p_index[i+1]:    
        obj1 = e_p_obj_relation_array[i]
        obj2 = e_p_obj_relation_array[i+1]
        for p in range(e_p_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(2)
                break
        
e_p_obj_relation_array = temp
e_p_index = temp2
del temp
del temp2

temp = []
temp2 = []
print('   c_f')
for i in range(len(c_f_index)-1):
    if c_f_index[i] == c_f_index[i+1]:    
        obj1 = c_f_obj_relation_array[i]
        obj2 = c_f_obj_relation_array[i+1]
        for p in range(c_f_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(c_f_index[i])
                break
            
c_f_obj_relation_array = temp
c_f_index = temp2
del temp
del temp2

temp = []
temp2 = []
print('   f_e')
for i in range(len(f_e_index)-1):
    if f_e_index[i] == f_e_index[i+1]:    
        obj1 = f_e_obj_relation_array[i]
        obj2 = f_e_obj_relation_array[i+1]
        for p in range(f_e_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(f_e_index[i])
                break
            
f_e_obj_relation_array = temp
f_e_index = temp2
del temp
del temp2


gc.collect()

print('indexifying index arrays')
# indexifying the index arrays
index_last = 0
for index in range(len(c_p_index)):
    c_p_index[index] = index_last + c_p_index[index]
    index_last = c_p_index[index]

c_p_sum = index_last

index_last = 0
for index in range(len(f_p_index)):
    f_p_index[index] = index_last + f_p_index[index]
    index_last = f_p_index[index]
    
f_p_sum = index_last

index_last = 0
for index in range(len(e_p_index)):
    e_p_index[index] = index_last + e_p_index[index]
    index_last = e_p_index[index]
    
e_p_sum = index_last

index_last = 0
for index in range(len(c_f_index)):
    c_f_index[index] = index_last + c_f_index[index]
    index_last = c_f_index[index]
    
c_f_sum = index_last

index_last = 0
for index in range(len(f_e_index)):
    f_e_index[index] = index_last + f_e_index[index]
    index_last = f_e_index[index]

f_e_sum = index_last

print('updating counts')
# update counts with non duplicate objects
nele  = len(c_p_index)
nface = len(f_p_index)
nedge = len(e_p_index)


print('beginning writing to file')
# outputting


print('creating raw file')
# opening/creating file
file_name = '../preprocessed_mesh_folder/raw_mesh_data.preprocessed_mesh_file'
file = open(file_name, "wb")


print('writing header')
# writing header
file.write(struct.pack('<4i' ,npoin,nedge,nface,nelem))
print(npoin,nedge,nface,nelem)


print('writing coordinates')
# writing coordinate data
for coord in x_coords:
    file.write(struct.pack('<d' ,coord))
for coord in y_coords:
    file.write(struct.pack('<d' ,coord))
for coord in z_coords:
    file.write(struct.pack('<d' ,coord))
    
    
print('writing connectivities')
# writing object connectivities
file.write(struct.pack('<2i' ,c_p_sum, 0))
for entry in c_p_index:
    file.write(struct.pack('<i' ,entry))
    
file.write(struct.pack('<2i' ,f_p_sum, 0))
for entry in f_p_index:
    file.write(struct.pack('<i' ,entry))


for obj in c_p_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry))

for obj in f_p_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry))
    
for obj in e_p_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry))
    
    
print('writing cell > face, face > edge mappings')
# writing object relation mappings
file.write(struct.pack('<2i' ,c_f_sum, 0))
for entry in c_f_index:
    file.write(struct.pack('<i' ,entry))
    
file.write(struct.pack('<2i' ,f_e_sum, 0))
for entry in f_e_index:
    file.write(struct.pack('<i' ,entry))
    
    
for obj in c_f_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry))

for obj in f_e_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry))
    
    
print('finished writing to file')
print('finished preprocessing')


