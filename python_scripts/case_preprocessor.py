from vtk import *
import pyvista
import os
import struct
import gc
import time

from vtkmodules.vtkIOEnSight import *
from vtkmodules.util.misc import *
from vtkmodules.vtkFiltersPoints import *

from vtkmodules.all import *

VTK_DATA_ROOT = vtkGetDataRoot()

start = time.time()



print('reading raw case file',time.time()-start)
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


print('cleaning case file',time.time()-start)
# general cleanup
algo = vtkStaticCleanUnstructuredGrid()
algo.SetInputData(polyblock)
algo.RemoveUnusedPointsOn()
algo.Update()
polyblock = algo.GetOutput()
del algo


print('converting element types to polyhedra',time.time()-start)
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


print('extracting point and element counts',time.time()-start)
# Pull general mesh info
npoin = polyblock.GetNumberOfPoints()
nelem = polyblock.GetNumberOfCells()
print(npoin,' nodes ', nelem,' cells')

gc.collect()

print('initialising output arrays',time.time()-start)
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


print('extracting coordinate data',time.time()-start)
# extract coordinate data

for p in range(npoin):
    x_coords.append(polyblock.GetPoint(p)[0])
    y_coords.append(polyblock.GetPoint(p)[1])
    z_coords.append(polyblock.GetPoint(p)[2])
    # node: these will include the duplicated points from the boundary data



# extract cell, face and edge datas from the primary data block
faceid = 0
edgeid = 0
print('extracting connectivity and mapping data',time.time()-start)
for c in range(nelem):
    
    cell = polyblock.GetCell(c)
    
    c_p_obj_relation_array.append([cell.GetPointId(p) for p in range(cell.GetNumberOfPoints())])
    
    
    nfaces = cell.GetNumberOfFaces()
    
    c_f_obj_relation_array.append( list(range(faceid,(faceid+nfaces-1))))
    
    for f in range(cell.GetNumberOfFaces()): 
        
        face = cell.GetFace(f)
        
        f_p_obj_relation_array.append([face.GetPointId(p) for p in range(face.GetNumberOfPoints())])
        
        nedges = face.GetNumberOfEdges()
        f_e_obj_relation_array.append( list(range(edgeid,(edgeid+nedges-1))))
    
        
        for e in range(nedges):
            
            edge = face.GetEdge(e)
            e_p_obj_relation_array.append([edge.GetPointId(0),edge.GetPointId(1)])
            
        edgeid = edgeid + nedges - 1
    
    faceid = faceid + nfaces - 1

del polyblock
del cell
del face
del edge

nface = faceid
nedge = edgeid

gc.collect()

print('sorting relation arrays',time.time()-start)

for i in range(len(c_p_obj_relation_array)):
    c_p_obj_relation_array[i]=sorted([c_p_obj_relation_array[i]])
    
for i in range(len(f_p_obj_relation_array)):
    f_p_obj_relation_array[i]=sorted([f_p_obj_relation_array[i]])
    
for i in range(len(e_p_obj_relation_array)):
    e_p_obj_relation_array[i]=sorted([e_p_obj_relation_array[i]])
    
for i in range(len(c_f_obj_relation_array)):
    c_f_obj_relation_array[i]=sorted([c_f_obj_relation_array[i]])
    
for i in range(len(f_e_obj_relation_array)):
    f_e_obj_relation_array[i]=sorted([f_e_obj_relation_array[i]])

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

print('generating index arrays',time.time()-start)
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


print('removing duplicate connectivites',time.time()-start)

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

print('indexifying index arrays',time.time()-start)
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

print('updating counts',time.time()-start)
# update counts with non duplicate objects
nele  = len(c_p_index)
nface = len(f_p_index)
nedge = len(e_p_index)


print('beginning writing to file',time.time()-start)
# outputting


print('creating raw file',time.time()-start)
# opening/creating file
file_name = '../preprocessed_mesh_folder/raw_mesh_data.preprocessed_mesh_file'
file = open(file_name, "wb")


print('writing header',time.time()-start)
# writing header
file.write(struct.pack('<4i' ,npoin,nedge,nface,nelem))
print(npoin,nedge,nface,nelem)


print('writing coordinates',time.time()-start)
# writing coordinate data
for coord in x_coords:
    file.write(struct.pack('<d' ,coord))
for coord in y_coords:
    file.write(struct.pack('<d' ,coord))
for coord in z_coords:
    file.write(struct.pack('<d' ,coord))
    
    
print('writing connectivities',time.time()-start)
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
    
    
print('writing cell > face, face > edge mappings',time.time()-start)
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
    
    
print('finished writing to file',time.time()-start)
print('finished preprocessing')


