from vtk import *
import pyvista
import os
import struct
import gc
import time
import sys
from vtkmodules.vtkIOEnSight import *
from vtkmodules.util.misc import *
from vtkmodules.vtkFiltersPoints import *

from vtkmodules.all import *

VTK_DATA_ROOT = vtkGetDataRoot()

start = time.time()



print('reading raw case file',time.time()-start); sys.stdout.flush()
# read the raw data
algo = vtkGenericEnSightReader()
cdp = vtkCompositeDataPipeline()
# Make sure all algorithms use the composite data pipeline
# might not be needed
algo.SetDefaultExecutivePrototype(cdp)
del cdp
# algo.SetCaseFileName("../case/2.case") 
algo.SetCaseFileName("../case/star1.case") 
# algo.SetCaseFileName("../case/star.case") 
# algo.SetCaseFileName("../case/2.case") 
# will want to change to somthing more generic
algo.Update()
raw_data = algo.GetOutput()
del algo


# block 0 is always the primary datablock here
polyblock = raw_data.GetBlock(0)  # vtkUnstructuredGrid format


del raw_data
# freeing up the multiblock view


print('cleaning case file',time.time()-start); sys.stdout.flush()
# general cleanup
algo = vtkStaticCleanUnstructuredGrid()
algo.SetInputData(polyblock)
algo.RemoveUnusedPointsOn()
algo.Update()
polyblock = algo.GetOutput()
del algo


print('converting element types to polyhedra',time.time()-start); sys.stdout.flush()
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


print('extracting point and element counts',time.time()-start); sys.stdout.flush()
# Pull general mesh info
npoin = polyblock.GetNumberOfPoints()
nele = polyblock.GetNumberOfCells()
print(npoin,' nodes ', nele,' cells')

gc.collect()

print('initialising output arrays',time.time()-start); sys.stdout.flush()
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


print('extracting coordinate data',time.time()-start); sys.stdout.flush()
# extract coordinate data

for p in range(npoin):
    x_coords.append(polyblock.GetPoint(p)[0])
    y_coords.append(polyblock.GetPoint(p)[1])
    z_coords.append(polyblock.GetPoint(p)[2])
    
    # node: these will include the duplicated points from the boundary data


# extract cell, face and edge datas from the primary data block
faceid = 0
edgeid = 0
print('extracting connectivity and mapping data',time.time()-start); sys.stdout.flush()
for c in range(nele):
    
    cell = polyblock.GetCell(c)
    
    c_p_obj_relation_array.append(sorted([cell.GetPointId(p) for p in range(cell.GetNumberOfPoints())]))
    
    
    nfaces = cell.GetNumberOfFaces()
    
    f_array = []
    
    for f in range(cell.GetNumberOfFaces()): 
        
        face = cell.GetFace(f)
        
        f_p_obj_relation_array.append(sorted([face.GetPointId(p) for p in range(face.GetNumberOfPoints())]))
        
        nedges = face.GetNumberOfEdges()
        
        e_array = []
        
        for e in range(nedges):
            
            edge = face.GetEdge(e)
            
            e_p_obj_relation_array.append(sorted([edge.GetPointId(0),edge.GetPointId(1)]))
            
            e_array.append(edgeid)
            edgeid = edgeid + 1
    
        f_array.append(faceid)
        faceid = faceid + 1
        
        f_e_obj_relation_array.append(e_array)
        
    c_f_obj_relation_array.append(f_array)
    

del polyblock
del cell
del face
del edge
    
nface = faceid
nedge = edgeid

gc.collect()

print('sorting connectivities',time.time()-start); sys.stdout.flush()

print('   c_p',time.time()-start); sys.stdout.flush()
# temp = list(range(len(c_p_obj_relation_array)))
# c_sort = temp

# coupled = zip(c_p_obj_relation_array, c_sort, c_f_obj_relation_array)
coupled = zip(c_p_obj_relation_array, c_f_obj_relation_array)

coupled = sorted(coupled, key=lambda element: element[0][0])

i=0
for x, z in coupled:
    c_p_obj_relation_array[i] = x
    # c_sort[i] = y
    c_f_obj_relation_array[i] = z
    i=i+1

del coupled

# del temp
# del c_sort

print('   f_p',time.time()-start); sys.stdout.flush()

temp = list(range(len(f_p_obj_relation_array)))
f_sort = temp

coupled = sorted(zip(f_p_obj_relation_array, f_sort, f_e_obj_relation_array), key=lambda element: element[0][0])

i=0
for x, y, z in coupled:
    f_p_obj_relation_array[i] = x
    f_sort[i] = y
    f_e_obj_relation_array[i] = z
    i=i+1

del coupled

coupled = sorted(zip(f_sort, temp), key=lambda element: element[0])

del temp

i=0
for x, y in coupled:
    f_sort[i] = x
    # temp[i] = y
    i=i+1
    
del coupled

for i in range(len(c_f_obj_relation_array)):
    obj = c_f_obj_relation_array[i]
    for j in range(len(obj)):
        obj[j] = f_sort[obj[j]] 
    c_f_obj_relation_array[i] = sorted(obj)

del f_sort


print('   e_p',time.time()-start); sys.stdout.flush()

temp = list(range(len(e_p_obj_relation_array)))
e_sort = temp

coupled = sorted(zip(e_p_obj_relation_array, e_sort), key=lambda element: element[0][0])

i=0
for x, y in coupled:
    e_p_obj_relation_array[i] = x
    e_sort[i] = y
    i=i+1

del coupled

coupled = sorted(zip(e_sort, temp), key=lambda element: element[0])

del temp

i=0
for x, y in coupled:
    e_sort[i] = x
    # temp[i] = y
    i=i+1
    
del coupled

for i in range(len(f_e_obj_relation_array)):
    obj = f_e_obj_relation_array[i]
    for j in range(len(obj)):
        obj[j] = e_sort[obj[j]] 
    f_e_obj_relation_array[i] = sorted(obj)

del e_sort

# print('   f_p',time.time()-start); sys.stdout.flush()
# temp = list(range(len(f_p_obj_relation_array)))
# f_sort = temp
# coupled = sorted(zip(f_p_obj_relation_array,f_sort,f_e_obj_relation_array))
# f_p_obj_relation_array, f_sort, f_e_obj_relation_array = list(zip(*coupled))
# coupled = sorted(zip(f_sort,temp))
# temp, f_sort = list(zip(*coupled))
# f_p_obj_relation_array = list(f_p_obj_relation_array)
# f_e_obj_relation_array= list(f_e_obj_relation_array)
# 
# for i in range(len(c_f_obj_relation_array)):
#     obj = c_f_obj_relation_array[i]
#     for j in range(len(obj)):
#         obj[j] = f_sort[obj[j]] 
#     c_f_obj_relation_array[i] = sorted(obj)
# 
# del temp
# del f_sort
# del coupled

# print('   e_p'); sys.stdout.flush()
# temp = list(range(len(e_p_obj_relation_array)))
# e_sort = temp
# coupled = sorted(zip(e_p_obj_relation_array,e_sort))
# e_p_obj_relation_array, e_sort = list(zip(*coupled))
# coupled = sorted(zip(e_sort,temp))
# temp, e_sort = list(zip(*coupled))
# e_p_obj_relation_array = list(e_p_obj_relation_array)
# 
# for i in range(len(f_e_obj_relation_array)):
#     obj = f_e_obj_relation_array[i]
#     for j in range(len(obj)):
#         obj[j] = e_sort[obj[j]] 
#     f_e_obj_relation_array[i] = sorted(obj)
# 
# del temp
# del e_sort
# del coupled

print('generating index arrays',time.time()-start); sys.stdout.flush()
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

print('removing duplicate connectivites',time.time()-start); sys.stdout.flush()

edge_to_rep=[]
face_to_rep=[]
cell_to_rep=[]

temp = []
temp2 = []
j=0
print('   c_p'); sys.stdout.flush()
for i in range(len(c_p_index)-1):
    cell_to_rep.append(j)
    if c_p_index[i] == c_p_index[i+1]:    
        obj1 = c_p_obj_relation_array[i]
        obj2 = c_p_obj_relation_array[i+1]
        for p in range(c_p_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(c_p_index[i])
                j=j+1
                break
                

temp.append(c_p_obj_relation_array[i+1])
temp2.append(c_p_index[i+1])
cell_to_rep.append(j+1)
c_p_obj_relation_array = temp
c_p_index = temp2
del temp
del temp2

temp = []
temp2 = []
j=0
print('   f_p'); sys.stdout.flush()
for i in range(len(f_p_index)-1):
    face_to_rep.append(j)
    if f_p_index[i] == f_p_index[i+1]:    
        obj1 = f_p_obj_relation_array[i]
        obj2 = f_p_obj_relation_array[i+1]
        for p in range(f_p_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(f_p_index[i])
                j=j+1
                break
            
temp.append(f_p_obj_relation_array[i+1])
temp2.append(f_p_index[i+1])
face_to_rep.append(j+1)
f_p_obj_relation_array = temp 
f_p_index = temp2
del temp
del temp2


temp = []
temp2 = []
j=0
print('   e_p'); sys.stdout.flush()
for i in range(len(e_p_index)-1):
    edge_to_rep.append(j)
    
    obj1 = e_p_obj_relation_array[i]
    obj2 = e_p_obj_relation_array[i+1]
    
    if obj1 !=obj2:
        temp.append(obj1)
        temp2.append(2)
        j=j+1
        
    
temp.append(e_p_obj_relation_array[i+1])
temp2.append(2)
edge_to_rep.append(j+1)
e_p_obj_relation_array = temp
e_p_index = temp2
del temp
del temp2
    
print('updating mappings',time.time()-start); sys.stdout.flush()

for i in range(len(f_e_obj_relation_array)):
    obj = f_e_obj_relation_array[i]
    for j in range(len(obj)):
        obj[j] = edge_to_rep[obj[j]] 
    f_e_obj_relation_array[i] = sorted(obj)
    
for i in range(len(c_f_obj_relation_array)):
    obj = c_f_obj_relation_array[i]
    for j in range(len(obj)):
        obj[j] = face_to_rep[obj[j]] 
    c_f_obj_relation_array[i] = sorted(obj)

print('removing duplicate mappings',time.time()-start); sys.stdout.flush()

temp = []
temp2 = []
print('   c_f'); sys.stdout.flush()
for i in range(len(c_f_index)-1):
    if c_f_index[i] == c_f_index[i+1]:    
        obj1 = c_f_obj_relation_array[i]
        obj2 = c_f_obj_relation_array[i+1]
        for p in range(c_f_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(c_f_index[i])
                break
            
temp.append(c_f_obj_relation_array[i+1])
temp2.append(c_f_index[i+1])
c_f_obj_relation_array = temp
c_f_index = temp2
del temp
del temp2

temp = []
temp2 = []
print('   f_e'); sys.stdout.flush()
for i in range(len(f_e_index)-1):
    if f_e_index[i] == f_e_index[i+1]:    
        obj1 = f_e_obj_relation_array[i]
        obj2 = f_e_obj_relation_array[i+1]
        for p in range(f_e_index[i]-1, 0, -1):
            if obj1[p]!=obj2[p]:
                temp.append(obj1)
                temp2.append(f_e_index[i])
                break
            
temp.append(f_e_obj_relation_array[i+1])
temp2.append(f_e_index[i+1])
f_e_obj_relation_array = temp
f_e_index = temp2
del temp
del temp2

gc.collect()

print('indexifying index arrays',time.time()-start); sys.stdout.flush()
# indexifying the index arrays
index_last = 0
for index in range(len(c_p_index)):
    c_p_index[index] = index_last + c_p_index[index]
    index_last = c_p_index[index]


index_last = 0
for index in range(len(f_p_index)):
    f_p_index[index] = index_last + f_p_index[index]
    index_last = f_p_index[index]
    

index_last = 0
for index in range(len(e_p_index)):
    e_p_index[index] = index_last + e_p_index[index]
    index_last = e_p_index[index]
    

index_last = 0
for index in range(len(c_f_index)):
    c_f_index[index] = index_last + c_f_index[index]
    index_last = c_f_index[index]
    

index_last = 0
for index in range(len(f_e_index)):
    f_e_index[index] = index_last + f_e_index[index]
    index_last = f_e_index[index]

print('updating counts',time.time()-start); sys.stdout.flush()
# update counts with non duplicate objects
nele  = len(c_p_index)
nface = len(f_p_index)
nedge = len(e_p_index)

c_p_sum=c_p_index[-1]
f_p_sum=f_p_index[-1]
e_p_sum=e_p_index[-1]
c_f_sum=c_f_index[-1]
f_e_sum=f_e_index[-1]

print()

print('beginning writing to file',time.time()-start); sys.stdout.flush()
# outputting


print('creating raw file',time.time()-start); sys.stdout.flush()
# opening/creating file
file_name = '../preprocessed_mesh_folder/raw_mesh_data.preprocessed_mesh_file'
file = open(file_name, "wb")


print('writing header',time.time()-start); sys.stdout.flush()
# writing header
file.write(struct.pack('<4i' ,npoin,nedge,nface,nele))
print(npoin,nedge,nface,nele)


print('writing coordinates',time.time()-start); sys.stdout.flush()
# writing coordinate data
for coord in x_coords:
    file.write(struct.pack('<f' ,coord))
for coord in y_coords:
    file.write(struct.pack('<f' ,coord))
for coord in z_coords:
    file.write(struct.pack('<f' ,coord))
    
    
print('writing connectivities',time.time()-start); sys.stdout.flush()
# writing object connectivities
file.write(struct.pack('<2i' ,c_p_sum, 0))
for entry in c_p_index:
    file.write(struct.pack('<i' ,entry))
    
file.write(struct.pack('<2i' ,f_p_sum, 0))
for entry in f_p_index:
    file.write(struct.pack('<i' ,entry))


for obj in c_p_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry + 1))

for obj in f_p_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry + 1))
    
for obj in e_p_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry + 1))
    
print('writing cell > face, face > edge mappings',time.time()-start); sys.stdout.flush()
# writing object relation mappings
file.write(struct.pack('<2i' ,c_f_sum, 0))
for entry in c_f_index:
    file.write(struct.pack('<i' ,entry))
    
file.write(struct.pack('<2i' ,f_e_sum, 0))
for entry in f_e_index:
    file.write(struct.pack('<i' ,entry))
    
for obj in c_f_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry + 1))

for obj in f_e_obj_relation_array:
    for entry in obj:
        file.write(struct.pack('<i' ,entry + 1))
    
print('finished writing to file',time.time()-start); sys.stdout.flush()
print('finished preprocessing'); sys.stdout.flush()


