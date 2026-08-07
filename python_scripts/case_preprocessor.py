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
# algo.SetCaseFileName("./case/tet_sphere.case") #tet sphere
# algo.SetCaseFileName("./case/666cube.case") #666 cube
# algo.SetCaseFileName("./case/2.case") #tiny cube
# algo.SetCaseFileName("./case/tet_cube.case") #TET cube
# algo.SetCaseFileName("./case/2_rad_sphere_poly.case") #4 diameter poly sphere
algo.SetCaseFileName("./case/star1.case") # 0.75 mil mesh
# algo.SetCaseFileName("./case/star.case") # 24 mil mesh
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
x_coords = [0]*npoin
y_coords = [0]*npoin
z_coords = [0]*npoin

c_p_obj_relation_array = [[]]*nele
f_p_obj_relation_array = []
e_p_obj_relation_array = []

c_f_obj_relation_array = [[]]*nele
f_e_obj_relation_array = []


print('extracting coordinate data',time.time()-start); sys.stdout.flush()
# extract coordinate data

for p in range(npoin):
    # x_coords.append(struct.pack('<f' ,polyblock.GetPoint(p)[0]))
    # y_coords.append(struct.pack('<f' ,polyblock.GetPoint(p)[1]))
    # z_coords.append(struct.pack('<f' ,polyblock.GetPoint(p)[2]))
    
    x_coords[p] = struct.pack('<d' ,polyblock.GetPoint(p)[0])
    y_coords[p] = struct.pack('<d' ,polyblock.GetPoint(p)[1])
    z_coords[p] = struct.pack('<d' ,polyblock.GetPoint(p)[2])
    
    # node: these will include the duplicated points from the boundary data

# extract cell, face and edge datas from the primary data block
faceid = 0
edgeid = 0
print('extracting connectivity and mapping data',time.time()-start); sys.stdout.flush()
for c in range(nele):
    
    cell = polyblock.GetCell(c)
    
    #c_p_obj_relation_array [c] = sorted([cell.GetPointId(p) for p in range(cell.GetNumberOfPoints())])
    c_p_obj_relation_array [c] = [cell.GetPointId(p) for p in range(cell.GetNumberOfPoints())]
    
    
    nfaces = cell.GetNumberOfFaces()
    
    f_array = []
    
    for f in range(cell.GetNumberOfFaces()): 
        
        face = cell.GetFace(f)
        
        # f_p_obj_relation_array.append(sorted([face.GetPointId(p) for p in range(face.GetNumberOfPoints())]))
        f_p_obj_relation_array.append([face.GetPointId(p) for p in range(face.GetNumberOfPoints())])
        
        nedges = face.GetNumberOfEdges()
        
        e_array = []
        
        for e in range(nedges):
            
            edge = face.GetEdge(e)
            
            e_p_obj_relation_array.append([edge.GetPointId(0),edge.GetPointId(1)])
            
            e_array.append(edgeid)
            edgeid = edgeid + 1
    
        f_array.append(faceid)
        faceid = faceid + 1
        
        f_e_obj_relation_array.append(e_array)
        
    c_f_obj_relation_array[c] = f_array
    
del polyblock
del cell
del face
del edge
    
nface = faceid
nedge = edgeid

gc.collect()

print('sorting connectivities',time.time()-start); sys.stdout.flush()

print('   c_p',time.time()-start); sys.stdout.flush()

sorted_c_p    = [sorted(obj) for obj in c_p_obj_relation_array]

coupled = zip(sorted_c_p, c_p_obj_relation_array, c_f_obj_relation_array)

coupled = sorted(coupled, key=lambda element: element[0])

i=0
for p, y, z in coupled:
    sorted_c_p[i]             = p
    c_p_obj_relation_array[i] = y
    c_f_obj_relation_array[i] = z
    i=i+1

del coupled; gc.collect()

print('   f_p',time.time()-start); sys.stdout.flush()

f_sort = list(range(len(f_p_obj_relation_array)))

sorted_f_p    = [sorted(obj) for obj in f_p_obj_relation_array]

coupled = sorted(zip(sorted_f_p, f_p_obj_relation_array, f_sort, f_e_obj_relation_array), key=lambda element: element[0])

i=0
for p, x, y, z in coupled:
    sorted_f_p[i]             = p
    f_p_obj_relation_array[i] = x
    f_sort[i]                 = y
    f_e_obj_relation_array[i] = z
    i=i+1



coupled = sorted(zip(f_sort, list(range(len(f_p_obj_relation_array)))), key=lambda element: element[0])

i=0
for x, y in coupled:
    f_sort[i] = y
    i=i+1

del coupled; gc.collect()

for i in range(len(c_f_obj_relation_array)):
    obj = c_f_obj_relation_array[i]
    for j in range(len(obj)):
        obj[j] = f_sort[obj[j]]
    c_f_obj_relation_array[i] = sorted(obj)

del f_sort

print('   e_p',time.time()-start); sys.stdout.flush()

e_sort = list(range(len(e_p_obj_relation_array)))

sorted_e_p    = [sorted(obj) for obj in e_p_obj_relation_array]

coupled = sorted(zip(sorted_e_p, e_p_obj_relation_array, e_sort), key=lambda element: element[0])

i=0
for p, x, y in coupled:
    sorted_e_p[i] = p
    e_p_obj_relation_array[i] = x
    e_sort[i] = y
    i=i+1
    
    
coupled = sorted(zip(e_sort, list(range(len(e_p_obj_relation_array)))), key=lambda element: element[0])

i=0
for x, y in coupled:
    e_sort[i] = y
    i=i+1
    
del coupled; gc.collect()

for i in range(len(f_e_obj_relation_array)):
    obj = f_e_obj_relation_array[i]
    for j in range(len(obj)):
        obj[j] = e_sort[obj[j]] 
    f_e_obj_relation_array[i] = sorted(obj)

del e_sort

print('generating index arrays',time.time()-start); sys.stdout.flush()

nele  = len(c_p_obj_relation_array)
nface = len(f_p_obj_relation_array)
nedge = len(e_p_obj_relation_array)

c_p_sum = 0
f_p_sum = 0
e_p_sum = 0

c_p_index = [0]*nele
f_p_index = [0]*nface
e_p_index = [0]*nedge

c_f_sum = 0
f_e_sum = 0

c_f_index = [0]*nele
f_e_index = [0]*nface

i=0
for obj in c_p_obj_relation_array:
    c_p_index[i] = len(obj)
    i=i+1
    
i=0
for obj in f_p_obj_relation_array:
    f_p_index[i] = len(obj)
    i=i+1
    
i=0
for obj in e_p_obj_relation_array:
    e_p_index[i] = len(obj)
    i=i+1
    
i=0
for obj in c_f_obj_relation_array:
    c_f_index[i] = len(obj)
    i=i+1
    
i=0
for obj in f_e_obj_relation_array:
    f_e_index[i] = len(obj)
    i=i+1
    
# ok I should Have all data I need to output now.
# and indexify the index arrays
# and sort the relation arrays and delete duplicated stuff
# I just need to filter out boundary duplicated nodes


print(len(c_p_index),len(c_f_index),len(c_f_obj_relation_array))
print(len(f_p_index),len(f_e_index),len(f_e_obj_relation_array))

print('removing duplicate connectivites',time.time()-start); sys.stdout.flush()

c_uniqe     = [False]*nele
cell_to_rep = [-1]*nele

k=0
print('   c arrays cp cf',time.time()-start); sys.stdout.flush()
for i in range(nele-1):
    
    cell_to_rep[i] = k
    
    if c_p_index[i] == c_p_index[i+1]:
        
        obj1 = sorted_c_p[i]
        obj2 = sorted_c_p[i+1]
        
        for p in range(c_p_index[i]):
            
            if obj1[p] != obj2[p]: 
                c_uniqe[i] = True
                k=k+1
                break
    else:
        c_uniqe[i] = True
        k=k+1

c_uniqe[-1] = True
cell_to_rep[-1] = k
temp  = []
temp2 = []
temp3 = []
temp4 = []

for i in range(nele):
    if c_uniqe[i] == True:
        temp.append( c_p_obj_relation_array[i])
        temp2.append(c_p_index[i])
        temp3.append(c_f_obj_relation_array[i])
        temp4.append(c_f_index[i])

c_p_obj_relation_array = temp 
c_p_index = temp2
c_f_obj_relation_array = temp3
c_f_index = temp4
del temp
del temp2
del temp3
del temp4; gc.collect()



f_uniqe     = [False]*nface
face_to_rep = [-1]*nface

k=0
print('   f arrays fp fe',time.time()-start); sys.stdout.flush()
for i in range(nface-1):
    
    face_to_rep[i] = k
    
    if f_p_index[i] == f_p_index[i+1]:
        
        obj1 = sorted_f_p[i]
        obj2 = sorted_f_p[i+1]
        
        for p in range(f_p_index[i]):
            
            if obj1[p] != obj2[p]: 
                f_uniqe[i] = True
                k=k+1
                break
    else:
        f_uniqe[i] = True
        k=k+1

f_uniqe[-1] = True
face_to_rep[-1] = k

temp  = []
temp2 = []
temp3 = []
temp4 = []

for i in range(nface):
    if f_uniqe[i] == True:
        temp.append( f_p_obj_relation_array[i])
        temp2.append(f_p_index[i])
        temp3.append(f_e_obj_relation_array[i])
        temp4.append(f_e_index[i])

f_p_obj_relation_array = temp 
f_p_index = temp2
f_e_obj_relation_array = temp3
f_e_index = temp4
del temp
del temp2
del temp3
del temp4; gc.collect()

e_uniqe     = [False]*nedge
edge_to_rep = [-1]*nedge

k=0
print('   e arrays ep',time.time()-start); sys.stdout.flush()
for i in range(nedge-1):
    
    edge_to_rep[i] = k
    
    obj1 = sorted_e_p[i]
    obj2 = sorted_e_p[i+1]
    
    if obj1 != obj2: 
        e_uniqe[i] = True
        k=k+1

e_uniqe[-1] = True
edge_to_rep[-1] = k
temp  = []
temp2 = []

for i in range(nedge):
    if e_uniqe[i] == True:
        temp.append( e_p_obj_relation_array[i])
        temp2.append(e_p_index[i])

e_p_obj_relation_array = temp 
e_p_index = temp2
del temp
del temp2; gc.collect()


print(len(c_p_index),len(c_f_index),len(c_f_obj_relation_array))
print(len(f_p_index),len(f_e_index),len(f_e_obj_relation_array))

print('   updating c f mapping',time.time()-start); sys.stdout.flush()

for i in range(len(c_f_obj_relation_array)):
    obj = c_f_obj_relation_array[i]
    for j in range(len(obj)):
        obj[j] = face_to_rep[obj[j]] 
    c_f_obj_relation_array[i] = sorted(obj)

print('   updating f e mapping',time.time()-start); sys.stdout.flush()
for i in range(len(f_e_obj_relation_array)):
    obj = f_e_obj_relation_array[i]
    for j in range(len(obj)):
        obj[j] = edge_to_rep[obj[j]] 
    f_e_obj_relation_array[i] = sorted(obj)

# temp = []
# temp2 = []
# print('   c_f'); sys.stdout.flush()
# for i in range(len(c_f_index)-1):
#     if c_f_index[i] == c_f_index[i+1]:    
#         obj1 = c_f_obj_relation_array[i]
#         obj2 = c_f_obj_relation_array[i+1]
# 
#         for p in range(c_f_index[i]):
#             if obj1[p]!=obj2[p]:
#                 temp.append(obj1)
#                 temp2.append(c_f_index[i])
#                 break
# 
# 
# temp.append(c_f_obj_relation_array[i+1])
# temp2.append(c_f_index[i+1])
# c_f_obj_relation_array = temp
# c_f_index = temp2
# del temp
# del temp2; gc.collect()
# 
# temp = []
# temp2 = []
# print('   f_e'); sys.stdout.flush()
# for i in range(len(f_e_index)-1):
#     if f_e_index[i] == f_e_index[i+1]:    
#         obj1 = f_e_obj_relation_array[i]
#         obj2 = f_e_obj_relation_array[i+1]
# 
#         for p in range(f_e_index[i]):
#             if obj1[p]!=obj2[p]:
#                 temp.append(obj1)
#                 temp2.append(f_e_index[i])
#                 break
#             
# temp.append(f_e_obj_relation_array[i+1])
# temp2.append(f_e_index[i+1])
# f_e_obj_relation_array = temp
# f_e_index = temp2
# del temp
# del temp2; gc.collect()
# 
# gc.collect()

print(len(c_p_index),len(c_f_index),len(c_f_obj_relation_array))
print(len(f_p_index),len(f_e_index),len(f_e_obj_relation_array))

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
file_name = './preprocessed_mesh_folder/raw_mesh_data.preprocessed_mesh_file'
file = open(file_name, "wb")


print('writing header',time.time()-start); sys.stdout.flush()
# writing header
file.write(struct.pack('<4i' ,npoin,nedge,nface,nele))
print('mesh contains : ')
print('points, edges, faces, elements ')
print(npoin,nedge,nface,nele)        
print( 'Relation array sums: ')
print( 'c_p_sum,  f_p_sum,  e_p_sum,  c_f_sum,  f_e_sum')
print(c_p_sum,f_p_sum,e_p_sum,c_f_sum,f_e_sum)
print(len(c_f_index),len(f_e_index))

print(len(c_p_index),len(c_f_index),len(c_f_obj_relation_array))
print(len(f_p_index),len(f_e_index),len(f_e_obj_relation_array))

print('writing coordinates',time.time()-start); sys.stdout.flush()
# writing coordinate data
for coord in x_coords:
    file.write(coord)
for coord in y_coords:
    file.write(coord)
for coord in z_coords:
    file.write(coord)
    
    
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


