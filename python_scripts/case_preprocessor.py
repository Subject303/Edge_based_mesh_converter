from vtk import *
import pyvista
import os
import struct
import numpy as np

from vtkmodules.vtkIOEnSight import *
from vtkmodules.util.misc import *
from vtkmodules.vtkFiltersPoints import *

from vtkmodules.all import *

VTK_DATA_ROOT = vtkGetDataRoot()

# set angle to split boundary edges
feature_angle = 15

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
fluidBlock = raw_data.GetBlock(0)  # vtkUnstructuredGrid format


del raw_data
# freeing up the multiblock view


# general cleanup
algo = vtkStaticCleanUnstructuredGrid()
algo.SetInputData(fluidBlock)
algo.RemoveUnusedPointsOn()
algo.Update()
fluidBlock = algo.GetOutput()
del algo


# convert all cells to polyhedra
algo = vtkConvertToPolyhedra()
algo.SetInputData(fluidBlock)
algo.OutputAllCellsOn()
algo.Update()
polyblock = algo.GetOutput()
del algo
del fluidBlock 
# I think I can get away deleting the non poly data


# this splits the points on feature edges
# we do this before assigning global IDs so the IDs are seperate
algo = vtkSplitSharpEdgesPolyData()
algo.SetInputData(polyblock)
algo.SetFeatureAngle(feature_angle)
algo.Update()
polyblock = algo.GetOutput()
del algo


# create global point IDs that I can reference down the line
algo = vtkGenerateGlobalIds()
algo.SetInputData(polyblock)
algo.Update()
polyblock = algo.GetOutput()
del algo


# Pull general mesh info
npoin = polyblock.GetNumberOfPoints()
nelem = polyblock.GetNumberOfCells()
print(npoin,' nodes ', nelem,' cells')


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

cells = []
faces = []
edges = []


# extract coordinate data

for p in range(npoin):
    x_coords.append(polyblock.GetPoint(p)[0])
    y_coords.append(polyblock.GetPoint(p)[1])
    z_coords.append(polyblock.GetPoint(p)[2])
    # node: these will include the duplicated points from the boundary data

# extract cell, face and edge datas from the primary data block

for c in range(nelem):
    cells.append(polyblock.GetCell(c))

del polyblock

faceid = 0

for cell in cells:
    nfaces_in_cell = cell.GetNumberOfFaces()
    
    c_p_index.append(cell.GetNumberOfPoints())
    c_f_index.append(nfaces_in_cell)
    
    for f in range(nfaces_in_cell): 
        faceid = faceid + 1
        c_f_obj_relation_array.append(faceid)
        faces.append(cell.GetFace(f))
    
edgeid = 0
    
for face in faces:
    nedge_in_face = face.GetNumberOfPoints()
    
    f_p_index.append(cell.GetNumberOfPoints())
    f_e_index.append(nedge_in_face)
    
    for e in range(nedge_in_face):
        edgeid = edgeid + 1
        f_e_obj_relation_array.append(edgeid)
        edges.append(face.GetEdge(e))


# extract pure connectivity arrays

cellid = 0
for cell in cells:
    cellid = cellid + 1
    
    for p in range(c_p_index(cellid)):
        c_p_obj_relation_array.append(cell.GetPointId(p))

faceid = 0
for face in faces:
    faceid = faceid + 1
    
    for p in range(f_p_index(faceid)):
        f_p_obj_relation_array.append(face.GetPointId(p))
        
nface = faceid

edgeid = 0
for edge in edges:
    edgeid = edgeid + 1
    
    e_p_obj_relation_array.append(edge.GetPointId(1))
    e_p_obj_relation_array.append(edge.GetPointId(2))
        
nedge = edgeid




# need to get the obj relation sums

c_p_sum = len(c_p_obj_relation_array)
f_p_sum = len(f_p_obj_relation_array)
e_p_sum = len(c_p_obj_relation_array)

c_f_sum = len(c_f_obj_relation_array)
f_e_sum = len(f_e_obj_relation_array)


# ok I should Have all data I need to output now.
# and indexify the index arrays
# and sort the relation arrays and delete duplicated stuff
# I just need to filter out boundary duplicated nodes


# indexifying the index arrays
index_last = 0
for index in c_p_index:
    index = index_last + index
    index_last = index

index_last = 0
for index in f_p_index:
    index = index_last + index
    index_last = index

e_p_index = 2 * (range(len(edges)) + 1)
# can be implicitly known

index_last = 0
for index in c_f_index:
    index = index_last + index
    index_last = index

index_last = 0
for index in f_e_index:
    index = index_last + index
    index_last = index
    
    
# sorting relation arrays
index_end = 0
for index in c_p_index:
    index_start = index_end + 1
    index_end = index
    
    c_p_obj_relation_array[index_start:index_end] = sorted(c_p_obj_relation_array[index_start:index_end])
    
index_end = 0
for index in f_p_index:
    index_start = index_end + 1
    index_end = index
    
    f_p_obj_relation_array[index_start:index_end] = sorted(f_p_obj_relation_array[index_start:index_end])
    
index_end = 0
for index in e_p_index:
    index_start = index_end + 1
    index_end = index
    
    e_p_obj_relation_array[index_start:index_end] = sorted(e_p_obj_relation_array[index_start:index_end])
    
index_end = 0
for index in c_f_index:
    index_start = index_end + 1
    index_end = index
    
    c_f_obj_relation_array[index_start:index_end] = sorted(c_f_obj_relation_array[index_start:index_end])
    
index_end = 0
for index in f_e_index:
    index_start = index_end + 1
    index_end = index
    
    f_e_obj_relation_array[index_start:index_end] = sorted(f_e_obj_relation_array[index_start:index_end])


# killing duplicate relations 
nele = len(c_p_index)
index_1_start = c_p_index[nele-1]
index_1_end   = c_p_index[nele] - 1
for i in range((nele-1),0,-1): # going backwards here means we can kill duplicates in one loop
    index_2_start = index_1_start
    index_2_end   = index_1_end
    index_1_start = c_p_index[i-1] + 1
    index_1_end   = c_p_index[i]
    if c_p_obj_relation_array[index_1_start:index_1_end]==c_p_obj_relation_array[index_2_start:index_2_end]: 
        del c_p_obj_relation_array[index_2_start:index_2_end]
        del c_p_index[i]

nface = len(f_p_index)
index_1_start = f_p_index[nface-1]
index_1_end   = f_p_index[nface] - 1
for i in range((nface-1),0,-1): # going backwards here means we can kill duplicates in one loop
    index_2_start = index_1_start
    index_2_end   = index_1_end
    index_1_start = f_p_index[i-1] + 1
    index_1_end   = f_p_index[i]
    if f_p_obj_relation_array[index_1_start:index_1_end]==f_p_obj_relation_array[index_2_start:index_2_end]: 
        del f_p_obj_relation_array[index_2_start:index_2_end]
        del f_p_index[i]

nedge = len(e_p_index)
index_1_start = e_p_index[nedge-1]
index_1_end   = e_p_index[nedge] - 1
for i in range((nedge-1),0,-1): # going backwards here means we can kill duplicates in one loop
    index_2_start = index_1_start
    index_2_end   = index_1_end
    index_1_start = e_p_index[i-1] + 1
    index_1_end   = e_p_index[i]
    if e_p_obj_relation_array[index_1_start:index_1_end]==e_p_obj_relation_array[index_2_start:index_2_end]: 
        del e_p_obj_relation_array[index_2_start:index_2_end]
        del e_p_index[i]

nele = len(c_f_index)
index_1_start = c_f_index[nele-1]
index_1_end   = c_f_index[nele] - 1
for i in range((nele-1),0,-1): # going backwards here means we can kill duplicates in one loop
    index_2_start = index_1_start
    index_2_end   = index_1_end
    index_1_start = c_f_index[i-1] + 1
    index_1_end   = c_f_index[i]
    if c_f_obj_relation_array[index_1_start:index_1_end]==c_f_obj_relation_array[index_2_start:index_2_end]: 
        del c_f_obj_relation_array[index_2_start:index_2_end]
        del c_f_index[i]

nface = len(f_e_index)
index_1_start = f_e_index[nface-1]
index_1_end   = f_e_index[nface] - 1
for i in range((nface-1),0,-1): # going backwards here means we can kill duplicates in one loop
    index_2_start = index_1_start
    index_2_end   = index_1_end
    index_1_start = f_e_index[i-1] + 1
    index_1_end   = f_e_index[i]
    if f_e_obj_relation_array[index_1_start:index_1_end]==f_e_obj_relation_array[index_2_start:index_2_end]: 
        del f_e_obj_relation_array[index_2_start:index_2_end]
        del f_e_index[i]


# update counts with non duplicate objects
nele  = len(c_p_index)
nface = len(f_p_index)
nedge = len(e_p_index)

c_p_sum = len(c_p_obj_relation_array)
f_p_sum = len(f_p_obj_relation_array)
e_p_sum = len(e_p_obj_relation_array)

c_f_sum = len(c_f_obj_relation_array)
f_e_sum = len(f_e_obj_relation_array)





# outputting


# opening/creating file
file_name = 'Raw_mesh_data'
file = open(file_name, "wb")


# writing header
file.write(struct.pack('<4d' ,npoin,nedge,nface,nelem))


# writing coordinate data
for coord in x_coords:
    file.write(struct.pack('<d' ,coord))
for coord in y_coords:
    file.write(struct.pack('<d' ,coord))
for coord in z_coords:
    file.write(struct.pack('<d' ,coord))


# writing object connectivities
file.write(struct.pack('<2d' ,c_p_sum, 0.0))
for entry in c_p_index:
    file.write(struct.pack('<d' ,entry))
    
file.write(struct.pack('<2d' ,f_p_sum, 0.0))
for entry in f_p_index:
    file.write(struct.pack('<d' ,entry))


for entry in c_p_obj_relation_array:
    file.write(struct.pack('<d' ,entry))

for entry in f_p_obj_relation_array:
    file.write(struct.pack('<d' ,entry))
    
for entry in e_p_obj_relation_array:
    file.write(struct.pack('<d' ,entry))
    
    
# writing object relation mappings
file.write(struct.pack('<2d' ,c_f_sum, 0.0))
for entry in c_f_index:
    file.write(struct.pack('<d' ,entry))
    
file.write(struct.pack('<2d' ,f_e_sum, 0.0))
for entry in f_e_index:
    file.write(struct.pack('<d' ,entry))
    
    
for entry in c_f_obj_relation_array:
    file.write(struct.pack('<d' ,entry))

for entry in f_e_obj_relation_array:
    file.write(struct.pack('<d' ,entry))
    
    
    
    
    
    
    
    
    
