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


# read the raw data
algo = vtkGenericEnSightReader()
cdp = vtkCompositeDataPipeline()
# Make sure all algorithms use the composite data pipeline
# might not be needed
algo.SetDefaultExecutivePrototype(cdp)
del cdp
algo.SetCaseFileName("/mnt/gpfs01/home/ws/wsjt10/python/case/star.case") 
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
algo.SetFeatureAngle(15)
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
#e_p_index = [] # implicit
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


#











