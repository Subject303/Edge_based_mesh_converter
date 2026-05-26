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

algo = vtkGenericEnSightReader()
# Make sure all algorithms use the composite data pipeline
# might not be needed

cdp = vtkCompositeDataPipeline()
algo.SetDefaultExecutivePrototype(cdp)
del cdp
algo.SetCaseFileName("/mnt/gpfs01/home/ws/wsjt10/python/case/star.case") 
# will want to change to somthing more generic
algo.Update()
raw_data = algo.GetOutput()
del algo
fluidBlock = raw_data.GetBlock(0)  # vtkUnstructuredGrid format
del raw_data

# freeing up the multiblock view

algo = vtkStaticCleanUnstructuredGrid()
algo.SetInputData(fluidBlock)
algo.RemoveUnusedPointsOn()
algo.Update()
fluidBlock = algo.GetOutput()
del algo


algo = vtkGenerateGlobalIds()
algo.SetInputData(fluidBlock)
algo.Update()
fluidBlock = algo.GetOutput()
del algo


# Overall mesh information
npoin = fluidBlock.GetNumberOfPoints()
nelem = fluidBlock.GetNumberOfCells()
print(npoin,' nodes ', nelem,' elements')
