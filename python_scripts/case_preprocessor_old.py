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


reader = vtkGenericEnSightReader()
# Make sure all algorithms use the composite data pipeline
cdp = vtkCompositeDataPipeline()
reader.SetDefaultExecutivePrototype(cdp)
del cdp
reader.SetCaseFileName("/mnt/gpfs01/home/ws/wsjt10/python/case/star.case")
reader.Update()
its_data = reader.GetOutput()

print(its_data.GetNumberOfBlocks (),' blocks in case')

print('Assuming Block = 0 is the fluid block')
fluidBlock = its_data.GetBlock(0)  # vtkUnstructuredGrid format

outer = vtkStaticCleanUnstructuredGrid()
outer.SetInputData(fluidBlock)
outer.RemoveUnusedPointsOn()
outer.Update()
fluidBlock = outer.GetOutput()

globalids = vtkGenerateGlobalIds()
globalids.SetInputData(fluidBlock)
globalids.Update()
fluidBlock = globalids.GetOutput()

nBlocks = its_data.GetNumberOfBlocks()
# Overall mesh information
npoin = fluidBlock.GetNumberOfPoints()
nelem = fluidBlock.GetNumberOfCells()
print(npoin,' nodes ', nelem,' elements')

#fluidBlock.InitializeFacesRepresentation(nelem)

###################################################################################################
            #finding all unique edges faces ect
###################################################################################################

convert = vtkConvertToPolyhedra()
convert.SetInputData(fluidBlock)
convert.OutputAllCellsOn()
convert.Update()
polyblock = convert.GetOutput()
listoffcons = []
listofecons = []

point_data = polyblock.GetPointData()

globalpointids = point_data.GetArray('GlobalPointIds')

ielem = []
ipoly = []
for i in range(nelem):
    ele = polyblock.GetCell(i)
    nfaces = ele.GetNumberOfFaces()
    elecon = []
    polyface = []
    for ip in range(ele.GetNumberOfPoints()):
        elecon.append(ele.GetPointId(ip))
        
    ielem.append(elecon)
    
    polyface.append(nfaces)
    
    for j in range(nfaces):
        unique = 0
        facecon =[]
        polycon =[]
        face = ele.GetFace(j)
        faceconvtk = face.GetPointIds()
        nfpoints = face.GetNumberOfPoints()
        polycon.append(nfpoints)
        
        for k in range(nfpoints):
            facecon.append(faceconvtk.GetId(k))
            polycon.append(faceconvtk.GetId(k)) # if I dont do this polycon magically becomes sorted
        
        polyface.append(polycon)

        
        facecon.sort()

        
        listoffcons.append(facecon)

        nfedge = face.GetNumberOfEdges()
        
        for q in range(nfedge):
            edgecon = []
            edge = face.GetEdge(q)
            edgecon.append(edge.GetPointId(0))
            edgecon.append(edge.GetPointId(1))
            
            edgecon.sort()
            listofecons.append(edgecon)

    ipoly.append(polyface)


listoffcons.sort()
listofecons.sort()
print(len(listofecons),' non - unique edges ',len(listoffcons), ' non - unique faces')

iea = []
ifa = []
if2 = []
ie2 = []
iface = []
iedge = []

iedge.append(listofecons[0])
ie2 = listofecons[0]
for iea in listofecons:
    for i in range(len(iea)):
        if iea[i] != ie2[i]:
            iedge.append(iea)
            break
    ie2 = iea

iface.append(listoffcons[0])
if2 = listoffcons[0]
for ifa in listoffcons:
    for i in range(len(ifa)):
        if ifa[i] != if2[i]:
            iface.append(ifa)
            break
    if2 = ifa

print(len(iedge),' unique edges ',len(iface), ' unique faces')

nedge = len(iedge)
nface = len(iface)

###################################################################################################
            #identifying boundaries
###################################################################################################

print('identifing boundary faces and edges')

angle = 15

bleh = vtkPolyData()
outer = vtkGeometryFilter()
outer.SetInputData(fluidBlock)
outer.Update()
bleh = outer.GetOutput()

outer = vtkCleanPolyData ()
outer.SetInputData(bleh)
outer.PointMergingOn ()
outer.Update()
bleh = outer.GetOutput()

outer = vtkRemoveDuplicatePolys()
outer.SetInputData(bleh)
outer.Update()
bleh = outer.GetOutput()

featedges = vtkPolyDataNormals()
featedges.SetInputData(bleh)
featedges.SetFeatureAngle(angle)
featedges.ComputePointNormalsOn ()
featedges.ComputeCellNormalsOn ()
featedges.Update()
bleh1 = featedges.GetOutput()

# I lose global node IDs here
# but not global cell ids
# so I can scan both the split and non split datasets,
# and relate local point ids in the split dataset to non split global node IDs

##regions = vtkConnectedPointsFilter()
##regions.SetInputData(bleh1)
##regions.InitializeSeedList()
##regions.SetExtractionModeToAllRegions()
##regions.AlignedNormalsOn ()
##regions.SetNormalAngle(angle)
##regions.ScalarConnectivityOn()
##regions.SetRadius (0.1)
##regions.Update()

regions = vtkPolyDataEdgeConnectivityFilter  ()
regions.SetInputData(bleh1)
regions.InitializeSeedList()
regions.CellRegionAreasOn ()
regions.SetExtractionModeToAllRegions()
regions.Update()

bleh2 = regions.GetOutput()
bleh3 = regions.GetNumberOfExtractedRegions ()

point_data = bleh.GetPointData()
point_data3 = bleh1.GetPointData()
point_data2 = bleh2.GetPointData()
Cell_data = bleh2.GetCellData()

CellNormals = Cell_data.GetArray('Normals')

globalpointids = point_data.GetArray('GlobalPointIds')
Normals = point_data2.GetArray('Normals')
Normals1 = point_data3.GetArray('Normals')
#RegionLabels = point_data2.GetArray('RegionLabels')
RegionLabels = Cell_data.GetArray('RegionId')

print(angle, ' degree angles used to identify ',bleh3, ' boundaries')

listFbounds = []
listNUEbounds = []
ibpoinandNormals = []
listSBBbounds = []

for i in range(bleh1.GetNumberOfCells()):
    boundcon = []
    boundcon2 = []
    ele = bleh1.GetCell(i) # I will get the point regions off this element for boundaries
    #cellid = globalcellids.GetValue(i) # gets the global cell ID off split dataset
    ele1 = bleh.GetCell(i)     # the cells are identical so I dont need a global cell id list

    for ip in range(ele1.GetNumberOfPoints()):
        bpoin = []
        boundcon.append(globalpointids.GetValue(ele1.GetPointId(ip))) # point id reletive to the original polydata
        
        #boundcon2.append(ele.GetPointId(ip))
        #boundcon2.append(RegionLabels.GetValue(ele.GetPointId(ip)))

        boundcon2.append(RegionLabels.GetValue(i))
        
        bpoin.append(globalpointids.GetValue(ele1.GetPointId(ip)))
        bpoin.append(RegionLabels.GetValue(i))
        #bpoin.append(i)
        bpoin.append(Normals1.GetComponent(ele.GetPointId(ip),0))
        bpoin.append(Normals1.GetComponent(ele.GetPointId(ip),1))
        bpoin.append(Normals1.GetComponent(ele.GetPointId(ip),2))
    
        ibpoinandNormals.append(bpoin)

    boundcon.sort()
    listFbounds.append(boundcon)
    listSBBbounds.append(boundcon2)

    nfedge = ele1.GetNumberOfEdges()
    for q in range(nfedge):
        edgecon = []
        edge = ele1.GetEdge(q)
        
        edgecon.append(globalpointids.GetValue(edge.GetPointId(0)))
        edgecon.append(globalpointids.GetValue(edge.GetPointId(1)))
        
        edgecon.sort()
        listNUEbounds.append(edgecon)

listNUEbounds.sort() # not needed but i want to be safe

ibpoinandNormals.sort()
print('Establishing the IBPOIN and ANP arrays')

listbpoin = []
listbpoin.append(ibpoinandNormals[0])
ie2 = ibpoinandNormals[0]  # set the first bound in the unique edge array
for iea in ibpoinandNormals: # for every edge in the non unique array
    for i in range(len(iea)): # for every point in the edge
        if iea[i] != ie2[i]: # if a point is not a duplicate with the one prior
            listbpoin.append(iea) # add it to the unique list
            break
    ie2 = iea

print(len(listbpoin),' boundary points (including duplicate corners nodes)')
listUEbounds = listNUEbounds
listUEbounds = []
listUEbounds.append(listNUEbounds[0])
ie2 = listNUEbounds[0]  # set the first bound in the unique edge array
for iea in listNUEbounds: # for every edge in the non unique array
    for i in range(len(iea)): # for every point in the edge
        if iea[i] != ie2[i]: # if a point is not a duplicate with the one prior
            listUEbounds.append(iea) # add it to the unique list
            break
    ie2 = iea


print('Removing Duplicate Boundary edges')
    
n4 = len(listUEbounds)

iedge2 = []
iface2 = []
n = 0
for i in range(len(iedge)): # for every edge

    ie = iedge[i]
    
    j = listUEbounds[n]
    if ie != j:
        iedge2.append(ie)
        continue
    n = n + 1
    if n == n4:
        if i < len(iedge):
            iedge2 = iedge2 + iedge[i+1:]
        break

print('edges complete moving to Boundary faces')

n4 = len(listFbounds)
n = 0
for i in range(len(iface)):

    ie = iface[i]
    j = listFbounds[n]
    if ie != j:
        iface2.append(ie)
        continue
    n = n + 1
    if n == n4:
        if i < len(iface):
            iface2 = iface2 + iface[i+1:]
        break

iedge = iedge2
iface = iface2
  
#for i in listFbounds: # it doesnt matter that half of listFbounds is region labels in euclidean geometry
#    iface.remove(i)

print(len(listUEbounds),' unique boundary edges ',len(listFbounds), ' unique boundary faces')
print(len(iedge),' unique internal edges ',len(iface), ' unique internal faces')

nedge = len(iedge)
nface = len(iface)
nbedge = len(listUEbounds)
nbface = len(listFbounds)



###################################################################################################
            #outputting
###################################################################################################

print('writing IO')

# Output file
file = 'FiniteElementData.txt'
thefile = open(file, "wb")
nbpoin = len(listbpoin)

thefile.write(struct.pack('<7id' ,npoin,nedge,nbedge,nface,nbface,nbpoin,nelem,angle))

#bin_array = struct.pack('<7if' ,npoin,nedge,nbedge,nface,nbface,nbpoin,nelem,angle)


for i in range(npoin):
    x = fluidBlock.GetPoint(i)[0]
    y = fluidBlock.GetPoint(i)[1]
    z = fluidBlock.GetPoint(i)[2]

    thefile.write(struct.pack('<3d' ,x, y, z))

    #thefile.writelines("%15.16f %15.16f %15.16f \n" %(x, y, z))

for i in range(nedge): #internal edges and barycentres
    ie = iedge[i]
    i1 = ie[0] + 1
    i2 = ie[1] + 1

    
    x1 = fluidBlock.GetPoint(ie[0])[0]
    y1 = fluidBlock.GetPoint(ie[0])[1]
    z1 = fluidBlock.GetPoint(ie[0])[2]

    x2 = fluidBlock.GetPoint(ie[1])[0]
    y2 = fluidBlock.GetPoint(ie[1])[1]
    z2 = fluidBlock.GetPoint(ie[1])[2]

    x = (x1+x2)/2
    y = (y1+y2)/2
    z = (z1+z2)/2
    
    thefile.write(struct.pack('<2i3d' ,i1,i2,x,y,z))
    
    #thefile.writelines("%8d %8d %15.16f %15.16f %15.16f \n" %(i1,i2,x,y,z))

for i in range(len(listUEbounds)): #Boundray edges and barycentres
    ie = listUEbounds[i]         # kbface
    i1 = ie[0] + 1
    i2 = ie[1] + 1
    x1 = fluidBlock.GetPoint(ie[0])[0]
    y1 = fluidBlock.GetPoint(ie[0])[1]
    z1 = fluidBlock.GetPoint(ie[0])[2]

    x2 = fluidBlock.GetPoint(ie[1])[0]
    y2 = fluidBlock.GetPoint(ie[1])[1]
    z2 = fluidBlock.GetPoint(ie[1])[2]

    x = (x1+x2)/2
    y = (y1+y2)/2
    z = (z1+z2)/2
    
    thefile.write(struct.pack('<2i3d' ,i1,i2,x,y,z))
    
    #thefile.writelines("%8d %8d %15.16f %15.16f %15.16f \n" %(i1,i2,x,y,z))
    

for i in range(nface): #internal face barycentres
    x = 0
    y = 0
    z = 0
    
    ele = iface[i]
    for ip in ele:
        x = x + fluidBlock.GetPoint(ip)[0]
        y = y + fluidBlock.GetPoint(ip)[1]
        z = z + fluidBlock.GetPoint(ip)[2]
        
    x = x/len(ele)
    y = y/len(ele)
    z = z/len(ele)

    thefile.write(struct.pack('<3d' ,x,y,z))

    #thefile.writelines("%15.16f %15.16f %15.16f \n" %(x,y,z))

for i in range(len(listFbounds)): #Boundary face barycentres
    x = 0
    y = 0
    z = 0
    
    ele = listFbounds[i]
    for ip in ele:
        x = x + fluidBlock.GetPoint(ip)[0]
        y = y + fluidBlock.GetPoint(ip)[1]
        z = z + fluidBlock.GetPoint(ip)[2]
        
    x = x/len(ele)
    y = y/len(ele)
    z = z/len(ele)

    thefile.write(struct.pack('<3d' ,x,y,z))
    
    #thefile.writelines("%15.16f %15.16f %15.16f \n" %(x,y,z))

for i in range(nelem): # elemental barycentres
    x = 0
    y = 0
    z = 0
    
    ele = ielem[i]
    for ip in ele:
        x = x + fluidBlock.GetPoint(ip)[0]
        y = y + fluidBlock.GetPoint(ip)[1]
        z = z + fluidBlock.GetPoint(ip)[2]
        
    x = x/len(ele)
    y = y/len(ele)
    z = z/len(ele)

    thefile.write(struct.pack('<3d' ,x,y,z))
    
    #thefile.writelines("%15.16f %15.16f %15.16f \n" %(x,y,z))

for i in range(len(listbpoin)): # ibpoin, boundary face cons and normals
    ibpoin = []
    ibpoin = listbpoin[i]

    
    thefile.write(struct.pack('<2i3d' ,(ibpoin[0]+1),ibpoin[1],ibpoin[2],ibpoin[3],ibpoin[4]))
    
    #thefile.writelines("%8d %8d %15.16f %15.16f %15.16f \n" %((ibpoin[0]+1),ibpoin[1],(ibpoin[2]),ibpoin[3],ibpoin[4]))

pointface = []
for i in range(nface):   # internal pointfaces
    ele = iface[i]
    for ip in ele:
        con = []
        con.append(ip+1)
        con.append(i+1)
        pointface.append(con)

pointface.sort()
thefile.write(struct.pack('<i' ,len(pointface)))
#thefile.writelines("%8d \n" %(len(pointface)))

for ie in pointface:
    thefile.write(struct.pack('<2i' ,ie[0],ie[1]))
    
    #thefile.writelines("%8d %8d \n" %(ie[0],ie[1]))

i1 = 0
pointface = []
for i in range(len(listFbounds)):   # boundary pointfaces
    ele = listFbounds[i]
    flag = listSBBbounds[i1]
    i1 = i1 + 1
    ib = 0 
    for ip in ele:
        con = []
        con.append(ip+1)
        con.append(i+1)
        con.append(flag[ib])
        pointface.append(con)
        ib = ib + 1
 
pointface.sort()
##
##ifa = []
##if2 = []
##
##pface = []
##if2.append(pointface[0])
##for ifa in pointface:
##    for i in range(len(ifa)):
##        if ifa[i] != if2[i]:
##            pface.append(ifa)
##            break
##    if2 = ifa

#thefile.writelines("%8d \n" %(len(pointface)))
thefile.write(struct.pack('<i' ,len(pointface)))
for ie in pointface:
    thefile.write(struct.pack('<3i' ,ie[0],ie[1],ie[2]))
    #thefile.writelines("%8d %8d %8d \n" %(ie[0],ie[1],ie[2]))

##pointface = []
##for i in range(len(listSBBbounds)):
##    ele = listSBBbounds[i]
##    for ip in ele:
##        con = []
##        con.append(ip+1)
##        con.append(i+1)
##        pointface.append(con)
##        
##pointface.sort()
##
##thefile.writelines("%8d \n" %(len(pointface)))
##for ie in pointface:
##    thefile.writelines("%8d %8d \n" %(ie[0],ie[1]))

pointelem = []
for i in range(nelem): # point elements
    
    ele = ielem[i]
    for ip in ele:
        con = []
        con.append(ip+1)
        con.append(i+1)
        pointelem.append(con)
        
pointelem.sort()

#thefile.writelines("%8d \n" %(len(pointelem)))

thefile.write(struct.pack('<i' ,len(pointelem)))

for ie in pointelem:
    thefile.write(struct.pack('<2i' ,ie[0],ie[1]))
    #thefile.writelines("%8d %8d \n" %(ie[0],ie[1]))

for ele in ielem:
    thefile.write(struct.pack('<i' ,len(ele)))
    #thefile.writelines("%8d \n" %(len(ele)))

for ele in ielem:
    for ip in ele:
        thefile.write(struct.pack('<i' ,ip))
        #thefile.writelines("%8d \n" %(ip))
    #thefile.writelines("\n" %())

i = 0
facecon = []
for ie in ipoly:# every element
    thefile.write(struct.pack('<i' ,ie[0]))
    #thefile.writelines("%8d  \n" %(ie[0]))

##    for f in ie[1:]:
##        thefile.writelines("%8d " %(f[0]))
##        thefile.writelines("\n" %())
##        for ip in range(f[0]):
##            thefile.writelines("%8d " %(f[ip+1]))
##        thefile.writelines("\n" %())
    #thefile.writelines("\n" %())


i = 0
facecon = []
for ie in ipoly:# every element
    #thefile.writelines("%8d  \n" %(ie[0]))

    for f in ie[1:]:
        thefile.write(struct.pack('<i' ,f[0]))
        #thefile.writelines("%8d \n" %(f[0]))
        #thefile.writelines("\n" %())
##        for ip in range(f[0]):
##            thefile.writelines("%8d " %(f[ip+1]))
##        thefile.writelines("\n" %())
    #thefile.writelines("\n" %())


i = 0
facecon = []
for ie in ipoly:# every element
    #thefile.writelines("%8d  \n" %(ie[0]))

    for f in ie[1:]:
##        thefile.writelines("%8d " %(f[0]))
##        thefile.writelines("\n" %())
        for ip in range(f[0]):
            thefile.write(struct.pack('<i' ,f[ip+1]))
            #thefile.writelines("%8d \n" %(f[ip+1]))
        #thefile.writelines("\n" %())
    #thefile.writelines("\n" %())






print('IO complete')
thefile.close()
