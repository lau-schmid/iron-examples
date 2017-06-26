#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import sys
import numpy as np

NumberGlobalXElements = 6
NumberGlobalYElements = 6
NumberGlobalZElements = 1

NumberOfElementsFE = NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements

NumberOfElementsInAtomX = 1
NumberOfElementsInAtomY = 1
NumberOfElementsInAtomZ = 1

NumberOfDomains = 10
 
print "input: x,y,z={},{},{}, ax,ay,az={},{},{}, NumberOfDomains={}".format(NumberGlobalXElements, NumberGlobalYElements, NumberGlobalZElements, NumberOfElementsInAtomX, NumberOfElementsInAtomY, NumberOfElementsInAtomZ, NumberOfDomains)
print "nElements: {}, n elments per domain: {}".format(NumberOfElementsFE, float(NumberOfElementsFE)/NumberOfDomains)

# An atom is a cuboid of NumberOfElementsInAtomX x NumberOfElementsInAtomY x NumberOfElementsInAtomZ 3D finite elasticity elements that will not be distributed to multiple processes.
# So this is an undividable unit for domain decomposition.
# Determine how many atoms are possible in each direction, fractions at the end count as full atom
nAtomsX = NumberGlobalXElements / NumberOfElementsInAtomX + np.sign(NumberGlobalXElements % NumberOfElementsInAtomX)
nAtomsY = NumberGlobalYElements / NumberOfElementsInAtomY + np.sign(NumberGlobalYElements % NumberOfElementsInAtomY)
nAtomsZ = NumberGlobalZElements / NumberOfElementsInAtomZ + np.sign(NumberGlobalZElements % NumberOfElementsInAtomZ)
nAtoms = nAtomsX*nAtomsY*nAtomsZ

print "nAtoms: {},{},{}={}".format(nAtomsX,nAtomsY,nAtomsZ,nAtoms)

# subdomain = all the atoms that belong to one process
# determine number of subdomains in each direction, such that domains are equally sized

# nSubdomainsXFloat = NumberGlobalXElements / OptimalSideLength
# nSubdomainsYFloat = NumberGlobalYElements / OptimalSideLength
# nSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
# nSubdomains = nSubdomainsXFloat * nSubdomainsYFloat * nSubdomainsZFloat 
#         = NumberGlobalXElements * NumberGlobalYElements * NumberGlobalZElements / (OptimalSideLength)**3
#         = NumberOfElementsFE / OptimalSideLength**3
# => OptimalSideLength = (NumberOfElementsFE / nSubdomains)**(1./3)

OptimalSideLength = (float(NumberOfElementsFE) / NumberOfDomains)**(1./3)

# if number of atoms is limited in z direction, begin with partioning in z direction
if float(NumberGlobalZElements) / OptimalSideLength > nAtomsZ:
  nSubdomainsZFloat = min(float(NumberGlobalZElements) / OptimalSideLength, nAtomsZ)
  print "nAtomsZ={}, OptimalSideLength={}, z={}, z/Opt={}, nSubdomainsZFloat={}".format(nAtomsZ, OptimalSideLength, NumberGlobalZElements, float(NumberGlobalZElements) / OptimalSideLength, nSubdomainsZFloat)

  OptimalSideLength = (nSubdomainsZFloat * NumberGlobalXElements * NumberGlobalYElements / NumberOfDomains)**(1./2)
  print "new OptimalSideLength={}".format(OptimalSideLength)
    
  nSubdomainsYFloat = min(float(NumberGlobalYElements) / OptimalSideLength, nAtomsY)
  print "nAtomsY={}, OptimalSideLength={}, y={}, y/Opt={}, nSubdomainsyFloat={}".format(nAtomsY, OptimalSideLength, NumberGlobalYElements, float(NumberGlobalYElements) / OptimalSideLength, nSubdomainsYFloat)

  nSubdomainsXFloat = min(float(NumberOfDomains) / (nSubdomainsZFloat * nSubdomainsYFloat), nAtomsX)
  print "nAtomsX={}, nSubdomainsXFloat={}, final: {}".format(nAtomsX, float(NumberOfDomains) / (nSubdomainsZFloat * nSubdomainsYFloat), nSubdomainsXFloat)

# begin partioning in  x direction
else:
  nSubdomainsXFloat = min(float(NumberGlobalXElements) / OptimalSideLength, nAtomsX)
  # now if nAtomX is smaller than the optimal nSubdomainsXFloat, we must take nAtomX instead
  # nSubdomainsXFloat is given
  # nSubdomainsYFloat = NumberGlobalYElements / OptimalSideLength
  # nSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
  # nSubdomains = nSubdomainsXFloat * nSubdomainsYFloat * nSubdomainsZFloat 
  #         = nSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / (OptimalSideLength)**2
  # => OptimalSideLength = (nSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / nSubdomains)**(1./2)
  # 

  print "nAtomsX={}, OptimalSideLength={}, x={}, x/Opt={}, nSubdomainsXFloat={}".format(nAtomsX, OptimalSideLength, NumberGlobalXElements, float(NumberGlobalXElements) / OptimalSideLength, nSubdomainsXFloat)

  OptimalSideLength = (nSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / NumberOfDomains)**(1./2)
  print "new OptimalSideLength={}".format(OptimalSideLength)

  nSubdomainsYFloat = min(float(NumberGlobalYElements) / OptimalSideLength, nAtomsY)
  print "nAtomsY={}, OptimalSideLength={}, y={}, y/Opt={}, nSubdomainsyFloat={}".format(nAtomsY, OptimalSideLength, NumberGlobalYElements, float(NumberGlobalYElements) / OptimalSideLength, nSubdomainsYFloat)

  # now if nAtomY is smaller than the optimal nSubdomainsYFloat, take nAtomsY
  # nSubdomainsXFloat and nSubdomainsYFloat are given
  # nSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
  # nSubdomains = nSubdomainsXFloat * nSubdomainsYFloat * nSubdomainsZFloat
  #         = nSubdomainsXFloat * nSubdomainsYFloat * NumberGlobalZElements / OptimalSideLength
  # => OptimalSideLength = nSubdomainsXFloat * nSubdomainsYFloat * NumberGlobalZElements / nSubdomains
  # nSubdomainsZFloat = NumberGlobalZElements / (nSubdomainsXFloat * nSubdomainsYFloat * NumberGlobalZElements / nSubdomains)
  #               = (NumberGlobalZElements * nSubdomains) / (nSubdomainsXFloat * nSubdomainsYFloat * NumberGlobalZElements)
  #               = nSubdomains / (nSubdomainsXFloat * nSubdomainsYFloat)
  #
   
  nSubdomainsZFloat = min(float(NumberOfDomains) / (nSubdomainsXFloat * nSubdomainsYFloat), nAtomsZ)

  print "nAtomsZ={}, nSubdomainsZFloat={}, final: {}".format(nAtomsZ, float(NumberOfDomains) / (nSubdomainsXFloat * nSubdomainsYFloat), nSubdomainsZFloat)
  
  
  
print "nSubdomainsFloat={},{},{}={}".format(nSubdomainsXFloat,nSubdomainsYFloat,nSubdomainsZFloat,nSubdomainsXFloat*nSubdomainsYFloat*nSubdomainsZFloat)

nSubdomainsX = max(1, int(np.round(nSubdomainsXFloat)))
nSubdomainsY = max(1, int(np.round(nSubdomainsYFloat)))
nSubdomainsZ = max(1, int(np.round(nSubdomainsZFloat)))

print "nSubdomains: {},{},{}={} ".format(nSubdomainsX, nSubdomainsY, nSubdomainsZ, nSubdomainsX*nSubdomainsY*nSubdomainsZ)

# adjust number of subdomains such that total number is <= number of domains (ideally '=')
while (nSubdomainsX*nSubdomainsY*nSubdomainsZ > NumberOfDomains):
  diffX = nSubdomainsX - nSubdomainsXFloat
  diffY = nSubdomainsY - nSubdomainsYFloat
  diffZ = nSubdomainsZ - nSubdomainsZFloat
  
  if diffX >= diffY and diffX >= diffZ:
    if nSubdomainsX != 1:
      nSubdomainsX -= 1
    elif diffY >= diffZ:
      if nSubdomainsY != 1:
        nSubdomainsY -= 1
      else:
        nSubdomainsZ -= 1
    else:
      if nSubdomainsZ != 1:
        nSubdomainsZ -= 1
      else:
        nSubdomainsY -= 1
        
  elif diffY >= diffZ:    # diffY >= diffZ, diffY > diffX
    if nSubdomainsY != 1:
      nSubdomainsY -= 1
    else:
      if diffX >= diffZ:
        if nSubdomainsX != 1:
          nSubdomainsX -= 1
        else:
          nSubdomainsZ -= 1
      else:
        if nSubdomainsZ != 1:
          nSubdomainsZ -= 1
        else:
          nSubdomainsX -= 1
      
  else:       # diffZ > diffY, diffZ >= diffX
    if nSubdomainsZ != 1:
      nSubdomainsZ -= 1
    else:
      if diffX >= diffY:
        if nSubdomainsX != 1:
          nSubdomainsX -= 1
        else:
          nSubdomainsY -= 1
      else:
        if nSubdomainsY != 1:
          nSubdomainsY -= 1
        else:
          nSubdomainsX -= 1
    
  print "diff:{},{},{}, nSubdomains: {},{},{}={} ".format(diffX,diffY,diffZ, nSubdomainsX, nSubdomainsY, nSubdomainsZ, nSubdomainsX*nSubdomainsY*nSubdomainsZ)
  
print "nSubdomains: {},{},{}={} ".format(nSubdomainsX, nSubdomainsY, nSubdomainsZ, nSubdomainsX*nSubdomainsY*nSubdomainsZ)

# determine shape of subdomain, i.e. nAtomsPerSubdomain
# nAtomsPerSubdomainX * nSubdomainsX = nAtomsX  =>  nAtomsPerSubdomainX = nAtomsX / nSubdomainsX

nAtomsPerSubdomainX = int(np.ceil(float(nAtomsX) / nSubdomainsX))
nAtomsPerSubdomainY = int(np.ceil(float(nAtomsY) / nSubdomainsY))
nAtomsPerSubdomainZ = int(np.ceil(float(nAtomsZ) / nSubdomainsZ))
nAtomsPerSubdomain = nAtomsPerSubdomainX*nAtomsPerSubdomainY*nAtomsPerSubdomainZ

# decrease number of subdomains to exclude now empty subdomains
nEmptySubdomainsX = int(float(nAtomsPerSubdomainX*nSubdomainsX - nAtomsX) / nAtomsPerSubdomainX)
nEmptySubdomainsY = int(float(nAtomsPerSubdomainY*nSubdomainsY - nAtomsY) / nAtomsPerSubdomainY)
nEmptySubdomainsZ = int(float(nAtomsPerSubdomainZ*nSubdomainsZ - nAtomsZ) / nAtomsPerSubdomainZ)

nEmptySubdomains = nSubdomainsX*nSubdomainsY*nSubdomainsZ - (nSubdomainsX-nEmptySubdomainsX)*(nSubdomainsY-nEmptySubdomainsY)*(nSubdomainsZ-nEmptySubdomainsZ)

print "nEmptySubdomains:{},{},{}, total:{}".format(nEmptySubdomainsX,nEmptySubdomainsY,nEmptySubdomainsZ, nEmptySubdomains)

nSubdomainsX -= nEmptySubdomainsX
nSubdomainsY -= nEmptySubdomainsY
nSubdomainsZ -= nEmptySubdomainsZ
nSubdomains = nSubdomainsX*nSubdomainsY*nSubdomainsZ

nUnusedSubdomains = NumberOfDomains - nSubdomains 

if nUnusedSubdomains != 0:
  print "Warning! {} Process(es) will not be used".format(nUnusedSubdomains)

print "nAtomsPerSubdomain: {},{},{}={}".format(nAtomsPerSubdomainX, nAtomsPerSubdomainY, nAtomsPerSubdomainZ, nAtomsPerSubdomain)
print "nSubdomains: {},{},{}={}, nAtomsPerSubdomain: {},{},{}={},  NumberOfElementsInAtom: {},{},{}={}, NumberGlobalElements: {},{},{}={}".\
  format(nSubdomainsX, nSubdomainsY, nSubdomainsZ, nSubdomainsX*nSubdomainsY*nSubdomainsZ,\
  nAtomsPerSubdomainX,nAtomsPerSubdomainY,nAtomsPerSubdomainZ,nAtomsPerSubdomain, \
  NumberOfElementsInAtomX, NumberOfElementsInAtomY, NumberOfElementsInAtomZ, NumberOfElementsInAtomX*NumberOfElementsInAtomY*NumberOfElementsInAtomX,\
  NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements,NumberGlobalXElements*NumberGlobalYElements,NumberGlobalZElements)

print "subdomain size: {}x{}x{} elements".format(NumberOfElementsInAtomX*nAtomsPerSubdomainX, NumberOfElementsInAtomY*nAtomsPerSubdomainY, NumberOfElementsInAtomZ*nAtomsPerSubdomainZ)
print "total size without considering partially filled subdomains {}x{}x{} elements >=   available {}x{}x{} elements".\
format(nSubdomainsX*nAtomsPerSubdomainX*NumberOfElementsInAtomX, nSubdomainsY*nAtomsPerSubdomainY*NumberOfElementsInAtomY, nSubdomainsZ*nAtomsPerSubdomainZ*NumberOfElementsInAtomZ, \
NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements)


# determine last subdomain sizes
# nAtomsLastSubdomainX = nAtomsPerSubdomainX - (nAtomsPerSubdomainX*nSubdomainsX - nAtomsX)
nAtomsLastSubdomainX = nAtomsPerSubdomainX*(1-nSubdomainsX) + nAtomsX
nAtomsLastSubdomainY = nAtomsPerSubdomainY*(1-nSubdomainsY) + nAtomsY
nAtomsLastSubdomainZ = nAtomsPerSubdomainZ*(1-nSubdomainsZ) + nAtomsZ

# NumberOfElementsLastAtomX = NumberOfElementsInAtomX - (nAtomsX*NumberOfElementsInAtomX - NumberGlobalXElements)
# NumberOfElementsLastAtomX = NumberOfElementsInAtomX*(1 - nAtomsX) + NumberGlobalXElements
NumberOfElementsLastAtomX = NumberOfElementsInAtomX*(1-nAtomsX) + NumberGlobalXElements
NumberOfElementsLastAtomY = NumberOfElementsInAtomY*(1-nAtomsY) + NumberGlobalYElements
NumberOfElementsLastAtomZ = NumberOfElementsInAtomZ*(1-nAtomsZ) + NumberGlobalZElements


print "nAtomsLastSubdomain: {},{},{}, NumberOfElementsLastAtom: {},{},{}".format(nAtomsLastSubdomainX, nAtomsLastSubdomainY, nAtomsLastSubdomainZ, NumberOfElementsLastAtomX, NumberOfElementsLastAtomY, NumberOfElementsLastAtomZ)

# for special subdomains (at x+,y+,z+ border) output the number of missing elements

normalNumberOfElements = \
  nAtomsPerSubdomainX*NumberOfElementsInAtomX \
  *nAtomsPerSubdomainY*NumberOfElementsInAtomY\
  *nAtomsPerSubdomainZ*NumberOfElementsInAtomZ

print "{} processes, normal number of elements per process: {}".format(nSubdomains, normalNumberOfElements)

# x+ face
actualNumberOfElements = \
  ((nAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) \
  *nAtomsPerSubdomainY*NumberOfElementsInAtomY\
  *nAtomsPerSubdomainZ*NumberOfElementsInAtomZ
      
if (nSubdomainsY-1)*(nSubdomainsZ-1) > 0 and actualNumberOfElements < normalNumberOfElements:
  print "{} process(es) on x+ boundary have {} elements less (only {}%)".format((nSubdomainsY-1)*(nSubdomainsZ-1), normalNumberOfElements-actualNumberOfElements, 100.*actualNumberOfElements/normalNumberOfElements)
      
# y+ face
actualNumberOfElements = \
  nAtomsPerSubdomainX*NumberOfElementsInAtomX \
  *((nAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) \
  *nAtomsPerSubdomainZ*NumberOfElementsInAtomZ
      
if (nSubdomainsX-1)*(nSubdomainsZ-1) > 0 and actualNumberOfElements < normalNumberOfElements:
  print "{} process(es) on y+ boundary have {} elements less (only {}%)".format((nSubdomainsX-1)*(nSubdomainsZ-1), normalNumberOfElements-actualNumberOfElements, 100.*actualNumberOfElements/normalNumberOfElements)
    
# z+ face
actualNumberOfElements = \
  nAtomsPerSubdomainX*NumberOfElementsInAtomX \
  *nAtomsPerSubdomainY*NumberOfElementsInAtomY \
  *((nAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
      
if (nSubdomainsX-1)*(nSubdomainsY-1) > 0 and actualNumberOfElements < normalNumberOfElements:
  print "{} process(es) on z+ boundary have {} elements less (only {}%)".format((nSubdomainsX-1)*(nSubdomainsY-1), normalNumberOfElements-actualNumberOfElements, 100.*actualNumberOfElements/normalNumberOfElements)
      
# x+ y+ edge
actualNumberOfElements = \
  ((nAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) \
  *((nAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) \
  *nAtomsPerSubdomainZ*NumberOfElementsInAtomZ
      
if nSubdomainsZ-1 > 0 and actualNumberOfElements < normalNumberOfElements:
  print "{} process(es) on x+,y+ boundary have {} elements less (only {}%)".format(nSubdomainsZ-1, normalNumberOfElements-actualNumberOfElements, 100.*actualNumberOfElements/normalNumberOfElements)
      
# x+ z+ edge
actualNumberOfElements = \
  ((nAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) \
  *nAtomsPerSubdomainY*NumberOfElementsInAtomY \
  *((nAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
      
if nSubdomainsY-1 > 0 and actualNumberOfElements < normalNumberOfElements:
  print "{} process(es) on x+,z+ boundary have {} elements less (only {}%)".format(nSubdomainsY-1, normalNumberOfElements-actualNumberOfElements, 100.*actualNumberOfElements/normalNumberOfElements)
      
# y+ z+ edge
actualNumberOfElements = \
  nAtomsPerSubdomainX*NumberOfElementsInAtomX \
  *((nAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) \
  *((nAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
      
if nSubdomainsX-1 > 0 and actualNumberOfElements < normalNumberOfElements:
  print "{} process(es) on y+,z+ boundary have {} elements less (only {}%)".format(nSubdomainsX-1, normalNumberOfElements-actualNumberOfElements, 100.*actualNumberOfElements/normalNumberOfElements)
    
# x+ y+ z+ subdomain
actualNumberOfElements = \
  ((nAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) \
  *((nAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) \
  *((nAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
      
if actualNumberOfElements < normalNumberOfElements:
  print "1 process on x+,y+,z+ boundary has {} elements less (only {}%)".format(normalNumberOfElements-actualNumberOfElements, 100.*actualNumberOfElements/normalNumberOfElements)

sys.exit(0)

# iterate over subdomains
for subdomainZ in range(nSubdomainsZ):
  for subdomainY in range(nSubdomainsY):
    for subdomainX in range(nSubdomainsX):
      
      DomainNo = subdomainZ*nSubdomainsY*nSubdomainsX + subdomainY*nSubdomainsX + subdomainX
      
      nAtomsCurrentSubdomainX = nAtomsPerSubdomainX
      if subdomainX == nSubdomainsX-1:  # last subdomain
        nAtomsCurrentSubdomainX = nAtomsLastSubdomainX
        
      nAtomsCurrentSubdomainY = nAtomsPerSubdomainY
      if subdomainY == nSubdomainsY-1:  # last subdomain
        nAtomsCurrentSubdomainY = nAtomsLastSubdomainY
        
      nAtomsCurrentSubdomainZ = nAtomsPerSubdomainZ
      if subdomainZ == nSubdomainsZ-1:  # last subdomain
        nAtomsCurrentSubdomainZ = nAtomsLastSubdomainZ
        
      print "subdomain ({},{},{}) has {},{},{} atoms".format(subdomainX,subdomainY,subdomainZ,nAtomsCurrentSubdomainX,nAtomsCurrentSubdomainY,nAtomsCurrentSubdomainZ)
      for atomZ in range(nAtomsCurrentSubdomainZ):
        for atomY in range(nAtomsCurrentSubdomainY):
          for atomX in range(nAtomsCurrentSubdomainX):
              
            nElementsCurrentAtomX = NumberOfElementsInAtomX
            if subdomainX == nSubdomainsX-1 and atomX == nAtomsCurrentSubdomainX-1:  # last atom
              nElementsCurrentAtomX = NumberOfElementsLastAtomX
              
            nElementsCurrentAtomY = NumberOfElementsInAtomY
            if subdomainY == nSubdomainsY-1 and atomY == nAtomsCurrentSubdomainY-1:  # last atom
              nElementsCurrentAtomY = NumberOfElementsLastAtomY
              
            nElementsCurrentAtomZ = NumberOfElementsInAtomZ
            if subdomainZ == nSubdomainsZ-1 and atomZ == nAtomsCurrentSubdomainZ-1:  # last atom
              nElementsCurrentAtomZ = NumberOfElementsLastAtomZ
              
            print "  atom ({},{},{}) has {},{},{} elements".format(atomX,atomY,atomZ,nElementsCurrentAtomX,nElementsCurrentAtomY,nElementsCurrentAtomZ)
              
            for z in range(nElementsCurrentAtomZ):
              for y in range(nElementsCurrentAtomY):
                for x in range(nElementsCurrentAtomX):
                  
                  zIndex = subdomainZ*nAtomsPerSubdomainZ*NumberOfElementsInAtomZ + atomZ*NumberOfElementsInAtomZ + z
                  yIndex = subdomainY*nAtomsPerSubdomainY*NumberOfElementsInAtomY + atomY*NumberOfElementsInAtomY + y
                  xIndex = subdomainX*nAtomsPerSubdomainX*NumberOfElementsInAtomX + atomX*NumberOfElementsInAtomX + x
                  
                  Index = zIndex * NumberGlobalXElements * NumberGlobalYElements + yIndex * NumberGlobalXElements + xIndex
                  
                  
                  print "      (x,y,z)=({},{},{}) i={} domain {}".format(xIndex,yIndex,zIndex,Index,DomainNo)
   
