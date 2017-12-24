
MODULE DECOMPOSITION

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE OpenCMISS_Variables
  USE REGION_MESH
  
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

CONTAINS

! return ceil(X/Y) for integers
FUNCTION CeilDiv(X,Y)
  INTEGER(CMISSIntg), INTENT(IN) :: X, Y
  REAL(CMISSRP) :: CeilDiv
  IF (MOD(X,Y) == 0) THEN
    CeilDiv = X/Y
  ELSE
    CeilDiv = X/Y + 1
  ENDIF
END FUNCTION CeilDiv

FUNCTION GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, &
  & ParameterNumberOfSubdomainsX, ParameterNumberOfSubdomainsY, ParameterNumberOfSubdomainsZ)
  INTEGER(CMISSIntg), INTENT(IN) :: NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ
  INTEGER(CMISSIntg), INTENT(IN) :: ParameterNumberOfSubdomainsX, ParameterNumberOfSubdomainsY, ParameterNumberOfSubdomainsZ
  INTEGER(CMISSIntg) :: NumberOfSubdomainsX, NumberOfSubdomainsY, NumberOfSubdomainsZ
  INTEGER(CMISSIntg) :: nEmptySubdomainsX, nEmptySubdomainsY, nEmptySubdomainsZ, nEmptySubdomains
  INTEGER(CMISSIntg) :: nSubdomains
  INTEGER(CMISSIntg) :: GetNumberOfUsedSubdomains, NumberOfAtomsPerSubdomain
  
  IF (ParameterNumberOfSubdomainsX <= 0 .OR. ParameterNumberOfSubdomainsY <= 0 .OR. ParameterNumberOfSubdomainsZ <= 0) THEN
    GetNumberOfUsedSubdomains = HUGE(GetNumberOfUsedSubdomains)
    RETURN
  ENDIF
  
  NumberOfAtomsPerSubdomainX = CEILING(DBLE(NumberOfAtomsX) / ParameterNumberOfSubdomainsX)
  NumberOfAtomsPerSubdomainY = CEILING(DBLE(NumberOfAtomsY) / ParameterNumberOfSubdomainsY)
  NumberOfAtomsPerSubdomainZ = CEILING(DBLE(NumberOfAtomsZ) / ParameterNumberOfSubdomainsZ)
  NumberOfAtomsPerSubdomain = NumberOfAtomsPerSubdomainX*NumberOfAtomsPerSubdomainY*NumberOfAtomsPerSubdomainX

  ! decrease number of subdomains to exclude now empty subdomains
  nEmptySubdomainsX = FLOOR(DBLE(NumberOfAtomsPerSubdomainX*ParameterNumberOfSubdomainsX - NumberOfAtomsX) &
    & / NumberOfAtomsPerSubdomainX)
  nEmptySubdomainsY = FLOOR(DBLE(NumberOfAtomsPerSubdomainY*ParameterNumberOfSubdomainsY - NumberOfAtomsY) &
    & / NumberOfAtomsPerSubdomainY)
  nEmptySubdomainsZ = FLOOR(DBLE(NumberOfAtomsPerSubdomainZ*ParameterNumberOfSubdomainsZ - NumberOfAtomsZ) &
    & / NumberOfAtomsPerSubdomainZ)

  nEmptySubdomains = ParameterNumberOfSubdomainsX*ParameterNumberOfSubdomainsY*ParameterNumberOfSubdomainsZ &
    & - (ParameterNumberOfSubdomainsX-nEmptySubdomainsX)*(ParameterNumberOfSubdomainsY-nEmptySubdomainsY) &
    & * (ParameterNumberOfSubdomainsZ-nEmptySubdomainsZ)

  NumberOfSubdomainsX = ParameterNumberOfSubdomainsX - nEmptySubdomainsX
  NumberOfSubdomainsY = ParameterNumberOfSubdomainsY - nEmptySubdomainsY
  NumberOfSubdomainsZ = ParameterNumberOfSubdomainsZ - nEmptySubdomainsZ
  nSubdomains = NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ

  GetNumberOfUsedSubdomains = nSubdomains
END FUNCTION GetNumberOfUsedSubdomains
  
SUBROUTINE ComputeSubdomainsWithAtoms()
  REAL(CMISSDP) :: OptimalSideLength
  INTEGER(CMISSIntg) :: NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfAtoms
  REAL(CMISSDP) :: NumberOfSubdomainsXFloat, NumberOfSubdomainsYFloat, NumberOfSubdomainsZFloat
  REAL(CMISSDP) :: DiffX, DiffY, DiffZ
  INTEGER(CMISSIntg) :: NumberOfAtomsPerSubdomain, nUnusedSubdomains
  INTEGER(CMISSIntg) :: nEmptySubdomainsX, nEmptySubdomainsY, nEmptySubdomainsZ, nEmptySubdomains
  INTEGER(CMISSIntg) :: nSubdomains
  INTEGER(CMISSIntg) :: NormalNumberOfElements, actualNumberOfElements
  INTEGER(CMISSIntg) :: DiffNumberOfDomainsXDecreased, DiffNumberOfDomainsYDecreased, DiffNumberOfDomainsZDecreased
  INTEGER(CMISSIntg) :: DiffNumberOfDomainsXYDecreased, DiffNumberOfDomainsXZDecreased, DiffNumberOfDomainsYZDecreased
  INTEGER(CMISSIntg) :: DiffNumberOfDomainsXYZDecreased, MinDiffNumberOfDomains
  INTEGER(CMISSIntg) :: PretendedNumberOfDomains     !< this is a value for the number of domains to be used in domain decomposition. It can be higher than the actual number of processes, because sometimes domain decomposition output produces less domains than requested.
  LOGICAL :: DEBUGGING = .FALSE.
  
  IF (DEBUGGING) DEBUGGING = (ComputationalNodeNumber == 0)  
  
  ! Test if a valid domain decomposition was specified via parameters NumberOfElementsInSubdomainXYZ and NumberOfSubdomainsXYZ
  ! if yes, use this domain decomposition and return, if not, control proceeds to computing the domain decomposition from paraameters of "atom" sizes
  IF (NumberOfElementsInSubdomainX /= 0 .AND. NumberOfElementsInSubdomainY /= 0 .AND. NumberOfElementsInSubdomainZ /= 0) THEN
    
    ! if specified number of subdomains does not match number of processes
    IF (NumberOfSubdomainsX * NumberOfSubdomainsY * NumberOfSubdomainsZ /= NumberOfComputationalNodes) THEN
      
      IF (ComputationalNodeNumber == 0) THEN
        PRINT *, "The specified number of subdomains ",NumberOfSubdomainsX, "x", NumberOfSubdomainsY, "x", NumberOfSubdomainsZ, & 
          & " = ", NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ, " does not match the number of processes (", &
          & NumberOfComputationalNodes,"). Now the domain decomposition is computed from atom sizes."
      ENDIF 
    ELSE
      ! number of subdomains does match number of processes, take the specified domain decomposition

      IF (NumberOfElementsInSubdomainX*NumberOfElementsInSubdomainY*NumberOfElementsInSubdomainZ*NumberOfComputationalNodes &
        & /= NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements) THEN

        IF (ComputationalNodeNumber == 0) THEN
          PRINT *, "The specified number of 3D elements per subdomain e=",NumberOfElementsInSubdomainX, "x", &
            & NumberOfElementsInSubdomainY, "x", NumberOfElementsInSubdomainZ, " * number of subdomains (", &
            & NumberOfComputationalNodes,"), p=",NumberOfSubdomainsX, "x", NumberOfSubdomainsY, "x", NumberOfSubdomainsZ, &
            & " does not match the total number of 3D elements (", &
            & NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements, &
            & "). Now the domain decomposition is computed from atom sizes."
        ENDIF 
      ELSE
        ! The following global variables are used later in the code to assign elements to domains and therefore need to be set here:
        ! NumberOfSubdomainsXYZ, NumberOfAtomsPerSubdomainXYZ, NumberOfAtomsLastSubdomainXYZ, NumberOfElementsInAtomXYZ, NumberOfElementsLastAtomXYZ

        ! For this approach without given "atom" sizes, set the atom size to be of size 1 element

        ! NumberOfSubdomainsXYZ: already set as parameter and checked, that it is valid

        ! NumberOfAtomsPerSubdomainXYZ
        NumberOfAtomsPerSubdomainX = NumberOfElementsInSubdomainX
        NumberOfAtomsPerSubdomainY = NumberOfElementsInSubdomainY
        NumberOfAtomsPerSubdomainZ = NumberOfElementsInSubdomainZ

        ! NumberOfAtomsLastSubdomainXYZ
        NumberOfAtomsLastSubdomainX = NumberOfAtomsPerSubdomainX
        NumberOfAtomsLastSubdomainY = NumberOfAtomsPerSubdomainY
        NumberOfAtomsLastSubdomainZ = NumberOfAtomsPerSubdomainZ

        ! NumberOfElementsInAtomXYZ
        ! set atom sizes to 1
        NumberOfElementsInAtomX = 1
        NumberOfElementsInAtomY = 1
        NumberOfElementsInAtomZ = 1

        ! NumberOfElementsLastAtomXYZ
        NumberOfElementsLastAtomX = 1
        NumberOfElementsLastAtomY = 1
        NumberOfElementsLastAtomZ = 1

        IF (ComputationalNodeNumber == 0) THEN
          PRINT *, "Using given domain decomposition."
        ENDIF

        RETURN
      ENDIF
    ENDIF

  ENDIF


  ! The following global variables are used later in the code to assign elements to domains and therefore need to be set here:
  ! NumberOfSubdomainsXYZ, NumberOfAtomsPerSubdomainXYZ, NumberOfAtomsLastSubdomainXYZ, NumberOfElementsInAtomXYZ, NumberOfElementsLastAtomXYZ

  ! Domain decomposition for the 3D finite elasticity mesh is computed as follows: 
  ! The user can specify the size of an 'atomic' cuboid, e.g. with ax*ay*az elements. These 'atoms' then must not be split to multiple processors.
  ! In that way it is possibile to forbid the subdivision of a 1D fibre mesh to multiple processes, or to restrict it to only 2 processes per fibre.
  ! The decomposition is then done such that the cuts are 2D plains, in a way that the subdomains are as cube-shaped as possible 
  ! and as evenly distributed as possible while fulfilling the constraint.
  ! This means that not all processes may get the same amount of elements and that some processes will be left idle.
  ! This is intended and makes sense for scaling measurements:
  ! E.g. the case of 6x6x1 elements total with 10 processes. It will give 2x2x1 portions to each of the 9 processes, not using 1 process. This is better than using all 10 processes.

  PretendedNumberOfDomains = PretendedNumberOfDomainsForDomainDecomposition
  
  ! ----------  
  ! An atom is a cuboid of NumberOfElementsInAtomX x NumberOfElementsInAtomY x NumberOfElementsInAtomZ 3D finite elasticity elements that will not be distributed to multiple processes.
  ! So this is an undividable unit for domain decomposition.
  ! Determine how many atoms are possible in each direction, fractions at the end count as full atom
  NumberOfAtomsX = CeilDiv(NumberGlobalXElements, NumberOfElementsInAtomX)
  NumberOfAtomsY = CeilDiv(NumberGlobalYElements, NumberOfElementsInAtomY)
  NumberOfAtomsZ = CeilDiv(NumberGlobalZElements, NumberOfElementsInAtomZ)
  NumberOfAtoms = NumberOfAtomsX * NumberOfAtomsY * NumberOfAtomsZ

  IF (DEBUGGING) PRINT *, "PretendedNumberOfDomains=",PretendedNumberOfDomains
  IF (DEBUGGING) PRINT *, "NumberOfAtoms: ",NumberOfAtomsX,NumberOfAtomsY,NumberOfAtomsZ,"=",NumberOfAtoms

  ! subdomain = all the atoms that belong to one process
  ! determine number of subdomains in each direction, such that domains are equally sized

  ! NumberOfSubdomainsXFloat = NumberGlobalXElements / OptimalSideLength
  ! NumberOfSubdomainsYFloat = NumberGlobalYElements / OptimalSideLength
  ! NumberOfSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
  ! nSubdomains = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberOfSubdomainsZFloat 
  !         = NumberGlobalXElements * NumberGlobalYElements * NumberGlobalZElements / (OptimalSideLength)**3
  !         = NumberOfElementsFE / OptimalSideLength**3
  ! => OptimalSideLength = (NumberOfElementsFE / nSubdomains)**(1./3)

  OptimalSideLength = (DBLE(NumberOfElementsFE) / PretendedNumberOfDomains)**(1./3)

  IF (DEBUGGING) PRINT *, "OptimalSideLength:", OptimalSideLength,", avg. number elements z: ", &
    & DBLE(NumberGlobalZElements) / OptimalSideLength, ">", NumberOfAtomsZ,"?"
  
  IF (DBLE(NumberGlobalZElements) / OptimalSideLength > NumberOfAtomsZ - 1e-3) THEN
      
    NumberOfSubdomainsZFloat = MIN(DBLE(NumberGlobalZElements) / OptimalSideLength, DBLE(NumberOfAtomsZ))
    
    IF (DEBUGGING) PRINT *, "YES, NumberOfAtomsZ =", NumberOfAtomsZ, ", OptimalSideLength =",OptimalSideLength,", z=",&
       & NumberGlobalZElements, ", z/Opt=", DBLE(NumberGlobalZElements) / OptimalSideLength, ", NumberOfSubdomainsZFloat=", &
       & NumberOfSubdomainsZFloat

    OptimalSideLength = SQRT(NumberOfSubdomainsZFloat * NumberGlobalXElements * NumberGlobalYElements / PretendedNumberOfDomains)
    IF (DEBUGGING) PRINT *, "new OptimalSideLength=", OptimalSideLength
    
    NumberOfSubdomainsYFloat = MIN(DBLE(NumberGlobalYElements) / OptimalSideLength, DBLE(NumberOfAtomsY))
    IF (DEBUGGING) PRINT *, "NumberOfAtomsY=", NumberOfAtomsY, ", OptimalSideLength=",OptimalSideLength,", y=",&
      & NumberGlobalYElements, ", y/Opt=", DBLE(NumberGlobalYElements) / OptimalSideLength, ", NumberOfSubdomainsYFloat=", &
      & NumberOfSubdomainsYFloat

    NumberOfSubdomainsXFloat = MIN(DBLE(PretendedNumberOfDomains) / (NumberOfSubdomainsZFloat * NumberOfSubdomainsYFloat), &
      & DBLE(NumberOfAtomsX))
    IF (DEBUGGING) PRINT *, "NumberOfAtomsX=", NumberOfAtomsX, ", NumberOfSubdomainsXFloat=", &
       & DBLE(PretendedNumberOfDomains) / (NumberOfSubdomainsZFloat * NumberOfSubdomainsYFloat), ", final: ", &
       & NumberOfSubdomainsXFloat
    
  ! begin partioning in x direction
  ELSE
    
    NumberOfSubdomainsXFloat = MIN(DBLE(NumberGlobalXElements) / OptimalSideLength, DBLE(NumberOfAtomsX))
    ! now IF nAtomX is smaller than the optimal NumberOfSubdomainsXFloat, we must take nAtomX instead
    ! NumberOfSubdomainsXFloat is given
    ! NumberOfSubdomainsYFloat = NumberGlobalYElements / OptimalSideLength
    ! NumberOfSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
    ! nSubdomains = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberOfSubdomainsZFloat 
    !         = NumberOfSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / (OptimalSideLength)**2
    ! => OptimalSideLength = (NumberOfSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / nSubdomains)**(1./2)
    ! 

    IF (DEBUGGING) PRINT *, "NO, NumberOfAtomsX =", NumberOfAtomsX, ", OptimalSideLength =",OptimalSideLength,", x=",&
       & NumberGlobalXElements, ", x/Opt=", DBLE(NumberGlobalXElements) / OptimalSideLength, ", NumberOfSubdomainsXFloat=", &
       & NumberOfSubdomainsXFloat

    OptimalSideLength = SQRT(NumberOfSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / PretendedNumberOfDomains)
    IF (DEBUGGING) PRINT *, "new OptimalSideLength=", OptimalSideLength, "=sqrt(",NumberOfSubdomainsXFloat," * ",&
      & NumberGlobalYElements," * ",NumberGlobalZElements,")/",PretendedNumberOfDomains

    NumberOfSubdomainsYFloat = MIN(DBLE(NumberGlobalYElements) / OptimalSideLength, DBLE(NumberOfAtomsY))
    IF (DEBUGGING) PRINT *, "NumberOfAtomsY=", NumberOfAtomsY, ", OptimalSideLength=",OptimalSideLength,", y=",&
      & NumberGlobalYElements, ", y/Opt=", DBLE(NumberGlobalYElements) / OptimalSideLength, ", NumberOfSubdomainsYFloat=", &
      & NumberOfSubdomainsYFloat

    ! now IF nAtomY is smaller than the optimal NumberOfSubdomainsYFloat, take NumberOfAtomsY
    ! NumberOfSubdomainsXFloat .AND. NumberOfSubdomainsYFloat are given
    ! NumberOfSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
    ! nSubdomains = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberOfSubdomainsZFloat
    !         = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberGlobalZElements / OptimalSideLength
    ! => OptimalSideLength = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberGlobalZElements / nSubdomains
    ! NumberOfSubdomainsZFloat = NumberGlobalZElements / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberGlobalZElements / nSubdomains)
    !               = (NumberGlobalZElements * nSubdomains) / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberGlobalZElements)
    !               = nSubdomains / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat)
    !
     
    NumberOfSubdomainsZFloat = MIN(DBLE(PretendedNumberOfDomains) / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat), &
      & DBLE(NumberOfAtomsZ))

    IF (DEBUGGING) PRINT *, "NumberOfAtomsZ=", NumberOfAtomsZ, ", NumberOfSubdomainsZFloat=", &
       & DBLE(PretendedNumberOfDomains) / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat), ", final: ", &
       & NumberOfSubdomainsZFloat
  ENDIF
       
  IF (DEBUGGING) PRINT *, "nSubdomainsFloat=", NumberOfSubdomainsXFloat,NumberOfSubdomainsYFloat,NumberOfSubdomainsZFloat,"=",&
    & NumberOfSubdomainsXFloat*NumberOfSubdomainsYFloat*NumberOfSubdomainsZFloat

  NumberOfSubdomainsX = MAX(1, CEILING(NumberOfSubdomainsXFloat))     ! try to round up to use the maximum number of processes possible, if value gets to high, it will be decreased anyway
  NumberOfSubdomainsY = MAX(1, CEILING(NumberOfSubdomainsYFloat))
  NumberOfSubdomainsZ = MAX(1, CEILING(NumberOfSubdomainsZFloat))

  IF (DEBUGGING) PRINT *, "nSubdomains: ",NumberOfSubdomainsX, NumberOfSubdomainsY, NumberOfSubdomainsZ, "=", &
    & NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ

  ! adjust number of subdomains such that total number is <= number of domains (ideally '=')
  
  IF (GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, NumberOfSubdomainsY, &
     & NumberOfSubdomainsZ) > PretendedNumberOfDomains) THEN
    DiffNumberOfDomainsXDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX-1, NumberOfSubdomainsY, &
      & NumberOfSubdomainsZ)
    DiffNumberOfDomainsYDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, NumberOfSubdomainsY-1, &
      & NumberOfSubdomainsZ)
    DiffNumberOfDomainsZDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, NumberOfSubdomainsY, &
      & NumberOfSubdomainsZ-1)
    DiffNumberOfDomainsXYDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX-1, NumberOfSubdomainsY-1, &
      & NumberOfSubdomainsZ)
    DiffNumberOfDomainsXZDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX-1, NumberOfSubdomainsY, &
      & NumberOfSubdomainsZ-1)
    DiffNumberOfDomainsYZDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, NumberOfSubdomainsY-1, &
      & NumberOfSubdomainsZ-1)
    DiffNumberOfDomainsXYZDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX-1, NumberOfSubdomainsY-1, &
      & NumberOfSubdomainsZ-1)

    IF (DiffNumberOfDomainsXDecreased < 0) DiffNumberOfDomainsXDecreased = HUGE(DiffNumberOfDomainsXDecreased)
    IF (DiffNumberOfDomainsYDecreased < 0) DiffNumberOfDomainsYDecreased = HUGE(DiffNumberOfDomainsYDecreased)
    IF (DiffNumberOfDomainsZDecreased < 0) DiffNumberOfDomainsZDecreased = HUGE(DiffNumberOfDomainsZDecreased)
    IF (DiffNumberOfDomainsXYDecreased < 0) DiffNumberOfDomainsXYDecreased = HUGE(DiffNumberOfDomainsXYDecreased)
    IF (DiffNumberOfDomainsXZDecreased < 0) DiffNumberOfDomainsXZDecreased = HUGE(DiffNumberOfDomainsXZDecreased)
    IF (DiffNumberOfDomainsYZDecreased < 0) DiffNumberOfDomainsYZDecreased = HUGE(DiffNumberOfDomainsYZDecreased)
    IF (DiffNumberOfDomainsXYZDecreased < 0) DiffNumberOfDomainsXYZDecreased = HUGE(DiffNumberOfDomainsXYZDecreased)
    
    MinDiffNumberOfDomains = MIN(DiffNumberOfDomainsXDecreased, MIN(DiffNumberOfDomainsYDecreased, &
      & MIN(DiffNumberOfDomainsZDecreased, MIN(DiffNumberOfDomainsXYDecreased, MIN(DiffNumberOfDomainsXZDecreased, &
      & MIN(DiffNumberOfDomainsYZDecreased, DiffNumberOfDomainsXYZDecreased))))))
      
    IF (DEBUGGING) PRINT *, "diffs: x:",DiffNumberOfDomainsXDecreased, ", y:",DiffNumberOfDomainsYDecreased, &
      & ", z:",DiffNumberOfDomainsZDecreased, ", xy:",DiffNumberOfDomainsXYDecreased, ", xz:",DiffNumberOfDomainsXZDecreased, &
      & ", yz:", DiffNumberOfDomainsYZDecreased, "xyz:",DiffNumberOfDomainsXYZDecreased

      
    IF (MinDiffNumberOfDomains == DiffNumberOfDomainsXDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) THEN
      IF (DEBUGGING) PRINT *, "best to decrease X by 1"
      NumberOfSubdomainsX = NumberOfSubdomainsX-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsYDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease Y by 1"
      NumberOfSubdomainsY = NumberOfSubdomainsY-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsZDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease Z by 1"
      NumberOfSubdomainsZ = NumberOfSubdomainsZ-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsXYDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease X and Y by 1"
      NumberOfSubdomainsX = NumberOfSubdomainsX-1
      NumberOfSubdomainsY = NumberOfSubdomainsY-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsXZDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease X and Z by 1"
      NumberOfSubdomainsX = NumberOfSubdomainsX-1
      NumberOfSubdomainsZ = NumberOfSubdomainsZ-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsYZDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease Y and Z by 1"
      NumberOfSubdomainsY = NumberOfSubdomainsY-1
      NumberOfSubdomainsZ = NumberOfSubdomainsZ-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsXYZDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains & 
    )) THEN
      IF (DEBUGGING) PRINT *, "best to decrease X, Y and Z by 1"
      NumberOfSubdomainsX = NumberOfSubdomainsX-1
      NumberOfSubdomainsY = NumberOfSubdomainsY-1
      NumberOfSubdomainsZ = NumberOfSubdomainsZ-1
    ELSE
      IF (DEBUGGING) PRINT *, "it does not help to decrease X,Y or Z by 1, start iterative procedure"
      DO WHILE (GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, &
        & NumberOfSubdomainsY, NumberOfSubdomainsZ) > PretendedNumberOfDomains)
        DiffX = NumberOfSubdomainsX - NumberOfSubdomainsXFloat
        DiffY = NumberOfSubdomainsY - NumberOfSubdomainsYFloat
        DiffZ = NumberOfSubdomainsZ - NumberOfSubdomainsZFloat
        
        IF (DiffX >= DiffY .AND. DiffX >= DiffZ) THEN
          IF (NumberOfSubdomainsX /= 1) THEN
            NumberOfSubdomainsX = NumberOfSubdomainsX - 1
          ELSEIF (DiffY >= DiffZ) THEN
            IF (NumberOfSubdomainsY /= 1) THEN
              NumberOfSubdomainsY = NumberOfSubdomainsY - 1
            ELSE
              NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
            ENDIF
          ELSE
            IF (NumberOfSubdomainsZ /= 1) THEN
              NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
            ELSE
              NumberOfSubdomainsY = NumberOfSubdomainsY - 1
            ENDIF
          ENDIF
              
        ELSEIF (DiffY >= DiffZ) THEN    ! DiffY >= DiffZ, DiffY > DiffX
          IF (NumberOfSubdomainsY /= 1) THEN
            NumberOfSubdomainsY = NumberOfSubdomainsY - 1
          ELSE
            IF (DiffX >= DiffZ) THEN
              IF (NumberOfSubdomainsX /= 1) THEN
                NumberOfSubdomainsX = NumberOfSubdomainsX - 1
              ELSE
                NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
              ENDIF
            ELSE
              IF (NumberOfSubdomainsZ /= 1) THEN
                NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
              ELSE
                NumberOfSubdomainsX = NumberOfSubdomainsX - 1
              ENDIF
            ENDIF
          ENDIF
      
        ELSE       ! DiffZ > DiffY, DiffZ >= DiffX
          IF (NumberOfSubdomainsZ /= 1) THEN
            NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
          ELSE
            IF (DiffX >= DiffY) THEN
              IF (NumberOfSubdomainsX /= 1) THEN
                NumberOfSubdomainsX = NumberOfSubdomainsX - 1
              ELSE
                NumberOfSubdomainsY = NumberOfSubdomainsY - 1
              ENDIF
            ELSE
              IF (NumberOfSubdomainsY /= 1) THEN
                NumberOfSubdomainsY = NumberOfSubdomainsY - 1
              ELSE
                NumberOfSubdomainsX = NumberOfSubdomainsX - 1
              ENDIF
            ENDIF
          ENDIF
        ENDIF
          
        IF (DEBUGGING) PRINT *, "Diff: ", DiffX,DiffY,DiffZ, ", nSubdomains:",NumberOfSubdomainsX, NumberOfSubdomainsY, &
          & NumberOfSubdomainsZ, "=", NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ
      ENDDO
    ENDIF  
  ENDIF
  
  
  IF (DEBUGGING) PRINT *, "nSubdomains: ",NumberOfSubdomainsX, NumberOfSubdomainsY, NumberOfSubdomainsZ, "=", &
    &  NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ

  ! determine shape of subdomain, i.e. NumberOfAtomsPerSubdomain
  ! NumberOfAtomsPerSubdomainX * NumberOfSubdomainsX = NumberOfAtomsX  =>  NumberOfAtomsPerSubdomainX = NumberOfAtomsX / NumberOfSubdomainsX

  NumberOfAtomsPerSubdomainX = CEILING(DBLE(NumberOfAtomsX) / NumberOfSubdomainsX)
  NumberOfAtomsPerSubdomainY = CEILING(DBLE(NumberOfAtomsY) / NumberOfSubdomainsY)
  NumberOfAtomsPerSubdomainZ = CEILING(DBLE(NumberOfAtomsZ) / NumberOfSubdomainsZ)
  NumberOfAtomsPerSubdomain = NumberOfAtomsPerSubdomainX*NumberOfAtomsPerSubdomainY*NumberOfAtomsPerSubdomainZ

  ! decrease number of subdomains to exclude now empty subdomains
  nEmptySubdomainsX = FLOOR(DBLE(NumberOfAtomsPerSubdomainX*NumberOfSubdomainsX - NumberOfAtomsX) / NumberOfAtomsPerSubdomainX)
  nEmptySubdomainsY = FLOOR(DBLE(NumberOfAtomsPerSubdomainY*NumberOfSubdomainsY - NumberOfAtomsY) / NumberOfAtomsPerSubdomainY)
  nEmptySubdomainsZ = FLOOR(DBLE(NumberOfAtomsPerSubdomainZ*NumberOfSubdomainsZ - NumberOfAtomsZ) / NumberOfAtomsPerSubdomainZ)

  IF (DEBUGGING) PRINT *, "NumberOfAtomsPerSubdomain: ",NumberOfAtomsPerSubdomainX,NumberOfAtomsPerSubdomainY, &
    & NumberOfAtomsPerSubdomainZ
  
  nEmptySubdomains = NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ &
    & - (NumberOfSubdomainsX-nEmptySubdomainsX)*(NumberOfSubdomainsY-nEmptySubdomainsY)*(NumberOfSubdomainsZ-nEmptySubdomainsZ)

  IF (DEBUGGING) PRINT *, "nEmptySubdomains: ",nEmptySubdomainsX,nEmptySubdomainsY,nEmptySubdomainsZ,&
    & ", total: ",nEmptySubdomains

  NumberOfSubdomainsX = NumberOfSubdomainsX - nEmptySubdomainsX
  NumberOfSubdomainsY = NumberOfSubdomainsY - nEmptySubdomainsY
  NumberOfSubdomainsZ = NumberOfSubdomainsZ - nEmptySubdomainsZ
  nSubdomains = NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ

  nUnusedSubdomains = NumberOfDomains - nSubdomains
  IF (DEBUGGING) PRINT *, "PretendedNumberOfDomains: ", PretendedNumberOfDomains,", nSubdomains:", nSubdomains, &
    & ", nUnusedSubdomains:",nUnusedSubdomains
  
  IF (nUnusedSubdomains /= 0) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT *, "The computed decomposition contains ", nUnusedSubdomains, " subdomains less than given processes. " // &
        & NEW_LINE('A') // " This is intended and no error. Restart program with ", nSubdomains, &
        & " processes instead of ", NumberOfDomains, " or use PretendedNumberOfDomains."
    ENDIF
    
    CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
    CALL MPI_FINALIZE(Err)
    STOP
  ENDIF

  NumberOfElementsInSubdomainX = NumberOfElementsInAtomX*NumberOfAtomsPerSubdomainX
  NumberOfElementsInSubdomainY = NumberOfElementsInAtomY*NumberOfAtomsPerSubdomainY
  NumberOfElementsInSubdomainZ = NumberOfElementsInAtomZ*NumberOfAtomsPerSubdomainZ

  IF (DEBUGGING) THEN
    PRINT *, "NumberOfAtomsPerSubdomain: ",NumberOfAtomsPerSubdomainX, NumberOfAtomsPerSubdomainY, NumberOfAtomsPerSubdomainZ, &
      & "=", NumberOfAtomsPerSubdomain
    PRINT *, "nSubdomains: ",NumberOfSubdomainsX, NumberOfSubdomainsY, NumberOfSubdomainsZ, "=", &
      & NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ, &
      & ", NumberOfAtomsPerSubdomain: ", NumberOfAtomsPerSubdomainX,NumberOfAtomsPerSubdomainY,NumberOfAtomsPerSubdomainZ, "=", &
      & NumberOfAtomsPerSubdomain, ", NumberOfElementsInAtom: ", NumberOfElementsInAtomX, NumberOfElementsInAtomY, &
      & NumberOfElementsInAtomZ, "=", NumberOfElementsInAtomX*NumberOfElementsInAtomY*NumberOfElementsInAtomX, & 
      ", NumberGlobalElements ", NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements,"=",&
      &  NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements

    PRINT *, "subdomain size: ",NumberOfElementsInSubdomainX,"x",&
       & NumberOfElementsInSubdomainY, "x", NumberOfElementsInSubdomainZ
    PRINT *, "total size: ", NumberOfSubdomainsX*NumberOfElementsInSubdomainX, "x",&
       & NumberOfSubdomainsY*NumberOfElementsInSubdomainY, "x", &
       & NumberOfSubdomainsZ*NumberOfElementsInSubdomainZ, & 
       " >= ", NumberGlobalXElements,"x",NumberGlobalYElements,"x",NumberGlobalZElements
  ENDIF

  ! determine last subdomain sizes
  ! NumberOfAtomsLastSubdomainX = NumberOfAtomsPerSubdomainX - (NumberOfAtomsPerSubdomainX*NumberOfSubdomainsX - NumberOfAtomsX)
  NumberOfAtomsLastSubdomainX = NumberOfAtomsPerSubdomainX*(1-NumberOfSubdomainsX) + NumberOfAtomsX
  NumberOfAtomsLastSubdomainY = NumberOfAtomsPerSubdomainY*(1-NumberOfSubdomainsY) + NumberOfAtomsY
  NumberOfAtomsLastSubdomainZ = NumberOfAtomsPerSubdomainZ*(1-NumberOfSubdomainsZ) + NumberOfAtomsZ

  ! NumberOfElementsLastAtomX = NumberOfElementsInAtomX - (NumberOfAtomsX*NumberOfElementsInAtomX - NumberGlobalXElements)
  ! NumberOfElementsLastAtomX = NumberOfElementsInAtomX*(1 - NumberOfAtomsX) + NumberGlobalXElements
  NumberOfElementsLastAtomX = NumberOfElementsInAtomX*(1-NumberOfAtomsX) + NumberGlobalXElements
  NumberOfElementsLastAtomY = NumberOfElementsInAtomY*(1-NumberOfAtomsY) + NumberGlobalYElements
  NumberOfElementsLastAtomZ = NumberOfElementsInAtomZ*(1-NumberOfAtomsZ) + NumberGlobalZElements


  IF (DEBUGGING) PRINT *, "NumberOfAtomsLastSubdomain: ", NumberOfAtomsLastSubdomainX, NumberOfAtomsLastSubdomainY, &
    & NumberOfAtomsLastSubdomainZ, ", NumberOfElementsLastAtom: ",NumberOfElementsLastAtomX, NumberOfElementsLastAtomY, &
    & NumberOfElementsLastAtomZ

  ! for special subdomains (at x+,y+,z+ border) output the number of missing elements
  
  NormalNumberOfElements = &
    NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
    *NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY &
    *NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ

  IF (ComputationalNodeNumber == 0) PRINT *, nSubdomains, " process(es), normal number of elements per process: ", &
    & normalNumberOfElements

  ! x+ face
  actualNumberOfElements = &
    ((NumberOfAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) &
    *NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY&
    *NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ
        
  IF ((NumberOfSubdomainsY-1)*(NumberOfSubdomainsZ-1) > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, (NumberOfSubdomainsY-1)*(NumberOfSubdomainsZ-1), &
      & " process(es) on    x+ boundary have ", normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
        
  ! y+ face
  actualNumberOfElements = &
    NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
    *((NumberOfAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) &
    *NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ
        
  IF ((NumberOfSubdomainsX-1)*(NumberOfSubdomainsZ-1) > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, (NumberOfSubdomainsX-1)*(NumberOfSubdomainsZ-1), &
      & " process(es) on    y+ boundary have ", normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
    
  ! z+ face
  actualNumberOfElements = &
    NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
    *NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY &
    *((NumberOfAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
        
  IF ((NumberOfSubdomainsX-1)*(NumberOfSubdomainsY-1) > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, (NumberOfSubdomainsX-1)*(NumberOfSubdomainsY-1), &
      & " process(es) on    z+ boundary have ", normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
  ! x+ y+ edge
  actualNumberOfElements = &
    ((NumberOfAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) &
    *((NumberOfAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) &
    *NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ
        
  IF (NumberOfSubdomainsZ-1 > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, NumberOfSubdomainsZ-1, " process(es) on x+,y+ boundary have ", &
      & normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
  ! x+ z+ edge
  actualNumberOfElements = &
    ((NumberOfAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) &
    *NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY &
    *((NumberOfAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
        
  IF (NumberOfSubdomainsY-1 > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, NumberOfSubdomainsY-1, " process(es) on x+,z+ boundary have ", &
      & normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
  ! y+ z+ edge
  actualNumberOfElements = &
    NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
    *((NumberOfAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) &
    *((NumberOfAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
        
  IF (NumberOfSubdomainsX-1 > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, NumberOfSubdomainsX-1, " process(es) on y+,z+ boundary have ", &
      & normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
  ! x+ y+ z+ subdomain
  actualNumberOfElements = &
    ((NumberOfAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) &
    *((NumberOfAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) &
    *((NumberOfAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
        
  IF (actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, "          1  process on  x+,y+,z+ boundary has  ", &
      & normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
END SUBROUTINE ComputeSubdomainsWithAtoms

SUBROUTINE CreateDecompositionFiniteElasticity()
  INTEGER(CMISSIntg) :: DomainNo
  INTEGER(CMISSIntg) :: ElementFENo
  INTEGER(CMISSIntg) :: SubdomainZ, SubdomainY, SubdomainX
  INTEGER(CMISSIntg) :: AtomX, AtomY, AtomZ
  INTEGER(CMISSIntg) :: X, Y, Z, XIndex, YIndex, ZIndex
  INTEGER(CMISSIntg) :: NumberOfAtomsCurrentSubdomainX, NumberOfAtomsCurrentSubdomainY, NumberOfAtomsCurrentSubdomainZ
  INTEGER(CMISSIntg) :: NumberOfElementsCurrentAtomX, NumberOfElementsCurrentAtomY, NumberOfElementsCurrentAtomZ
  LOGICAL :: DEBUGGING = .FALSE.
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Create a decomposition for global FE elements
  ! The following global variables are read: NumberOfSubdomainsXYZ, NumberOfAtomsPerSubdomainXYZ, NumberOfAtomsLastSubdomainXYZ, NumberOfElementsInAtomXYZ, NumberOfElementsLastAtomXYZ
  CALL cmfe_Decomposition_Initialise(DecompositionFE,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumberFE,MeshFE,DecompositionFE,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(DecompositionFE,.TRUE.,Err)

  IF (DEBUGGING) DEBUGGING = (ComputationalNodeNumber == 2)
  
  IF(NumberOfDomains > 1) THEN
    CALL cmfe_Decomposition_TypeSet(DecompositionFE, CMFE_DECOMPOSITION_USER_DEFINED_TYPE, Err)
    
    ! assign domain values
        
    ! iterate over subdomains
    DO SubdomainZ = 1,NumberOfSubdomainsZ
      DO SubdomainY = 1,NumberOfSubdomainsY
        DO SubdomainX = 1,NumberOfSubdomainsX
        
          DomainNo = (SubdomainZ-1)*NumberOfSubdomainsY*NumberOfSubdomainsX + (SubdomainY-1)*NumberOfSubdomainsX + SubdomainX-1
          
          NumberOfAtomsCurrentSubdomainX = NumberOfAtomsPerSubdomainX
          IF (SubdomainX == NumberOfSubdomainsX) THEN  ! last subdomain
            NumberOfAtomsCurrentSubdomainX = NumberOfAtomsLastSubdomainX
          ENDIF
            
          NumberOfAtomsCurrentSubdomainY = NumberOfAtomsPerSubdomainY          
          IF (SubdomainY == NumberOfSubdomainsY) THEN  ! last subdomain
            NumberOfAtomsCurrentSubdomainY = NumberOfAtomsLastSubdomainY
          ENDIF
            
          NumberOfAtomsCurrentSubdomainZ = NumberOfAtomsPerSubdomainZ
          IF (SubdomainZ == NumberOfSubdomainsZ) THEN  ! last subdomain
            NumberOfAtomsCurrentSubdomainZ = NumberOfAtomsLastSubdomainZ
          ENDIF
            
          IF (DEBUGGING) THEN
            PRINT *, "subdomain (",SubdomainX,",",SubdomainY,",",SubdomainZ,") has ", NumberOfAtomsCurrentSubdomainX,",", &
              & NumberOfAtomsCurrentSubdomainY,",",NumberOfAtomsCurrentSubdomainZ," atoms"
          ENDIF
          
          DO AtomZ = 1,NumberOfAtomsCurrentSubdomainZ
            DO AtomY = 1,NumberOfAtomsCurrentSubdomainY
              DO AtomX = 1,NumberOfAtomsCurrentSubdomainX
                  
                NumberOfElementsCurrentAtomX = NumberOfElementsInAtomX
                IF (SubdomainX == NumberOfSubdomainsX .AND. AtomX == NumberOfAtomsCurrentSubdomainX) THEN  ! last atom
                  NumberOfElementsCurrentAtomX = NumberOfElementsLastAtomX
                ENDIF
                  
                NumberOfElementsCurrentAtomY = NumberOfElementsInAtomY
                IF (SubdomainY == NumberOfSubdomainsY .AND. AtomY == NumberOfAtomsCurrentSubdomainY) THEN  ! last atom
                  NumberOfElementsCurrentAtomY = NumberOfElementsLastAtomY
                ENDIF
                  
                NumberOfElementsCurrentAtomZ = NumberOfElementsInAtomZ
                IF (SubdomainZ == NumberOfSubdomainsZ .AND. AtomZ == NumberOfAtomsCurrentSubdomainZ) THEN  ! last atom
                  NumberOfElementsCurrentAtomZ = NumberOfElementsLastAtomZ
                ENDIF
                  
                IF (DEBUGGING) PRINT *, "  atom (",AtomX,",",AtomY,",",AtomZ,") has ",NumberOfElementsCurrentAtomX,",",&
                  & NumberOfElementsCurrentAtomY,",",NumberOfElementsCurrentAtomZ," elements"
                  
                DO Z = 1,NumberOfElementsCurrentAtomZ
                  DO Y = 1,NumberOfElementsCurrentAtomY
                    DO X = 1,NumberOfElementsCurrentAtomX
                      
                      ZIndex = (SubdomainZ-1)*NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ &
                        & + (AtomZ-1)*NumberOfElementsInAtomZ + Z
                      YIndex = (SubdomainY-1)*NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY &
                        & + (AtomY-1)*NumberOfElementsInAtomY + Y
                      XIndex = (SubdomainX-1)*NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
                        & + (AtomX-1)*NumberOfElementsInAtomX + X
                      
                      ElementFENo = (ZIndex-1)*NumberGlobalXElements*NumberGlobalYElements + (YIndex-1)*NumberGlobalXElements &
                        & + XIndex
                      
                      IF (DEBUGGING) PRINT *,"      (x,y,z)=(",XIndex,",",YIndex,",",ZIndex,") i=",ElementFENo," domain ",DomainNo
     
                      !                                        DECOMPOSITION,   GLOBAL_ELEMENT_NUMBER, DOMAIN_NUMBER
                      CALL cmfe_Decomposition_ElementDomainSet(DecompositionFE, ElementFENo,           DomainNo, Err)
                    
                    ENDDO ! X
                  ENDDO   ! Y
                ENDDO     ! Z
              ENDDO ! AtomX
            ENDDO   ! AtomY
          ENDDO     ! AtomZ
        ENDDO ! SubdomainX
      ENDDO   ! SubdomainY
    ENDDO     ! SubdomainZ
                    
        
    IF (DEBUGGING) THEN
      CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
      CALL MPI_ABORT(MPI_COMM_WORLD, Err, Err)
      PRINT *, "stop in FortranExample.f90:1733"
      STOP
    ENDIF
    
  ELSE
    ! single process
    CALL cmfe_Decomposition_TypeSet(DecompositionFE,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  ENDIF

  ! finish decomposition
  CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionFE,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(DecompositionFE,Err)
END SUBROUTINE CreateDecompositionFiniteElasticity

SUBROUTINE CreateDecompositionMonodomain()

  INTEGER(CMISSIntg) :: DomainNo
  INTEGER(CMISSIntg) :: ElementMGlobalNumber, ElementMInFibreNo
  INTEGER(CMISSIntg) :: FibreNo
  INTEGER(CMISSIntg) :: FEElementGlobalNumber
  INTEGER(CMISSIntg) :: FEElementXIdx, FEElementYIdx, FEElementZIdx
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Create a decompositions for monodomain elements on fibres
  !PRINT "(I3.3,A)", ComputationalNodeNumber, ": Create decomposition for monodomain elements"
  
  CALL cmfe_Decomposition_Initialise(DecompositionM,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumberM,MeshM,DecompositionM,Err)
         
  IF(NumberOfDomains > 1) THEN
    ! set decomposition to be of user defined type
    CALL cmfe_Decomposition_TypeSet(DecompositionM, CMFE_DECOMPOSITION_USER_DEFINED_TYPE, Err)
    
    ! assign the same domains to bioelectric nodes as the containing FE elements
    ! loop over elementsM
    ElementMGlobalNumber = 1
    DO FibreNo = 1, NumberOfFibres
      !PRINT *, "Fibre ",FibreNo
    
      DO ElementMInFibreNo = 1,NumberOfElementsMPerFibre
        
        ! compute global ElementFE that contains the current elementM
        FEElementZIdx = INT(INT((FibreNo-1) / NumberGlobalYFibres) / NumberOfNodesInXi3)
        FEElementYIdx = INT(MOD((FibreNo-1), NumberGlobalYFibres)  / NumberOfNodesInXi2)
        FEElementXIdx = (ElementMInFibreNo-1) / NumberOfNodesInXi1
        
        FEElementGlobalNumber = FEElementZIdx * NumberGlobalYElements * NumberGlobalXElements &
         & + FEElementYIdx * NumberGlobalXElements + FEElementXIdx + 1
        
        !PRINT "(I1.1,2(A,I2.2),4(A,I3.3))", ComputationalNodeNumber,": Fibre", FibreNo, ", ElementM ", ElementMGlobalNumber, &
        !  & ", Element (",FEElementXIdx,",",FEElementYIdx,",",FEElementZIdx, &
        !  & ")=", FEElementGlobalNumber
      
        ! get the domain of the global ElementFE
        !                                        DECOMPOSITION,   USER_ELEMENT_NUMBER,   DOMAIN_NUMBER
        CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE, FEElementGlobalNumber, DomainNo, Err)

        ! set the domain of the ElementM to the same domain as the containing global element
        !                                        DECOMPOSITION,   GLOBAL_ELEMENT_NUMBER, DOMAIN_NUMBER
        CALL cmfe_Decomposition_ElementDomainSet(DecompositionM,  ElementMGlobalNumber,  DomainNo, Err)
        !PRINT "(I1.1,A,I5.5,A,I2,A,I2)", ComputationalNodeNumber, ": fibre ", FibreNo, ", 1D el. no. ", ElementMGlobalNumber, &
        !  & " to domain no. ", DomainNo
          
        ElementMGlobalNumber = ElementMGlobalNumber + 1
      ENDDO
    ENDDO
  ELSE
    ! single process
    CALL cmfe_Decomposition_TypeSet(DecompositionM,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  ENDIF

  ! finish decomposition
  CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionM,NumberOfDomains,Err)
    
  ! In this call it is checked whether all processes have at least one node
  CALL cmfe_Decomposition_CreateFinish(DecompositionM,Err)

  IF (.FALSE.) THEN
    PRINT*, "Print elements mapping in FortranExample.f90:1565"
    PRINT*, "--------------- ElementsMapping --------------------------"
    CALL cmfe_PrintElementsMapping(DecompositionM,Err)
  ENDIF
  
  !PRINT*, "Print nodes mapping in FortranExample.f90:1549"
  !CALL cmfe_PrintNodesMapping(DecompositionM,Err)
  
END SUBROUTINE CreateDecompositionMonodomain


END MODULE DECOMPOSITION