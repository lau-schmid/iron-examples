
MODULE REGION_MESH

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE OPENCMISS_VARIABLES
  
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

SUBROUTINE CreateRegionMesh()
  INTEGER(CMISSIntg) :: ElementInFibreNo
  INTEGER(CMISSIntg) :: FibreNo
  INTEGER(CMISSIntg) :: NodeUserNumber
  INTEGER(CMISSIntg) :: ElementMGlobalNumber
  !-------------------------------------------------------------------------------------------------------------------------------
  
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  
  IF (ComputationalNodeNumber == 0) THEN
    IF (NumberOfComputationalNodes == 1) THEN
      PRINT "(A)", "Running with 1 process."
    ELSE 
      PRINT "(A,I6,A)", "Running with ",NumberOfComputationalNodes," processes."
    ENDIF
  ENDIF
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystemFE,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumberFE,CoordinateSystemFE,Err)
  !Set the coordinate system to be 3D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemFE,3,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemFE,Err)


  ! CREATE A 1D COORDINATE SYSTEM
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystemM,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumberM,CoordinateSystemM,Err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemM,1,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL cmfe_Region_Initialise(RegionFE,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumberFE,WorldRegion,RegionFE,Err)
  CALL cmfe_Region_CoordinateSystemSet(RegionFE,CoordinateSystemFE,Err)
  CALL cmfe_Region_LabelSet(RegionFE,"Region3D",Err)
  CALL cmfe_Region_CreateFinish(RegionFE,Err)


  ! CREATE A SECOND REGION
  CALL cmfe_Region_Initialise(RegionM,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumberM,WorldRegion,RegionM,Err)
  CALL cmfe_Region_CoordinateSystemSet(RegionM,CoordinateSystemM,Err)
  CALL cmfe_Region_LabelSet(RegionM,"Region1D",Err)
  CALL cmfe_Region_CreateFinish(RegionM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the bases
  !Define basis functions - tri-Quadratic Lagrange
  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_TypeSet(QuadraticBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  !Define basis functions - tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_TypeSet(LinearBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  ! CREATE A SECOND LINEAR BASIS FOR THE 1D GRID
  !Define basis functions - tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(LinearBasisM,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumberM,LinearBasisM,Err)
  CALL cmfe_Basis_TypeSet(LinearBasisM,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasisM,1,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearBasisM,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasisM,[2],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasisM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a mesh with 8 three-dimensional elements
  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMeshFE,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,RegionFE,GeneratedMeshFE,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMeshFE,[QuadraticBasis,LinearBasis],Err)
  !Define the mesh on the region

  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMeshFE,[PhysicalLength,PhysicalWidth,PhysicalHeight],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMeshFE,[NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements],Err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(MeshFE,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMeshFE,MeshUserNumberFE,MeshFE,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! create the monodomain mesh
  !Create a mesh in the region
    
  CALL cmfe_Mesh_Initialise(MeshM,Err)
  
  !                          user_number,     region,  n_dim, mesh(out)
  CALL cmfe_Mesh_CreateStart(MeshUserNumberM, RegionM, 1,     MeshM, Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(MeshM,1,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(MeshM,NumberOfElementsM,Err)

  CALL cmfe_MeshElements_Initialise(ElementsM,Err)   
  !                                  mesh,  meshComponentNumber,           basis,        meshElements
  CALL cmfe_MeshElements_CreateStart(MeshM, MonodomainMeshComponentNumber, LinearBasisM, ElementsM, Err)
  
  ! Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(RegionM, NumberOfNodesM, Nodes, Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  
  
  ! Set adjacent nodes for each element
  NodeUserNumber = 1
  ElementMGlobalNumber = 1
  DO FibreNo = 1, NumberOfFibres
    !PRINT *, "Fibre ", FibreNo
    DO ElementInFibreNo = 1,NumberOfElementsMPerFibre
      
      !                               meshElements, globalElementNumber,  elementUserNodes (user numbers)
      CALL cmfe_MeshElements_NodesSet(ElementsM,    ElementMGlobalNumber, [NodeUserNumber,NodeUserNumber+1], Err)

      !PRINT "(A,I3.3,A,I3.3,A,I5.5,A,I7.7,A,I7.7)", &
      !  & "a ", ComputationalNodeNumber, ": fibre ", FibreNo, " 1D el. no. ", ElementInFibreNo, " has nodes ", &
      !  & NodeUserNumber, ", ", NodeUserNumber+1

      ! If at the end of a fibre line, increment node to starting node of next fibre line
      IF (MOD(ElementInFibreNo, NumberOfElementsMPerFibreLine) == 0) THEN
        NodeUserNumber = NodeUserNumber+1
      ENDIF
     
      NodeUserNumber = NodeUserNumber+1
      ElementMGlobalNumber = ElementMGlobalNumber+1
    ENDDO
  ENDDO
  !write(*,*) "Finished setting up 1D elements"

  CALL cmfe_MeshElements_CreateFinish(ElementsM, Err)
  CALL cmfe_Mesh_CreateFinish(MeshM, Err)
                 
  !PRINT*, "internal elements ", ELEMENTS_MAPPING%INTERNAL_START,"to",ELEMENTS_MAPPING%INTERNAL_FINISH
  !PRINT*, "element parameters:"
  
  !PRINT*, "index        element_no      interpolation_parameters"
  !DO element_idx=ELEMENTS_MAPPING%INTERNAL_START, ELEMENTS_MAPPING%INTERNAL_FINISH
  
  !  ne = ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
    
    ! get interpolation parameters of element
    ! version which is used with real preallocated variable names:
  !  CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne,&
  !    & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

  !  INTERPOLATION_PARAMETERS=> &
  !    & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%INTERPOLATION_PARAMETERS                
    
    ! direct version
  !  CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
    
  !  DO component_idx = 1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
  !    PRINT*, element_idx, ne, INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx)
  !  ENDDO
  !ENDDO
  


END SUBROUTINE CreateRegionMesh

END MODULE REGION_MESH