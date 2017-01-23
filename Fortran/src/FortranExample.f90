!> \file
!> \author Adam Reeve
!> \brief This is an example program to solve a finite elasticity equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>

!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example FiniteElasticity/LargeUniAxialExtension/src/LargeUniAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-gnu'>Linux GNU Build</a>
!<


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING MODULE INCLUDES THE DATA STRUCTURES AND  "CONTAIN SUBROUTINES AND FUNCTIONS" USED IN PARSING ALGORITHM  !!!!!!!!!!!!!!!!1
include "parsing_module.f90"

!> Main program
PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE




  USE OpenCMISS
  USE OpenCMISS_Iron
  USE parsing
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

 ! REAL(CMISSRP), PARAMETER :: Gravity(3)=[0.0_CMISSRP,0.0_CMISSRP,-9.8_CMISSRP] !in m s^-2
 ! REAL(CMISSRP), PARAMETER :: Density=   9.0E-4_CMISSRP
  REAL(CMISSRP) :: HEIGHT
  REAL(CMISSRP) :: WIDTH
  REAL(CMISSRP) :: LENGTH
  INTEGER(CMISSIntg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=300) :: CommandArgument,inputFile

  
!  INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg) 		 :: InterpolationType
  INTEGER(CMISSIntg) 		 :: PressureInterpolationType
  LOGICAL			 :: UsePressureBasis 
  INTEGER(CMISSIntg)             :: ScalingType
!  LOGICAL, PARAMETER :: UsePressureBasis=.FALSE.
  INTEGER(CMISSIntg)             :: NumberOfGaussXi

!!!!!!!!!!!!!!!!!!!!		COUTER SIZE		!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! declare 
  INTEGER(CMISSIntg)    :: num_of_PressureBasis
  INTEGER(CMISSIntg)    :: num_of_Basis
  INTEGER(CMISSIntg)    :: num_of_BoundaryCondition
  INTEGER(CMISSIntg)    :: num_of_WorldCoordinateSystem
  INTEGER(CMISSIntg)    :: num_of_CoordinateSystem
  INTEGER(CMISSIntg)    :: num_of_Mesh
  INTEGER(CMISSIntg)    :: num_of_Decomposition
  INTEGER(CMISSIntg)    :: num_of_Equation
  INTEGER(CMISSIntg)    :: num_of_EquationsSet
  INTEGER(CMISSIntg)    :: num_of_GeometricField
  INTEGER(CMISSIntg)    :: num_of_FiberField
  INTEGER(CMISSIntg)    :: num_of_MaterialField
  INTEGER(CMISSIntg)    :: num_of_DependentField
  INTEGER(CMISSIntg)    :: num_of_EquationSetField
  INTEGER(CMISSIntg)    :: num_of_Field
  INTEGER(CMISSIntg)    :: num_of_Problem
  INTEGER(CMISSIntg)    :: num_of_Region
  INTEGER(CMISSIntg)    :: num_of_WorldRegion
  INTEGER(CMISSIntg)    :: num_of_Solver
  INTEGER(CMISSIntg)    :: num_of_LinearSolver
  INTEGER(CMISSIntg)    :: num_of_SolverEquations
  INTEGER(CMISSIntg)    :: num_of_ControlLoop
  INTEGER(CMISSIntg)    :: num_of_GeneratedMesh
  INTEGER(CMISSIntg)    :: j ,component_idx,num_of_dirichelet,num_of_traction_neumann,num_of_pressure_neumann
  character(len = 100)  :: constraint 

 
  
 
!!!!!!!!!!!!!!!!!!!!!!  	END COUNTER SIZE 		!!!!!!!!!!!!!!!!!!!

  INTEGER(CMISSIntg), allocatable :: CoordinateSystemUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: RegionUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: BasisUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: PressureBasisUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: GeneratedMeshUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: MeshUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: DecompositionUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: FieldGeometryUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: FieldFibreUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: FieldMaterialUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: FieldDependentUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: EquationSetUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: EquationsSetFieldUserNumber(:)
  INTEGER(CMISSIntg), allocatable :: ProblemUserNumber(:)
  INTEGER(CMISSIntg) ,parameter   :: FieldSourceUserNumber = 20
  !Program types

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err
  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: SurfaceNodes(:),LeftSurfaceNodes(:),FrontSurfaceNodes(:)

  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi

! INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 3 !nearly incompressible
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 4 !fully incompressible

!!!!!!!!!!!!!!!!!!!!!!!!!!! 		          MY CONTRIBUTION 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  include "DerivedTypes.f90" !!DECLARING THE DERIVED TYPES HERE      

  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif



  NumberOfArguments = COMMAND_ARGUMENT_COUNT()
  
  IF(NumberOfArguments .NE. 1) THEN
    CALL HANDLE_ERROR("Please provide only the input file")
  ELSE
    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgStatus)
    InputFile = trim(CommandArgument)
  ENDIF
 fileplace = InputFile
 print *, fileplace

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 		INITIALIZE THE DERIVED DATA STRUCTURES 			 !!!!!!!!!!!!!!!!!!!!!
 include "AllocatingDerivedDataStructures.f90"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              ALLOCATE  "USERNUMBER"  DATA STRUCTURE SIZE 		 !!!!!!!!!!!!!!!!!!!!!!
 include "UserNumberDataStructures.f90"

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        INCLUDES I/0 STATEMENT THAT PARSE THROUGH THE INPUT FILE        !!!!!!!!!!!!!!!!!!!!!!
 include "parsing_algorithm.f90"





#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif
  !Intialise cmiss
  CALL cmfe_Initialise(all_WorldCoordinateSystem%WorldCoordinateSystem(1),all_WorldRegion%WorldRegion(1),Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
   !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  NumberOfDomains=NumberOfComputationalNodes
  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

DO  i = 1,num_of_CoordinateSystem 			!!!!!!!!!!!!! introducing a LOOP for assigning  coordinate systems to thier respective regions !!!!!!!!!!
 
  !Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(all_CoordinateSystem%CoordinateSystem(i),Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber(i),all_CoordinateSystem%CoordinateSystem(i),Err)
  Call cmfe_CoordinateSystem_TypeSet(CoordinateSystemUserNumber(i), &
       & match_coordinate_system(CoordinateSystem_arguments(2,i)), err)
  CALL cmfe_CoordinateSystem_CreateFinish(all_CoordinateSystem%CoordinateSystem(i),Err)
  



  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(all_Region%Region(i),Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber(i),all_WorldRegion%WorldRegion(i),all_Region%Region(i),Err)
  CALL cmfe_Region_LabelSet(all_Region%Region(i),"Region",Err)
  CALL cmfe_Region_CoordinateSystemSet(all_Region%Region(i),(all_CoordinateSystem%CoordinateSystem(i)),Err)
  CALL cmfe_Region_CreateFinish(all_Region%Region(i),Err)

END DO 


Basis: DO i = 1, num_of_basis  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Assigning basis to thier respective edges       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  InterpolationType= match_basis(Basis_arguments(4,i))
  PressureInterpolationType= match_basis(Basis_arguments(3,i))      !!!!! DATA READ FROM THE INPUT FILE !!!!!!!
  UsePressureBasis  = (str2int(Basis_arguments(2,i)) == 4) 
  NumberGlobalXElements=str2int(Mesh_arg4(1,i))
  NumberGlobalYElements=str2int(Mesh_arg4(2,i))
  NumberGlobalZElements=str2int(Mesh_arg4(3,i))
  NumberOfGaussXi = str2int(Basis_arguments(5,i))

  !Define geometric basis
  CALL cmfe_Basis_Initialise(all_Basis%Basis(i),Err)
  CALL cmfe_Basis_CreateStart(1,all_Basis%Basis(i),Err)
  SELECT CASE(InterpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(all_Basis%Basis(i),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(all_Basis%Basis(i),CMFE_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  IF(NumberGlobalZElements==0) THEN
    CALL cmfe_Basis_NumberOfXiSet(all_Basis%Basis(i),2,Err)
    CALL cmfe_Basis_InterpolationXiSet(all_Basis%Basis(i),[InterpolationType,InterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis%Basis(i),[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    CALL cmfe_Basis_NumberOfXiSet(all_Basis%Basis(i),3,Err)
    CALL cmfe_Basis_InterpolationXiSet(all_Basis%Basis(i),[InterpolationType,InterpolationType,InterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis%Basis(i),[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ENDIF
  CALL cmfe_Basis_CreateFinish(all_Basis%Basis(i),Err)




  !Define pressure basis
  IF(UsePressureBasis) THEN
    CALL cmfe_Basis_Initialise(all_PressureBasis%PressureBasis(i),Err)
    CALL cmfe_Basis_CreateStart(PressureBasisUserNumber(i),all_PressureBasis%PressureBasis(i),Err)
    SELECT CASE(PressureInterpolationType)
    CASE(1,2,3,4)
      CALL cmfe_Basis_TypeSet(all_PressureBasis%PressureBasis(i),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CASE(7,8,9)
      CALL cmfe_Basis_TypeSet(all_PressureBasis%PressureBasis(i),CMFE_BASIS_SIMPLEX_TYPE,Err)
    END SELECT
    IF(NumberGlobalZElements==0) THEN
      CALL cmfe_Basis_NumberOfXiSet(all_PressureBasis%PressureBasis(i),2,Err)
      CALL cmfe_Basis_InterpolationXiSet(all_PressureBasis%PressureBasis(i),[PressureInterpolationType,PressureInterpolationType] &
                                        & ,Err)
      IF(NumberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_PressureBasis%PressureBasis(i),[NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ELSE
      CALL cmfe_Basis_NumberOfXiSet(all_PressureBasis%PressureBasis(i),3,Err)
      CALL cmfe_Basis_InterpolationXiSet(all_PressureBasis%PressureBasis(i), &
        & [PressureInterpolationType,PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_PressureBasis%PressureBasis(i),[NumberOfGaussXi,NumberOfGaussXi, & 
                                                     & NumberOfGaussXi],Err)
      ENDIF
    ENDIF
    CALL cmfe_Basis_CreateFinish(all_PressureBasis%PressureBasis(i),Err)
  ENDIF
end do Basis


Mesh: do i = 1,num_of_Mesh
   HEIGHT=str2int(Mesh_arg3(1,i))
   WIDTH=str2int(Mesh_arg3(2,i))
   LENGTH=str2int(Mesh_arg3(3,i))

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(all_GeneratedMesh%GeneratedMesh(i),Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber(i),all_Region%Region(i),all_GeneratedMesh%GeneratedMesh(i),Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(all_GeneratedMesh%GeneratedMesh(i),CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  IF(UsePressureBasis) THEN
    CALL cmfe_GeneratedMesh_BasisSet(all_GeneratedMesh%GeneratedMesh(i),[all_Basis%Basis(i),all_PressureBasis%PressureBasis(i)],Err)
  ELSE
    CALL cmfe_GeneratedMesh_BasisSet(all_GeneratedMesh%GeneratedMesh(i),[all_Basis%Basis(i)],Err)
  ENDIF
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh%GeneratedMesh(i),[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh%GeneratedMesh(i), & 
                                                & [NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh%GeneratedMesh(i),[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh%GeneratedMesh(i), &
                                    & [NumberGlobalXElements,NumberGlobalYElements, &
                                    & NumberGlobalZElements],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(all_Mesh%Mesh(i),Err)
  CALL cmfe_GeneratedMesh_CreateFinish(all_GeneratedMesh%GeneratedMesh(i),MeshUserNumber(i),all_Mesh%Mesh(i),Err)

END do Mesh

Decomposition: do i = 1,num_of_Decomposition
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(all_Decomposition%Decomposition(i),Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber(i),all_Mesh%Mesh(i),all_Decomposition%Decomposition(i),Err)
  CALL cmfe_Decomposition_TypeSet(all_Decomposition%Decomposition(i),CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(all_Decomposition%Decomposition(i),NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(all_Decomposition%Decomposition(i),Err)
end do Decomposition

GeometricField: DO i = 1, num_of_GeometricField

  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(all_GeometricField%GeometricField(i),Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber(i),all_Region%Region(i),all_GeometricField%GeometricField(i),Err)
  CALL cmfe_Field_MeshDecompositionSet(all_GeometricField%GeometricField(i),all_Decomposition%Decomposition(i),Err)
  CALL cmfe_Field_VariableLabelSet(all_GeometricField%GeometricField(i),CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ScalingTypeSet(all_GeometricField%GeometricField(i),CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_Field_CreateFinish(all_GeometricField%GeometricField(i),Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(all_GeneratedMesh%GeneratedMesh(i),all_GeometricField%GeometricField(i),Err)
end do GeometricField  !!!!! END LOOP FOR GEOMETRIC FIELD

FiberField: DO i = 1, num_of_FiberField
  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(all_FibreField%FibreField(i),Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber(i),all_Region%Region(i),all_FibreField%FibreField(i),Err)
  CALL cmfe_Field_TypeSet(all_FibreField%FibreField(i),CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(all_FibreField%FibreField(i),all_Decomposition%Decomposition(i),Err)
  CALL cmfe_Field_GeometricFieldSet(all_FibreField%FibreField(i),all_GeometricField%GeometricField(i),Err)
  CALL cmfe_Field_VariableLabelSet(all_FibreField%FibreField(i),CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_ScalingTypeSet(all_FibreField%FibreField(i),CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_Field_CreateFinish(all_FibreField%FibreField(i),Err)
    !Set the fibre field
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    					    & 1,str2real(FiberField_arg2(1,i)),Err)
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    					    & 2,str2real(FiberField_arg2(2,i)),Err)
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    					    & 3,str2real(FiberField_arg2(3,i)),Err)
end do FiberField  !!!!! END LOOP FOR Fiber FIELD

  !Create the dependent field
DependentField: DO i = 1, num_of_DependentField
  CALL cmfe_Field_Initialise(all_DependentField%DependentField(i),Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber(i),all_Region%Region(i),all_DependentField%DependentField(i),Err)

!  DO component_idx=1,3
!    CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
!    CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(i),CMFE_FIELD_DELUDELN_VARIABLE_TYPE,& 
!                                              & component_idx,1,Err)
!  ENDDO

  CALL cmfe_Field_TypeSet(all_DependentField%DependentField(i),CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(all_DependentField%DependentField(1),all_Decomposition%Decomposition(i),Err)
  CALL cmfe_Field_GeometricFieldSet(all_DependentField%DependentField(i),all_GeometricField%GeometricField(i),Err)
  CALL cmfe_Field_DependentTypeSet(all_DependentField%DependentField(i),CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(all_DependentField%DependentField(i),2,Err)
  CALL cmfe_Field_VariableTypesSet(all_DependentField%DependentField(i), & 
                                   & [CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_VariableLabelSet(all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_NumberOfComponentsSet(all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE, & 
       & str2int(DependentField_arguments(3,1)),Err)
  CALL cmfe_Field_NumberOfComponentsSet(all_DependentField%DependentField(i), & 
       & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,str2int(DependentField_arguments(3,1)),Err)
  IF(UsePressureBasis) THEN
    !Set the pressure to be nodally based and use the second mesh component if required
    CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(i),&
    & CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(i),CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(i),CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  ELSE
    CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE,4, &
                                            & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
                                             & Err)
    CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(i),CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
                                              & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  END IF
  CALL cmfe_Field_CreateFinish(all_DependentField%DependentField(i),Err)


 !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  if  (trim(DependentField_arguments(4,i)) == "UNDEFORMED") then
   
  	CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(i), & 
       	match_dependent_field(DependentField_arguments(2,i)),CMFE_FIELD_VALUES_SET_TYPE, &
 	& 1,all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)

  	CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(i), & 
                                match_dependent_field(DependentField_arguments(2,i)),CMFE_FIELD_VALUES_SET_TYPE, &
 	& 2,all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)

  	if (NumberGlobalZElements .NE. 0) then
  	CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(i), & 
                                match_dependent_field(DependentField_arguments(2,i)),CMFE_FIELD_VALUES_SET_TYPE, &
 	     & 3,all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
        end if



  else 
         CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(i), & 
              &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
              & str2real(DependentField_arguments(5,i)),Err)

         CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(i), & 
              &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
              & str2real(DependentField_arguments(5,i)),Err)

  	if (NumberGlobalZElements .NE. 0) then
         	CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(i), & 
              		&  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
              		& str2real(DependentField_arguments(5,i)),Err)
        end if
   
  end if 
 
  IF(UsePressureBasis) THEN
  	CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(i), & 
                      &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
                      & str2real(DependentField_arguments(5,i)),Err)
  end if 

  CALL cmfe_Field_ParameterSetUpdateStart(all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE & 
                                          & ,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(all_DependentField%DependentField(i),CMFE_FIELD_U_VARIABLE_TYPE, & 
                                          & CMFE_FIELD_VALUES_SET_TYPE,Err)

end do DependentField


  !Create the material field
MaterialField: do i = 1,num_of_MaterialField
  CALL cmfe_Field_Initialise(all_MaterialField%MaterialField(i),Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber(i),all_Region%Region(i),all_MaterialField%MaterialField(i),Err)
  CALL cmfe_Field_TypeSet(all_MaterialField%MaterialField(i),CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(all_MaterialField%MaterialField(i),all_Decomposition%Decomposition(i),Err)
  CALL cmfe_Field_GeometricFieldSet(all_MaterialField%MaterialField(i),all_GeometricField%GeometricField(i),Err)
  CALL cmfe_Field_NumberOfVariablesSet(all_MaterialField%MaterialField(i),1,Err)
  CALL cmfe_Field_VariableLabelSet(all_MaterialField%MaterialField(i),CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
!  CALL cmfe_Field_VariableLabelSet(all_MaterialField%MaterialField(i),CMFE_FIELD_V_VARIABLE_TYPE,"Density",Err)
  CALL cmfe_Field_NumberOfComponentsSet(all_MaterialField%MaterialField(i),CMFE_FIELD_U_VARIABLE_TYPE, & 
                                        & material_parameters(EquationsSet_arg4(1,i)),Err)
  CALL cmfe_Field_CreateFinish(all_MaterialField%MaterialField(i),Err)

  do j = 1,material_parameters(EquationsSet_arg4(1,i))

   CALL  cmfe_Field_ComponentValuesInitialise(all_MaterialField%MaterialField(i), & 
                                              &CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,& 
                                              & j,str2real(MaterialField_arg2(j,i)),Err)
  end do 
!  CALL cmfe_Field_ComponentValuesInitialise(all_MaterialField%MaterialField(i),CMFE_FIELD_V_VARIABLE_TYPE, & 
!                                            & CMFE_FIELD_VALUES_SET_TYPE,1,1,Err)


end do MaterialField

  !Create the equations_set


EquationsSet: DO i = 1,Num_of_EquationsSet            !!!! loop for eqaution set

  CALL cmfe_Field_Initialise(all_EquationsSetField%EquationsSetField(i),Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber(i),all_Region%Region(i),all_FibreField%FibreField(i),& 
       & [match_equations_set(EquationsSet_arg2(1,i)), &
       & match_equations_set(EquationsSet_arg3(1,i)),match_equations_set(EquationsSet_arg4(1,i))], &
       & EquationsSetFieldUserNumber(i),all_EquationsSetField%EquationsSetField(i),all_EquationsSet%EquationsSet(i),Err)
  CALL cmfe_EquationsSet_CreateFinish(all_EquationsSet%EquationsSet(i),Err)

  !Create the equations set dependent field
  CALL cmfe_EquationsSet_DependentCreateStart(all_EquationsSet%EquationsSet(i),FieldDependentUserNumber(i), & 
                                              & all_DependentField%DependentField(i),Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(all_EquationsSet%EquationsSet(i),Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(all_EquationsSet%EquationsSet(i),FieldMaterialUserNumber(i) & 
                                              & ,all_MaterialField%MaterialField(i),Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(all_EquationsSet%EquationsSet(i),Err)


  !Create the equations set equations
  CALL cmfe_Equations_Initialise(all_Equations%Equations(i),Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(all_EquationsSet%EquationsSet(i),all_Equations%Equations(i),Err)
  CALL cmfe_Equations_SparsityTypeSet(all_Equations%Equations(i),output_type(Output_arguments(1,1)),Err)
  CALL cmfe_Equations_OutputTypeSet(all_Equations%Equations(i),output_type(Output_arguments(2,1)),Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(all_EquationsSet%EquationsSet(i),Err)

  !Create the source field with the gravity vector
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(all_EquationsSet%EquationsSet(i),FieldSourceUserNumber,SourceField,Err)
  CALL cmfe_Field_ScalingTypeSet(SourceField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_EquationsSet_SourceCreateFinish(all_EquationsSet%EquationsSet(i),Err)
  DO component_idx=1,3
    	CALL cmfe_Field_ComponentValuesInitialise(SourceField, & 
                       & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
                       & component_idx,Gravity(component_idx),Err)

  ENDDO

END DO EquationsSet 


  

DO i = 1, num_of_Problem

  !Define the problem
  CALL cmfe_Problem_Initialise(all_Problem%Problem(i),Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber(i),[match_problem(Problem_arg2(1,i)),match_problem(Problem_arg3(1,i)), &
    & match_problem(Problem_arg4(1,i))],all_Problem%Problem(i),Err)
  CALL cmfe_Problem_CreateFinish(all_Problem%Problem(i),Err)

END DO

do i = 1,num_of_ControlLoop 
  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(all_Problem%Problem(i),Err)
  CALL cmfe_ControlLoop_Initialise(all_ControlLoop%ControlLoop(i),Err)
  CALL cmfe_Problem_ControlLoopGet(all_Problem%Problem(i),CMFE_CONTROL_LOOP_NODE,all_ControlLoop%ControlLoop(i),Err)
  CALL cmfe_ControlLoop_TypeSet(all_ControlLoop%ControlLoop(i),control_loop_def(ControlLoop_arguments(2,i)),Err)
  CALL cmfe_ControlLoop_LoadOutputSet(all_ControlLoop%ControlLoop(i),1,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(all_ControlLoop%ControlLoop(i),str2int(ControlLoop_arguments(3,i)),Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(all_Problem%Problem(i),Err)

end do 

Solver: do i = 1,num_of_solver

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(all_Solver%Solver(i),Err)
  CALL cmfe_Solver_Initialise(all_LinearSolver%LinearSolver(i),Err)
  CALL cmfe_Problem_SolversCreateStart(all_Problem%Problem(i),Err)
  CALL cmfe_Problem_SolverGet(all_Problem%Problem(i),CMFE_CONTROL_LOOP_NODE,1,all_Solver%Solver(i),Err)
  CALL cmfe_Solver_OutputTypeSet(all_Solver%Solver(i),CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(all_Solver%Solver(i),CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(all_Solver%Solver(i),all_LinearSolver%LinearSolver(i),Err)
  CALL cmfe_Solver_LinearTypeSet(all_LinearSolver%LinearSolver(i),solver_def(Solvers_arguments(2,i)),Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(all_Solver%Solver(i),str2real(Solvers_arguments(5,i)),Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(all_Solver%Solver(i),str2real(Solvers_arguments(7,i)),Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(all_Solver%Solver(i),str2int(Solvers_arguments(6,i)),Err)
  CALL cmfe_Problem_SolversCreateFinish(all_Problem%Problem(i),Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(all_Solver%Solver(i),Err)
  CALL cmfe_SolverEquations_Initialise(all_SolverEquations%SolverEquations(i),Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(all_Problem%Problem(i),Err)
  CALL cmfe_Problem_SolverGet(all_Problem%Problem(i),CMFE_CONTROL_LOOP_NODE,1,all_Solver%Solver(i),Err)
  CALL cmfe_Solver_SolverEquationsGet(all_Solver%Solver(i),all_SolverEquations%SolverEquations(i),Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(all_SolverEquations%SolverEquations(i),all_EquationsSet%EquationsSet(i), & 
                                           & EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(all_Problem%Problem(i),Err)

end do Solver 

 

 
do j = 1,num_of_BOundaryCondition

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(all_BoundaryConditions%BoundaryConditions(j),Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(all_SolverEquations%SolverEquations(j), & 
       & all_BoundaryConditions%BoundaryConditions(j),Err)

  do i = 1,num_of_dirichelet

  	CALL cmfe_GeneratedMesh_SurfaceGet(all_GeneratedMesh%GeneratedMesh(j),&
      					   & bc_def(BC_arg2(2+4*(i-1),j)),SurfaceNodes,LeftNormalXi,Err)

        DO node_idx=1,SIZE(SurfaceNodes,1)
                NodeNumber=SurfaceNodes(node_idx)

		CALL cmfe_Decomposition_NodeDomainGet(all_Decomposition%Decomposition(j),NodeNumber,1,NodeDomain,Err)

                IF(NodeDomain==ComputationalNodeNumber) THEN

                        DO component_idx=1,3

                                constraint =  BC_arg2(3+4*(i-1),j)
                                
                                if (trim(constraint(component_idx:component_idx)) == "1") then 

 					CALL cmfe_BoundaryConditions_SetNode(all_BoundaryConditions%BoundaryConditions(j),all_DependentField%DependentField(j), & 
 					CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, component_idx,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, & 
 				 	& str2real(BC_arg2(4+4*(i-1),j)),Err)

                                 end if 
                        ENDDO
                 ENDIF
        ENDDO

 deallocate(SurfaceNodes)

 end do



! Dear Andreas : Please vaerify is the following way of applying the pressure boundary condition looks alright 
 
! Applying Pressure Neuman Boundary condition

! CALL cmfe_GeneratedMesh_SurfaceGet(all_GeneratedMesh%GeneratedMesh(1),&
 !                                  & CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,SurfaceNodes,LeftNormalXi,Err)


!DO node_idx=1,SIZE(SurfaceNodes,1)

! NodeNumber=SurfaceNodes(node_idx)

! CALL cmfe_Decomposition_NodeDomainGet(all_Decomposition%Decomposition(1),NodeNumber,1,NodeDomain,Err)


! IF(NodeDomain==ComputationalNodeNumber) THEN

!CALL cmfe_BoundaryConditions_SetNode(all_BoundaryConditions%BoundaryConditions(1),all_DependentField%DependentField(1), & 
! 		CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, abs(LeftNormalXi),CMFE_BOUNDARY_CONDITION_PRESSURE, & 
! 		& 2.1_CMISSRP,Err)



 ! ENDIF
 CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(all_SolverEquations%SolverEquations(j),Err)
ENDDO


!CMFE_BOUNDARY_CONDITION_PRESSURE


Problems: do i = 1, num_of_problem
  
  !Solve problem
  CALL cmfe_Problem_Solve(all_Problem%Problem(1),Err)

end do Problems

  !Output solution
FiberFields: do i = 1,1

  CALL cmfe_Fields_Initialise(all_Fields%Fields(i),Err)
  CALL cmfe_Fields_Create(all_Region%Region(i),all_Fields%Fields(i),Err)
  CALL cmfe_Fields_NodesExport(all_Fields%Fields(i),"LargeUniaxialExtension","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(all_Fields%Fields(i),"LargeUniaxialExtension","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(all_Fields%Fields(i),Err)

end do FiberFields
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

 contains 

 include "problems_match_library.f90"

END PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE






















