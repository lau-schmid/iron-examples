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
!> Contributor(s): Adam Reeve, Thomas Heidlauf
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


  REAL(CMISSRP) :: HEIGHT
  REAL(CMISSRP) :: WIDTH
  REAL(CMISSRP) :: LENGTH
!  INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg) 		 :: InterpolationType
  INTEGER(CMISSIntg) 		 :: PressureInterpolationType
  LOGICAL			 :: UsePressureBasis 
!  LOGICAL, PARAMETER :: UsePressureBasis=.FALSE.
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3

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
  INTEGER(CMISSIntg)    :: j ,component_idx
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

  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: SurfaceNodes(:)

  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,BackNormalXi

!  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 3 !nearly incompressible
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 4 !fully incompressible

!!!!!!!!!!!!!!!!!!!!!!!!!!! 		          MY CONTRIBUTION 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  include "DerivedTypes.f90" !!DECLARING THE DERIVED TYPES HERE      

  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif


  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err


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



DO  i = 1,num_of_CoordinateSystem 			!!!!!!!!!!!!! introducing a LOOP for assigning  coordinate systems to thier respective regions !!!!!!!!!!
  !Intialise cmiss

  CALL cmfe_Initialise(all_WorldCoordinateSystem%WorldCoordinateSystem(i),all_WorldRegion%WorldRegion(i),Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=str2int(Mesh_arg4(1,i))
  NumberGlobalYElements=str2int(Mesh_arg4(2,i))
  NumberGlobalZElements=str2int(Mesh_arg4(3,i))
  NumberOfDomains=NumberOfComputationalNodes
  

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


DO i = 1, num_of_basis  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Assigning basis to thier respective edges       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  InterpolationType= match_basis(Basis_arguments(4,i))
  PressureInterpolationType= match_basis(Basis_arguments(3,i))      !!!!! DATA READ FROM THE INPUT FILE !!!!!!!
  UsePressureBasis  = (str2int(Basis_arguments(2,i)) == 4) 

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


END DO

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(all_Decomposition%Decomposition(1),Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber(1),all_Mesh%Mesh(1),all_Decomposition%Decomposition(1),Err)
  CALL cmfe_Decomposition_TypeSet(all_Decomposition%Decomposition(1),CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(all_Decomposition%Decomposition(1),NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(all_Decomposition%Decomposition(1),Err)

DO i = 1, num_of_GeometricField



  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(all_GeometricField%GeometricField(i),Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber(i),all_Region%Region(i),all_GeometricField%GeometricField(i),Err)
  CALL cmfe_Field_MeshDecompositionSet(all_GeometricField%GeometricField(i),all_Decomposition%Decomposition(i),Err)
  CALL cmfe_Field_VariableLabelSet(all_GeometricField%GeometricField(i),CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ScalingTypeSet(all_GeometricField%GeometricField(i),CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_Field_CreateFinish(all_GeometricField%GeometricField(i),Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(all_GeneratedMesh%GeneratedMesh(i),all_GeometricField%GeometricField(i),Err)


  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(all_FibreField%FibreField(i),Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber(i),all_Region%Region(i),all_FibreField%FibreField(i),Err)
  CALL cmfe_Field_TypeSet(all_FibreField%FibreField(i),CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(all_FibreField%FibreField(i),all_Decomposition%Decomposition(i),Err)
  CALL cmfe_Field_GeometricFieldSet(all_FibreField%FibreField(i),all_GeometricField%GeometricField(i),Err)
  CALL cmfe_Field_VariableLabelSet(all_FibreField%FibreField(i),CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(all_FibreField%FibreField(i),Err)
end do  !!!!! END LOOP FOR GEOMETRIC FIELD

  !Create the dependent field
  CALL cmfe_Field_Initialise(all_DependentField%DependentField(1),Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber(1),all_Region%Region(1),all_DependentField%DependentField(1),Err)
  CALL cmfe_Field_TypeSet(all_DependentField%DependentField(1),CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(all_DependentField%DependentField(1),all_Decomposition%Decomposition(1),Err)
  CALL cmfe_Field_GeometricFieldSet(all_DependentField%DependentField(1),all_GeometricField%GeometricField(1),Err)
  CALL cmfe_Field_DependentTypeSet(all_DependentField%DependentField(1),CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(all_DependentField%DependentField(1),2,Err)
  CALL cmfe_Field_VariableTypesSet(all_DependentField%DependentField(1), & 
                                   & [CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_VariableLabelSet(all_DependentField%DependentField(1),CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_NumberOfComponentsSet(all_DependentField%DependentField(1),CMFE_FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,Err)
  CALL cmfe_Field_NumberOfComponentsSet(all_DependentField%DependentField(1), & 
       & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,Err)
  IF(UsePressureBasis) THEN
    !Set the pressure to be nodally based and use the second mesh component if required
    CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(1),&
    & CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(1),CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(1),CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(1),CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  END IF
  CALL cmfe_Field_CreateFinish(all_DependentField%DependentField(1),Err)

 

  !Create the material field
  CALL cmfe_Field_Initialise(all_MaterialField%MaterialField(1),Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber(1),all_Region%Region(1),all_MaterialField%MaterialField(1),Err)
  CALL cmfe_Field_TypeSet(all_MaterialField%MaterialField(1),CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(all_MaterialField%MaterialField(1),all_Decomposition%Decomposition(1),Err)
  CALL cmfe_Field_GeometricFieldSet(all_MaterialField%MaterialField(1),all_GeometricField%GeometricField(1),Err)
  CALL cmfe_Field_NumberOfVariablesSet(all_MaterialField%MaterialField(1),1,Err)
  CALL cmfe_Field_VariableLabelSet(all_MaterialField%MaterialField(1),CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_Field_NumberOfComponentsSet(all_MaterialField%MaterialField(1),CMFE_FIELD_U_VARIABLE_TYPE,5,Err)
  CALL cmfe_Field_CreateFinish(all_MaterialField%MaterialField(1),Err)


  !Create the equations_set


DO i = 1,Num_of_EquationsSet            !!!! loop for eqaution set

  CALL cmfe_Field_Initialise(all_EquationsSetField%EquationsSetField(i),Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber(i),all_Region%Region(i),all_FibreField%FibreField(i),& 
       & [match_equations_set(EquationsSet_arg2(1,i)), &
       & match_equations_set(EquationsSet_arg3(1,i)),match_equations_set(EquationsSet_arg4(1,i))], &
       & EquationsSetFieldUserNumber(i),all_EquationsSetField%EquationsSetField(i),all_EquationsSet%EquationsSet(i),Err)
  CALL cmfe_EquationsSet_CreateFinish(all_EquationsSet%EquationsSet(i),Err)

                            !!! end loop for equation set 


 
  !Create the equations set dependent field
  CALL cmfe_EquationsSet_DependentCreateStart(all_EquationsSet%EquationsSet(i),FieldDependentUserNumber(i), & 
                                              & all_DependentField%DependentField(1),Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(all_EquationsSet%EquationsSet(i),Err)


if (NUMBER_OF_COMPONENTS > 1)  then 		!!!! SO  in case of the laplacian problem the following bit does not work 
  !Create the equations set material field 


  CALL cmfe_EquationsSet_MaterialsCreateStart(all_EquationsSet%EquationsSet(i),FieldMaterialUserNumber(i) & 
                                              & ,all_MaterialField%MaterialField(i),Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(all_EquationsSet%EquationsSet(i),Err)

do j = 1,material_parameters(EquationsSet_arg4(1,i))

   CALL cmfe_Field_ComponentValuesInitialise(all_MaterialField%MaterialField(i), & 
       & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,j,str2real(MaterialField_arg2(j,i)),Err)
  end do 

end if 


!Set the fibre field
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,str2real(FiberField_arg2(1,i)),Err)
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,str2real(FiberField_arg2(2,i)),Err)
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(i),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,str2real(FiberField_arg2(3,i)),Err)


  !Create the equations set equations
  CALL cmfe_Equations_Initialise(all_Equations%Equations(i),Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(all_EquationsSet%EquationsSet(i),all_Equations%Equations(i),Err)
  CALL cmfe_Equations_SparsityTypeSet(all_Equations%Equations(i),CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(all_Equations%Equations(i),CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(all_EquationsSet%EquationsSet(i),Err)

END DO  




DO i = 1,num_of_DependentField
  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(1), & 
                                        match_dependent_field(DependentField_arguments(2,i)),CMFE_FIELD_VALUES_SET_TYPE, &
 & 1,all_DependentField%DependentField(1),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)

  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(1), & 
                                        match_dependent_field(DependentField_arguments(2,i)),CMFE_FIELD_VALUES_SET_TYPE, &
 & 2,all_DependentField%DependentField(1),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)

  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(1)& 
                                       ,match_dependent_field(DependentField_arguments(2,i)),CMFE_FIELD_VALUES_SET_TYPE, &
 & 3,all_DependentField%DependentField(1),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
END DO



DO i = 1, num_of_Problem

  !Define the problem
  CALL cmfe_Problem_Initialise(all_Problem%Problem(1),Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber(1),[match_problem(Problem_arg2(1,i)),match_problem(Problem_arg3(1,i)), &
    & match_problem(Problem_arg4(1,i))],all_Problem%Problem(1),Err)
  CALL cmfe_Problem_CreateFinish(all_Problem%Problem(1),Err)

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



  !Create the problem solvers
  CALL cmfe_Solver_Initialise(all_Solver%Solver(1),Err)
  CALL cmfe_Solver_Initialise(all_LinearSolver%LinearSolver(1),Err)
  CALL cmfe_Problem_SolversCreateStart(all_Problem%Problem(1),Err)
  CALL cmfe_Problem_SolverGet(all_Problem%Problem(1),CMFE_CONTROL_LOOP_NODE,1,all_Solver%Solver(1),Err)
  CALL cmfe_Solver_OutputTypeSet(all_Solver%Solver(1),CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(all_Solver%Solver(1),CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(all_Solver%Solver(1),all_LinearSolver%LinearSolver(1),Err)
  CALL cmfe_Solver_LinearTypeSet(all_LinearSolver%LinearSolver(1),solver_def(Solvers_arguments(2,1)),Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(all_Solver%Solver(1),str2real(Solvers_arguments(5,1)),Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(all_Solver%Solver(1),str2real(Solvers_arguments(7,1)),Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(all_Solver%Solver(1),str2int(Solvers_arguments(6,1)),Err)
  CALL cmfe_Problem_SolversCreateFinish(all_Problem%Problem(1),Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(all_Solver%Solver(1),Err)
  
  CALL cmfe_SolverEquations_Initialise(all_SolverEquations%SolverEquations(1),Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(all_Problem%Problem(1),Err)
  CALL cmfe_Problem_SolverGet(all_Problem%Problem(1),CMFE_CONTROL_LOOP_NODE,1,all_Solver%Solver(1),Err)
  CALL cmfe_Solver_SolverEquationsGet(all_Solver%Solver(1),all_SolverEquations%SolverEquations(1),Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(all_SolverEquations%SolverEquations(1),all_EquationsSet%EquationsSet(1), & 
                                           & EquationsSetIndex,Err)


!Create the problem solver equations
  
 
! CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)

  CALL cmfe_Problem_SolverEquationsCreateFinish(all_Problem%Problem(1),Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(all_BoundaryConditions%BoundaryConditions(1),Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(all_SolverEquations%SolverEquations(1), & 
       & all_BoundaryConditions%BoundaryConditions(1),Err)

!!!!!! 				Generic  form of Assining dirichelet BCs 			!!!!!!!!!!!!
do j = 1,num_of_Mesh

  do i = 1,4

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

 end do 
 
! Dear Andreas : Please vaerify is the following way of applying the pressure boundary condition looks alright 
 
! Applying Pressure Neuman Boundary condition

! DO NN=1,SIZE(SurfaceNodes,1)

! NODE=SurfaceNodes(NN)

! CALL cmfe_Decomposition_NodeDomainGet(all_Decomposition%Decomposition(j),NODE,1,NodeDomain,Err)

! IF(NodeDomain==ComputationalNodeNumber) THEN

! CALL !cmfe_BoundaryConditions_SetNode(all_BoundaryConditions%BoundaryConditions(j),DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NODE,&

!ABS(NormalXi), CMFE_BOUNDARY_CONDITION_PRESSURE,2,Err)


!ENDIF

!ENDDO





  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(all_SolverEquations%SolverEquations(1),Err)

  !Solve problem
  CALL cmfe_Problem_Solve(all_Problem%Problem(1),Err)

  !Output solution
  CALL cmfe_Fields_Initialise(all_Fields%Fields(1),Err)
  CALL cmfe_Fields_Create(all_Region%Region(1),all_Fields%Fields(1),Err)
  CALL cmfe_Fields_NodesExport(all_Fields%Fields(1),"LargeUniaxialExtension","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(all_Fields%Fields(1),"LargeUniaxialExtension","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(all_Fields%Fields(1),Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

 contains 

 include "problems_match_library.f90"

END PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE






















