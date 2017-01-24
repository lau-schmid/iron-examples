!> \file
!> \author: Waleed Mirza
!> \brief This is a generic OpenCMISS-iron example program to solve various different types of equations using OpenCMISS calls.
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
!> Contributor(s): Waleed Mirza, Andreas Hessenthaler, Thomas Heidlauf, Thomas Klotz
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING MODULE INCLUDES THE DATA STRUCTURES AND  "CONTAIN SUBROUTINES AND FUNCTIONS" USED IN PARSING ALGORITHM  !!!!!!!!!!!!!!!!1
include "parsing_module.f90"

!> Main program
PROGRAM GENERICEXAMPLE

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
  INTEGER(CMISSIntg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=300) :: CommandArgument,inputFile
  INTEGER(CMISSIntg) 		 :: InterpolationType
  INTEGER(CMISSIntg) 		 :: PressureInterpolationType
  LOGICAL			 :: UsePressureBasis 
  INTEGER(CMISSIntg)             :: ScalingType
  INTEGER(CMISSIntg)             :: NumberOfGaussXi

!!!!!!!!!!!!!!!!!!!!		COUTER SIZE		!!!!!!!!!!!!!!!!!!!!!!!!

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
  INTEGER(CMISSIntg)    :: j,k ,component_idx
  INTEGER(CMISSIntg)    :: Csys_idx, Region_idx,Basis_idx,Decomposition_idx,BoundaryCondition_idx, & 
                         & Problem_idx,Solver_idx,EquationSet_idx,MaterialField_idx,DependentField_idx, & 
 			   FiberField_idx,GeometricField_idx,mesh_idx,ControlLoop_idx	
  INTEGER(CMISSIntg)    :: num_of_dirichelet,num_of_traction_neumann,num_of_pressure_neumann
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
  INTEGER(CMISSIntg) ,allocatable  :: FieldUserNumber(:)
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


  INTEGER(CMISSIntg)		   NUMBER_OF_COMPONENTS  

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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 		INITIALIZE THE DERIVED DATA STRUCTURES 			 !!!!!!!!!!!!!!!!!!!!!
  include "AllocatingDerivedDataStructures.f90"



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              ALLOCATE  "USERNUMBER"  DATA STRUCTURE SIZE 		 !!!!!!!!!!!!!!!!!!!!!!
  include "UserNumberDataStructures.f90"

 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        INCLUDES I/0 STATEMENT THAT PARSE THROUGH THE INPUT FILE        !!!!!!!!!!!!!!!!!!!!!!
  include "parsing_algorithm.f90"




  NUMBER_OF_COMPONENTS = str2int(DependentField_arguments(3,1))

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

  !!!!!!!!!!!!! introducing a LOOP for assigning  coordinate systems to their respective regions !!!!!!!!!!
  DO  Csys_idx = 1,num_of_CoordinateSystem
    !Create a 3D rectangular cartesian coordinate system
    CALL cmfe_CoordinateSystem_Initialise(all_CoordinateSystem%CoordinateSystem(Csys_idx),Err)
    CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber(Csys_idx),all_CoordinateSystem%CoordinateSystem(Csys_idx),Err)
    Call cmfe_CoordinateSystem_TypeSet(CoordinateSystemUserNumber(Csys_idx), &
      & match_coordinate_system(CoordinateSystem_arguments(2,Csys_idx)), err)
    CALL cmfe_CoordinateSystem_CreateFinish(all_CoordinateSystem%CoordinateSystem(Csys_idx),Err)
  END DO ! Csys_idx

  DO Region_idx = 1, num_of_region
    !Create a region and assign the coordinate system to the region
    CALL cmfe_Region_Initialise(all_Region%Region(Region_idx),Err)
    CALL cmfe_Region_CreateStart(RegionUserNumber(Region_idx),all_WorldRegion%WorldRegion(Region_idx), & 
      & all_Region%Region(Region_idx),Err)
    CALL cmfe_Region_LabelSet(all_Region%Region(Region_idx),"Region",Err)
    CALL cmfe_Region_CoordinateSystemSet(all_Region%Region(Region_idx), & 
      & (all_CoordinateSystem%CoordinateSystem(Region_idx)),Err)
    CALL cmfe_Region_CreateFinish(all_Region%Region(Region_idx),Err)
  END DO  ! region_idx


DO basis_idx = 1, num_of_basis  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Assigning basis to their respective edges       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  UsePressureBasis  = (str2int(Basis_arguments(2,basis_idx)) == 4)

  if (UsePressureBasis)  then
 	j=2
  else 
        j=1
  end if 

  do k = 1,j
  
  	InterpolationType= match_basis(Basis_arguments(4-k+1,basis_idx))
  
  

!  PressureInterpolationType= match_basis(Basis_arguments(3,i))      !!!!! DATA READ FROM THE INPUT FILE !!!!!!!
   
  	NumberGlobalXElements=str2int(Mesh_arg4(1,basis_idx))
  	NumberGlobalYElements=str2int(Mesh_arg4(2,basis_idx))
  	NumberGlobalZElements=str2int(Mesh_arg4(3,basis_idx))
  	NumberOfGaussXi = str2int(Basis_arguments(5,basis_idx))

  !Define geometric basis
  	CALL cmfe_Basis_Initialise(all_Basis%Basis(k),Err)
  	CALL cmfe_Basis_CreateStart(BasisUserNumber(k),all_Basis%Basis(k),Err)
  	SELECT CASE(InterpolationType)
  		CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, & 
		    CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION )
    			CALL cmfe_Basis_TypeSet(all_Basis%Basis(k),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  		CASE(CMFE_BASIS_QUADRATIC1_HERMITE_INTERPOLATION ,CMFE_BASIS_QUADRATIC2_HERMITE_INTERPOLATION , & 
		     CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION )
    			CALL cmfe_Basis_TypeSet(all_Basis%Basis(k),CMFE_BASIS_SIMPLEX_TYPE,Err)
  	END SELECT
  	IF(NumberGlobalZElements==0) THEN
   		CALL cmfe_Basis_NumberOfXiSet(all_Basis%Basis(k),2,Err)
    		CALL cmfe_Basis_InterpolationXiSet(all_Basis%Basis(k),[InterpolationType,InterpolationType],Err)
    		IF(NumberOfGaussXi>0) THEN
      			CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis%Basis(k),[NumberOfGaussXi,NumberOfGaussXi],Err)
    		ENDIF
  	ELSE
    		CALL cmfe_Basis_NumberOfXiSet(all_Basis%Basis(k),3,Err)
    		CALL cmfe_Basis_InterpolationXiSet(all_Basis%Basis(k),[InterpolationType,InterpolationType,InterpolationType],Err)
    		IF(NumberOfGaussXi>0) THEN
      			CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis%Basis(k),[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    		ENDIF
  	ENDIF ! end_k
  	CALL cmfe_Basis_CreateFinish(all_Basis%Basis(k),Err)
end do ! basis_idx 


end do 


do mesh_idx = 1,num_of_Mesh

   HEIGHT=str2int(Mesh_arg3(1,mesh_idx))
   WIDTH=str2int(Mesh_arg3(2,mesh_idx))
   LENGTH=str2int(Mesh_arg3(3,mesh_idx))

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(all_GeneratedMesh%GeneratedMesh(mesh_idx),Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber(mesh_idx),all_Region%Region(mesh_idx),& 
 				      all_GeneratedMesh%GeneratedMesh(mesh_idx),Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(all_GeneratedMesh%GeneratedMesh(mesh_idx),CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  IF(UsePressureBasis) THEN
    CALL cmfe_GeneratedMesh_BasisSet(all_GeneratedMesh%GeneratedMesh(mesh_idx),[all_Basis%Basis(mesh_idx), & 
                                     & all_Basis%Basis(mesh_idx+1)],Err)
  ELSE
    CALL cmfe_GeneratedMesh_BasisSet(all_GeneratedMesh%GeneratedMesh(mesh_idx),[all_Basis%Basis(mesh_idx)],Err)
  ENDIF
  !Define the mesh on the region
  IF(NumberGlobalZElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh%GeneratedMesh(mesh_idx),[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh%GeneratedMesh(mesh_idx), & 
                                                & [NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh%GeneratedMesh(mesh_idx),[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh%GeneratedMesh(mesh_idx), &
                                    & [NumberGlobalXElements,NumberGlobalYElements, &
                                    & NumberGlobalZElements],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(all_Mesh%Mesh(mesh_idx),Err)
  CALL cmfe_GeneratedMesh_CreateFinish(all_GeneratedMesh%GeneratedMesh(mesh_idx),MeshUserNumber(mesh_idx), & 
					& all_Mesh%Mesh(mesh_idx),Err)

END do !mesh_idx

do Decomposition_idx = 1,num_of_Decomposition
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(all_Decomposition%Decomposition(Decomposition_idx),Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber(Decomposition_idx), & 
				   all_Mesh%Mesh(Decomposition_idx),all_Decomposition%Decomposition(Decomposition_idx),Err)
  CALL cmfe_Decomposition_TypeSet(all_Decomposition%Decomposition(Decomposition_idx),CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(all_Decomposition%Decomposition(Decomposition_idx),NumberOfDomains,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(all_Decomposition%Decomposition(Decomposition_idx),.TRUE.,Err) 
  CALL cmfe_Decomposition_CreateFinish(all_Decomposition%Decomposition(Decomposition_idx),Err)
end do !decomp_idx

DO GeometricField_idx  = 1, num_of_GeometricField

  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(all_GeometricField%GeometricField(GeometricField_idx),Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber(GeometricField_idx),all_Region%Region(GeometricField_idx) & 
						 	,all_GeometricField%GeometricField(GeometricField_idx),Err)
  CALL cmfe_Field_MeshDecompositionSet(all_GeometricField%GeometricField(GeometricField_idx), & 
					all_Decomposition%Decomposition(GeometricField_idx),Err)
  CALL cmfe_Field_VariableLabelSet(all_GeometricField%GeometricField(GeometricField_idx), & 
						& CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ScalingTypeSet(all_GeometricField%GeometricField(GeometricField_idx), & 
						& CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_Field_CreateFinish(all_GeometricField%GeometricField(GeometricField_idx),Err)

  !Update the geometric field parameters
   CALL cmfe_GeneratedMesh_GeometricParametersCalculate(all_GeneratedMesh%GeneratedMesh(GeometricField_idx), & 
 						all_GeometricField%GeometricField(GeometricField_idx),Err)

end do ! GeometricField_idx  !!!!! END LOOP FOR GEOMETRIC FIELD


  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(all_FibreField%FibreField(1),Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber(1),all_Region%Region(1), & 
				all_FibreField%FibreField(1),Err)
  CALL cmfe_Field_TypeSet(all_FibreField%FibreField(1),CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(all_FibreField%FibreField(1), & 
				all_Decomposition%Decomposition(1),Err)
  CALL cmfe_Field_GeometricFieldSet(all_FibreField%FibreField(1), & 
				all_GeometricField%GeometricField(1),Err)
  CALL cmfe_Field_VariableLabelSet(all_FibreField%FibreField(1), & 
					CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_ScalingTypeSet(all_FibreField%FibreField(1), & 
 						CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_Field_CreateFinish(all_FibreField%FibreField(1),Err)

DO FiberField_idx = 1, num_of_FiberField
    !Set the fibre field
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(FiberField_idx), & 
                                CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    					    & 1,str2real(FiberField_arg2(1,FiberField_idx)),Err)
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(FiberField_idx), & 
				CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    					    & 2,str2real(FiberField_arg2(2,FiberField_idx)),Err)
  CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(FiberField_idx),& 
				CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    					    & 3,str2real(FiberField_arg2(3,FiberField_idx)),Err)
end do   !!!!! END LOOP FOR Fiber FIELD

DO EquationSet_idx  = 1,Num_of_EquationsSet            !!!! loop for eqaution set


  CALL cmfe_EquationsSet_Initialise(all_EquationsSet%EquationsSet(1),Err)
  CALL cmfe_Field_Initialise(all_EquationsSetField%EquationsSetField(1),Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber(1),all_Region%Region(1) & 
	 ,all_FibreField%FibreField(1),& 
       & [match_equations_set(EquationsSet_arg2(1,EquationSet_idx)), &
       & match_equations_set(EquationsSet_arg3(1,EquationSet_idx)), & 
						    match_equations_set(EquationsSet_arg4(1,EquationSet_idx))], &
       & EquationsSetFieldUserNumber(1),all_EquationsSetField%EquationsSetField(1), & 
         all_EquationsSet%EquationsSet(1),Err)
  CALL cmfe_EquationsSet_CreateFinish(all_EquationsSet%EquationsSet(1),Err)


end do ! EquationSet_idx


  !Create the dependent field
  DO DependentField_idx = 1, num_of_DependentField
  	  !Create the dependent field
  	CALL cmfe_Field_Initialise(all_DependentField%DependentField(1),Err)
  	CALL cmfe_EquationsSet_DependentCreateStart(all_EquationsSet%EquationsSet(1),FieldDependentUserNumber(1), & 
									& all_DependentField%DependentField(1),Err)
  	CALL cmfe_Field_VariableLabelSet(all_DependentField%DependentField(1),CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
 	 DO component_idx=1,NUMBER_OF_COMPONENTS
    		CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(1), & 
		CMFE_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    		CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(1), & 
		CMFE_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  	ENDDO


  IF(UsePressureBasis) THEN
     CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(1), & 
		CMFE_FIELD_U_VARIABLE_TYPE,4, 1 ,Err)
     CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(1), &
		CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,1,Err)
    !Set the pressure to be nodally based and use the second mesh component if required
     CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(DependentField_idx),&
     & CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
     CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(DependentField_idx), & 
								CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
                                                             & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
     CALL cmfe_Field_ComponentMeshComponentSet( & 
		all_DependentField%DependentField(DependentField_idx),CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
     CALL cmfe_Field_ComponentMeshComponentSet( & 
	 all_DependentField%DependentField(DependentField_idx),CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)

  END IF
  CALL cmfe_EquationsSet_DependentCreateFinish(all_EquationsSet%EquationsSet(1),Err)


  
!!!
  if  (trim(DependentField_arguments(4,DependentField_idx)) == "UNDEFORMED") then
   
  	CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(DependentField_idx), & 
       	match_dependent_field(DependentField_arguments(2,DependentField_idx)),CMFE_FIELD_VALUES_SET_TYPE, &
 	& 1,all_DependentField%DependentField(DependentField_idx),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
	if (NUMBER_OF_COMPONENTS .ge. 2) then
  		CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(DependentField_idx), & 
                                match_dependent_field(DependentField_arguments(2,DependentField_idx)),CMFE_FIELD_VALUES_SET_TYPE, &
 		& 2,all_DependentField%DependentField(DependentField_idx),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
	end if 

  	if (NumberGlobalZElements .NE. 0) then
  	CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(DependentField_idx), & 
                                match_dependent_field(DependentField_arguments(2,DependentField_idx)),CMFE_FIELD_VALUES_SET_TYPE, &
 	     & 3,all_DependentField%DependentField(DependentField_idx),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
        end if



  else 

         CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(DependentField_idx), & 
              &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
              & str2real(DependentField_arguments(4,DependentField_idx)),Err)
	if (NUMBER_OF_COMPONENTS .ge. 2) then
         	CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(DependentField_idx), & 
              	&  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
              	& str2real(DependentField_arguments(4,DependentField_idx)),Err)
	end if 
  	if (NumberGlobalZElements .NE. 0 .AND. NUMBER_OF_COMPONENTS .ge. 3) then
         	CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(DependentField_idx), & 
              		&  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
              		& str2real(DependentField_arguments(4,DependentField_idx)),Err)
        end if
   
  end if 
 
  IF(UsePressureBasis) THEN
  	CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(DependentField_idx), & 
                      &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
                      & str2real(DependentField_arguments(5,DependentField_idx)),Err)
  end if 

  CALL cmfe_Field_ParameterSetUpdateStart(all_DependentField%DependentField(DependentField_idx),CMFE_FIELD_U_VARIABLE_TYPE & 
                                          & ,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(all_DependentField%DependentField(DependentField_idx),CMFE_FIELD_U_VARIABLE_TYPE, & 
                                          & CMFE_FIELD_VALUES_SET_TYPE,Err)


end do ! dependent fields


  


  !Create the material field
do MaterialField_idx  = 1,num_of_MaterialField
  CALL cmfe_Field_Initialise(all_MaterialField%MaterialField(MaterialField_idx),Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber(MaterialField_idx),all_Region%Region(MaterialField_idx), & 
							all_MaterialField%MaterialField(MaterialField_idx),Err)
  CALL cmfe_Field_TypeSet(all_MaterialField%MaterialField(MaterialField_idx),CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(all_MaterialField%MaterialField(MaterialField_idx), & 
							all_Decomposition%Decomposition(MaterialField_idx),Err)
  CALL cmfe_Field_GeometricFieldSet(all_MaterialField%MaterialField(MaterialField_idx), & 
							all_GeometricField%GeometricField(MaterialField_idx),Err)
  CALL cmfe_Field_NumberOfVariablesSet(all_MaterialField%MaterialField(MaterialField_idx),1,Err)
  CALL cmfe_Field_VariableLabelSet(all_MaterialField%MaterialField(MaterialField_idx), & 
				   CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)

  CALL cmfe_Field_NumberOfComponentsSet(all_MaterialField%MaterialField(MaterialField_idx),CMFE_FIELD_U_VARIABLE_TYPE, & 
                                        & material_parameters(EquationsSet_arg4(1,MaterialField_idx)),Err)
  CALL cmfe_Field_CreateFinish(all_MaterialField%MaterialField(MaterialField_idx),Err)

  do j = 1,material_parameters(EquationsSet_arg4(1,MaterialField_idx))

   	CALL  cmfe_Field_ComponentValuesInitialise(all_MaterialField%MaterialField(MaterialField_idx), & 
                                              &CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,& 
                                              & j,str2real(MaterialField_arg2(j,MaterialField_idx)),Err)
  end do 


end do  ! MaterialField_idx

  !Create the equations_set


DO EquationSet_idx  = 1,Num_of_EquationsSet            !!!! loop for eqaution set


if (num_of_MaterialField .ge. 1) then
  CALL cmfe_EquationsSet_MaterialsCreateStart(all_EquationsSet%EquationsSet(EquationSet_idx), & 
						     FieldMaterialUserNumber(EquationSet_idx) & 
                                              & ,all_MaterialField%MaterialField(EquationSet_idx),Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(all_EquationsSet%EquationsSet(EquationSet_idx),Err)
end if 

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(all_Equations%Equations(EquationSet_idx),Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(all_EquationsSet%EquationsSet(EquationSet_idx), & 
						all_Equations%Equations(EquationSet_idx),Err)
  CALL cmfe_Equations_SparsityTypeSet(all_Equations%Equations(EquationSet_idx),output_type(Output_arguments(1,1)),Err)
  CALL cmfe_Equations_OutputTypeSet(all_Equations%Equations(EquationSet_idx),output_type(Output_arguments(2,1)),Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(all_EquationsSet%EquationsSet(EquationSet_idx),Err)


if (Field_arg2(1,1) == "GRAVITY" ) then
  !Create the source field with the gravity vector
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(all_EquationsSet%EquationsSet(EquationSet_idx),FieldUserNumber(num_of_field) & 
       ,SourceField,Err)
  CALL cmfe_Field_ScalingTypeSet(SourceField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_EquationsSet_SourceCreateFinish(all_EquationsSet%EquationsSet(EquationSet_idx),Err)
  DO component_idx=1,3
	
    	CALL cmfe_Field_ComponentValuesInitialise(SourceField, & 
                       & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
                       & component_idx,str2real(Field_arg3(component_idx,1)),Err)

  ENDDO
end if 
END DO ! EquationsSet_idx 

DO Problem_idx = 1, num_of_Problem

  !Define the problem
  CALL cmfe_Problem_Initialise(all_Problem%Problem(Problem_idx),Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber(Problem_idx),[match_problem(Problem_arg2(1,Problem_idx)), & 
						 &    match_problem(Problem_arg3(1,Problem_idx)), &
    				     & match_problem(Problem_arg4(1,Problem_idx))],all_Problem%Problem(Problem_idx),Err)
  CALL cmfe_Problem_CreateFinish(all_Problem%Problem(Problem_idx),Err)

END DO ! problem_idx

do ControlLoop_idx = 1,num_of_ControlLoop 
  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(all_Problem%Problem(ControlLoop_idx),Err)
  CALL cmfe_ControlLoop_Initialise(all_ControlLoop%ControlLoop(ControlLoop_idx),Err)
  CALL cmfe_Problem_ControlLoopGet(all_Problem%Problem(ControlLoop_idx),CMFE_CONTROL_LOOP_NODE, & 
						all_ControlLoop%ControlLoop(ControlLoop_idx),Err)
  CALL cmfe_ControlLoop_TypeSet(all_ControlLoop%ControlLoop(ControlLoop_idx), & 
				  control_loop_def(ControlLoop_arguments(2,ControlLoop_idx)),Err)
  CALL cmfe_ControlLoop_LoadOutputSet(all_ControlLoop%ControlLoop(ControlLoop_idx),1,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(all_ControlLoop%ControlLoop(ControlLoop_idx), & 
					str2int(ControlLoop_arguments(3,ControlLoop_idx)),Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(all_Problem%Problem(ControlLoop_idx),Err)

end do ! ControlLoop_idx

do solver_idx = 1,num_of_solver

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(all_Solver%Solver(solver_idx),Err)
  CALL cmfe_Solver_Initialise(all_LinearSolver%LinearSolver(solver_idx),Err)
  CALL cmfe_Problem_SolversCreateStart(all_Problem%Problem(solver_idx),Err)
  CALL cmfe_Problem_SolverGet(all_Problem%Problem(solver_idx),CMFE_CONTROL_LOOP_NODE,1,& 
	 		      all_Solver%Solver(solver_idx),Err)
  CALL cmfe_Solver_OutputTypeSet(all_Solver%Solver(solver_idx),CMFE_SOLVER_PROGRESS_OUTPUT,Err)

 if (trim(Solvers_arguments(8,solver_idx)) == "ACTIVE") then 
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(all_Solver%Solver(solver_idx),CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)

  CALL cmfe_Solver_NewtonLinearSolverGet(all_Solver%Solver(solver_idx),all_LinearSolver%LinearSolver(solver_idx),Err)
  CALL cmfe_Solver_LinearTypeSet(all_LinearSolver%LinearSolver(solver_idx),solver_def(Solvers_arguments(2,solver_idx)),Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(all_Solver%Solver(solver_idx),str2real(Solvers_arguments(5,solver_idx)),Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(all_Solver%Solver(solver_idx),str2real(Solvers_arguments(7,solver_idx)),Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(all_Solver%Solver(solver_idx),str2int(Solvers_arguments(6,solver_idx)),Err)
 end if 
  CALL cmfe_Problem_SolversCreateFinish(all_Problem%Problem(solver_idx),Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(all_Solver%Solver(solver_idx),Err)
  CALL cmfe_SolverEquations_Initialise(all_SolverEquations%SolverEquations(solver_idx),Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(all_Problem%Problem(solver_idx),Err)
  CALL cmfe_Problem_SolverGet(all_Problem%Problem(solver_idx),CMFE_CONTROL_LOOP_NODE,1,all_Solver%Solver(solver_idx),Err)
  CALL cmfe_Solver_SolverEquationsGet(all_Solver%Solver(solver_idx),all_SolverEquations%SolverEquations(solver_idx),Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(all_SolverEquations%SolverEquations(solver_idx), & 
						  all_EquationsSet%EquationsSet(solver_idx), & 
                                           	   & EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(all_Problem%Problem(solver_idx),Err)

end do ! solver_idx 


do BOundaryCondition_idx = 1,num_of_BOundaryCondition

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(all_BoundaryConditions%BoundaryConditions(BOundaryCondition_idx),Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(all_SolverEquations%SolverEquations(BOundaryCondition_idx),& 
       & all_BoundaryConditions%BoundaryConditions(BOundaryCondition_idx),Err)

  do i = 1,num_of_dirichelet

  	CALL cmfe_GeneratedMesh_SurfaceGet(all_GeneratedMesh%GeneratedMesh(BOundaryCondition_idx),&
      				    & bc_def(BC_arg2(2+4*(i-1),BOundaryCondition_idx)),SurfaceNodes,LeftNormalXi,Err)

        DO node_idx=1,SIZE(SurfaceNodes,1)
                NodeNumber=SurfaceNodes(node_idx)

		CALL cmfe_Decomposition_NodeDomainGet(all_Decomposition%Decomposition(BOundaryCondition_idx), & 
										 & NodeNumber,1,NodeDomain,Err)

                IF(NodeDomain==ComputationalNodeNumber) THEN

                        DO component_idx=1,3

                                constraint =  BC_arg2(3+4*(i-1),BOundaryCondition_idx)
                                
                                if (trim(constraint(component_idx:component_idx)) == "1") then 

 				 CALL cmfe_BoundaryConditions_SetNode(all_BoundaryConditions%BoundaryConditions(BOundaryCondition_idx), & 
				        all_DependentField%DependentField(BOundaryCondition_idx), & 
 					CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, component_idx,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, & 
 				 	& str2real(BC_arg2(4+4*(i-1),BOundaryCondition_idx)),Err)

                                 end if 
                        ENDDO
                 ENDIF
        ENDDO

 deallocate(SurfaceNodes)

 end do ! BOundaryCondition_idx



!#DO node_idx=1,SIZE(SurfaceNodes,1)

! NodeNumber=SurfaceNodes(node_idx)

! CALL cmfe_Decomposition_NodeDomainGet(all_Decomposition%Decomposition(1),NodeNumber,1,NodeDomain,Err) 


! IF(NodeDomain==ComputationalNodeNumber) THEN


  !    CALL cmfe_BoundaryConditions_SetNode(all_BoundaryConditions%BoundaryConditions(1),all_DependentField%DependentField(1), &
  !      & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber, &
  !      & ABS(LeftNormalXi), &
  !      & CMFE_BOUNDARY_CONDITION_PRESSURE,2.1_CMISSRP,Err)


!  ENDIF
!ENDDO

 CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(all_SolverEquations%SolverEquations(BOundaryCondition_idx),Err)
ENDDO


!CMFE_BOUNDARY_CONDITION_PRESSURE


 do problem_idx = 1, num_of_problem
  
  !Solve problem
  CALL cmfe_Problem_Solve(all_Problem%Problem(problem_idx),Err)

end do  ! problem_idx

  !Output solution
 do i = 1,1

  CALL cmfe_Fields_Initialise(all_Fields%Fields(i),Err)
  CALL cmfe_Fields_Create(all_Region%Region(i),all_Fields%Fields(i),Err)
  CALL cmfe_Fields_NodesExport(all_Fields%Fields(i),"LargeUniaxialExtension","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(all_Fields%Fields(i),"LargeUniaxialExtension","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(all_Fields%Fields(i),Err)

end do 
#
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

 contains 

 include "problems_match_library.f90"

END PROGRAM GENERICEXAMPLE
