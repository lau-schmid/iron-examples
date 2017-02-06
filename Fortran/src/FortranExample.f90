!> \file
!> \author: Waleed Mirza
!> \brief This is a generic OpenCMISS-iron example program to solve various different types of equations using OpenCMISS CALLs.
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
!> KingDOm. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Waleed Mirza, Andreas Hessenthaler, Thomas Heidlauf, Thomas Klotz
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. IF you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. IF you DO not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE INCLUDES THE DATA STRUCTURES  SUBROUTINES AND FUNCTIONS EMPLOYED IN PARSING ALGORITHM.  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INCLUDE "parsing_module.f90"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAIN PROGRAM STARTS FROM HERE. !!!!!!!!!!!!!!!!!!
PROGRAM GENERICEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE Parsing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOGICAL PREPROCESSOR  DIRECTIVES. !!!!!!!!!!!!!!!
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#INCLUDE "mpif.h"
#endif


  !!!!!!!!!!!!!!!!!!!!VARIABLES DEFINING SIZES OF DERIVED DATA STRUCTURES. !!!!!!!!!!!!!!!

  INTEGER(CMISSIntg)              :: NumberOfPressureBasis
  INTEGER(CMISSIntg)              :: NumberOfBasis
  INTEGER(CMISSIntg)              :: NumberOfBoundaryCondition
  INTEGER(CMISSIntg),PARAMETER    :: NumberOfWorldCoordinateSystem = 1 ! ALways equal to 1.
  INTEGER(CMISSIntg),PARAMETER    :: NumberOfWorldRegion = 1           ! Always equal to 1.
  INTEGER(CMISSIntg)              :: NumberOfCoordinateSystem
  INTEGER(CMISSIntg)              :: NumberOfMesh
  INTEGER(CMISSIntg)              :: NumberOfDecomposition
  INTEGER(CMISSIntg)              :: NumberOfEquation
  INTEGER(CMISSIntg)              :: NumberOfEquationsSet
  INTEGER(CMISSIntg)              :: NumberOfGeometricField
  INTEGER(CMISSIntg)              :: NumberOfFiberField
  INTEGER(CMISSIntg)              :: NumberOfMaterialField
  INTEGER(CMISSIntg)              :: NumberOfDependentField
  INTEGER(CMISSIntg)              :: NumberOfEquationSetField
  INTEGER(CMISSIntg)              :: NumberOfSourceField
  INTEGER(CMISSIntg)              :: NumberOfProblem
  INTEGER(CMISSIntg)              :: NumberOfRegion
  INTEGER(CMISSIntg)              :: NumberOfSolver
  INTEGER(CMISSIntg)              :: NumberOfLinearSolver
  INTEGER(CMISSIntg)              :: NumberOfSolverEquations
  INTEGER(CMISSIntg)              :: NumberOfControlLoop
  INTEGER(CMISSIntg)              :: NumberOfGeneratedMesh
  INTEGER(CMISSIntg)              :: NumberOfOutput

  !!  These store the number of boundary conditions defined by user in the input file.
  INTEGER(CMISSIntg)                :: NumberOfDirichelet,NumberOfTractionNeumann,NumberOfPressureNeumann

  !!!!!!!!!!!!!!!!!!!!        COUNTERS USED IN FORTRAN EXAMPLE FILE.                    !!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!   FOR INSTANCE "SolverIdx" IS USED IN THE BLOCK ........     !!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!  ......WHERE SOLVER IS DEFINED.                                  !!! !!!!!!!!!!!!

  INTEGER(CMISSIntg)            :: CsysIdx, RegionIdx,BasisIdx,DecompositionIdx,BoundaryConditionIdx
  INTEGER(CMISSIntg)            :: ProblemIdx,SolverIdx,EquationSetIdx,MaterialFieldIdx,DependentFieldIdx
  INTEGER(CMISSIntg)            :: FiberFieldIdx,GeometricFieldIdx,MeshIdx,ControlLoopIdx, FieldIdx
  !! Counters used in different sections of FortranExample.f90.
  INTEGER(CMISSIntg)            :: j,k , UserNumberLabel, ComponentIdx,MaterialParametersIdx
  INTEGER(CMISSIntg)            :: ProblemId, BasisId , BoundaryId, OutputId, RegionId,EquationId
  INTEGER(CMISSIntg)            :: ControlLoopId, SolverSettingsId , DependentFieldId, CoordinateSystemId, DecompositionId
  INTEGER(CMISSIntg)            :: MaterialFieldId, FibreFieldId , MeshSettingsId, SourceFieldId
  !! Parameter to be used in Boundary Condition bit.
  CHARACTER(LEN = 100)          :: Constraint
  
  !!!!!!!!!!!!!!!! FOLLOWING DATA STRUCTURE STORE THE ID´s ASSINGED TO REGIONS, BASIS etc. !!!!!!!!

  INTEGER(CMISSIntg), ALLOCATABLE :: CoordinateSystemUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: RegionUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: BasisUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: PressureBasisUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: GeneratedMeshUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: MeshUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: DecompositionUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FieldGeometryUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FieldFibreUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FieldMaterialUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FieldDependentUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: EquationSetUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: EquationsSetFieldUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: ProblemUserNumber(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FieldUserNumber(:)


  !!!!!!!!!!!!!!! FOLLOWING ARE THE PROGRAM VARIABLES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER(CMISSIntg)              :: Err
  REAL(CMISSRP)                   :: Height,Width,Length                                                        !! These variables stores dimensions of the DOmain.
  INTEGER(CMISSIntg)              :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements          !! These variables stores  number of elements in each dimension.
  INTEGER(CMISSIntg)              :: NumberOfArguments,ArgumentLength,ArgStatus                                 !! These variables are used in reading input arguments "path + file name"" ..
  CHARACTER(LEN=300)              :: CommandArgument,InputFile                                                  !! ... from the command line. 
  INTEGER(CMISSIntg)              :: InterpolationType,NumberOfGaussXi(3),NumberOfComponents,DependentFieldComponent !! These variables store information to be used in "BASIS BLOCK"  defined by ..
  LOGICAL                         :: UsePressureBasis                                                           !! ... the user in the input file.
  INTEGER(CMISSIntg)              :: ScalingType                                                                !! As of now, i have no idea what this is for , i will look into it later.
  INTEGER(CMISSIntg)              :: EquationsSetIndex                                                          !! As of now, i have no idea what this is for , i will look into it later.
  INTEGER(CMISSIntg)              :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber         !! THese variables are used in "DECOMPOSITION BLOCK".
  INTEGER(CMISSIntg)              :: NodeNumber,NodeDOmain,nodeIdx                                              !! Same as the last comment
  INTEGER(CMISSIntg),ALLOCATABLE  :: SurfaceNodes(:)                                                            !! Possible surface where Dirichelet  or Natural BCs are imposed.
  INTEGER(CMISSIntg)              :: LeftNormalXi                                                               !! Normal directive of possible surfaces where Dirichelet  or Natural BCs are imposed.

  !!!!!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING HEADER FILE INITIALIZE  THE DATA STRUCTURES WITH DERIVED TYPES.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INCLUDE "DerivedTypes.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOGICAL PREPROCESSOR  DIRECTIVES. !!!!!!!!!!!!!!!
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING BLOCK READS THE PATH AND NAME OF THE INPUT FILE.   !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IT MAKES USE OF THE INTRINSIC FUNCTIONS. : ...................   !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ... 1- COMMAND_ARGUMENT_COUNT(). 2- GET_COMMAND_ARGUMENT().      !!!!!!!!!!!!!!!!!!!!!!!!

  NumberOfArguments = COMMAND_ARGUMENT_COUNT()

  IF(NumberOfArguments .NE. 1) THEN

    CALL HANDLE_ERROR("Please provide only the input file")                          !! Throws warning IF the path is wrong.

  ELSE

    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgStatus)
    InputFile = TRIM(CommandArgument)

  END IF



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  FOLLLOWING HEADER FILE INITIALIZE THE SIZE OF DERIVED DATA STRUCTURES. !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INCLUDE "AllocatingDerivedDataStructures.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING HEADER FILE ALLOCATES  "USERNUMBER"  DATA STRUCTURE SIZE           !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INCLUDE "UserNumberDataStructures.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        INCLUDES I/0 STATEMENT THAT PARSE THROUGH THE INPUT FILE        !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INCLUDE "parsing_algorithm.f90"


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


  !!!!!!!!!!!!! INTRODUCING  LOOPS FOR ASSIGNING COORDINATE SYSTEM TO THIER RESPECTIVE REGIONS !!!!!!!!!!

  DO  CsysIdx = 1,NumberOfCoordinateSystem
    !Create a 3D rectangular cartesian coordinate system

    CALL cmfe_CoordinateSystem_Initialise(all_CoordinateSystem%CoordinateSystem(CsysIdx),Err)
    CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber(CsysIdx),all_CoordinateSystem%CoordinateSystem(CsysIdx),Err)
    CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystemUserNumber(CsysIdx), &
     & match_coordinate_system(CoordinateSystem%CoordinateSystemType(1,CsysIdx)), err)
    CALL cmfe_CoordinateSystem_CreateFinish(all_CoordinateSystem%CoordinateSystem(CsysIdx),Err)

  END DO ! CsysIdx

  DO RegionIdx = 1, NumberOfRegion
    !Create a region and assign the coordinate system to the region

    CALL cmfe_Region_Initialise(all_Region%Region(RegionIdx),Err)
    CALL cmfe_Region_CreateStart(RegionUserNumber(RegionIdx),all_WorldRegion%WorldRegion(RegionIdx), &
      & all_Region%Region(RegionIdx),Err)
    CALL cmfe_Region_LabelSet(all_Region%Region(RegionIdx),"Region",Err)
    CALL cmfe_Region_CoordinateSystemSet(all_Region%Region(RegionIdx), &
      & (all_CoordinateSystem%CoordinateSystem(RegionIdx)),Err)
    CALL cmfe_Region_CreateFinish(all_Region%Region(RegionIdx),Err)

  END DO  ! RegionIdx
 !!!!!!!!!!!!!!!!!!!!!!!!!  ASSIGNING BASIS TO RESPECTIVE REGIONS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  UsePressureBasis = (NumberOfPressureBasis .GE. 1)

  DO BasisIdx = 1, NumberOfBasis

    InterpolationType= match_basis(Basis%BasisType(1,BasisIdx))

    NumberGlobalXElements=str2int(Mesh%MeshNumberOfElements(1,1))
    NumberGlobalYElements=str2int(Mesh%MeshNumberOfElements(2,1))
    NumberGlobalZElements=str2int(Mesh%MeshNumberOfElements(3,1))
    NumberOfGaussXi(1) = str2int(Basis%BasisNumberOfGaussPoints(1,BasisIdx))
    NumberOfGaussXi(2) = str2int(Basis%BasisNumberOfGaussPoints(2,BasisIdx))
    NumberOfGaussXi(3) = str2int(Basis%BasisNumberOfGaussPoints(3,BasisIdx))

      !Define geometric basis

    CALL cmfe_Basis_Initialise(all_Basis%Basis(BasisIdx),Err)
    CALL cmfe_Basis_CreateStart(BasisUserNumber(BasisIdx),all_Basis%Basis(BasisIdx),Err)

    SELECT CASE(InterpolationType)

      CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
        & CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION )

        CALL cmfe_Basis_TypeSet(all_Basis%Basis(BasisIdx),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)

      CASE(CMFE_BASIS_QUADRATIC1_HERMITE_INTERPOLATION ,CMFE_BASIS_QUADRATIC2_HERMITE_INTERPOLATION, &
        & CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION )

        CALL cmfe_Basis_TypeSet(all_Basis%Basis(BasisIdx),CMFE_BASIS_SIMPLEX_TYPE,Err)

    END SELECT

    ! 1D case
    IF(NumberGlobalYElements==0 .AND. NumberGlobalZElements==0) THEN

      CALL cmfe_Basis_NumberOfXiSet(all_Basis%Basis(BasisIdx),1,Err)
      CALL cmfe_Basis_InterpolationXiSet(all_Basis%Basis(BasisIdx),[InterpolationType],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis%Basis(BasisIdx), &
        & [NumberOfGaussXi(1),NumberOfGaussXi(2),NumberOfGaussXi(3)],Err)

    END IF
      ! 2D case
    IF(NumberGlobalZElements==0) THEN

      CALL cmfe_Basis_NumberOfXiSet(all_Basis%Basis(BasisIdx),2,Err)
      CALL cmfe_Basis_InterpolationXiSet(all_Basis%Basis(BasisIdx),[InterpolationType,InterpolationType],Err)

      IF(NumberOfGaussXi(1)>0 .AND. NumberOfGaussXi(2)>0  ) THEN

        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis%Basis(BasisIdx),[NumberOfGaussXi(1),NumberOfGaussXi(2)],Err)

      END IF

    ELSE

      CALL cmfe_Basis_NumberOfXiSet(all_Basis%Basis(BasisIdx),3,Err)
      CALL cmfe_Basis_InterpolationXiSet(all_Basis%Basis(BasisIdx),[InterpolationType,InterpolationType,InterpolationType],Err)

      IF(NumberOfGaussXi(1)>0 .AND. NumberOfGaussXi(2)>0 .AND. NumberOfGaussXi(3)>0 ) THEN 
 
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis%Basis(BasisIdx), &
          & [NumberOfGaussXi(1),NumberOfGaussXi(2),NumberOfGaussXi(3)],Err)

      END IF
    END IF

    CALL cmfe_Basis_CreateFinish(all_Basis%Basis(BasisIdx),Err)


  END DO  ! BasisIdx 
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DEFINING MESH AND ASSIGNING RESPECTIVE MESHES TO THIER  REGIONS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO MeshIdx = 1,NumberOfMesh

    !! number of state variables ( for instance 3 ( displacements) +1(Pressure) for 3D cantiliverbeam study).
    NumberOfComponents                       = str2int(DependentField%DependentFieldNumberOfComponents(1,MeshIdx))
    Width                                    = str2int(Mesh%MeshGeometricParameters(1,1))
    Height                                   = str2int(Mesh%MeshGeometricParameters(2,1))
    Length                                   = str2int(Mesh%MeshGeometricParameters(3,1))


    !Start the creation of a generated mesh in the region.

    CALL cmfe_GeneratedMesh_Initialise(all_GeneratedMesh%GeneratedMesh(MeshIdx),Err)
    
    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber(MeshIdx),all_Region%Region(MeshIdx), &
     & all_GeneratedMesh%GeneratedMesh(MeshIdx),Err)

    !Set up a regular x*y*z mesh

    CALL cmfe_GeneratedMesh_TypeSet(all_GeneratedMesh%GeneratedMesh(MeshIdx),CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)

    !Set the default basis

    IF(UsePressureBasis) THEN

      CALL cmfe_GeneratedMesh_BasisSet(all_GeneratedMesh%GeneratedMesh(MeshIdx),[all_Basis%Basis(MeshIdx), &
        & all_Basis%Basis(MeshIdx+1)],Err)

    ELSE

      CALL cmfe_GeneratedMesh_BasisSet(all_GeneratedMesh%GeneratedMesh(MeshIdx),[all_Basis%Basis(MeshIdx)],Err)

    END IF

    ! 1D case
    IF(NumberGlobalYElements==0 .AND. NumberGlobalZElements==0) THEN

      CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh%GeneratedMesh(MeshIdx),[Width],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh%GeneratedMesh(MeshIdx), &
       & [NumberGlobalXElements],Err)

    ELSEIF (NumberGlobalZElements==0) THEN

      CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh%GeneratedMesh(MeshIdx),[Width,Height],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh%GeneratedMesh(MeshIdx), &
        & [NumberGlobalXElements,NumberGlobalYElements],Err)

    ELSE

      CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh%GeneratedMesh(MeshIdx),[Width,Height,Length],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh%GeneratedMesh(MeshIdx), &
       & [NumberGlobalXElements,NumberGlobalYElements, NumberGlobalZElements],Err)

    END IF
    !Finish the creation of a generated mesh in the region
    CALL cmfe_Mesh_Initialise(all_Mesh%Mesh(MeshIdx),Err)

    CALL cmfe_GeneratedMesh_CreateFinish(all_GeneratedMesh%GeneratedMesh(MeshIdx),MeshUserNumber(MeshIdx), &
      & all_Mesh%Mesh(MeshIdx),Err)

  END DO !MeshIdx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DECOMPOSING BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO DecompositionIdx = 1,NumberOfDecomposition
    !Create a decomposition

    CALL cmfe_Decomposition_Initialise(all_Decomposition%Decomposition(DecompositionIdx),Err)
    CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber(DecompositionIdx), &
      & all_Mesh%Mesh(DecompositionIdx),all_Decomposition%Decomposition(DecompositionIdx),Err)
    CALL cmfe_Decomposition_TypeSet(all_Decomposition%Decomposition(DecompositionIdx),CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(all_Decomposition%Decomposition(DecompositionIdx),NumberOfDomains,Err)

    IF (Decomposition%DecompositionFaceActive(1,1) == "ACTIVE") then 
      CALL cmfe_Decomposition_CalculateFacesSet(all_Decomposition%Decomposition(DecompositionIdx),.TRUE.,Err)
    END IF

    CALL cmfe_Decomposition_CreateFinish(all_Decomposition%Decomposition(DecompositionIdx),Err)

  END DO !DecompositionIdx



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   GEOMETRIC  BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO GeometricFieldIdx  = 1, NumberOfGeometricField
    !Create a field to put the geometry (default is geometry)

    CALL cmfe_Field_Initialise(all_GeometricField%GeometricField(GeometricFieldIdx),Err)
    CALL cmfe_Field_CreateStart(FieldGeometryUserNumber(GeometricFieldIdx),all_Region%Region(GeometricFieldIdx), &
      & all_GeometricField%GeometricField(GeometricFieldIdx),Err)
    CALL cmfe_Field_MeshDecompositionSet(all_GeometricField%GeometricField(GeometricFieldIdx), &
      & all_Decomposition%Decomposition(GeometricFieldIdx),Err)
    CALL cmfe_Field_VariableLabelSet(all_GeometricField%GeometricField(GeometricFieldIdx), &
      & CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
    CALL cmfe_Field_ScalingTypeSet(all_GeometricField%GeometricField(GeometricFieldIdx), &
      & CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    CALL cmfe_Field_CreateFinish(all_GeometricField%GeometricField(GeometricFieldIdx),Err)

    !Update the geometric field parameters
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(all_GeneratedMesh%GeneratedMesh(GeometricFieldIdx), &
      & all_GeometricField%GeometricField(GeometricFieldIdx),Err)

  END DO ! GeometricFieldIdx  !!!!! END LOOP FOR GEOMETRIC FIELD

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FIBER FIELD  BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO FiberFieldIdx = 1, NumberOfFiberField

    CALL cmfe_Field_Initialise(all_FibreField%FibreField(FiberFieldIdx),Err)
    CALL cmfe_Field_CreateStart(FieldFibreUserNumber(FiberFieldIdx),all_Region%Region(FiberFieldIdx), &
      & all_FibreField%FibreField(FiberFieldIdx),Err)
    CALL cmfe_Field_TypeSet(all_FibreField%FibreField(FiberFieldIdx),CMFE_FIELD_FIBRE_TYPE,Err)
    CALL cmfe_Field_MeshDecompositionSet(all_FibreField%FibreField(FiberFieldIdx),  &
      all_Decomposition%Decomposition(FiberFieldIdx),Err)
    CALL cmfe_Field_GeometricFieldSet(all_FibreField%FibreField(FiberFieldIdx), &
      all_GeometricField%GeometricField(FiberFieldIdx),Err)
    CALL cmfe_Field_VariableLabelSet(all_FibreField%FibreField(FiberFieldIdx),CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
    CALL cmfe_Field_ScalingTypeSet(all_FibreField%FibreField(FiberFieldIdx),CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    CALL cmfe_Field_CreateFinish(all_FibreField%FibreField(FiberFieldIdx),Err)

    !Set the fibre field
    CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(FiberFieldIdx), &
      & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
        & str2real(FibreField%FibreFieldParameters(1,NumberOfFiberField)),Err)

    IF  (NumberGlobalYElements .NE. 0) THEN
      CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(FiberFieldIdx), &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
          & str2real(FibreField%FibreFieldParameters(2,NumberOfFiberField)),Err)
    END IF

    IF ( NumberGlobalZElements .NE. 0) THEN
      CALL cmfe_Field_ComponentValuesInitialise(all_FibreField%FibreField(FiberFieldIdx), &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
          & str2real(FibreField%FibreFieldParameters(3,NumberOfFiberField)),Err)
    END IF

  END DO   !!!!! END LOOP FOR Fiber FIELD


  DO EquationSetIdx  = 1,NumberOfEquationsSet            !!!! loop for eqaution set

    CALL cmfe_EquationsSet_Initialise(all_EquationsSet%EquationsSet(EquationSetIdx),Err)
    CALL cmfe_Field_Initialise(all_EquationsSetField%EquationsSetField(EquationSetIdx),Err)
    IF (NumberOfFiberField .GE. 1) then
      CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber(EquationSetIdx),all_Region%Region(EquationSetIdx), &
        & all_FibreField%FibreField(EquationSetIdx), &
          & [match_equations_set(EquationsSet%EquationSetClass(1,EquationSetIdx)), &
            & match_equations_set(EquationsSet%EquationSetType(1,EquationSetIdx)), &
              & match_equations_set(EquationsSet%EquationSetSubType(1,EquationSetIdx))], &
                & EquationsSetFieldUserNumber(EquationSetIdx),all_EquationsSetField%EquationsSetField(EquationSetIdx), &
                  & all_EquationsSet%EquationsSet(EquationSetIdx),Err)

    ELSE


      CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber(EquationSetIdx),all_Region%Region(EquationSetIdx), &
        & all_GeometricField%GeometricField(EquationSetIdx), &
          & [match_equations_set(EquationsSet%EquationSetClass(1,EquationSetIdx)), &
            & match_equations_set(EquationsSet%EquationSetType(1,EquationSetIdx)), &
              & match_equations_set(EquationsSet%EquationSetSubType(1,EquationSetIdx))], &
                & EquationsSetFieldUserNumber(EquationSetIdx),all_EquationsSetField%EquationsSetField(EquationSetIdx), &
                  & all_EquationsSet%EquationsSet(EquationSetIdx),Err)
    END IF 
    CALL cmfe_EquationsSet_CreateFinish(all_EquationsSet%EquationsSet(EquationSetIdx),Err)

  END DO ! EquationSetIdx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Dependent  FIELD  BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Create the dependent field
  DO DependentFieldIdx = 1, NumberOfDependentField
    !Create the Dependent field

    CALL cmfe_Field_Initialise(all_DependentField%DependentField(DependentFieldIdx),Err)
    CALL cmfe_EquationsSet_DependentCreateStart(all_EquationsSet%EquationsSet(DependentFieldIdx), &
      & FieldDependentUserNumber(DependentFieldIdx), all_DependentField%DependentField(DependentFieldIdx),Err)
    CALL cmfe_Field_VariableLabelSet(all_DependentField%DependentField(DependentFieldIdx), &
      & CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)

    
    DO ComponentIdx=1,NumberOfComponents

      CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(DependentFieldIdx), &
        & CMFE_FIELD_U_VARIABLE_TYPE,ComponentIdx,1,Err)

      CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField%DependentField(DependentFieldIdx), &
        & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,ComponentIdx,1,Err)

    END DO


    IF(UsePressureBasis) THEN

      CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(DependentFieldIdx), &
        & CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
      CALL cmfe_Field_ComponentInterpolationSet(all_DependentField%DependentField(DependentFieldIdx), &
        & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
      CALL cmfe_Field_ComponentMeshComponentSet( &
        & all_DependentField%DependentField(DependentFieldIdx),CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
      CALL cmfe_Field_ComponentMeshComponentSet( &
        & all_DependentField%DependentField(DependentFieldIdx),CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)

    END IF

    CALL cmfe_EquationsSet_DependentCreateFinish(all_EquationsSet%EquationsSet(DependentFieldIdx),Err)

    IF  (TRIM(DependentField%DependentFieldInitialValueOfStateVector(1,DependentFieldIdx)) == "UNDEFORMED") THEN

      DO ComponentIdx=1,NumberOfComponents
        CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField%GeometricField(DependentFieldIdx), &
          & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
            & ComponentIdx,all_DependentField%DependentField(DependentFieldIdx), &
              & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,ComponentIdx,Err)
      END DO


    ELSE
      DO ComponentIdx=1,NumberOfComponents

        CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(DependentFieldIdx), &
          &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,ComponentIdx, &
            & str2real(DependentField%DependentFieldInitialValueOfStateVector(1,DependentFieldIdx)),Err)

       END DO 
    END IF  ! TRIM(DependentField_arguments(4,DependentFieldIdx)) == "UNDEFORMED")
 
    IF(UsePressureBasis) THEN

      CALL cmfe_Field_ComponentValuesInitialise(all_DependentField%DependentField(DependentFieldIdx), &
        &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
          & str2real(DependentField%DependentFieldInitialValueOfStateScalar(1,DependentFieldIdx)),Err)

    END IF

    CALL cmfe_Field_ParameterSetUpdateStart(all_DependentField%DependentField(DependentFieldIdx),CMFE_FIELD_U_VARIABLE_TYPE &
      & ,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(all_DependentField%DependentField(DependentFieldIdx),CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,Err)

  END DO ! DependentFieldIdxx


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   MATERIAL FIELD  BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO MaterialFieldIdx  = 1,NumberOfMaterialField

    CALL cmfe_Field_Initialise(all_MaterialField%MaterialField(MaterialFieldIdx),Err)


    CALL cmfe_EquationsSet_MaterialsCreateStart(all_EquationsSet%EquationsSet(MaterialFieldIdx), &
        & FieldMaterialUserNumber(EquationSetIdx),all_MaterialField%MaterialField(MaterialFieldIdx),Err)
    CALL cmfe_Field_VariableLabelSet(all_MaterialField%MaterialField(MaterialFieldIdx), &
      & CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
    CALL cmfe_EquationsSet_MaterialsCreateFinish(all_EquationsSet%EquationsSet(MaterialFieldIdx),Err)

    DO  MaterialParametersIdx = 1,material_parameters(EquationsSet%EquationSetSubType(1,MaterialFieldIdx))

      CALL  cmfe_Field_ComponentValuesInitialise(all_MaterialField%MaterialField(MaterialFieldIdx), &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, MaterialParametersIdx , &
          & str2real(MaterialField%MaterialFieldParameters(MaterialParametersIdx,NumberOfMaterialField)),Err)

    END DO

  END DO  ! MaterialFieldIdx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   EQUATION SET BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

  DO EquationSetIdx  = 1,NumberOfEquationsSet            !!!! loop for eqaution set

    !Create the equations set equations
    CALL cmfe_Equations_Initialise(all_Equations%Equations(EquationSetIdx),Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(all_EquationsSet%EquationsSet(EquationSetIdx), &
      & all_Equations%Equations(EquationSetIdx),Err)
    CALL cmfe_Equations_SparsityTypeSet(all_Equations%Equations(EquationSetIdx),CMFE_EQUATIONS_SPARSE_MATRICES,Err)
    CALL cmfe_Equations_OutputTypeSet(all_Equations%Equations(EquationSetIdx),output_type(Output%OutputType(1,EquationSetIdx)),Err)
    CALL cmfe_EquationsSet_EquationsCreateFinish(all_EquationsSet%EquationsSet(EquationSetIdx),Err)


    IF (SourceField%SourceFieldType(1,EquationSetIdx) == "GRAVITY" ) THEN   !! look input file IF the given field is gravitational
      !Create the source field with the gravity vector
      CALL cmfe_Field_Initialise(all_SourceField%SourceField(EquationSetIdx),Err)
      CALL cmfe_EquationsSet_SourceCreateStart(all_EquationsSet%EquationsSet(EquationSetIdx), &
            & FieldUserNumber(NumberOfSourceField) ,all_SourceField%SourceField(NumberOfSourceField),Err)
      CALL cmfe_Field_ScalingTypeSet(all_SourceField%SourceField(NumberOfSourceField),CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
      CALL cmfe_EquationsSet_SourceCreateFinish(all_EquationsSet%EquationsSet(EquationSetIdx),Err)

      DO ComponentIdx=1,NumberOfComponents

        CALL cmfe_Field_ComponentValuesInitialise(all_SourceField%SourceField(NumberOfSourceField),  &
          & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,ComponentIdx, &
            & str2real(SourceField%SourceFieldComponents(ComponentIdx,EquationSetIdx)),Err)
            
      END DO

    END IF
  END DO ! EquationsSetIdx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   PROBLEM BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  DO ProblemIdx = 1, NumberOfProblem

    !Define the problem
    CALL cmfe_Problem_Initialise(all_Problem%Problem(ProblemIdx),Err)
    CALL cmfe_Problem_CreateStart(ProblemUserNumber(ProblemIdx),[match_problem(Problem%ProblemClass(1,ProblemIdx)), &
      & match_problem(Problem%ProblemType(1,ProblemIdx)),match_problem(Problem%ProblemSubType(1,ProblemIdx))], &
        & all_Problem%Problem(ProblemIdx),Err)
    CALL cmfe_Problem_CreateFinish(all_Problem%Problem(ProblemIdx),Err)

  END DO ! ProblemIdx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   CONTROL LOOP BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO ControlLoopIdx = 1,NumberOfControlLoop 
    !Create the problem control loop
    CALL cmfe_Problem_ControlLoopCreateStart(all_Problem%Problem(ControlLoopIdx),Err)
    CALL cmfe_ControlLoop_Initialise(all_ControlLoop%ControlLoop(ControlLoopIdx),Err)
    CALL cmfe_Problem_ControlLoopGet(all_Problem%Problem(ControlLoopIdx),CMFE_CONTROL_LOOP_NODE, &
      & all_ControlLoop%ControlLoop(ControlLoopIdx),Err)
    CALL cmfe_ControlLoop_TypeSet(all_ControlLoop%ControlLoop(ControlLoopIdx), &
      & control_loop_def(ControlLoop%ControlLoopType(1,ControlLoopIdx)),Err)
   !CALL cmfe_ControlLoop_TimesSet(ControlLoop,StartTime,StopTime,TimeIncrement,Err)
   !CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,1,Err)
    CALL cmfe_ControlLoop_LoaDOutputSet(all_ControlLoop%ControlLoop(ControlLoopIdx),1,Err)
    CALL cmfe_ControlLoop_MaximumIterationsSet(all_ControlLoop%ControlLoop(ControlLoopIdx), &
      & str2int(ControlLoop%ControlLoopLoadIncrement(1,ControlLoopIdx)),Err)
    CALL cmfe_Problem_ControlLoopCreateFinish(all_Problem%Problem(ControlLoopIdx),Err)

  END DO ! ControlLoopIdx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SOLVER BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO solverIdx = 1,NumberOfSolver

    !Create the problem solvers
    CALL cmfe_Solver_Initialise(all_Solver%Solver(solverIdx),Err)
    CALL cmfe_Solver_Initialise(all_LinearSolver%LinearSolver(solverIdx),Err)
    CALL cmfe_Problem_SolversCreateStart(all_Problem%Problem(solverIdx),Err)
    CALL cmfe_Problem_SolverGet(all_Problem%Problem(solverIdx),CMFE_CONTROL_LOOP_NODE,1,all_Solver%Solver(solverIdx),Err)
    CALL cmfe_Solver_OutputTypeSet(all_Solver%Solver(solverIdx),CMFE_SOLVER_PROGRESS_OUTPUT,Err)
    IF (TRIM(Solver%NonlinearSolver(1,solverIdx)) == "ACTIVE") THEN
      CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(all_Solver%Solver(solverIdx), &
                                                         & solver_def(Solver%NewtonJacobianType(1,solverIdx)),Err)
      CALL cmfe_Solver_NewtonLinearSolverGet(all_Solver%Solver(solverIdx),all_LinearSolver%LinearSolver(solverIdx),Err)
      CALL cmfe_Solver_LinearTypeSet(all_LinearSolver%LinearSolver(solverIdx),solver_def(Solver%SolverType(1,solverIdx)),Err)
      CALL cmfe_Solver_NewtonRelativeToleranceSet(all_Solver%Solver(solverIdx), &
                                                   & str2real(Solver%NewtonTolerance(1,solverIdx)),Err)
      CALL cmfe_Solver_NewtonAbsoluteToleranceSet(all_Solver%Solver(solverIdx), &
                                                    & str2real(Solver%NewtonTolerance(1,solverIdx)),Err)
      CALL cmfe_Solver_NewtonMaximumIterationsSet(all_Solver%Solver(solverIdx),str2int(Solver%NewtonMaxIterations(1,solverIdx)),Err)
    ENDIF

    CALL cmfe_Problem_SolversCreateFinish(all_Problem%Problem(solverIdx),Err)

    !Create the problem solver equations
    CALL cmfe_Solver_Initialise(all_Solver%Solver(solverIdx),Err)
    CALL cmfe_SolverEquations_Initialise(all_SolverEquations%SolverEquations(solverIdx),Err)
    CALL cmfe_Problem_SolverEquationsCreateStart(all_Problem%Problem(solverIdx),Err)
    CALL cmfe_Problem_SolverGet(all_Problem%Problem(solverIdx),CMFE_CONTROL_LOOP_NODE,1,all_Solver%Solver(solverIdx),Err)
    CALL cmfe_Solver_SolverEquationsGet(all_Solver%Solver(solverIdx),all_SolverEquations%SolverEquations(solverIdx),Err)
    CALL cmfe_SolverEquations_EquationsSetAdd(all_SolverEquations%SolverEquations(solverIdx), &
                                               & all_EquationsSet%EquationsSet(solverIdx), EquationsSetIndex,Err)
    CALL cmfe_Problem_SolverEquationsCreateFinish(all_Problem%Problem(solverIdx),Err)

  END DO ! solverIdx 


  DO BOundaryConditionIdx = 1,NumberOfBoundaryCondition

    !Prescribe boundary conditions (absolute nodal parameters)
    CALL cmfe_BoundaryConditions_Initialise(all_BoundaryConditions%BoundaryConditions(BOundaryConditionIdx),Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(all_SolverEquations%SolverEquations(BOundaryConditionIdx), &
       & all_BoundaryConditions%BoundaryConditions(BOundaryConditionIdx),Err)

    DO i = 1,NumberOfDirichelet

      CALL cmfe_GeneratedMesh_SurfaceGet(all_GeneratedMesh%GeneratedMesh(BOundaryConditionIdx), &
        & bc_def(BoundaryDirichelet(2+4*(i-1),BOundaryConditionIdx)),SurfaceNodes,LeftNormalXi,Err)

      DO nodeIdx=1,SIZE(SurfaceNodes,1)
        NodeNumber=SurfaceNodes(nodeIdx)

        CALL cmfe_Decomposition_NodeDOmainGet(all_Decomposition%Decomposition(BOundaryConditionIdx),NodeNumber,1,NodeDOmain,Err)

        IF(NodeDOmain==ComputationalNodeNumber) THEN


          DO ComponentIdx=1,NumberOfComponents

            Constraint =  BoundaryDirichelet(3+4*(i-1),BOundaryConditionIdx)

            IF (TRIM(Constraint(ComponentIdx:ComponentIdx)) == "1") THEN

              CALL cmfe_BoundaryConditions_SetNode(all_BoundaryConditions%BoundaryConditions(BOundaryConditionIdx), &
                & all_DependentField%DependentField(BOundaryConditionIdx), &
                  & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, ComponentIdx,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
                    & str2real(BoundaryDirichelet(4+4*(i-1),BOundaryConditionIdx)),Err)

            END IF  ! (TRIM(Constraint(ComponentIdx:ComponentIdx)) == "1")
          END DO  ! ComponentIdx=1,3
        END IF ! NodeDOmain==ComputationalNodeNumber
      END DO  !! nodeIdx=1,SIZE(SurfaceNodes,1)

      DEALLOCATE(SurfaceNodes)

    END DO  !! i = 1,NumberOfDirichelet

  !!__________________________________________________________________________________________________________________________!!!
  ! Will INCLUDE the pressure Neumann later.

  !#DO nodeIdx=1,SIZE(SurfaceNodes,1)

  ! NodeNumber=SurfaceNodes(nodeIdx)

  ! CALL cmfe_Decomposition_NodeDOmainGet(all_Decomposition%Decomposition(1),NodeNumber,1,NodeDOmain,Err)


  ! IF(NodeDOmain==ComputationalNodeNumber) THEN


    !    CALL cmfe_BoundaryConditions_SetNode(all_BoundaryConditions%BoundaryConditions(1),all_DependentField%DependentField(1), &
    !      & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber, &
    !      & ABS(LeftNormalXi), &
    !      & CMFE_BOUNDARY_CONDITION_PRESSURE,2.1_CMISSRP,Err)


  !  END IF
  !END DO
  !!__________________________________________________________________________________________________________________________!!!
    CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(all_SolverEquations%SolverEquations(BOundaryConditionIdx),Err)
  END DO ! BOundaryConditionIdx


  DO ProblemIdx = 1, NumberOfProblem

    !Solve problem
    CALL cmfe_Problem_Solve(all_Problem%Problem(ProblemIdx),Err)

  END DO  ! ProblemIdx

  !Output solution
  DO FieldIdx = 1,1

    CALL cmfe_Fields_Initialise(all_Fields%Fields(FieldIdx),Err)
    CALL cmfe_Fields_Create(all_Region%Region(FieldIdx),all_Fields%Fields(FieldIdx),Err)
    CALL cmfe_Fields_NodesExport(all_Fields%Fields(FieldIdx),"LargeUniaxialExtension","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(all_Fields%Fields(FieldIdx),"LargeUniaxialExtension","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(all_Fields%Fields(FieldIdx),Err)

  END DO

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

  CONTAINS

  INCLUDE "problems_match_library.f90"

END PROGRAM GENERICEXAMPLE
