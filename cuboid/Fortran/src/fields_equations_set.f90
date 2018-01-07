
MODULE FIELDS_EQUATIONS_SET

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE OpenCMISS_Variables
  USE REGION_MESH
  USE DECOMPOSITION 
  
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

SUBROUTINE CreateFieldFiniteElasticity()

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for finite elasticity - quadratic interpolation
  CALL cmfe_Field_Initialise(GeometricFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberFE,RegionFE,GeometricFieldFE,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldFE,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldFE,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldFE,Err)

!  CALL cmfe_Field_ParameterSetUpdateStart(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
!  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMeshFE,GeometricFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a fibre field and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,RegionFE,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a material field for Finite Elasticity and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(MaterialFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumberFE,RegionFE,MaterialFieldFE,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldFE,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldFE,FieldMaterialNumberOfVariablesFE,Err)
  CALL cmfe_Field_VariableTypesSet(MaterialFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE],Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,8,Err)
  ELSEIF (ModelType == 2) THEN  ! 4, "Titin"
    CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,7,Err)
  ENDIF

  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1,Err)
   
  IF (ModelType == 0 .OR. ModelType == 1) THEN ! "MultiPhysStrain"
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,7,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  ENDIF

  IF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,7,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,8,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  ENDIF

  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"MaterialFE",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,"Gravity",Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldFE,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    !Set Mooney-Rivlin constants c10 and c01.
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,MAT_FE(1),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,MAT_FE(2),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,MAT_FE(3),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,MAT_FE(4),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0.0_CMISSRP, &
     & Err)
     
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable

    !Set Material-Parameters [mu(1) mu(2) mu(3) alpha(1) alpha(2) alpha(3) mu_0]
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,C(1),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,C(2),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,C(3),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,C(4),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,C(5),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,C(6),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7,C(7),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,8,C(8),Err)
 
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    !Set Mooney-Rivlin constants c10 and c01.
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,MAT_FE(1),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,MAT_FE(2),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,MAT_FE(3),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,MAT_FE(4),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0.0_CMISSRP, &
     & Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,PMax,Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7,TKLinParam,&
    & Err)
    
  ENDIF



  !Create the dependent field for FE with 2 variables and * components
  !  3-d: 3 displacement (quad interpol), 1 pressure (lin interpol) --> * = 4
  !  2-d: 2 displacement (quad interpol), 1 pressure (lin interpol) --> * = 3
  !  1-d: 1 displacement (quad interpol), 1 pressure (lin interpol) --> * = 2
  CALL cmfe_Field_Initialise(DependentFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumberFE,RegionFE,DependentFieldFE,Err)
  CALL cmfe_Field_TypeSet(DependentFieldFE,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldFE,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldFE,FieldDependentNumberOfVariablesFE,Err)

  CALL cmfe_Field_VariableTypesSet(DependentFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
!  CALL cmfe_Field_ScalingTypeSet(DependentFieldFE,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"DependentFE",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"Reaction_Force",Err)
  CALL cmfe_Field_CreateFinish(DependentFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the mechanics mesh
  CALL cmfe_Field_Initialise(IndependentFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldIndependentUserNumberFE,RegionFE,IndependentFieldFE,Err)
  CALL cmfe_Field_TypeSet(IndependentFieldFE,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(IndependentFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(IndependentFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_DependentTypeSet(IndependentFieldFE,CMFE_FIELD_INDEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(IndependentFieldFE,2,Err)
  CALL cmfe_Field_VariableTypesSet(IndependentFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE],Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_DimensionSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_DimensionSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6,Err)
  ENDIF

  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,5,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1, &
       & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !A_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !A_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !x_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !x_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !lambda_a
  
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force (unbound)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force (bound)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force in XF-direction (unbound)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force in XF-direction (bound)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! activation for titin
    CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  ENDIF

  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) ! number of nodes in XI1
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,2, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) ! number of nodes in XI2
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,3, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) ! number of nodes in XI3
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,4, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) ! if fibre starts in current FE element
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,5, &
    & CMFE_FIELD_CONSTANT_INTERPOLATION,Err) ! number of in series fibres
  CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_DataTypeSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DataTypeSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_INTG_TYPE,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"Active_Stress_FE",Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"XB_state_variables_FE",Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"Active_Stress_FE",Err)
  ENDIF

  CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,"subgrid_info",Err)
  CALL cmfe_Field_CreateFinish(IndependentFieldFE,Err)

END SUBROUTINE CreateFieldFiniteElasticity

SUBROUTINE CreateFieldMonodomain()

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for monodomain - quadratic interpolation
  CALL cmfe_Field_Initialise(GeometricFieldM,Err)
  ! FIELD_CREATE_START(fieldUserNumber,REGION,FIELD)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberM,RegionM,GeometricFieldM,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldM,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldM,DecompositionM,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldM,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldM,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"GeometryM",Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldM,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a materials field for monodomain and attach it to the geometric field - constant interpolation
  CALL cmfe_Field_Initialise(MaterialFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumberM,RegionM,MaterialFieldM,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldM,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldM,DecompositionM,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldM,FieldMaterialNumberOfVariablesM,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponentsM,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"MaterialM",Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldM,Err)
  !Set Am
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Am,Err)
  !Set Cm
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,CmFast,Err)
  !Set Conductivity
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & Conductivity,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for monodomain with 3 variables. [U: 1 component, DelUdelN: 1 component, V: 3]
  !DependentFieldM: FIELD_U_VARIABLE_TYPE: 1) Vm, FIELD_DELUDELN_VARIABLE_TYPE: 1)dVm/dn, FIELD_V_VARIABLE_TYPE: 1),2),3): GeometryM3D, 3D-position of geometry
  CALL cmfe_Field_Initialise(DependentFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumberM,RegionM,DependentFieldM,Err)
  CALL cmfe_Field_TypeSet(DependentFieldM,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldM,DecompositionM,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldM,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldM,FieldDependentNumberOfVariablesM,Err)
  CALL cmfe_Field_VariableTypesSet(DependentFieldM,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
   & CMFE_FIELD_V_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,Err)
  !additional the V_Var_Type with 3 components for the 3-d position of the geometry
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,3,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM, & 
    & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
!  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
  CALL cmfe_Field_DimensionSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_DimensionSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"Vm",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dn",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,"GeometryM3D",Err)
  CALL cmfe_Field_CreateFinish(DependentFieldM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the bioelectrics mesh
  CALL cmfe_Field_Initialise(IndependentFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldIndependentUserNumberM,RegionM,IndependentFieldM,Err)
  CALL cmfe_Field_TypeSet(IndependentFieldM,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(IndependentFieldM,DecompositionM,Err)
  CALL cmfe_Field_GeometricFieldSet(IndependentFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_DependentTypeSet(IndependentFieldM,CMFE_FIELD_INDEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(IndependentFieldM,4,Err)
  CALL cmfe_Field_VariableTypesSet(IndependentFieldM,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE, &
   & CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_U2_VARIABLE_TYPE],Err)

  !first variable:   CMFE_FIELD_U_VARIABLE_TYPE
  CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    !first variable:   CMFE_FIELD_U_VARIABLE_TYPE -- 1) active stress
    CALL cmfe_Field_DimensionSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"Active_Stress_M",Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !A_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !A_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !x_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !x_2
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"XB_state_variables_M",Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_DimensionSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,6,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !normalised sarcomere-based active stress
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the unbound titin stress
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the bound titin stress
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the unbound titin stress in the XF-direction
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,5, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the bound titin stress in the XF-direction
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,6, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the titin activation (stress without force-length relation)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"Active_Stress_M",Err)
  ENDIF

  !second variable:   CMFE_FIELD_V_VARIABLE_TYPE -- 1) motor unit number   2) fibre type   3) fibre number   4) nearest Gauss point   5) containing element local number 6) in-fibre contiguous node number
  CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_INTG_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM, & 
    & CMFE_FIELD_V_VARIABLE_TYPE,FieldIndependentNumberOfComponentsM2,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,1, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,2, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,3, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,4, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,5, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,6, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,"fibre_info",Err)

  !third variable:   FIELD_U1_VARIABLE_TYPE -- 1) half-sarcomere length   2) inital half-sarcomere length   3) initial node distance
  CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)

IF (ModelType == 0 .OR. ModelType == 1) THEN ! 3,3a MultiPhysStrain
  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,3,Err)
ELSEIF (ModelType == 2) THEN ! 4, "Titin"
  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,4,Err)
ENDIF
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,1, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,2, &
   & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,3, &
   & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)

  IF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  ENDIF

  CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,"half-sarcomere_length",Err)

  !fourth variable:   FIELD_U2_VARIABLE_TYPE -- 1) old node distance   2) maximum contraction velocity   3) relative contraction velocity   4) velocity before 1 time step   5) velocity before 2 time step   6) velocity before 3 time steps   7) node distance to right node (only computed on some nodes)
  CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,7,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,1, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,2, &
   & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,3, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,4, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,5, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,6, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,7, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,"contraction_velocity",Err)

  CALL cmfe_Field_CreateFinish(IndependentFieldM,Err)

END SUBROUTINE CreateFieldMonodomain

SUBROUTINE CreateEquationsSet()
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for Finite Elasticity

  CALL cmfe_Field_Initialise(EquationsSetFieldFE,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetFE,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberFE,RegionFE,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
     & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE], &
     & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberFE,RegionFE,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
     & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE], &
     & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberFE,RegionFE,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
     & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE], & 
     & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  ENDIF
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetFE,Err)

  !Create the equations set dependent field variables for Finite Elasticity
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetFE,FieldDependentUserNumberFE,DependentFieldFE,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetFE,FieldIndependentUserNumberFE,IndependentFieldFE,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set materials field variables for Finite Elasticity
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetFE,FieldMaterialUserNumberFE,MaterialFieldFE,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for monodomain
  CALL cmfe_Field_Initialise(EquationsSetFieldM,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetM,Err)
  !Set the equations set to be a Monodomain equations set
  !> \todo solve the monodomain problem on the fibre field rather than on the geometric field: GeometricField <--> FibreField

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberM,RegionM,GeometricFieldM, &
     & [CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
     & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE], &
     & FieldEquationsSetUserNumberM,EquationsSetFieldM,EquationsSetM,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberM,RegionM,GeometricFieldM, &
     & [CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
     & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE], &
     & FieldEquationsSetUserNumberM,EquationsSetFieldM,EquationsSetM,Err)
  
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberM,RegionM,GeometricFieldM, &
     & [CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
     & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE], &
     & FieldEquationsSetUserNumberM,EquationsSetFieldM,EquationsSetM,Err)
  ENDIF
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetM,Err)

  !Create the equations set dependent field variables for monodomain
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetM,FieldDependentUserNumberM,DependentFieldM,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetM,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetM,FieldIndependentUserNumberM,IndependentFieldM,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetM,Err)

  !Create the equations set materials field variables for monodomain
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetM,FieldMaterialUserNumberM,MaterialFieldM,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetM,Err)
    
END SUBROUTINE CreateEquationsSet

SUBROUTINE InitializeFieldMonodomain()
  INTEGER(CMISSIntg) :: FibreLineNo
  INTEGER(CMISSIntg) :: FibreInLineNo
  INTEGER(CMISSIntg) :: FibreNo, NodeNo, NodeIdx
  INTEGER(CMISSIntg) :: MotorUnitRank
  INTEGER(CMISSIntg) :: FEElementGlobalNumber, FEElementLocalNumber
  INTEGER(CMISSIntg) :: FEElementXIdx, FEElementYIdx, FEElementZIdx
  INTEGER(CMISSIntg) :: NodeMUserNumber
  INTEGER(CMISSIntg) :: MElementUserNumber, ElementIdx
  INTEGER(CMISSIntg) :: ElementDomain, NodeDomain
  INTEGER(CMISSIntg), DIMENSION(2) :: ElementUserNodes
  REAL(CMISSDP) :: InitialBioelectricNodeDistance

  !UPDATE THE INDEPENDENT FIELD IndependentFieldM
  ! OldTomoMechanics (ModelType 0)
  !first variable (U)
  !  components:
  !    1) active stress
  !
  ! .NOT. OldTomoMechanics (ModelType 1)
  !first variable (U)
  !  components:
  !    1) A_1
  !    2) A_2
  !    3) x_1
  !    4) x_2
  !
  !second variable (V)
  !  components:
  !    1) motor unit number
  !    2) fibre type
  !    3) fibre number
  !    4) nearest Gauss point
  !    5) FE element local number that contains bioelectric node
  !
  !init the motor unit number to 101
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, & 
   & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,101,Err) !

  ! initialize values
  !init the fibre type to 1
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, & 
    & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,1,Err) !Ftype=1
  !init the fibre number, the nearest Gauss point info and the inElem info to 0
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, & 
    & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, &
    & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, & 
    & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0,Err)  ! local element number
  !third variable (U1):
  !  components:
  !    1) half-sarcomere length
  !    2) initial half-sarcomere length
  !    3) initial node distance
  InitialBioelectricNodeDistance = PhysicalLength / (NumberOfElementsMPerFibreLine)
  !PRINT *, "InitialBioelectricNodeDistance = ", InitialBioelectricNodeDistance
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & 1.0_CMISSRP,Err) ! lengths in the cell model are in /micro/meters!!!
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & InitialBioelectricNodeDistance, Err)  !lengths in the cell model are in /micro/meters!!!
  !fourth variable (U2):
  !  components:
  !    1) old node distance
  !    2) maximum contraction velocity
  !    3) relative contraction velocity
  !    4) node distance to previous node before 1 time step   
  !    5) node distance to previous node before 2 time steps
  !    6) node distance to previous node before 3 time steps
  !    7) node distance to right node (only computed on some nodes)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
   & InitialBioelectricNodeDistance,Err)  !lengths in the cell model are in /micro/meters!!!
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & Vmax / (NumberOfElementsMPerFibreLine),Err)    !velocity in the cell model is in micro/meters/millisecond!!!
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSRP,Err)    !velocity in the cell model is in /micro/meters/per/millisecond!!!
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
   & InitialBioelectricNodeDistance,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
   & InitialBioelectricNodeDistance,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6, &
   & InitialBioelectricNodeDistance,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7, &
   & 0.0_CMISSRP,Err)
     
  ! Store containing FE element local numbers for bioelectric nodes in IndependentFieldM V comp. 5
  ! loop over fibres and bioelectric nodes
  DO FibreNo = 1, NumberOfFibres
    
    ! get rank of fibre
    MotorUnitRank = MUDistribution(MOD(FibreNo-1, 4050)+1)
    
    ! loop over elements of fibre, for each element, the right node is considered, except for the first element (on ElementIdx=0), then the left node of the first element is considered
    DO ElementIdx = 0, NumberOfElementsMPerFibre
      
      ! compute the bioelectric element user number
      IF (ElementIdx == 0) THEN   ! index 0 and 1 are both for the first element, but first for the left node, then for the right node
        MElementUserNumber = (FibreNo-1) * NumberOfElementsMPerFibre + 1
      ELSE
        MElementUserNumber = (FibreNo-1) * NumberOfElementsMPerFibre + ElementIdx
      ENDIF
      
      ! get the nodes of the current element
      CALL cmfe_MeshElements_NodesGet(ElementsM, MElementUserNumber, ElementUserNodes, Err)
      ! ElementUserNodes(1:2): node user number
          
      IF (ElementIdx == 0) THEN
        NodeMUserNumber = ElementUserNodes(1)   ! take left node of first element in first iteration per fibre
      ELSE
        NodeMUserNumber = ElementUserNodes(2)
      ENDIF
      
      !PRINT "(I1.1,3(A,I3))", ComputationalNodeNumber,": Fibre", FibreNo, ", Element in fibre ", ElementIdx, &
      !  & ", NodeMUserNumber:", NodeMUserNumber        
      
      ! get domain number of element and node
      !                                        decomposition,  elementUserNumber, domain
      CALL cmfe_Decomposition_ElementDomainGet(DecompositionM, MElementUserNumber, ElementDomain, Err)
      
      !                                     decomposition,  nodeUserNumber,  meshComponentNumber, domain
      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM, NodeMUserNumber, 1,                   NodeDomain, Err)
      
      IF (NodeDomain == ComputationalNodeNumber) THEN
      
        ! compute number of containing FE element
        FEElementZIdx = INT(INT((FibreNo-1) / NumberGlobalYFibres) / NumberOfNodesInXi3)
        FEElementYIdx = INT(MOD((FibreNo-1), NumberGlobalYFibres)  / NumberOfNodesInXi2)
        FEElementXIdx = (ElementIdx-1) / NumberOfNodesInXi1
        IF (ElementIdx == 0) THEN
          FEElementXIdx = 0
        ENDIF
      
        FEElementGlobalNumber = FEElementZIdx * NumberGlobalYElements * NumberGlobalXElements &
         & + FEElementYIdx * NumberGlobalXElements + FEElementXIdx + 1
      
        ! get local FEElementLocalNumber from FEElementGlobalNumber
        CALL cmfe_BioelectricFiniteElasticity_GetLocalElementNumber(GeometricFieldFE, FEElementGlobalNumber, FEElementLocalNumber, &
          & Err)
      
        !PRINT "(I1.1,2(A,I2.2),6(A,I3.3))", ComputationalNodeNumber,": Fibre", FibreNo, ", Element in fibre ", ElementIdx, &
        !  & ", Element (",FEElementXIdx,",",FEElementYIdx,",",FEElementZIdx, &
        !  & ")=", FEElementGlobalNumber, ", MU ",MotorUnitRank, ", local FE ", FEElementLocalNumber
      
        ! set number of containing FE element
        !                                     FIELD,              VARIABLE_TYPE,              
        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM, CMFE_FIELD_V_VARIABLE_TYPE, &
        !   FIELD_SET_TYPE,             VERSION_NO, DERIVATIVE_NO, USER_NODE_NUMBER, COMPONENT_NUMBER, VALUE, ERR
          & CMFE_FIELD_VALUES_SET_TYPE, 1,          1,             NodeMUserNumber,  5,                FEElementLocalNumber, Err)
          
      
        ! set motor unit number
        !                                      FIELD,             VARIABLE_TYPE               
        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM, CMFE_FIELD_V_VARIABLE_TYPE, &
          ! FIELD_SET_TYPE              VERSION_NO, DERIVATIVE_NO, USER_NODE_NUMBER, COMPONENT_NUMBER, VALUE
          & CMFE_FIELD_VALUES_SET_TYPE, 1,          1,             NodeMUserNumber,  1,                MotorUnitRank, Err)
        
        ! set fibre number
        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM, CMFE_FIELD_V_VARIABLE_TYPE, &
          ! FIELD_SET_TYPE              VERSION_NO, DERIVATIVE_NO, USER_NODE_NUMBER, COMPONENT_NUMBER, VALUE
          & CMFE_FIELD_VALUES_SET_TYPE, 1,          1,             NodeMUserNumber,  3,                FibreNo, Err)
        
      ENDIF
    ENDDO
  ENDDO
  
END SUBROUTINE InitializeFieldMonodomain

SUBROUTINE InitializeFieldFiniteElasticity()
  INTEGER(CMISSIntg) :: ElementUserNumber
  !UPDATE THE INDEPENDENT FIELD IndependentFieldFE
  !second variable of IndependentFieldFE
  !  components:
  !    1) number of nodes in Xi(1) direction per element
  !    2) number of nodes in Xi(2) direction per element
  !    3) number of nodes in Xi(3) direction per element
  !    4) beginning of fibres in this FE element? 1=yes, 0=no
  !    5) number of in series fibres
  !
  !initialise as if the fibres would not start in any element, and adjust below
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
   & NumberOfNodesInXi1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & NumberOfNodesInXi2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & NumberOfNodesInXi3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
   & 0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
   & NumberOfInSeriesFibres,Err)

  ! Set V component 4) to 1 if fibre starts in element
  DO ElementUserNumber=1, NumberOfElementsFE, NumberGlobalXElements/NumberOfInSeriesFibres

    ! only if element is assigned to own domain
    !                                        DECOMPOSITION,   USER_ELEMENT_NUMBER, DOMAIN_NUMBER
    CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE, ElementUserNumber,   ElementDomain,Err)
    IF (ElementDomain == ComputationalNodeNumber) THEN
      !fibres begin in this element
      !                                        FIELD,              VARIABLE_TYPE,             FIELD_SET_TYPE,
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      !  USER_ELEMENT_NUMBER, COMPONENT_NUMBER, VALUE
       & ElementUserNumber,   4,                1,     Err)
      !PRINT "(A,I3.3,A,I5.5,A)", "x ", ComputationalNodeNumber, ": 3D el. ", ElementUserNumber, " has beginning fibre "
    ENDIF
  ENDDO

END SUBROUTINE InitializeFieldFiniteElasticity


END MODULE FIELDS_EQUATIONS_SET