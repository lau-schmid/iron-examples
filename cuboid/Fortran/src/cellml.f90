
MODULE CELLML

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE OpenCMISS_Variables
  USE REGION_MESH
  USE DECOMPOSITION 
  USE FIELDS_EQUATIONS_SET 
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

SUBROUTINE InitializeCellML()
  REAL(CMISSRP) :: InitialPressure
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellMLEnvironment,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,RegionM,CellMLEnvironment,Err)
  !Import the Shorten et al. 2007 model from a file
  
  !                            CellML env.    ,URI,                modelIndex (out)
  CALL cmfe_CellML_ModelImport(CellMLEnvironment,CellMLModelFilename,shortenModelIndex,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !,- to set from this side

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"wal_environment/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"razumova/L_S",Err) !(= l_{hs})
   ! muss theoretisch auch übergeben werden, ist aber noch nicht im subcell model enthalten:
   ! CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"razumova/velo",Err) ( = d{l_{hs}} / d{t})
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"Aliev_Panfilov/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"Razumova/l_hs",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"Razumova/velo",Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"Aliev_Panfilov/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"Razumova/l_hs",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"Razumova/rel_velo",Err)
  ENDIF
!  CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex,"razumova/rel_velo",Err)
!
!  CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex2,"wal_environment/I_HH",Err)
!  CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex2,"razumova/L_S",Err)
!  CALL cmfe_CellML_VariableSetAsKnown(CellMLEnvironment,shortenModelIndex2,"razumova/rel_velo",Err)
  !,- to get from the CellML side
!  CALL cmfe_CellML_VariableSetAsWanted(CellMLEnvironment,shortenModelIndex,"wal_environment/I_T",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(CellMLEnvironment,shortenModelIndex,"wal_environment/I_ionic_s",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(CellMLEnvironment,shortenModelIndex,"wal_environment/I_ionic_t",Err)
  !
  !NOTE: If an INTERMEDIATE (or ALGEBRAIC) variable should be used in a mapping, it has to be set as known or wanted first!
  !,  --> set "razumova/stress" as wanted!
  !,  --> no need to set "wal_environment/vS" since all STATE variables are automatically set as wanted!
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
  ! wird in CellML als Stress berechnet, dann aber als Rate ausgegeben. passt das so?: Achtung, Modell nicht implementiert in ModelType==0.
    CALL cmfe_CellML_VariableSetAsWanted(CellMLEnvironment,shortenModelIndex,"razumova/stress",Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_CellML_VariableSetAsWanted(CellMLEnvironment,shortenModelIndex,"Razumova/ActiveStress",Err)
    CALL cmfe_CellML_VariableSetAsWanted(CellMLEnvironment,shortenModelIndex,"Razumova/Activation",Err)
  ENDIF
!  CALL cmfe_CellML_VariableSetAsWanted(CellMLEnvironment,shortenModelIndex2,"razumova/stress",Err)
  !,- and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
  !Finish the CellML environment
  
  IF (ODESolverId==2) THEN      ! BDF solver
   ! To initialize the BDF solver, OpenCMISS requires to have set CELLML%MAX_NUMBER_OF_INTERMEDIATE.
   ! For now, we can manipulate it here:
    CALL cmfe_CellML_IntermediateMaxNumberSet(CellMLEnvironment,1,Err) ! todo: exact number required or just something > 0?
  ENDIF
  
  CALL cmfe_CellML_CreateFinish(CellMLEnvironment,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellMLEnvironment,Err)
  
  !Map the half-sarcomere length L_S
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    !Map the transmembrane voltage V_m
    !                                       CellML env., "from"-field
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,      DependentFieldM, & 
    !  variable type,          comp. no., "from"-field param. set
     & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    !  "to"-cellml model user no., "to"-variable,      "to"-cellml variable parameter set
     & shortenModelIndex,          "wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE,Err)
     
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"wal_environment/vS", &
     & CMFE_FIELD_VALUES_SET_TYPE,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    
    !Map the active stress
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"razumova/stress", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,IndependentFieldM, &
     & CMFE_FIELD_U1_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Razumova/l_hs",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the sarcomere relative contraction velocity
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,IndependentFieldM, & 
     & CMFE_FIELD_U2_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Razumova/velo",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the transmembrane voltage V_m
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,DependentFieldM, &
     & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Aliev_Panfilov/V_m",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"Aliev_Panfilov/V_m", &
     & CMFE_FIELD_VALUES_SET_TYPE,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the active stress
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"Razumova/A_1", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"Razumova/A_2", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"Razumova/x_1", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"Razumova/x_2", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_VALUES_SET_TYPE,Err)
  
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    !Map the sarcomere half length L_S
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,1,&
     & CMFE_FIELD_VALUES_SET_TYPE,shortenModelIndex,"Razumova/l_hs",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the sarcomere relative contraction velocity
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,3,&
     & CMFE_FIELD_VALUES_SET_TYPE,shortenModelIndex,"Razumova/rel_velo",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the transmembrane voltage V_m
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,&
     & CMFE_FIELD_VALUES_SET_TYPE,shortenModelIndex,"Aliev_Panfilov/V_m",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"Aliev_Panfilov/V_m",CMFE_FIELD_VALUES_SET_TYPE, &
     & DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the active stress
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"Razumova/ActiveStress",&
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex,"Razumova/Activation",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_VALUES_SET_TYPE,Err)

  ENDIF

!  CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"razumova/L_S",CMFE_FIELD_VALUES_SET_TYPE,Err)
!  !Map the sarcomere relative contraction velocity
!  CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"razumova/rel_velo",CMFE_FIELD_VALUES_SET_TYPE,Err)
!  !Map the transmembrane voltage V_m
!  CALL cmfe_CellML_CreateFieldToCellMLMap(CellMLEnvironment,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE,Err)
!  CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex2,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE, &
!   & DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
!  !Map the active stress
!  CALL cmfe_CellML_CreateCellMLToFieldMap(CellMLEnvironment,shortenModelIndex2,"razumova/stress",CMFE_FIELD_VALUES_SET_TYPE, &
!   & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_FieldMapsCreateFinish(CellMLEnvironment,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Initialise dependent field for monodomain
  !> \todo - get V_m initialial value.
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & -79.974_CMISSRP,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & 0.0_CMISSRP,Err)

  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & 0.0_CMISSRP,Err)

  ENDIF

  !Initialise dependent field for Finite Elasticity from undeformed geometry and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
   & CMFE_FIELD_VALUES_SET_TYPE,1,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
   & CMFE_FIELD_VALUES_SET_TYPE,2,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
   & CMFE_FIELD_VALUES_SET_TYPE,3,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)

  IF (ModelType == 0 .OR. ModelType == 1) THEN    ! 3,3a, "MultiPhysStrain"
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
     & 0.0_CMISSRP,Err)

  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    InitialPressure = -2.0_CMISSRP*MAT_FE(2)-MAT_FE(1)
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
     & InitialPressure,Err)

  ENDIF

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML models field
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellMLEnvironment,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellMLEnvironment,Err)

  !Set up the models field
  CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
   & shortenModelIndex,Err)

  !Create the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellMLEnvironment,CellMLStateFieldUserNumber,CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateFinish(CellMLEnvironment,Err)

  !Create the CellML intermediate field
  IF (ModelType==0 .OR. ModelType==2 .OR. ODESolverId==2) THEN ! ModelType==0: 3a, "MultiPhysStrain", ModelType==2: 4, "Titin", ODESolverId== 2: BDF
    CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
    CALL cmfe_CellML_IntermediateFieldCreateStart(CellMLEnvironment,CellMLIntermediateFieldUserNumber, & 
     & CellMLIntermediateField,Err)
    CALL cmfe_CellML_IntermediateFieldCreateFinish(CellMLEnvironment,Err)
  ENDIF

  !Create the CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellMLEnvironment,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellMLEnvironment,Err)

END SUBROUTINE InitializeCellML

END MODULE CELLML