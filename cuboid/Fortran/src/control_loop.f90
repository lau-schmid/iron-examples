
MODULE CONTROL_LOOP

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE OpenCMISS_Variables
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


SUBROUTINE CreateControlLoops()

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)

  !set the main control loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopMain,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoopMain,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  !Loop in time for StimDuration with the Stimulus applied.
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,0.0_CMISSRP,ElasticityTimeStep,ElasticityTimeStep,Err)

  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoopMain,OutputTimeStepStride,Err)
  IF (DebuggingOutput) THEN
    CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,cmfe_CONTROL_LOOP_TIMING_OUTPUT,Err)
  ELSE
    ! output types:
    ! CONTROL_LOOP_NO_OUTPUT = 0 !<No output from the control loop (no output of MainTime_* files)
    ! CONTROL_LOOP_PROGRESS_OUTPUT = 1 !<Progress output from control loop (also output MainTime_* files)
    ! CONTROL_LOOP_TIMING_OUTPUT = 2 !<Timing output from the control loop (also output MainTime_* files)
    ! CONTROL_LOOP_FILE_OUTPUT = -1 <Only MainTime_* files output
    CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,cmfe_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,cmfe_CONTROL_LOOP_TIMING_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,cmfe_CONTROL_LOOP_FILE_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,-1,Err)
  ENDIF


  !set the monodomain loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopM,'MONODOMAIN_TIME_LOOP',Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoopM,0.0_CMISSRP,ElasticityTimeStep,PDETimeStep,Err)
  ! question from Aaron: is the following necessary, or is the TimesSet-CALL before sufficient?
  IF (PdeNSteps /= -1) THEN
    CALL cmfe_ControlLoop_NumberOfIterationsSet(ControlLoopM,PdeNSteps,Err)
  ENDIF
  
  IF (DebuggingOutput) THEN
    CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)
  ELSE
    CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_FILE_OUTPUT,Err)
  ENDIF

  !set the finite elasticity loop (simple type)
  IF (.NOT. elasticityDisabled) THEN
    CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
    CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
    CALL cmfe_ControlLoop_TypeSet(ControlLoopFE,CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
    CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoopFE,ElasticityLoopMaximumNumberOfIterations,Err)
    CALL cmfe_ControlLoop_LabelSet(ControlLoopFE,'ELASTICITY_LOOP',Err)

    IF (DebuggingOutput) THEN
      CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopFE,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)
    ELSE
      CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopFE,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)
    ENDIF
  ENDIF
    
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

END SUBROUTINE CreateControlLoops


END MODULE CONTROL_LOOP