!> \file
!> \author Adam Reeve, Benjamin Maier, Nehzat Emamy, Aaron Krämer, Thomas Klotz
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

!> Main program
PROGRAM SKELETAL_MUSCLE_SIMULATION

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE SETUP_PARAMETERS
  USE OpenCMISS_Variables
  USE REGION_MESH
  USE CELLML
  USE CONTROL_LOOP 
  USE DECOMPOSITION 
  USE FIELDS_EQUATIONS_SET 
  USE PROBLEM_ROUTINES
  USE SOLVERS
  USE STIMULATION_BOUNDARY_CONDITIONS 
  USE UTILITY
  USE EQUATIONS
  
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


!##################################################################################################################################
!##################################################################################################################################
!##################################################################################################################################

  WRITE(*,'(A,A)') TRIM(GetTimeStamp()), " Program started."

  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeStart = DurationSystemUser(2)
    
  CALL ParseParameters()
    
  !================================================================================================================================
  !  G E N E R A L   F E A T U R E S
  !================================================================================================================================
  CALL CreateRegionMesh()
  CALL CreateDecompositionFiniteElasticity()
  CALL CreateDecompositionMonodomain()

  !================================================================================================================================
  !  F I N I T E   E L A S T I C I T Y
  !================================================================================================================================
  CALL CreateFieldFiniteElasticity()

  !================================================================================================================================
  !  M O N O D O M A I N
  !================================================================================================================================

  CALL CreateFieldMonodomain()

  !================================================================================================================================
  !  EQUATIONS SET
  CALL CreateEquationsSet()
  !UPDATE THE INDEPENDENT FIELD IndependentFieldM
  CALL InitializeFieldMonodomain()

  !UPDATE THE INDEPENDENT FIELD IndependentFieldFE
  CALL InitializeFieldFiniteElasticity()

  CALL CreateEquations()
  CALL InitializeCellML()
  CALL CreateProblem()

  CALL CreateControlLoops()
  CALL CreateSolvers()
  !--------------------------------------------------------------------------------------------------------------------------------
  !boundary conditions
  CALL SetBoundaryConditions()
  
  ! Output the data structure Problem
  !IF (DEBUGGING_PROBLEM_OUTPUT .AND. ComputationalNodeNumber == 0) THEN
  !  PRINT*, ""
  !  PRINT*, ""
  !  CALL cmfe_PrintProblem(Problem,6,30,Err)
  !  PRINT*, ""
  !  PRINT*, ""
    !PRINT*, "End the program after output of problem datastructure"
    !STOP
  !ENDIF

  !--------------------------------------------------------------------------------------------------------------------------------
  !Calculate the bioelectrics geometric field
  CALL CalculateBioelectrics()
  
 !PRINT*, "Size IndependentFieldFE:",cmfe_getFieldSize(IndependentFieldFE, Err),"Bytes"
  
  !PRINT*, "Abort program in FortranExample.f90:450"
  CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
  !STOP
  CALL ExportEMG()
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem

  !Solve the problem -- bring to new length before applying the stimulus
  IF (InitialStretch == 1.0_CMISSRP) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, "1.) Pre-stretch disabled, initial stretch =", InitialStretch
  ELSE
    IF (ComputationalNodeNumber == 0) PRINT *, "1.) Pre-stretch simulation, initial stretch =", InitialStretch
  ENDIF
  Temp = GetMemoryConsumption()
  CALL cmfe_CustomSolverInfoReset(Err)
  IF (DEBUGGING_PARALLEL_BARRIER) CALL gdbParallelDebuggingBarrier()

  ! store duration
  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeInitFinshed = DurationSystemUser(2)
  
  IF (InitialStretch /= 1.0_CMISSRP) THEN
    CALL cmfe_Problem_Solve(Problem,Err) ! ### PAPER SETTING: not necessary
  ENDIF
    
  ! store duration
  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeStretchSimFinished = DurationSystemUser(2)
  
  CALL HandleSolverInfo(-1.0_CMISSRP)

  !reset the relative contraction velocity to 0
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSRP,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_ParameterSetUpdateConstant(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,PMax, Err)
  ENDIF


  CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoopFE,1,Err)

! no change for BCs -- fix at this length!!!

  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  IF (ComputationalNodeNumber == 0) PRINT*, "2.) Simulate with stimulation"
  MemoryConsumptionBeforeSim = GetMemoryConsumption()

  CALL cmfe_CustomTimingGet(CustomTimingOdeSolverPreLoad, CustomTimingParabolicSolverPreLoad, &
    & CustomTimingFESolverPreLoad, CustomTimingFileOutputUserPreLoad, CustomTimingFileOutputSystemPreLoad, Err)
  CALL cmfe_CustomTimingReset(Err)
  CALL cmfe_CustomProfilingReset(Err)   ! reset and discard all custom profiling data that was recorded so far
  
  ! store duration
  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeMainSimulationStart = DurationSystemUser(2)

  time = 0.0_CMISSRP
  VALUE = 0.0_CMISSRP
  k = 1       ! row in firing_times input (time)
  ! main time loop
  DO WHILE(time < TimeStop-1e-10)
  
    CALL cmfe_CustomProfilingStart("level 0: stimulation handling",Err)

    IF (ComputationalNodeNumber == 0) PRINT "(A,F0.5,A)","t = ",time," s"
    !-------------------------------------------------------------------------------------------------------------------------------
    !Set the Stimulus for monodomain at the middle of the fibres

  !  VALUE = VALUE-ABS(Vmax)/20.0_CMISSRP*StimDuration
  !!  VALUE = VALUE+ABS(Vmax)/10.0_CMISSRP*StimDuration
  !  CALL cmfe_ControlLoop_BoundaryConditionUpdate(ControlLoopFE,1,1,VALUE,Err)
  !  DO NodeIdx=1,SIZE(LeftSurfaceNodes,1)
  !    NodeNumber=LeftSurfaceNodes(NodeIdx)
  !    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
  !    IF(NodeDomain==ComputationalNodeNumber) THEN
  !      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
  !        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
  !    ENDIF
  !  ENDDO

    CALL SetStimulationAtNodes(StimValuePerNode)

    IF (ComputationalNodeNumber == 0) THEN
      PRINT "(A,I4,A,I6)", "   Number of stimulated fibres on process 0: ", NumberFiringFibres," of",NumberOfFibres
    ENDIF

    !-------------------------------------------------------------------------------------------------------------------------------
    !Solve the problem for the stimulation time
    IF (ComputationalNodeNumber == 0) Print*, "  Solve with stimulation,    time span: ", time, " to ",time+StimDuration
    CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,time,time+StimDuration,ElasticityTimeStep,Err)

    CALL cmfe_CustomSolverInfoReset(Err)
    CALL cmfe_CustomProfilingStop("level 0: stimulation handling",Err)
    
    CALL cmfe_Problem_Solve(Problem,Err)
    
    CALL cmfe_CustomProfilingStart("level 0: stimulation handling",Err)
    CALL HandleSolverInfo(time)
    Temp = GetMemoryConsumption()
    IF (DebuggingOnlyRunShortPartOfSimulation) EXIT

    !-------------------------------------------------------------------------------------------------------------------------------
    !Now turn the stimulus off
    CALL SetStimulationAtNodes(0.0_CMISSRP)
    
    !-------------------------------------------------------------------------------------------------------------------------------
    !Solve the problem for the rest of the period
    IF (ComputationalNodeNumber == 0) PRINT*, "  Solve without stimulation, time span: ", time+StimDuration, " to ",time+StimPeriod
    CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,time+StimDuration,time+StimPeriod,ElasticityTimeStep,Err)

    CALL cmfe_CustomSolverInfoReset(Err)    
    CALL cmfe_CustomProfilingStop("level 0: stimulation handling",Err)
    
    CALL cmfe_Problem_Solve(Problem,Err)
    
    CALL cmfe_CustomProfilingStart("level 0: stimulation handling",Err)
    CALL HandleSolverInfo(time+StimDuration)
    !-------------------------------------------------------------------------------------------------------------------------------
    time = time + StimPeriod
    k = k+1

    IF (k == 2) THEN
      MemoryConsumption1StTimeStep = GetMemoryConsumption()
    ELSE
      Temp = GetMemoryConsumption()
    ENDIF
    CALL cmfe_CustomProfilingStop("level 0: stimulation handling",Err)

  ENDDO

  ! store duration
  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeMainSimulationFinished = DurationSystemUser(2)
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------------------------------------------------
  CALL ExportEMG()
  CALL cmfe_TimingSummaryOutput(Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  CALL cmfe_CustomTimingGet(CustomTimingOdeSolver, CustomTimingParabolicSolver, CustomTimingFESolver, CustomTimingFileOutputUser, &
    & CustomTimingFileOutputSystem, Err)

  CALL WriteTimingFile()

  CALL WriteCustomProfilingFile()

  IF (ComputationalNodeNumber == 0 .AND. CustomProfilingEnabled) THEN
    PRINT*, TRIM(cmfe_CustomProfilingGetInfo(Err))
  ENDIF

  PRINT*, ""
  PRINT*, "--------------------------------------------------"
  PRINT*, "Process ", ComputationalNodeNumber
  PRINT*, "Timing (user time):"
  PRINT*, "   Ode Solver:       preload: ", CustomTimingOdeSolverPreLoad, " s, main: ", CustomTimingOdeSolver, " s"
  PRINT*, "   Parabolic Solver: preload: ", CustomTimingParabolicSolverPreLoad, " s, main: ", CustomTimingParabolicSolver, " s"
  PRINT*, "   3D FE Solver:     preload: ", CustomTimingFESolverPreLoad, " s, main: ", CustomTimingFESolver, " s"
  PRINT*, "   Node File Output: preload: ", CustomTimingFileOutputUserPreLoad, " s, main: ", CustomTimingFileOutputUser, " s"
  PRINT*, "           (system): preload: ", CustomTimingFileOutputSystemPreLoad, " s, main: ", CustomTimingFileOutputSystem, " s"
  PRINT*, "   Total Simulation: preload: ", (TimeStretchSimFinished - TimeInitFinshed), " s, main: ", &
    & (TimeMainSimulationFinished-TimeMainSimulationStart), " s"
  PRINT*, "   EMG Output: user: ", TimingExportEMGUser, " s, system: ", TimingExportEMGSystem, " s"
  PRINT*, ""

  WRITE(*,'(A,A)') TRIM(GetTimeStamp()), " Program successfully completed."
  STOP
  
CONTAINS

END PROGRAM SKELETAL_MUSCLE_SIMULATION


