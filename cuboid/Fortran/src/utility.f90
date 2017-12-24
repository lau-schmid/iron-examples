
MODULE UTILITY

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


SUBROUTINE CalculateBioelectrics()
  !Calculate the bioelectrics geometric field
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)

  !PRINT*, "before cmfe_BioelectricsFiniteElasticity_UpdateGeometricField"
  !CALL cmfe_OutputInterpolationParameters(Problem,DependentFieldM, SolverParabolic,Err)
  
  CALL cmfe_BioelectricsFiniteElasticity_UpdateGeometricField(ControlLoopM,.TRUE.,Err)

  CALL SLEEP(2)
  
  !PRINT*, "after cmfe_BioelectricsFiniteElasticity_UpdateGeometricField"
  !CALL cmfe_OutputInterpolationParameters(Problem,DependentFieldM, SolverParabolic,Err)
  
  !reset the relative contraction velocity to 0
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSRP,Err)
END SUBROUTINE CalculateBioelectrics

SUBROUTINE ExportEMG()
  REAL(CMISSSP), DIMENSION(2) :: TimeStart, TimeStop
  REAL(CMISSSP) :: Total
 
  IF (.NOT. EnableExportEMG) RETURN
   
  IF (ComputationalNodeNumber == 0) WRITE(*,'(A)',advance='no') "Output EMG Data ... "
  CALL ETIME(TimeStart, Total)
  
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RegionM,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"EMGExample_M","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"EMGExample_M","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RegionFE,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"EMGExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"EMGExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF
  
  CALL ETIME(TimeStop, Total)
  TimingExportEMGUser = TimingExportEMGUser + (TimeStop(1) - TimeStart(1))
  TimingExportEMGSystem = TimingExportEMGSystem + (TimeStop(2) - TimeStart(2))
  
  IF (ComputationalNodeNumber == 0) PRINT*, " done"

END SUBROUTINE ExportEMG

FUNCTION GetTimeStamp()
  INTEGER(CMISSIntg), DIMENSION(3) :: Today, Now
  CHARACTER(LEN=100) :: TimeString
  CHARACTER(LEN=100) :: GetTimeStamp

  CALL IDATE(Today)  ! day, month, year
  CALL ITIME(Now)    ! h, min, s

  WRITE (TimeString, '(I2.2, A, I2.2, A, I4.4, A, I2.2, A, I2.2, A, I2.2)') &
  & Today(1), '.', Today(2), '.', Today(3), ' ', Now(1), ':', Now(2), ':', Now(3)

  GetTimeStamp = TimeString
END FUNCTION GetTimeStamp

FUNCTION GetMemoryConsumption()
  CHARACTER(LEN=100) :: GetMemoryConsumption
  CHARACTER(LEN=100) ::  &
   & pid, comm, state, ppid, pgrp, session, tty_nr, &
   & tpgid, flags, minflt, cminflt, majflt, cmajflt, &
   & utime, stime, cutime, cstime, priority, nice, &
   & O, itrealvalue, starttime, Description, Limit
  CHARACTER(LEN=10000) :: Debug
  INTEGER(CMISSIntg) :: I
  INTEGER(CMISSLintg) :: MemoryConsumption, &
   & VmSize, VmRSS, Shared, Text, Lib, Data, Dt, RssLimBytes, RssAnon, Pagesize, VSizeBytes

  ! Critical Section
  DO I=0,NumberOfComputationalNodes
    !CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
    IF (I == ComputationalNodeNumber) THEN

      ! read memory page size
      Stat = 0
      IF (I==0) CALL SYSTEM("getconf PAGESIZE > pagesize", Stat)
      IF (Stat == 0) THEN
        OPEN(UNIT=10, FILE="pagesize", ACTION="read", IOSTAT=Stat)
        IF (Stat == 0) THEN
          READ(10,*, IOSTAT=Stat) Pagesize
          CLOSE(UNIT=10)
        ELSE
          PRINT*, "Error opening pagesize."
        ENDIF
      ELSE
        PRINT*, "Error calling 'getconf PAGESIZE'"
        Pagesize = 4096
      ENDIF

      ! read from /proc/self/stat
      ! see http://man7.org/linux/man-pages/man5/proc.5.html for reference
      OPEN(UNIT=10, FILE="/proc/self/stat", ACTION="read", IOSTAT=stat)
      IF (STAT /= 0) THEN
        PRINT*, "Could not read memory consumption from /proc/self/stat."
      ELSE
        !READ(10,"(A)",IOSTAT=stat,advance='no') Debug
        !PRINT*, "proc/self/stat: "//TRIM(Debug)
        READ(10,*, IOSTAT=stat) pid, comm, state, ppid, pgrp, session, tty_nr, &
           & tpgid, flags, minflt, cminflt, majflt, cmajflt, &
           & utime, stime, cutime, cstime, priority, nice, &
           & O, itrealvalue, starttime, VSizeBytes, VmRSS, RssLimBytes
        CLOSE(UNIT=10)

        MemoryConsumption = VSizeBytes

        IF (ComputationalNodeNumber == 0) THEN
          WRITE(*, "(A,F7.3,A,I11,A,F7.3,A)") "     MemoryConsumption:", (MemoryConsumption/1e9), " GB (", &
            & MemoryConsumption, " B), Resident: ", (VmRSS*PageSize/1e9), " GB"
        ENDIF

        WRITE(GetMemoryConsumption, *) MemoryConsumption
      ENDIF

      ! read from /proc/self/limits
      IF (.TRUE.) THEN
        OPEN(UNIT=10, FILE="/proc/self/limits", ACTION="read", IOSTAT=stat)
        IF (STAT /= 0) THEN
          PRINT*, "Could not read limits from /proc/self/limits."
        ELSE
          DO
            READ(10, "(A26,A21)", IOSTAT=Stat) Description, Limit
            !PRINT*, "Description:["//TRIM(Description)//"], Limit:["//TRIM(Limit)//"]"
            IF (Stat /= 0) EXIT
            !IF (INDEX(Description, "Max resident set") /= 0) THEN
              !IF (TRIM(Limit) == "unlimited") THEN
              !  IF (ComputationalNodeNumber == 0) PRINT*, "    (Resident has no soft limit)"
              !ELSE
              !  IF (ComputationalNodeNumber == 0) PRINT*, "    (Resident is limited to ", Limit,")"
              !ENDIF
            !ENDIF
          ENDDO
          CLOSE(UNIT=10)
        ENDIF
      ENDIF

      ! read from /proc/self/statm
      OPEN(UNIT=10, FILE="/proc/self/statm", ACTION="read", IOSTAT=stat)
      IF (STAT /= 0) THEN
        PRINT*, "Could not read memory consumption from /proc/self/statm."
      ELSE
        READ(10,*, IOSTAT=stat) VmSize, VmRSS, Shared, Text, Lib, Data, Dt
        CLOSE(UNIT=10)

        ! VmRSS = RssAnon + Shared (all as number of pages)
        ! RssAnon is the resident set size of anonymous memory (real memory in RAM, not laid out in files)
        RssAnon = VmRSS - Shared

        IF (ComputationalNodeNumber == 0) THEN
          !PRINT*, "VmSize: ", VmSize, ", VmRSS: ", VmRSS, ", Shared: ", Shared, ", Text: ", Text, ", Lib:", Lib, ", Data:", Data
          !PRINT*, "RssAnon: ", RSSAnon, ", RssLimBytes: ", RssLimBytes, ", Pagesize: ", Pagesize

          ! Output Percentage
          !WRITE(*, "(3(A,F7.3),A,F5.1,A)") "     VmSize:", (VmSize*Pagesize/1.e9), " GB, RssAnon:", (RssAnon*Pagesize/1.e9), &
          !  & " GB, RssLimBytes: ", (RssLimBytes/1.e9), " GB (", (REAL(RssAnon*Pagesize) / RssLimBytes * 100.0), "%)"
        ENDIF
      ENDIF

      !PRINT*, TRIM(cmfe_CustomProfilingGetInfo(Err))
    ENDIF
  ENDDO

END FUNCTION GetMemoryConsumption

SUBROUTINE WriteTimingFile()

  CHARACTER(len=256) :: Filename = "duration."
  CHARACTER(len=20)  :: ComputationalNodeNumberStr
  CHARACTER(len=100) :: Hostname
  LOGICAL :: FileExists
  REAL(CMISSSP), DIMENSION(2) :: Elapsed     ! For receiving user and system time
  REAL(CMISSSP) :: DurationTotal
  REAL(CMISSRP) :: DurationInit, DurationStretchSim, DurationIntInit, DurationMainSim
  CHARACTER(LEN=100) :: TimeStampStr, MemoryConsumptionEnd

  ! create filename
  WRITE(ComputationalNodeNumberStr, '(I0.5)') ComputationalNodeNumber     ! convert ComputationalNodeNumber to string
  Filename = TRIM(Filename) // TRIM(ComputationalNodeNumberStr) // ".csv"

  ! Check if file exists
  INQUIRE(file=Filename, exist=FileExists)

  ! Write Comment in first line if file does not yet exist
  IF(.NOT. FileExists) THEN
    OPEN(unit=123, file=Filename, iostat=stat)
    IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'
    WRITE(123,'(A)') '# Stamp; Host; NProc; X; Y; Z; F; Total FE; Total M; End Time; ' // &
      & 'Dur. Init; Stretch Sim; Int. Init; Main Sim; Total; Total (User); Total (System); ' // &
      & 'ODE; Parabolic; FE; FE before Main Sim; Mem. Consumption after 1st timestep; Memory Consumption At End; ' // &
      & 'Parabolic reason; Newton reason; parabolic n. iter; min; max; newton n. iter; min; max; ' // &
      & 'level 0: stimulation handling; level 0: problem solve; level 1: MAIN_TIME_LOOP overhead; ' // &
      & 'level 1: MONODOMAIN_TIME_LOOP overhead; level 1: ELASTICITY_LOOP overhead; level 1: SolverDAE solve; ' // &
      & 'level 1: SolverParabolic solve; level 1: SolverFE solve; level 1: interpolate 1D->3D; level 1: interpolate 3D->1D; ' // &
      & 'level 1: file output; level 2: solver overhead; level 2: 0D solve; level 2: 1D solve; level 2: 3D solve; ' // &
      & 'level 3: 1D assembly; level 3: 1D solve; level 3: 1D other; level 3: 3D assembly; level 3: 3D solve; ' // &
      & 'level 3: 3D other;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;' // &
      & 'duration FESolverPreLoad; duration OdeSolverPreLoad; duration ParabolicSolverPreLoad; ' // &
      & 'duration FileOutputPreLoad (user); duration export EMG user; duration export EMG system; duration FileOutput (user); ' // &
      & 'duration FileOutput (system); duration FileOutputPreload (system); MonodomainSolverId; MonodomainPreconditionerId; ' // &
      & 'ODESolverId; NumberOfElementsInAtomX; NumberOfElementsInAtomY; NumberOfElementsInAtomZ; NumberOfSubdomainsX; ' // & 
      & 'NumberOfSubdomainsY; NumberOfSubdomainsZ; ModelType; ScenarioName'
      
    CLOSE(unit=123)
  ENDIF

  ! Write line in file
  OPEN(unit=123, file=Filename, iostat=stat, access='append')
  IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'

  DurationInit = TimeInitFinshed - TimeStart
  DurationStretchSim = TimeStretchSimFinished - TimeInitFinshed
  DurationIntInit = TimeMainSimulationStart - TimeStretchSimFinished
  DurationMainSim = TimeMainSimulationFinished - TimeMainSimulationStart

  CALL ETIME(Elapsed, DurationTotal)
  CALL HOSTNM(Hostname)
  MemoryConsumptionEnd = GetMemoryConsumption()

  TimeStampStr = GetTimeStamp()

  IF (CustomProfilingEnabled) THEN

    WRITE(123,"(4A,7(I11,A),(F8.3,A),11(F0.8,A),2(A,A),8(I7,A),21(F25.13,A),A,A,9(F8.3,A),3(I2,A),6(I8,A),I1,2A)") &
      & TRIM(TimeStampStr), ';', &
      & TRIM(Hostname(1:22)), ';', &          ! end of 4A
      & NumberOfComputationalNodes, ';', &
      & NumberGlobalXElements, ';', &
      & NumberGlobalYElements, ';', &
      & NumberGlobalZElements, ';', &
      & NumberOfInSeriesFibres, ';', &
      & NumberOfElementsFE, ';', &
      & NumberOfElementsM, ';', &             ! end of 7(I11,A)
      & TimeStop, ';', &                      ! (F8.3,A)
      & DurationInit, ';', &
      & DurationStretchSim, ';', &
      & DurationIntInit, ';', &
      & DurationMainSim,';',  &
      & DurationTotal, ';', &
      & Elapsed(1), ';', &
      & Elapsed(2), ';', &
      & CustomTimingOdeSolver, ';', &
      & CustomTimingParabolicSolver, ';', &
      & CustomTimingFESolver, ';', &
      & CustomTimingFESolverPreLoad, ';', &   ! end of 11(F0.8,A)
      & TRIM(ADJUSTL(MemoryConsumption1StTimeStep)), ';', &
      & TRIM(ADJUSTL(MemoryConsumptionEnd)), ';', &   ! end of 2(A,A)
      & CustomSolverConvergenceReasonParabolic, ';', &
      & CustomSolverConvergenceReasonNewton, ';', &
      & CustomSolverNumberIterationsParabolic, ';', &
      & CustomSolverNumberIterationsParabolicMin,';',  &
      & CustomSolverNumberIterationsParabolicMax, ';', &
      & CustomSolverNumberIterationsNewton, ';', &
      & CustomSolverNumberIterationsNewtonMin, ';', &
      & CustomSolverNumberIterationsNewtonMax, ';', &   ! end of 8(I7,A)
      ! Custom Profiling durations
      & cmfe_CustomProfilingGetDuration("level 0: stimulation handling", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 0: problem solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: MAIN_TIME_LOOP overhead", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: MONODOMAIN_TIME_LOOP overhead", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: ELASTICITY_LOOP overhead", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: SolverDAE solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: SolverParabolic solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: SolverFE solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: interpolate 1D->3D", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: interpolate 3D->1D", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: file output", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 2: solver overhead", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 2: 0D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 2: 1D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 2: 3D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 1D assembly", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 1D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 1D other", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 3D assembly", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 3D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 3D other", Err), ';', &                 ! end of 21(F25.13,A)
      & ';;;;;;;;;;;;;;', &             ! end of A (14 free columns)
      & ';;;;;;;;;;;;;;;;;;;;;;;;', &   ! end of A (3*8=24 free columns) (8(I17,A,I5,A,I7,A))
      ! custom profiling memory consumption
      !& cmfe_CustomProfilingGetMemory("distributed vector cmiss DP", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed vector cmiss DP", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed vector cmiss DP", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed vector cmiss INTG", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed vector cmiss INTG", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed vector cmiss INTG", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage diag", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage diag", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage diag", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage, offdiag", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage, offdiag", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage, offdiag", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage, row ind.", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage, row ind.", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage, row ind.", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage, col. ind.", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage, col. ind.", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage, col. ind.", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage (local to global mapping)", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage (local to global mapping)", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage (local to global mapping)", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed vector petsc", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed vector petsc", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed vector petsc", Err), ';', &       ! end of 8(I17,A,I5,A,I7,A)
      & CustomTimingFESolverPreLoad, ';', &
      & CustomTimingOdeSolverPreLoad, ';', &
      & CustomTimingParabolicSolverPreLoad, ';', &
      & CustomTimingFileOutputUserPreLoad, ';', &
      & TimingExportEMGUser, ';', &
      & TimingExportEMGSystem, ';', &
      & CustomTimingFileOutputUser, ';', &
      & CustomTimingFileOutputSystem, ';', &
      & CustomTimingFileOutputSystemPreLoad, ';', &     ! end of 9(F8.3,A)
      & MonodomainSolverId, ';', &
      & MonodomainPreconditionerId, ';', &
      & ODESolverId, ';', &                             ! end of 3(I2,A)
      & NumberOfElementsInAtomX, ';', &
      & NumberOfElementsInAtomY, ';', &
      & NumberOfElementsInAtomZ, ';', &                 
      & NumberOfSubdomainsX, ';', &
      & NumberOfSubdomainsY, ';', &
      & NumberOfSubdomainsZ, ';', &                            ! end of 6(I8,A)
      & ModelType, ';', &                                       ! I1
      & TRIM(ScenarioName)

  ELSE  ! custom profiling is disabled
    
    WRITE(123,"(4A,7(I11,A),(F8.3,A),11(F0.8,A),2(A,A),8(I7,A),9(F8.3,A),3(I2,A),6(I8,A),I1,2A)") &
      & TRIM(TimeStampStr), ';', &
      & TRIM(Hostname(1:22)), ';', &                    ! end of 4A
      & NumberOfComputationalNodes, ';', &
      & NumberGlobalXElements, ';', &
      & NumberGlobalYElements, ';', &
      & NumberGlobalZElements, ';', &
      & NumberOfInSeriesFibres, ';', &
      & NumberOfElementsFE, ';', &
      & NumberOfElementsM, ';', &                       ! end of 7(I11,A)
      & TimeStop, ';', &                                ! (F8.3,A)
      & DurationInit, ';', &
      & DurationStretchSim, ';', &
      & DurationIntInit, ';', &
      & DurationMainSim,';',  &
      & DurationTotal, ';', &
      & Elapsed(1), ';', &
      & Elapsed(2), ';', &
      & CustomTimingOdeSolver, ';', &
      & CustomTimingParabolicSolver, ';', &
      & CustomTimingFESolver, ';', &
      & CustomTimingFESolverPreLoad, ';', &             ! end of 11(F0.8,A)
      & TRIM(ADJUSTL(MemoryConsumption1StTimeStep)), ';', &
      & TRIM(ADJUSTL(MemoryConsumptionEnd)), ';', &     ! end of 2(A,A)
      & CustomSolverConvergenceReasonParabolic, ';', &
      & CustomSolverConvergenceReasonNewton, ';', &
      & CustomSolverNumberIterationsParabolic, ';', &
      & CustomSolverNumberIterationsParabolicMin,';',  &
      & CustomSolverNumberIterationsParabolicMax, ';', &
      & CustomSolverNumberIterationsNewton, ';', &
      & CustomSolverNumberIterationsNewtonMin, ';', &
      & CustomSolverNumberIterationsNewtonMax, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;', &    ! end of 8(I7,A)
      & CustomTimingFESolverPreLoad, ';', &
      & CustomTimingOdeSolverPreLoad, ';', &
      & CustomTimingParabolicSolverPreLoad, ';', &
      & CustomTimingFileOutputUserPreLoad, ';', &
      & TimingExportEMGUser, ';', &
      & TimingExportEMGSystem, ';', &
      & CustomTimingFileOutputUser, ';', &
      & CustomTimingFileOutputSystem, ';', &
      & CustomTimingFileOutputSystemPreLoad, ';', &   ! end of 9(F8.3,A)
      & MonodomainSolverId, ';', &
      & MonodomainPreconditionerId, ';', &
      & ODESolverId, ';', &                           ! end of 3(I2,A)
      & NumberOfElementsInAtomX, ';', &
      & NumberOfElementsInAtomY, ';', &
      & NumberOfElementsInAtomZ, ';', &
      & NumberOfSubdomainsX, ';', &
      & NumberOfSubdomainsY, ';', &
      & NumberOfSubdomainsZ, ';', &                          ! end of 6(I8,A)
      & ModelType, ';', &                                     ! I1
      & TRIM(ScenarioName)
  ENDIF

  CLOSE(unit=123)
END SUBROUTINE WriteTimingFile

SUBROUTINE WriteCustomProfilingFile()
  CHARACTER(len=256) :: Filename = "profiling."
  CHARACTER(len=20)  :: ComputationalNodeNumberStr
  CHARACTER(len=100) :: Hostname

  CALL HOSTNM(Hostname)

  ! create filename
  WRITE(ComputationalNodeNumberStr, '(I0.5)') ComputationalNodeNumber     ! convert ComputationalNodeNumber to string
  Filename = TRIM(Filename) // TRIM(ComputationalNodeNumberStr) // ".txt"

  OPEN(unit=123, file=Filename, iostat=stat, access='append')
  IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'

  ! Write header to file
  WRITE(123,'(3A,I5)') '# CustomProfiling at ',GetTimeStamp(), ", number of computational nodes: ", NumberOfComputationalNodes
  WRITE(123,'(A)') '# Hostname, NProc, X,Y,Z,F, nFE, nM, tEnd'
  WRITE(123,'(2A,7(I11,A),(F8.3,A))')  &
    & TRIM(Hostname(1:22)), ';', &
    & NumberOfComputationalNodes, ';', &
    & NumberGlobalXElements, ';', &
    & NumberGlobalYElements, ';', &
    & NumberGlobalZElements, ';', &
    & NumberOfInSeriesFibres, ';', &
    & NumberOfElementsFE, ';', &
    & NumberOfElementsM, ';', &
    & TimeStop, ';'

  ! write info
  WRITE(123,*) cmfe_CustomProfilingGetInfo(Err)

  CLOSE(unit=123)
END SUBROUTINE WriteCustomProfilingFile

SUBROUTINE HandleSolverInfo(TimeStep)
  REAL(CMISSRP), INTENT(IN) :: TimeStep
  CHARACTER(len=256) :: Filename = "iterations.csv"
  LOGICAL :: FileExists

  RETURN      ! disable solver info for the moment
  
  IF (ComputationalNodeNumber == 0) THEN

    CALL cmfe_CustomSolverInfoGet(CustomSolverConvergenceReasonParabolic, CustomSolverConvergenceReasonNewton, &
      & CustomSolverNumberIterationsParabolic, CustomSolverNumberIterationsParabolicMin, CustomSolverNumberIterationsParabolicMax, &
      & CustomSolverNumberIterationsNewton, CustomSolverNumberIterationsNewtonMin, CustomSolverNumberIterationsNewtonMax, Err)

    ! Check if file exists
    INQUIRE(file=Filename, exist=FileExists)

    ! Write Comment in first line if file does not yet exist
    IF(.NOT. FileExists) THEN
      OPEN(unit=123, file=Filename, iostat=stat)
      IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'
      WRITE(123,'(A)') '# Timestep; Parabolic reason; Newton reason; parabolic n. iter; min; max; newton n. iter; min; max '
      WRITE(123,'(A)') '# Parabolic reason: 1=RTOL_NORMAL, 2=RTOL, 3=ATOL, 4=ITS, 5=CG_NEG_CURVE, 6=CG_CONSTRAINED, ' // &
        & '7=STEP_LENGTH, 8=HAPPY_BREAKDOWN, 9=ATOL_NORMAL'
      WRITE(123,'(A)') '# Newton reason: 2=FNORM_ABS ||F||<atol, 3=FNORM_RELATIVE ||F|| < rtol*||F_initial||, ' // &
        & '4=SNORM_RELATIVE Newton computed step size small || delta x || < stol || x ||, ' // &
        & '5=ITS, 7=TR_DELTA'
      CLOSE(unit=123)
    ENDIF

    ! Write line in file
    OPEN(unit=123, file=Filename, iostat=stat, access='append')
    IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'

    WRITE(123,"(F0.8,A,8(I7,A))") &
      & TimeStep, ';', &
      & CustomSolverConvergenceReasonParabolic, ';', &
      & CustomSolverConvergenceReasonNewton, ';', &
      & CustomSolverNumberIterationsParabolic, ';', &
      & CustomSolverNumberIterationsParabolicMin,';',  &
      & CustomSolverNumberIterationsParabolicMax, ';', &
      & CustomSolverNumberIterationsNewton, ';', &
      & CustomSolverNumberIterationsNewtonMin, ';', &
      & CustomSolverNumberIterationsNewtonMax, ';'

    CLOSE(unit=123)
  ENDIF

END SUBROUTINE HandleSolverInfo

SUBROUTINE gdbParallelDebuggingBarrier()
  INTEGER(CMISSIntg) :: Gdb_Resume
  INTEGER(CMISSIntg) :: MPI_IERROR, MPI_COMM_WORLD
  INTEGER(CMISSIntg) :: ComputationalNodeNumber, NumberOfComputationalNodes
  Gdb_Resume = 0

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ComputationalNodeNumber,MPI_IERROR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NumberOfComputationalNodes,MPI_IERROR)

  IF (NumberOfComputationalNodes > 1) THEN
    PRINT*, "Node ", ComputationalNodeNumber, ", UID ",GETPID()," is waiting for Gdb_Resume=", Gdb_Resume &
      & , " to become 1 " // NEW_LINE('A') // "sudo gdb cuboid ",GETPID(), NEW_LINE('A') //"select-frame 2" // &
      & NEW_LINE('A') // "set var gdb_resume = 1" // NEW_LINE('A') // &
      & "info locals" // NEW_LINE('A') // "next"
    DO WHILE (Gdb_Resume == 0)
      CALL Sleep(1)
    ENDDO
    PRINT*, "Node ", ComputationalNodeNumber, " resumes because gdb_resume=", Gdb_Resume, "."
  ENDIF

END SUBROUTINE gdbParallelDebuggingBarrier


END MODULE UTILITY