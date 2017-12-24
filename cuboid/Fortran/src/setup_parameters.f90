
MODULE SETUP_PARAMETERS

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE OpenCMISS_Variables
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

! Test whether parametrization and number of processes matches and is valid
FUNCTION CheckGeometry()
  LOGICAL :: CheckGeometry
  CheckGeometry = .TRUE.

  IF (NumberGlobalXElements <= 0 .OR. NumberGlobalYElements <= 0 .OR. NumberGlobalZElements <= 0 &
    & .OR. NumberOfInSeriesFibres <= 0) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: Number of elements (", NumberGlobalXElements, ", ",NumberGlobalYElements, ", ", &
        & NumberGlobalZElements, ", ", NumberOfInSeriesFibres, ") is invalid!"
    ENDIF
    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsInAtomX <= 0 ) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: NumberOfElementsInAtomX=", NumberOfElementsInAtomX, " is invalid!"
    ENDIF
    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsInAtomY <= 0 ) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: NumberOfElementsInAtomY=", NumberOfElementsInAtomY, " is invalid!"
    ENDIF
    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsInAtomZ <= 0 ) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: NumberOfElementsInAtomZ=", NumberOfElementsInAtomZ, " is invalid!"
    ENDIF
    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsMPerFibre * NumberOfInSeriesFibres /= NumberOfElementsMPerFibreLine) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: NumberOfElementsMPerFibreLine=",NumberOfElementsMPerFibreLine, &
        & " is not an integer multiple of F=",NumberOfInSeriesFibres,"!"
    ENDIF

    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsInAtomX > NumberGlobalXElements) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Warning: Decreasing NumberOfElementsInAtomX=",NumberOfElementsInAtomX, &
        & " to value of NumberGlobalXElements=",NumberGlobalXElements,"."
    ENDIF
    NumberOfElementsInAtomX = NumberGlobalXElements
  ENDIF
  
  IF (NumberOfElementsInAtomY > NumberGlobalYElements) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Warning: Decreasing NumberOfElementsInAtomY=",NumberOfElementsInAtomY, &
        & " to value of NumberGlobalYElements=",NumberGlobalYElements,"."
    ENDIF
    NumberOfElementsInAtomY = NumberGlobalYElements
  ENDIF
  
  IF (NumberOfElementsInAtomZ > NumberGlobalZElements) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Warning: Decreasing NumberOfElementsInAtomZ=",NumberOfElementsInAtomZ, &
        & " to value of NumberGlobalZElements=",NumberGlobalZElements,"."
    ENDIF
    NumberOfElementsInAtomZ = NumberGlobalZElements
  ENDIF

END FUNCTION CheckGeometry

SUBROUTINE ReadInputFiles()
  LOGICAL :: FileExists
  INTEGER :: Status
  INTEGER :: I
  CHARACTER(len=100000) :: Buffer
  INTEGER :: NumberOfEntries

  ! Read in firing times file
  FiringTimesFile = TRIM(InputDirectory) // TRIM(FiringTimesFile)
  INQUIRE(file=FiringTimesFile, exist=FileExists)

  IF (.NOT. FileExists) THEN
    PRINT*, "Error: File """ // TRIM(FiringTimesFile) // """ does not exist!"
    STOP
  ENDIF

  IF (ComputationalNodeNumber == 0) PRINT*,  "Open file """ // TRIM(FiringTimesFile) // """."

  OPEN(unit=5, file=FiringTimesFile, action="read", iostat=Status)

  ! loop over maximum 10000 lines
  DO I = 1, 10000
    READ(5,*,iostat=Status) MotorUnitFiringTimes(I,:)
    IF (Status /= 0) THEN
      IF (ComputationalNodeNumber == 0) PRINT*,  "File """ // TRIM(FiringTimesFile) // """ contains ", I, " firing times."
      EXIT
    ENDIF
  ENDDO
  CLOSE(unit=5)


  !--------------------------------------------------------------------------------------------------------------------------------
  ! Read in InnervationZoneFile
  ! Gaussian distribution with mean 0 and std 2 (MATLAB: k = 2*randn(len,1), kk = int64(k);)
  InnervationZoneFile = TRIM(InputDirectory) // TRIM(InnervationZoneFile)
  INQUIRE(file=InnervationZoneFile, exist=FileExists)

  IF (.NOT. FileExists) THEN
    PRINT*, "Error: File """ // TRIM(InnervationZoneFile) // """ does not exist!"
    STOP
  ENDIF

  OPEN(unit=5, file=InnervationZoneFile, action="read", iostat=Status)
  READ(unit=5, fmt='(A)') Buffer

  ! Determine the number of entries in the first line
  NumberOfEntries = 1
  DO I = 0, LEN_TRIM(Buffer)
    IF (Buffer(I:I) == ',') THEN
      NumberOfEntries = NumberOfEntries + 1
    ENDIF
  ENDDO

  REWIND(UNIT=5)
  IF (ComputationalNodeNumber == 0) THEN
    PRINT*, "File  """ // TRIM(InnervationZoneFile) // """ contains ", NumberOfEntries, " Entries, NumberOfFibres=", &
      & NumberOfFibres, "."
  ENDIF

  ALLOCATE(InnervationZoneOffset(NumberOfFibres), Stat=stat)
  
  IF (stat /= 0) THEN
    PRINT *, "Could not allocate ", NumberOfFibres, " objects for InnervationZoneOffset!"
    STOP
  ENDIF

  IF (NumberOfFibres <= NumberOfEntries) THEN
    READ(unit=5, fmt=*, iostat=Status) InnervationZoneOffset(1:NumberOfFibres)
  ELSE
    READ(unit=5, fmt=*, iostat=Status) InnervationZoneOffset(1:NumberOfEntries)

    ! fill with already read data
    DO I = NumberOfEntries,NumberOfFibres
      InnervationZoneOffset(I) = InnervationZoneOffset(MOD(I, NumberOfEntries)+1)
      !PRINT*, "fill I=",I," with entry ", InnervationZoneOffset(MOD(I, NumberOfEntries)+1), " from index ", &
      !  & MOD(I, NumberOfEntries)+1
    ENDDO
  ENDIF
  CLOSE(unit=5)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Read in motor unit fibre distribution
  !the MU (=motor unit) number the fiber belongs to
  FibreDistributionFile = TRIM(InputDirectory) // TRIM(FibreDistributionFile)
  INQUIRE(file=FibreDistributionFile, exist=FileExists)

  IF (.NOT. FileExists) THEN
    PRINT*, "Error: File """ // TRIM(FibreDistributionFile) // """ does not exist!"
    STOP
  ENDIF

  OPEN(unit=5, file=FibreDistributionFile, action="read", iostat=Status)
  READ(unit=5, fmt='(A)') Buffer

  ! Determine the number of entries in the first line
  NumberOfEntries = 1
  DO I = 0, LEN_TRIM(Buffer)
    IF (Buffer(I:I) == ',') THEN
      NumberOfEntries = NumberOfEntries + 1
    ENDIF
  ENDDO

  REWIND(UNIT=5)
  IF (ComputationalNodeNumber == 0) THEN
    PRINT*, "File  """ // TRIM(FibreDistributionFile) // """ contains ", NumberOfEntries, " Entries, NumberOfFibres=", &
      & NumberOfFibres, "."
  ENDIF

  ALLOCATE(MUDistribution(NumberOfFibres))

  IF (NumberOfFibres <= NumberOfEntries) THEN
    READ(unit=5, fmt=*, iostat=Status) MUDistribution(1:NumberOfFibres)
  ELSE
    READ(unit=5, fmt=*, iostat=Status) MUDistribution(1:NumberOfEntries)

    ! fill with already read data
    DO I = NumberOfEntries,NumberOfFibres
      MUDistribution(I) = MUDistribution(MOD(I, NumberOfEntries)+1)
      !PRINT*, "fill I=",I," with entry ", MUDistribution(MOD(I, NumberOfEntries)+1), " from index ", &
      !  & MOD(I, NumberOfEntries)+1
    ENDDO
  ENDIF
  CLOSE(unit=5)
END SUBROUTINE ReadInputFiles

SUBROUTINE ToLower(Str)
  CHARACTER(*), INTENT(INOUT) :: Str
  INTEGER(CMISSIntg) :: I
  DO I = 1, LEN(Str)
    SELECT CASE(Str(I:I))
      CASE("A":"Z")
        STR(I:I) = ACHAR(IACHAR(Str(I:I))+32)
    END SELECT
  ENDDO
END SUBROUTINE ToLower

SUBROUTINE ParseAssignment(Line, LineNumber, ScenarioInputFile)
  CHARACTER(LEN=*), INTENT(IN) :: Line
  INTEGER(CMISSIntg), INTENT(IN) :: LineNumber
  CHARACTER(LEN=*), INTENT(IN) :: ScenarioInputFile
  CHARACTER(LEN=1024) :: VariableName, StrValue
  INTEGER(CMISSIntg) :: StringLength
  LOGICAL :: OldTomoMechanics   

  !PRINT *, "Read line """//TRIM(Line)//"""."
  
  ! Extract variable name and value
  IF (INDEX(Line, "=") /= 0) THEN
    VariableName = TRIM(ADJUSTL(Line(1:INDEX(Line, "=")-1)))
    CALL ToLower(VariableName)
    StrValue = TRIM(ADJUSTL(Line(INDEX(Line, "=")+1:)))
    
    !PRINT *, "VariableName=["//TRIM(VariableName)//"], Value=["//TRIM(StrValue)//"]."
  ELSE
    RETURN
  ENDIF
  
  ! Parse variables
  SELECT CASE(VariableName)
    CASE ("x","numberglobalxelements")
      READ(StrValue, *, IOSTAT=Stat) NumberGlobalXElements
    CASE ("y","numberglobalyelements")
      READ(StrValue, *, IOSTAT=Stat) NumberGlobalYElements
    CASE ("z","numberglobalzelements")
      READ(StrValue, *, IOSTAT=Stat) NumberGlobalZElements
    CASE ("f","numberofinseriesfibres")
      READ(StrValue, *, IOSTAT=Stat) NumberOfInSeriesFibres
    CASE ("a","numberofelementsinatomicportionperdomain","numberofelementsinatomx","ax")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInAtomX
    CASE ("numberofelementsinatomy","ay")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInAtomY
    CASE ("numberofelementsinatomz","az")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInAtomZ
    CASE ("numberofelementsinsubdomainx","ex")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInSubdomainX
    CASE ("numberofelementsinsubdomainy","ey")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInSubdomainY
    CASE ("numberofelementsinsubdomainz","ez")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInSubdomainZ
    CASE ("numberofsubdomainsx","px")
      READ(StrValue, *, IOSTAT=Stat) NumberOfSubdomainsX
    CASE ("numberofsubdomainsy","py")
      READ(StrValue, *, IOSTAT=Stat) NumberOfSubdomainsY
    CASE ("numberofsubdomainsz","pz")
      READ(StrValue, *, IOSTAT=Stat) NumberOfSubdomainsZ
    CASE ("odesolverid")
      READ(StrValue, *, IOSTAT=Stat) ODESolverId
    CASE ("monodomainsolverid")
      READ(StrValue, *, IOSTAT=Stat) MonodomainSolverId
    CASE ("monodomainpreconditionerid")
      READ(StrValue, *, IOSTAT=Stat) MonodomainPreconditionerId
    CASE ("t","tend","timestop")
      READ(StrValue, *, IOSTAT=Stat) TimeStop
    CASE ("odetimestep")
      READ(StrValue, *, IOSTAT=Stat) ODETimeStep
    CASE ("pdetimestep")
      READ(StrValue, *, IOSTAT=Stat) PDETimeStep
    CASE ("elasticitytimestep")
      READ(StrValue, *, IOSTAT=Stat) ElasticityTimeStep
    CASE ("stimvalue")
      READ(StrValue, *, IOSTAT=Stat) StimValue
    CASE ("stimduration", "stimstop")
      READ(StrValue, *, IOSTAT=Stat) StimDuration
    CASE ("stimperiod", "periodd")
      READ(StrValue, *, IOSTAT=Stat) StimPeriod
    CASE ("outputtimestepstride")
      READ(StrValue, *, IOSTAT=Stat) OutputTimeStepStride
    CASE ("xi1", "numberofnodesinxi1")
      READ(StrValue, *, IOSTAT=Stat) NumberOfNodesInXi1
    CASE ("xi2", "numberofnodesinxi2")
      READ(StrValue, *, IOSTAT=Stat) NumberOfNodesInXi2
    CASE ("xi3", "numberofnodesinxi3")
      READ(StrValue, *, IOSTAT=Stat) NumberOfNodesInXi3
    CASE ("pretendednumberofdomainsfordomaindecomposition", "pretendednumberofdomains")
      READ(StrValue, *, IOSTAT=Stat) PretendedNumberOfDomainsForDomainDecomposition
    CASE ("newtonmaximumnumberofiterations")
      READ(StrValue, *, IOSTAT=Stat) NewtonMaximumNumberOfIterations
    CASE ("daerelativetolerance")
      READ(StrValue, *, IOSTAT=Stat) DAERelativeTolerance
    CASE ("daeabsolutetolerance")
      READ(StrValue, *, IOSTAT=Stat) DAEAbsoluteTolerance
    CASE ("newtontolerance")
      READ(StrValue, *, IOSTAT=Stat) NewtonTolerance
    CASE ("elasticityloopmaximumnumberofiterations")
      READ(StrValue, *, IOSTAT=Stat) ElasticityLoopMaximumNumberOfIterations
    CASE ("enableexportemg")
      READ(StrValue, *, IOSTAT=Stat) EnableExportEMG
    CASE ("oldtomomechanics")                   ! deprecated variable, use ModelType instead
      READ(StrValue, *, IOSTAT=Stat) OldTomoMechanics
      IF (OldTomoMechanics) THEN
        ModelType = 0
      ELSE
        ! If OldTomoMechanics is false, it means that ModelType=0 should not be used. 
        ! But is then ModelType=1 or ModelType=2 intended? Not clear (assume ModelType=1), thus deprecated.
        
        ModelType = 1
        PRINT *, "Warning: Parameter ""OldTomoMechanics=F"" was given, which is deprecated. Use ModelType=1 or ModelType=2 instead!"
      ENDIF
    CASE ("modeltype")
      READ(StrValue, *, IOSTAT=Stat) ModelType
    CASE ("splittingtype")
      READ(StrValue, *, IOSTAT=Stat) SplittingType
    CASE ("tklinparam")
      READ(StrValue, *, IOSTAT=Stat) TkLinParam
    CASE ("debuggingoutput")
      READ(StrValue, *, IOSTAT=Stat) DebuggingOutput
    CASE ("debuggingonlyrunshortpartofsimulation")
      READ(StrValue, *, IOSTAT=Stat) DebuggingOnlyRunShortPartOfSimulation
    CASE ("elasticitydisabled")
      READ(StrValue, *, IOSTAT=Stat) ElasticityDisabled
    CASE ("physicallength")
      READ(StrValue, *, IOSTAT=Stat) PhysicalLength
    CASE ("physicalwidth")
      READ(StrValue, *, IOSTAT=Stat) PhysicalWidth
    CASE ("physicalheight")
      READ(StrValue, *, IOSTAT=Stat) PhysicalHeight
    CASE ("physicalstimulationlength", "physicalstimulationwidth")
      READ(StrValue, *, IOSTAT=Stat) PhysicalStimulationLength
    CASE ("pmax")
      READ(StrValue, *, IOSTAT=Stat) PMax
    CASE ("conductivity")
      READ(StrValue, *, IOSTAT=Stat) Conductivity
    CASE ("am")
      READ(StrValue, *, IOSTAT=Stat) Am
    CASE ("cmfast")
      READ(StrValue, *, IOSTAT=Stat) CmFast
    CASE ("cmslow")
      READ(StrValue, *, IOSTAT=Stat) CmSlow
    CASE ("vmax")
      READ(StrValue, *, IOSTAT=Stat) Vmax
    CASE ("initialstretch")
      READ(StrValue, *, IOSTAT=Stat) InitialStretch
    CASE ("odensteps","nodesteps","nstepsode","nsteps")
      READ(StrValue, *, IOSTAT=Stat) OdeNSteps
    CASE ("pdensteps","npdesteps","nstepspde")
      READ(StrValue, *, IOSTAT=Stat) PdeNSteps
    CASE ("inputdirectory")
      InputDirectory = TRIM(ADJUSTL(StrValue))
      
      ! Append slash to input directory if necessary
      StringLength = LEN_TRIM(InputDirectory)

      IF (.NOT. InputDirectory(StringLength:StringLength) == "/") THEN
        InputDirectory(StringLength+1:StringLength+1) = "/"
      ENDIF
      InputDirectory = TRIM(InputDirectory)
    CASE ("firingtimesfile")
      FiringTimesFile = TRIM(ADJUSTL(StrValue))
    CASE ("innervationzonefile")
      InnervationZoneFile = TRIM(ADJUSTL(StrValue))
    CASE ("fibredistributionfile")
      FibreDistributionFile = TRIM(ADJUSTL(StrValue))
    CASE ("cellmlmodelfilename")
      CellMLModelFilename = TRIM(ADJUSTL(StrValue))
    CASE ("scenarioname")
      ScenarioName = TRIM(ADJUSTL(StrValue))
    CASE DEFAULT
      Read(LineNumber, *, IOSTAT=Stat) StrValue
      WRITE(*,'(A)') TRIM(ScenarioInputFile) // ":" // TRIM(StrValue) // ": Unrecognized variable """ // TRIM(VariableName) // """."
  END SELECT
  IF (Stat /= 0) THEN
    Read(LineNumber, *, IOSTAT=Stat) StrValue
    WRITE(*,'(A)') TRIM(ScenarioInputFile) // ": " // &
     &  "Could not parse value """ // TRIM(StrValue) // """ for variable """//TRIM(VariableName)//"""."            
  ENDIF
  
END SUBROUTINE ParseAssignment

SUBROUTINE ParseScenarioInputFile(ScenarioInputFile)
  CHARACTER(LEN=1024), INTENT(IN) :: ScenarioInputFile
  CHARACTER(LEN=10000) :: Line
  INTEGER(CMISSIntg) :: LineNumber, StringLength
  LineNumber = 0  
  
  OPEN(UNIT=100, FILE=ScenarioInputFile, ACTION="read", IOSTAT=Stat)
  
  IF (Stat /= 0) THEN
    PRINT*, "Could not read from Input file """ // TRIM(ScenarioInputFile) // """."
    RETURN
  ENDIF
  
  DO
    READ(UNIT=100, FMT='(A)', IOSTAT=Stat) Line
    LineNumber = LineNumber + 1
    
    IF (Stat < 0) EXIT  ! end of file reached
    
    ! Remove comments starting with # or !
    IF (INDEX(Line, "#") /= 0) THEN
      Line = TRIM(Line(1:INDEX(Line, "#")-1))
    ENDIF
    IF (INDEX(Line, "!") /= 0) THEN
      Line = Line(1:INDEX(Line, "!")-1)
    ENDIF
    
    CALL ParseAssignment(Line, LineNumber, ScenarioInputFile)
    
  ENDDO
  
END SUBROUTINE ParseScenarioInputFile

! Parse command line parameters and set numbers of elements
SUBROUTINE ParseParameters()

  INTEGER(CMISSLINTg) :: Factor, NumberArguments, ValueArgumentCount
  INTEGER(CMISSINTg) :: StringLength
  CHARACTER(LEN=10000) :: Arg
  CHARACTER(LEN=1024) :: ScenarioInputFile
  LOGICAL :: GeometryIsValid, FileExists
  INTEGER(CMISSIntg) :: I

  ! Default values
  NumberGlobalXElements = 2 ! ### PAPER SETTING
  NumberGlobalYElements = 2 ! ### PAPER SETTING
  NumberGlobalZElements = 2 ! ### PAPER SETTING
  NumberOfInSeriesFibres = 1 !1
  NumberOfElementsInAtomX = 1
  NumberOfElementsInAtomY = 1
  NumberOfElementsInAtomZ = 1

  ! loop over arguments
  NumberArguments = IARGC()
  ValueArgumentCount = 1
  DO I = 1, NumberArguments
    CALL GETARG(I, Arg)   ! get argument
    
    ! if argument ends in .sce, consider this as a scenario input file
    IF (INDEX(Arg, ".sce", .TRUE.) /= 0 .AND. INDEX(Arg, ".sce", .TRUE.) == LEN_TRIM(Arg)-3) THEN
      ScenarioInputFile = TRIM(Arg)
      CALL ParseScenarioInputFile(ScenarioInputFile)
      
    ! if argument contains '=' consider it as assignment
    ElSEIF (INDEX(Arg, "=") /= 0) THEN
      CALL ParseAssignment(Arg, I, "command line argument")
    ELSE
      SELECT CASE (ValueArgumentCount)     ! Input directory
      CASE(1)
        InputDirectory = Arg
      
        ! Append slash to input directory if necessary
        StringLength = LEN_TRIM(InputDirectory)

        IF (.NOT. InputDirectory(StringLength:StringLength) == "/") THEN
          InputDirectory(StringLength+1:StringLength+1) = "/"
        ENDIF
        InputDirectory = TRIM(InputDirectory)
      
      CASE(2)
        READ(Arg,*,Iostat=Stat)  NumberGlobalXElements
      CASE(3)
        READ(Arg,*,Iostat=Stat)  NumberGlobalYElements
      CASE(4)
        READ(Arg,*,Iostat=Stat)  NumberGlobalZElements
      CASE(5)
        READ(Arg,*,Iostat=Stat)  NumberOfInSeriesFibres
      CASE(6)
        READ(Arg,*,Iostat=Stat)  NumberOfElementsInAtomX
      CASE(7)
        READ(Arg,*,Iostat=Stat)  ODESolverId
      CASE(8)
        READ(Arg,*,Iostat=Stat)  MonodomainSolverId
      CASE(9)
        READ(Arg,*,Iostat=Stat)  MonodomainPreconditionerId
      CASE(10)
        READ(Arg,*,Iostat=Stat)  OdeNSteps       
      CASE(11)
        READ(Arg,*,Iostat=Stat)  SplittingType 
      CASE(12)
        READ(Arg,*,Iostat=Stat)  PDETimeStep
      CASE(13)
        READ(Arg,*,Iostat=Stat)  ODETimeStep
        
      ENDSELECT
      
      ValueArgumentCount = ValueArgumentCount + 1
    
    ENDIF
  ENDDO
  
  IF (NumberArguments == 0) THEN
    PRINT*, "Using default values. " // NEW_LINE('A') &
     & // "Usage: " // NEW_LINE('A') // &
     & "1)   ./cuboid [<input folder> [<X> <Y> <Z> [<F> [<NumberOfElementsInAtomX> " // &
     & "[<ODEsolver> [<Msolver> [<Mprecond> ]]]]]]] " // NEW_LINE('A') // &
     & "2)   ./cuboid [<variable>=<value>] [<input.sce>] [<variable>=<value>] " // NEW_LINE('A') // &
     & "     See the example scenario file for file format and variable names. Variables will be set in order of the arguments."
  ENDIF
!#################### the following is only necessary to correct user output, i think. ############################################
  IF (PdeNSteps/=-1) THEN
    PDETimeStep = ElasticityTimeStep/PdeNSteps
  ENDIF
  IF (OdeNSteps/=-1) THEN
    ODETimeStep = PDETimeStep/OdeNSteps
  END IF
!##################################################################################################################################

  ! direction of fibres is in Xi1=Global X direction
  NumberOfElementsFE = NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements

  NumberGlobalYFibres = NumberOfNodesInXi2 * NumberGlobalYElements
  NumberOfFibreLinesPerGlobalElement = NumberOfNodesInXi2 * NumberOfNodesInXi3
  NumberOfGlobalElementLines = NumberGlobalYElements * NumberGlobalZElements
  NumberOfFibreLinesTotal = NumberOfFibreLinesPerGlobalElement * NumberOfGlobalElementLines
  NumberOfFibres = NumberOfFibreLinesTotal * NumberOfInSeriesFibres
  NumberOfElementsMInXi1 = NumberOfNodesInXi1 - 1
  NumberOfElementsMPerFibreLine = NumberOfElementsMInXi1 * NumberGlobalXElements
  NumberOfElementsMPerFibre = NumberOfElementsMPerFibreLine / NumberOfInSeriesFibres
  NumberOfNodesPerShortFibre = NumberOfElementsMPerFibre
  NumberOfNodesPerLongFibre = NumberOfElementsMPerFibre + 1
  ! the end point nodes of fibres can not be shared between different fibres, therefore some fibres have 1 more node (the ones at the left)

  ! total number of 1D nodes
  NumberOfNodesMPerFibreLine = NumberOfNodesPerShortFibre * (NumberOfInSeriesFibres-1) + NumberOfNodesPerLongFibre
  NumberOfNodesM = NumberOfNodesMPerFibreLine * NumberOfFibreLinesTotal
  NumberOfElementsM = NumberOfElementsMPerFibre * NumberOfFibres

  ! compute number of bioelectric nodes that will be stimulated
  NumberStimulatedNodesPerFibre = MAX(1, NINT(DBLE(PhysicalStimulationLength) * (NumberOfElementsMPerFibre / PhysicalWidth)))
  
  ! previous implementation changed:: StimValuePerNode = StimValue / NumberStimulatedNodesPerFibre to::
  StimValuePerNode = StimValue
  ! ^ This was done because we need to keep the stimulation in a cell constant, no matter, how many cells there are. the dynamic of a cell is independent of the spacial discretization. Checked by Nehzat, Thomas, Aaron
!##################################################################################################################################
!  fast_twitch=.true.
!  if(fast_twitch) then
!  pathname="/home/heidlauf/OpenCMISS/OpenCMISS/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/"
!  CellMLModelFilename=trim(pathname)//"fast_2014_03_25_no_Fl_no_Fv.xml" !FAST
 
  IF (TRIM(CellMLModelFilename) == "standard") THEN
    IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
      CellMLModelFilename = TRIM(inputDirectory) // "new_slow_TK_2014_12_08.cellml"
      !StimValue = 20000.0_CMISSRP !700.0_CMISSRP!700.0_CMISSRP

    ELSEIF (ModelType == 1) THEN ! 3, , "MultiPhysStrain", numerically more stable
      CellMLModelFilename = TRIM(inputDirectory) // "Aliev_Panfilov_Razumova_2016_08_22.cellml"
      !StimValue = 90.0_CMISSRP !90.0_CMISSRP

    ELSEIF (ModelType == 2) THEN  ! 4, "Titin"
      CellMLModelFilename = TRIM(inputDirectory) // "Aliev_Panfilov_Razumova_Titin_2016_10_10_noFv.cellml"
      
    ENDIF
  ENDIF

  ! check if file exists
  INQUIRE(file=CellMLModelFilename, exist=FileExists)
  IF (.NOT. FileExists) THEN
    PRINT*, "Error: CellML file """ // TRIM(CellMLModelFilename) // """ does not exist!"
    STOP
  ENDIF


!   &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/shorten_mod_2011_07_04.xml"
!    pathname="/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles"
!    CellMLModelFilename=trim(pathname)//"/fast_shortening_0.1vmax.xml"
!    StimValue=700.0_CMISSRP!2000.0_CMISSRP!700.0_CMISSRP
!  else !slow twitch
!    filename2= &
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/fast_stim_2012_07_23.xml"
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_2012_07_23.xml_0.401"
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_twitch_2012_01_27.xml"
!    StimValue=2000.0_CMISSRP
!  endif
!##################################################################################################################################

 !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  ! set diagnostics
!  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"dignostics",["FIELD_MAPPINGS_CALCULATE"],err)

  !                          type                level list,  output file  routine list
  !CALL cmfe_DiagnosticsSetOn(cmfe_FROM_DIAG_TYPE, [3], "", ["cmfe_Problem_Solve", "PROBLEM_SOLVE     "], Err)

  !CALL cmfe_DiagnosticsSetOn(cmfe_CALL_STACK_DIAG_TYPE, [5], "diagnostics", ["FIELD_MAPPINGS_CALCULATE"], Err)

  !CALL cmfe_OutputSetOn("output.txt", Err)
  
  !CMFE_ALL_DIAG_TYPE    !<Type for setting diagnostic output in all routines \see OPENCMISS_DiagnosticTypes,OPENCMISS
  !CMFE_IN_DIAG_TYPE     !<Type for setting diagnostic output in one routine \see OPENCMISS_DiagnosticTypes,OPENCMISS
  !CMFE_FROM_DIAG_TYPE   !<Type for setting diagnostic output in one routine downwards \see OPENCMISS_DiagnosticTypes,OPENCMISS
  !                          in which routine,   levelList(:), diagFilename
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE, [1],          "",&
  ! routine list
  !  & ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"], Err)
  !  & ["FIELD_MAPPINGS_CALCULATE"], Err)
  
  
  !                     type                  not output directly, filename
  ! cmfe_IN_TIMING_TYPE, cmfe_FROM_TIMING_TYPE, cmfe_ALL_TIMING_TYPE
  !CALL cmfe_TimingSetOn(cmfe_ALL_TIMING_TYPE, .TRUE.,             "", ["cmfe_Problem_Solve", "PROBLEM_SOLVE     "], Err)

  CALL cmfe_OutputSetOn("EMG",Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains = NumberOfComputationalNodes
  IF (PretendedNumberOfDomainsForDomainDecomposition == 0) THEN
    PretendedNumberOfDomainsForDomainDecomposition = NumberOfDomains
  ENDIF
  
  IF (NumberOfDomains > 0) THEN
    CALL ComputeSubdomainsWithAtoms()   ! compute domain decomposition for 3D mesh considering the atoms
  ENDIF
  GeometryIsValid = CheckGeometry()

  ! print the current directory from which the program was launched
  CALL PrintWorkingDirectory()

  ! Read in input files, stops execution if files do not exist
  CALL ReadInputFiles()

  ! output scenario information
  IF (ComputationalNodeNumber == 0) THEN
    PRINT *, "CellML file: """ // TRIM(CellMLModelFilename) // """"
    PRINT *, "Scenario name: """ // TRIM(ScenarioName) // """"
    PRINT *, ""
    PRINT *, "---------- Timing parameters -----------------------------------------------"
    PRINT *, "The time unit is 1 ms."
    PRINT "(A,F7.2,A,F5.2,A,F5.2)", "  Main loop, Δt = ", TimeStop, ", dt = ", ElasticityTimeStep
    PRINT "(A,F5.2)", "  - stimulation enabled:  Δt = ", StimDuration
    PRINT "(A,F5.2)", "  - stimulation disabled: Δt = ", (StimPeriod - StimDuration)
    PRINT *, ""

    PRINT "(A,F7.2,A,F0.5,A,I5)", "- MAIN_TIME_LOOP,         Δt =", TimeStop, ", dt = ", ElasticityTimeStep, &
      & ", N. Iter: ", CEILING(TimeStop/ElasticityTimeStep)
    
    
    IF (PdeNSteps /= -1) THEN
      PRINT "(A,F0.4,A,F0.5,2(A,I5),A)", "  - MONODOMAIN_TIME_LOOP, Δt = ", ElasticityTimeStep, ",  dt = ", PDETimeStep,&
        & ", N. Iter: ", CEILING(ElasticityTimeStep/PDETimeStep), " (set by PdeNSteps=", PdeNSteps, ")"
    ELSE
      PRINT "(A,F0.4,A,F0.5,A,I5)", "  - MONODOMAIN_TIME_LOOP, Δt = ", ElasticityTimeStep, ",  dt = ", PDETimeStep,&
        & ", N. Iter: ", CEILING(ElasticityTimeStep/PDETimeStep)
    ENDIF
      
    IF (OdeNSteps /= -1) THEN
      PRINT "(A,F0.5,A,I5,A,I5,A)", "    - SolverDAE,                       dt = ", ODETimeStep, &
        & ", N. Iter: ", CEILING(PDETimeStep/ODETimeStep), " (set by OdeNSteps=", OdeNSteps, ")"
    ELSE 
      PRINT "(A,F0.5,A,I5,A,I5,A)", "    - SolverDAE,                       dt = ", ODETimeStep, &
        & ", N. Iter: ", CEILING(PDETimeStep/ODETimeStep)
    ENDIF
    
    PRINT "(A,F0.4)", "    - SolverParabolic"
    
    IF (ElasticityDisabled) THEN
      PRINT "(A)",               "  - ELASTICITY_LOOP     (disabled) "
    ELSE
      PRINT "(A,I5)",               "  - ELASTICITY_LOOP,                                N. Iter: ",&
        & ElasticityLoopMaximumNumberOfIterations
      PRINT "(A,I4,A,E10.4)", "    - SolverFE,                 N. Iter (max): ", NewtonMaximumNumberOfIterations, &
        & ", Tol.: ",NewtonTolerance
      PRINT "(A,I4)", "      - LinearSolverFE"
      PRINT *, ""
    ENDIF
    
    IF (DebuggingOnlyRunShortPartOfSimulation) PRINT *, "Abort after first stimulation."
    PRINT *, "OutputTimeStepStride:", OutputTimeStepStride, ", EnableExportEMG:", EnableExportEMG

    ! ElasticityTimeStep should be a integer multiple of StimDuration

    ! Output problem size information
    PRINT *, ""
    PRINT *, "---------- Problem size parameters ------------------------------------------"

    PRINT "(A,3(I6,A),I12)", "# global FE-elements:      ", NumberGlobalXElements, ", ", NumberGlobalYElements, ", ", &
      & NumberGlobalZElements, &
      & ", Total: ", NumberOfElementsFE
    PRINT "(A,3(I6,A),I12)", "# local nodes per element: ", NumberOfNodesInXi1, ", ", NumberOfNodesInXi2, ", ", NumberOfNodesInXi3,&
      & ", Total: ", NumberOfNodesInXi1*NumberOfNodesInXi2*NumberOfNodesInXi3
    !PRINT "(A,I8)", "                  NumberOfInSeriesFibres: ", NumberOfInSeriesFibres
    PRINT "(A,I8)", "      NumberOfFibreLinesPerGlobalElement: ", NumberOfFibreLinesPerGlobalElement
    PRINT "(A,I8)", "              NumberOfGlobalElementLines: ", NumberOfGlobalElementLines
    !PRINT "(A,I8)", "                 NumberOfFibreLinesTotal: ", NumberOfFibreLinesTotal
    PRINT "(A,I8)", "                          NumberOfFibres: ", NumberOfFibres
    PRINT "(A,I8)", "               NumberOfElementsMPerFibre: ", NumberOfElementsMPerFibre
    PRINT "(A,I8)", "               NumberOfNodesPerLongFibre: ", NumberOfNodesPerLongFibre
    !PRINT "(A,I8)", "           NumberOfElementsMPerFibreLine: ", NumberOfElementsMPerFibreLine
    PRINT "(A,I8)", "                       NumberOfElementsM: ", NumberOfElementsM
    !PRINT "(A,I8)", "              NumberOfNodesPerShortFibre: ", NumberOfNodesPerShortFibre
    !PRINT "(A,I8)", "              NumberOfNodesMPerFibreLine: ", NumberOfNodesMPerFibreLine
    PRINT "(A,I8)", "                          NumberOfNodesM: ", NumberOfNodesM
    PRINT *,""
    PRINT "(3(A,I5),A)", "              Non-decomposable atom: ", NumberOfElementsInAtomX, &
      & "x", NumberOfElementsInAtomY, "x", NumberOfElementsInAtomZ, " elements"
    PRINT "(A,3(I5,A))", "                        Subdomain e: ", NumberOfElementsInSubdomainX, "x", &
      & NumberOfElementsInSubdomainY, "x", &
      & NumberOfElementsInSubdomainZ, " elements"
    PRINT "(A,4(I5,A))", "             Domain decomposition p: ", NumberOfSubdomainsX, "x", NumberOfSubdomainsY, "x", &
      & NumberOfSubdomainsZ, &
      & " subdomains"
    PRINT "(2(A,I5))", " Number of initial domains for domain decomposition: ", PretendedNumberOfDomainsForDomainDecomposition, &
      & ", number of processes: ", NumberOfDomains
    PRINT "(A,I5)", " Number of different processes for a fibre: ", NumberOfSubdomainsX
    
    PRINT *, ""
    PRINT *, "---------- Physical parameters -----------------------------------------------"
    PRINT "(4(A,F5.2))", "       Dimensions [cm]: ",PhysicalLength,"x",PhysicalWidth,"x",PhysicalHeight,&
     & ",          InitialStretch: ", InitialStretch
    IF (NumberStimulatedNodesPerFibre == 1) THEN
      PRINT "(A,F11.2,A,F8.5,A,I3,A,F9.2,A)", " Stimulation [uA/cm^2]: ",StimValue," on ",PhysicalStimulationLength," cm, i.e. ",&
        & NumberStimulatedNodesPerFibre," node"
    ELSE
      PRINT "(A,F11.2,A,F8.5,A,I3,A,F9.2,A)", " Stimulation [uA/cm^2]: ",StimValue," on ",PhysicalStimulationLength," cm, i.e. ",&
        & NumberStimulatedNodesPerFibre," nodes, ",StimValuePerNode, " per node"
    ENDIF
    PRINT "(3(A,F5.2))", " PMax:", PMax, ",        VMax: ", VMax, ",   Conductivity: ", Conductivity
    PRINT "(3(A,F7.2))", "   Am:", Am, ", Cm (fast):", CmFast, ",     Cm (slow):", CmSlow
    
    IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
      PRINT *, "ModelType: 0 (""3a"", MultiPhysStrain, Old mechanics formulation that works in parallel.)"
    ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
      PRINT *, "ModelType: 1 (""3"", MultiPhysStrain)"
    ELSEIF (ModelType == 2) THEN  ! 4, "Titin"
      PRINT *, "ModelType: 2 (""4"", Titin)"
    ENDIF

    PRINT *, ""
    PRINT *, "---------- Solvers -----------------------------------------------------------"
  
  
    SELECT CASE(SplittingType)
      CASE(0)
        PRINT *, "Splitting:       0 Godunov"
      CASE(1)
        PRINT *, "Splitting:       1 Strang"
      CASE DEFAULT
        PRINT *, "Splitting:       0 Godunov (default)"
    END SELECT
  
    SELECT CASE(ODESolverId)
      CASE(1)
        PRINT *, "(0D) ODE:        1 Explicit Euler"
      CASE(2)
        PRINT *, "(0D) ODE:        2 BDF"
        PRINT *, "                   Absolute Tolerance: ",DAEAbsoluteTolerance
        PRINT *, "                   Relative Tolerance: ",DAERelativeTolerance
      CASE(3)
        PRINT *, "(0D) ODE:        3 General Linear (GL)"
      CASE(4)
        PRINT *, "(0D) ODE:        4 Crank-Nicolson"
      CASE(5)
        PRINT *, "(0D) ODE:        5 Improved Euler (Heun's method)"
      CASE DEFAULT
        PRINT *, "(0D) ODE:        Explicit Euler (default)" 
    END SELECT

    SELECT CASE(MonodomainSolverId)
      CASE(1)
        PRINT *, "(1D) Monodomain: 1 SOLVER_DIRECT_LU"
      CASE(2)
        PRINT *, "(1D) Monodomain: 2 SOLVER_ITERATIVE_GMRES"
      CASE(3)
        PRINT *, "(1D) Monodomain: 3 SOLVER_ITERATIVE_CONJUGATE_GRADIENT"
      CASE(4)
        PRINT *, "(1D) Monodomain: 4 SOLVER_ITERATIVE_CONJGRAD_SQUARED"
      CASE DEFAULT
        PRINT *, "(1D) Monodomain: SOLVER_ITERATIVE_GMRES (default)" 
    END SELECT
    
    SELECT CASE(MonodomainPreconditionerId)
      CASE(1)
        PRINT *, "                 1 NO_PRECONDITIONER"
      CASE(2)
        PRINT *, "                 2 JACOBI_PRECONDITIONER"
      CASE(3)
        PRINT *, "                 3 BLOCK_JACOBI_PRECONDITIONER"
      CASE(4)
        PRINT *, "                 4 SOR_PRECONDITIONER"
      CASE(5)
        PRINT *, "                 5 INCOMPLETE_CHOLESKY_PRECONDITIONER"
      CASE(6)
        PRINT *, "                 6 INCOMPLETE_LU_PRECONDITIONER"
      CASE(7)
        PRINT *, "                 7 ADDITIVE_SCHWARZ_PRECONDITIONER"
      CASE DEFAULT
        PRINT *, "                 NO_PRECONDITIONER (default)"
    END SELECT
    PRINT *, "------------------------------------------------------------------------------"
    PRINT *, ""

    ! Output some static (compile-time) settings
    CALL cmfe_CustomProfilingGetEnabled(CustomProfilingEnabled, TAUProfilingEnabled, Err)

    IF (CustomProfilingEnabled) THEN
      PRINT*, "Custom Profiling is enabled."
    ELSE
      PRINT*, "Custom Profiling is disabled. (Enable with -DUSE_CUSTOM_PROFILING)"
    ENDIF

    IF (TAUProfilingEnabled) THEN
      PRINT*, "TAU Profiling is enabled."
    ELSE
      PRINT*, "TAU Profiling is disabled."
    ENDIF

    IF (.NOT. GeometryIsValid) THEN
      PRINT*, "Parametrization is invalid. Aborting."
    ENDIF
  ENDIF

  IF (.NOT. GeometryIsValid) THEN
    STOP
  ENDIF

  IF (ODESolverId/=5 .AND. SplittingType==1) THEN
    PRINT *, "Strang-Splitting must be used with Improved Euler method! &
     & Use ODESolverId=5 instead."
    STOP
  ENDIF
  
END SUBROUTINE ParseParameters

SUBROUTINE PrintWorkingDirectory()

  IF (ComputationalNodeNumber == 0) THEN
    CALL SYSTEM("pwd > pwd.txt", Stat)
    IF (Stat == 0) THEN
      OPEN(UNIT=100, FILE="pwd.txt", ACTION="read", IOSTAT=Stat)
      IF (Stat == 0) THEN
        READ(100,"(A)", IOSTAT=Stat) WorkingDirectory
        CLOSE(UNIT=100)
        IF (Stat == 0) THEN
          PRINT*, "Working Directory: """ // TRIM(WorkingDirectory) // """."
        ELSE
          PRINT*, "Error reading pwd.txt"
        ENDIF
      ELSE
        PRINT*, "Error opening pwd.txt"
      ENDIF
    ELSE
      PRINT*, "Error calling 'pwd'"
    ENDIF
  ENDIF
END SUBROUTINE PrintWorkingDirectory
  
  
END MODULE SETUP_PARAMETERS