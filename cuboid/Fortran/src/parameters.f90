
MODULE PARAMETERS

  USE OpenCMISS
  USE OpenCMISS_Iron
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
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Test program parameters
  LOGICAL :: DebuggingOnlyRunShortPartOfSimulation = .FALSE.    ! only run one timestep of MAIN_LOOP with stimulus
  LOGICAL, PARAMETER :: DEBUGGING_PARALLEL_BARRIER = .FALSE.   ! execute a barrier where all processes wait when running in parallel, this is helpful to only debug the execution of single processes in a parallel scenario with gdb
  LOGICAL, PARAMETER :: DEBUGGING_PROBLEM_OUTPUT = .FALSE.     ! output the 'solver' object after it is created
  LOGICAL :: DebuggingOutput = .FALSE.    ! enable information from solvers
  
  !--------------------------------------------------------------------------------------------------------------------------------

  INTEGER(CMISSIntg) :: ModelType = 0 ! ### PAPERBRANCH SETTING     ! type of the model (was OldTomoMechanics): 0 = "3a","MultiPhysStrain", old version of tomo that works in parallel, 1 = "3","MultiPhysStrain", new version of tomo that is more stable in numerical sense, 2 = "4","Titin"
  INTEGER(CMISSIntg) :: SplittingType = 0   ! 0 = godunov splitting, 1 = strang splitting
  LOGICAL :: EnableExportEMG = .FALSE.
  LOGICAL :: ElasticityDisabled = .FALSE. ! do create the elasticity control loop
  
  ! physical dimensions in [cm]
  REAL(CMISSRP) :: PhysicalLength=1.0_CMISSRP ! ### PAPERBRANCH SETTING !    X-direction
  REAL(CMISSRP) :: PhysicalWidth =1.0_CMISSRP ! ### PAPERBRANCH SETTING !    Y-direction
  REAL(CMISSRP) :: PhysicalHeight=1.0_CMISSRP ! ### PAPERBRANCH SETTING !    Z-direction
  REAL(CMISSRP) :: PhysicalStimulationLength = 0.03125_CMISSRP  ! X-direction   ### PAPERBRANCH SETTING: value 0.03125 is unphysical (roughly 20 times too large). !NMJ area: 200 (um)² -> NMJ diameter: 16 um = 0.0016cm. Based on Tse et al., 2014, The Neuromuscular Junction: Measuring Synapse Size, Fragmentation and Changes in Synaptic Protein Density Using Confocal Fluorescence Microscopy
  
  !all times in [ms]
  REAL(CMISSRP) :: time
  REAL(CMISSRP) :: StimPeriod=1.0_CMISSRP ! ### PAPERBRANCH SETTING
  REAL(CMISSRP) :: TimeStop=10.0_CMISSRP ! ### PAPERBRANCH SETTING

  ! time discretizatzion settings
  REAL(CMISSRP) :: ODETimeStep = 0.0001_CMISSRP ! ### PAPERBRANCH SETTING: 0.0001
  REAL(CMISSRP) :: PDETimeStep = 0.0005_CMISSRP ! ### PAPERBRANCH SETTING: 0.0005
  REAL(CMISSRP) :: ElasticityTimeStep = 0.1_CMISSRP ! ### PAPERBRANCH SETTING
  INTEGER(CMISSIntg) :: OdeNSteps = -1 ! can be used to set ODETimeStep implicitly.
  INTEGER(CMISSIntg) :: PdeNSteps = -1  ! overrides PDETimeStep

  REAL(CMISSRP) :: StimDuration=0.1_CMISSRP ! ### PAPERBRANCH SETTING ! should be the same as ElasticityTimeStep

  INTEGER(CMISSIntg)  :: OutputTimestepStride=10  ! (10)

  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------

  !stimulation current in [uA/cm^2]
  REAL(CMISSRP) :: StimValue = 1200.0_CMISSRP ! ### PAPERBRANCH SETTING ! will be applied to all nodes in PhysicalStimulationLength.

  REAL(CMISSRP) :: PMax=7.3_CMISSRP ! ### PAPERBRANCH SETTING ! N/cm^2

  !conductivity in [mS/cm] - This is sigma
  REAL(CMISSRP) :: Conductivity=3.828_CMISSRP ! ### PAPERBRANCH SETTING

  !surface area to volume ratio in [cm^-1]
  REAL(CMISSRP) :: Am=500.0_CMISSRP ! ### PAPERBRANCH SETTING

  !membrane capacitance in [uF/cm^2]
  REAL(CMISSRP) :: CmFast=0.58_CMISSRP ! ### PAPERBRANCH SETTING: Which one is needed? fast/slow? I guess slow.
  REAL(CMISSRP) :: CmSlow=0.58_CMISSRP ! ### PAPERBRANCH SETTING: Which one is needed? fast/slow? I guess slow.

  !Material-Parameters C=[mu_1, mu_2, mu_3, alpha_1, alpha_2, alpha_3, mu_0, XB_stiffness (= e or \eta in different papers)]
  REAL(CMISSRP), PARAMETER, DIMENSION(8) :: C = &
    & [0.0085_CMISSRP*5.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP, & ! [N/cm^2 = 10^4 J/m^3]            ! ### NOT RELEVANT FOR PAPER (other scenario)
    &  11.0_CMISSRP,1.0_CMISSRP,6.0_CMISSRP, &                                                  ! ### NOT RELEVANT FOR PAPER
    &  1.0_CMISSRP,2.2e-9_CMISSRP]                                                              ! ### NOT RELEVANT FOR PAPER

  !maximum contraction velocity in [cm/ms] not relevant for paper
  REAL(CMISSRP) :: VMax=-0.02_CMISSRP ! =0.2 m/s, rat GM                                        ! ### NOT RELEVANT FOR PAPER (isometric conditions)

  !CAUTION - what are the units???   ![N/cm^2]?
!  REAL(CMISSRP), PARAMETER, DIMENSION(4) :: MAT_FE= &
!    &[0.0000000000635201_CMISSRP,0.3626712895523322_CMISSRP,0.0000027562837093_CMISSRP,43.372873938671383_CMISSRP] ! ### SKIPPED due to citation differences

  REAL(CMISSRP), PARAMETER, DIMENSION(4) :: MAT_FE = [0.00000000006352_CMISSRP,0.3627_CMISSRP,0.000002756_CMISSRP,43.373_CMISSRP] ! ### PAPERBRANCH SETTING

  REAL(CMISSRP) :: TkLinParam=1.0_CMISSRP ! 1: With Actin-Tintin Interaction 0: No Actin-Titin Interactions    ! ### NOT RELEVANT FOR PAPER

  !Inital Conditions
  REAL(CMISSRP) :: InitialStretch=1.0_CMISSRP   ! previous value in new mechanical description: 1.2_CMISSRP
  
  INTEGER(CMISSIntg) :: ElasticityLoopMaximumNumberOfIterations = 5     ! only relevant if InitialStretch /= 1.0
  INTEGER(CMISSIntg) :: NewtonMaximumNumberOfIterations = 500           ! 
  REAL(CMISSRP) :: NewtonTolerance = 1.E-8_CMISSRP                      ! 
  REAL(CMISSRP) :: DAERelativeTolerance = 1.E-7_CMISSRP, DAEAbsoluteTolerance = 1.E-7_CMISSRP   ! ### Something we dont have to tell about.

  ! ### PAPERBRANCH SETTING DESCRIPTION:
  ! we have 36 fibers spread over 2X2X2 FEM elements. Each fiber has 32 elements, 33 nodes (cells). The 16nth cell is stimulated (counting from 1 to 32):
  INTEGER(CMISSIntg) :: NumberGlobalXElements      ! ### PAPERBRANCH SETTING, set later!
  INTEGER(CMISSIntg) :: NumberGlobalYElements      ! ### PAPERBRANCH SETTING, set later!
  INTEGER(CMISSIntg) :: NumberGlobalZElements      ! ### PAPERBRANCH SETTING, set later!
  INTEGER(CMISSIntg) :: NumberOfInSeriesFibres = 1 ! ### has no good effect. (set /=1 to break code) TO BE DISCUSSED? Only Benni will need this. Default 1 makes sense.
  INTEGER(CMISSIntg) :: NumberOfNodesInXi1 = 16    ! ### PAPERBRANCH SETTING, number of 1D elements per 3D element, i.e. number of nodes is +1 (but nodes on edges of 3d elements are shared)
  INTEGER(CMISSIntg) :: NumberOfNodesInXi2 = 3     ! ### PAPERBRANCH SETTING
  INTEGER(CMISSIntg) :: NumberOfNodesInXi3 = 3     ! ### PAPERBRANCH SETTING
    
  INTEGER(CMISSLIntg) :: NumberOfElementsFE
  INTEGER(CMISSIntg) :: NumberOfNodesM
  INTEGER(CMISSIntg) :: NumberOfElementsM
  INTEGER(CMISSIntg) :: NumberOfFibres
  INTEGER(CMISSIntg) :: NumberOfNodesPerLongFibre   ! fibre that touches right boundary has one additional electricity node
  INTEGER(CMISSIntg) :: NumberOfNodesPerShortFibre  ! the number of nodes on ordinary fibres not lying on the rightt boundary
  INTEGER(CMISSIntg) :: NumberOfElementsInAtomX
  INTEGER(CMISSIntg) :: NumberOfElementsInAtomY
  INTEGER(CMISSIntg) :: NumberOfElementsInAtomZ
  INTEGER(CMISSIntg) :: NumberOfElementsInSubdomainX = 0  ! if zero it will be computed from size of atoms
  INTEGER(CMISSIntg) :: NumberOfElementsInSubdomainY = 0
  INTEGER(CMISSIntg) :: NumberOfElementsInSubdomainZ = 0
  INTEGER(CMISSIntg) :: NumberOfElementsMInXi1
  INTEGER(CMISSIntg) :: NumberGlobalYFibres
  INTEGER(CMISSINTg) :: NumberOfFibreLinesPerGlobalElement
  INTEGER(CMISSIntg) :: NumberOfGlobalElementLines
  INTEGER(CMISSIntg) :: NumberOfFibreLinesTotal
  INTEGER(CMISSIntg) :: NumberOfElementsMPerFibre ! this is the total number of (0D sub-)cells (nodes) per global fibre.
  INTEGER(CMISSIntg) :: NumberOfElementsMPerFibreLine
  INTEGER(CMISSIntg) :: NumberOfNodesMPerFibreLine
  INTEGER(CMISSIntg) :: NumberOfSubdomainsX = 0, NumberOfSubdomainsY = 0, NumberOfSubdomainsZ = 0    ! if NumberOfElementsInSubdomain is set, NumberOfSubdomains will also be used. Otherwise is will be computed.
  INTEGER(CMISSIntg) :: NumberOfAtomsPerSubdomainX, NumberOfAtomsPerSubdomainY, NumberOfAtomsPerSubdomainZ
  INTEGER(CMISSIntg) :: NumberOfAtomsLastSubdomainX, NumberOfAtomsLastSubdomainY, NumberOfAtomsLastSubdomainZ
  INTEGER(CMISSIntg) :: NumberOfElementsLastAtomX, NumberOfElementsLastAtomY, NumberOfElementsLastAtomZ
  INTEGER(CMISSIntg) :: NumberStimulatedNodesPerFibre
  INTEGER(CMISSIntg) :: PretendedNumberOfDomainsForDomainDecomposition = 0
  
  INTEGER(CMISSIntg) :: Stat
  CHARACTER(len=256) :: CellMLModelFilename = "standard" ! standard will be replaced by the standard model file
  CHARACTER(len=1024) :: inputDirectory = "input/"
  CHARACTER(len=1024) :: FiringTimesFile = "MU_firing_times_10s.txt"
  CHARACTER(len=1024) :: InnervationZoneFile = "innervation_zone_18.txt"
  CHARACTER(len=1024) :: FibreDistributionFile = "MU_fibre_distribution_4050.txt"
  CHARACTER(len=256) :: MemoryConsumption1StTimeStep = "", MemoryConsumptionBeforeSim, Temp
  CHARACTER(len=10000) :: WorkingDirectory
  CHARACTER(len=1024) :: ScenarioName = ""

  
END MODULE PARAMETERS