# Scenario name: default values
# Date: 24.4.17, 6.7.17, 4.9.17
# Author: Benjamin Maier
# dihu-stuttgart/iron version: b54928c2f5f9e903fc71f0bf02a730727e670d47
# dihu-stuttgart/iron-examples version: 79c9113d74e8c4d64e26b89653e35be737e18997
#
# Designated number of Processors: 1
# Estimated run time (walltime):
#
# Scenario file for cuboid example
#
# Comment (sub-)lines start with # or !
# Specify assignments as Variable=Value, logical values as `T` or `F`.
# Variables are case-insensitive, whitespace between specifiers is ignored.
# Pass this file's name to cuboid as command argument. It must end in .sce and can be mixed with the usual arguments which are:
# ./cuboid <input folder> <X> <Y> <Z> <F> <A> <ODESolverId> <MonodomainSolverId> <MonodomainPreconditionerId>
# where X,Y,Z = NumberGlobal{X,Y,Z}Elements, A = NumberOfElementsInAtomicPortionPerDomain, F = NumberOfInSeriesFibres
# Command line arguments can also be an assignment, like `DebuggingOutput=T` or `X=1`.
# The command line arguments are processed in the order they occur, when this file appears, all its assignments are parsed.
# Therefore values can be overwritten when the filename is the first argument.

# ------------- discretization ----------------
# orientation: X = Xi1 = length, Y = Xi2 = width, Z = Xi3 = height
# fibres are parallel to X axis

NumberGlobalXElements = 2       # number of finite elasticty elements in x-direction, alias x
NumberGlobalYElements = 2       # number of finite elasticty elements in y-direction, alias y
NumberGlobalZElements = 2       # number of finite elasticty elements in z-direction, alias z

NumberOfNodesInXi1 = 16         # number of bioelectric elements (i.e. number nodes+1) per 3D FE element in direction of fibre (X direction), alias xi1
NumberOfNodesInXi2 = 3          # number of fibres per FE element in Y direction, alias xi2
NumberOfNodesInXi3 = 3          # number of fibres per FE element in Z direction, alias xi3

NumberOfInSeriesFibres = 1      # number of fibres that are in a series and mechanically connected. This is not completely tested, set to 1, alias f

NumberOfElementsInAtomX = 1     # x-size of an non-decomposable "atom" of finite elasticity elements, that are guaranteed to be on the same subdomain, alias ax
NumberOfElementsInAtomY = 1     # y-size of an non-decomposable "atom" of finite elasticity elements, that are guaranteed to be on the same subdomain, alias ay
NumberOfElementsInAtomZ = 1     # z-size of an non-decomposable "atom" of finite elasticity elements, that are guaranteed to be on the same subdomain, alias az

# ------------- debugging --------------------
OutputTimeStepStride = 10       # write output *.exnode files every nth timestep
EnableExportEMG = F             # if EMG should be output
DebuggingOutput = F
DebuggingOnlyRunShortPartOfSimulation = F    # abort simulation after first stimulation
ElasticityDisabled = F          # disable the elasticity control loop

# ------------- numerics -------------------
# time steps
TimeStop = 1.0                 # total simulated time, note that the real simulated time is always a multiple of StimPeriod
ODETimeStep = 0.0001            # timestep size for 0D problem
PDETimeStep = 0.0005            # timestep size for 1D problem
ElasticityTimeStep = 0.1        # timestep size for 3D problem
SplittingType = 1               # 0 = Godunov splitting, 1 = Strang splitting

# solvers
ODESolverId = 5                 # 0D problem, ODE solver type: 1=explicit Euler, 2=BDF, 5=Improved Euler (Heun's Method)
MonodomainSolverId = 2          # 1D problem, solver
MonodomainPreconditionerId = 1  # 1D problem, preconditioner

# ODESolverId
# 1 explicit Euler (default)
# 2 BDF
# 3 General Linear (GL)
# 4 Crank-Nicolson
# 5 Improved Euler (Heun's Method)

# MonodomainSolverId
# 1 SOLVER_DIRECT_LU
# 2 SOLVER_ITERATIVE_GMRES (default)
# 3 SOLVER_ITERATIVE_CONJUGATE_GRADIENT
# 4 SOLVER_ITERATIVE_CONJGRAD_SQUARED

# MonodomainPreconditionerId
# 1 NO_PRECONDITIONER (default)
# 2 JACOBI_PRECONDITIONER
# 3 BLOCK_JACOBI_PRECONDITIONER
# 4 SOR_PRECONDITIONER
# 5 INCOMPLETE_CHOLESKY_PRECONDITIONER
# 6 INCOMPLETE_LU_PRECONDITIONER
# 7 ADDITIVE_SCHWARZ_PRECONDITIONER

NewtonMaximumNumberOfIterations = 500         # 3D solver maximum number of iterations
ElasticityLoopMaximumNumberOfIterations = 5   # number of load increments in pre-stretch simulation
NewtonTolerance = 1.E-8                       # abs. and rel. tolerance of 3D problem newton solver

# ------------- physical parameters ------------------

# stimulation
StimValue = 1200.0                   # current value with which fibres are stimulated, this value will be set on all stimulated nodes
PhysicalStimulationLength = 0.03125  # length of neuromuscular junction, length of line segment where stimulus is applied [cm] (set to 0 to always use 1 node)

# stimulation is applied periodic with periodic time PERIODD, length of stimulation burst is stim_stop
StimDuration = 0.1                   # the duration for which the stimulation is present, former STIM_STOP, a reasonable value is the same as ElastcityTimeStep
StimPeriod = 1.0                    # the duration between subsequent stimulation spikes, former PERIODD

# physical dimension [cm]
PhysicalLength = 1.0
PhysicalWidth  = 1.0
PhysicalHeight = 1.0

OldTomoMechanics = T            # whether to use the old mechanical description of Thomas Heidlauf that works also in parallel (deprecated, T=ModelType 0, F=ModelType 1)
ModelType = 0                   # which physical model to use

# ModelType
# 0 "MultiPhysStrain", no. 3a, OldTomoMechanics, old model of Thomas Heidlauf, that worked in parallel from the beginning
# 1 "MultiPhysStrain", no. 3, Model of Thomas Heidlauf, his comment in a e-mail of 25 Oct 2016: "multiscale and multiphysics active strain"
# 2 "Titin", no. 4, with additional Titin term in stress tensor

PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
CmFast = 0.58           # membrane capacitance [uF/cm^2]
CmSlow = 0.58
Vmax = -0.02            # maximum contraction velocity [m/s]
InitialStretch = 1.0   ! 1.2, 1.6   # Pre-Stretch factor, if set to 1.0, pre-stretch is disabled
TkLinParam = 1.0                # parameter for Titin model (ModelType=2) 0: No Actin-Titin Interactions, 1: With Actin-Tintin Interaction

# -------------- input files --------------------------
InputDirectory = input                        # directory where all input files are found, relative to working directory
FiringTimesFile = MU_firing_times_10s.txt     # file that contains the time steps (rows) for each fibre (columns) when the fibre is stimulated (value 1 else 0)
InnervationZoneFile = innervation_zone_18.txt # the position of the neuromuscular junction for each fibre
FibreDistributionFile = MU_fibre_distribution_4050.txt    # the motor unit number distribution for the fibres
CellMLModelFilename = standard                # the CELLML model file (XML), "standard" gets replaced by the actual standard file, depending on ModelType

# CellMLModelFilename
# ModelType = 0 (MultiPhysStrain): {InputDirectory}/slow_TK_2014_12_08.xml
# ModelType = 1 (MultiPhysStrain): {InputDirectory}/Aliev_Panfilov_Razumova_2016_08_22.cellml
# ModelType = 2 (Titin): {InputDirectory}/Aliev_Panfilov_Razumova_Titin_2016_10_10_noFv.cellml
