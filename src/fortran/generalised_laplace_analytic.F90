PROGRAM GeneralisedLaplaceAnalytic

  USE OpenCMISS
  USE OpenCMISS_Iron
  
  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Program parameters
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
 
  REAL(CMISSRP), PARAMETER :: SIGMA_11=3.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: SIGMA_22=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: SIGMA_12=0.0_CMISSRP  
  
  REAL(CMISSRP), PARAMETER :: PI=3.141592653589793_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: MAX_NUMBER_OF_REFINEMENTS=5
  INTEGER(CMISSIntg), PARAMETER :: MAX_NUMBER_OF_INTERPOLATIONS=9
  
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=6
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=7
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=8
  INTEGER(CMISSIntg), PARAMETER :: FIBRE_FIELD_USER_NUMBER=9
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=12
  INTEGER(CMISSIntg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=13
  INTEGER(CMISSIntg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=14
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=15

  !Program types

  !Program variables
  INTEGER(CMISSIntg) :: argumentLength,numberOfArguments,status
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex
  INTEGER(CMISSIntg) :: baseNumberOfGlobalXElements,baseNumberOfGlobalYElements,fibreAngleDegrees,interpolationType, &
    & numberOfGaussXi,numberOfGlobalXElements,numberOfGlobalYElements
  INTEGER(CMISSIntg) :: err
  REAL(CMISSRP) :: fibreAngleRadians
  INTEGER(CMISSIntg) :: convergenceDataIdx,numberOfDOFs,numberOfInterpolations,numberOfRefinements,refinementIdx
  REAL(CMISSRP) :: integralValues(2),localValues(8),localGhostValues(8),ghostIntegralValues(2), globalValues(8)
  REAL(CMISSRP) :: convergenceData(4,MAX_NUMBER_OF_REFINEMENTS,MAX_NUMBER_OF_INTERPOLATIONS), &
    & logConvergenceData(4,MAX_NUMBER_OF_REFINEMENTS,MAX_NUMBER_OF_INTERPOLATIONS), &
    & meanLogConvergenceData(4,MAX_NUMBER_OF_INTERPOLATIONS),slopeConvergenceData(4,MAX_NUMBER_OF_INTERPOLATIONS), &
    & sumLogConvergenceData(6,MAX_NUMBER_OF_INTERPOLATIONS)
  CHARACTER(LEN=255) :: commandArgument,filename

  !CMISS variables
  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: analyticField,dependentField,geometricField,equationsSetField,fibreField,materialsField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_NodesType) :: nodes
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 4) THEN
    !If we have enough arguments then use the first three for setting up the problem.
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) baseNumberOfGlobalXElements
    IF(numberOfGlobalXElements<=0) CALL HandleError("Invalid number of base X elements.")
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) baseNumberOfGlobalYElements
    IF(numberOfGlobalYElements<=0) CALL HandleError("Invalid number of base Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) fibreAngleDegrees
    IF(fibreAngleDegrees<=-180_CMISSIntg.Or.fibreAngleDegrees>180_CMISSIntg) CALL HandleError("Invalid fibre angle in degrees.")
    CALL GET_COMMAND_ARGUMENT(4,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 4.")
    READ(commandArgument(1:argumentLength),*) numberOfRefinements
    IF(numberOfRefinements<=0) CALL HandleError("Invalid number of refinements.")
    IF(numberOfRefinements>MAX_NUMBER_OF_REFINEMENTS) CALL HandleError("Increase the maximum number of refinements.")
  ELSE
    !If there are not enough arguments default the problem specification
    baseNumberOfGlobalXElements=4
    baseNumberOfGlobalYElements=2
    fibreAngleDegrees=30_CMISSIntg
    numberOfRefinements=4
  ENDIF
  fibreAngleRadians=REAL(fibreAngleDegrees,CMISSRP)*2.0_CMISSRP*PI/360.0_CMISSRP

  numberOfInterpolations=3

  convergenceData=0.0_CMISSRP
  logConvergenceData=0.0_CMISSRP
  meanLogConvergenceData=0.0_CMISSRP
  slopeConvergenceData=0.0_CMISSRP
  sumLogConvergenceData=0.0_CMISSRP

  !Initialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  
  DO interpolationType=1,numberOfInterpolations
    
    DO refinementIdx=1,numberOfRefinements

      numberOfGlobalXElements=refinementIdx*baseNumberOfGlobalXElements
      numberOfGlobalYElements=refinementIdx*baseNumberOfGlobalYElements

      convergenceData(1,refinementIdx,interpolationType)=SQRT((WIDTH*HEIGHT)/(REAL(numberOfGlobalXElements,CMISSRP)* &
        & REAL(numberOfGlobalYElements,CMISSRP)))

      !Create the context
      CALL cmfe_Context_Initialise(context,err)
      CALL cmfe_Context_Create(refinementIdx,context,err)
      CALL cmfe_Region_Initialise(worldRegion,err)
      CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
      CALL cmfe_Context_RandomSeedsSet(context,9999,err)
      !CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["Laplace_FiniteElementCalculate"],err)

      WRITE(filename,'(A,"_",I0,"x",I0,"_",I0,"_",I0)') "GeneralisedLaplace",numberOfGlobalXElements,numberOfGlobalYElements, &
        & interpolationType,fibreAngleDegrees

      CALL cmfe_OutputSetOn(filename,err)

      !Get the computational nodes information
      CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
      CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)

      CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
      CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
      CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
      CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

      !-----------------------------------------------------------------------------------------------------------
      ! COORDINATE SYSTEM
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of a new RC coordinate system
      CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
      CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
      !Set the coordinate system to be 2D
      CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,2,err)
      !Finish the creation of the coordinate system
      CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! REGION
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the region
      CALL cmfe_Region_Initialise(region,err)
      CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
      !Set the regions coordinate system to the 2D RC coordinate system that we have created
      CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
      !Set the region label
      CALL cmfe_Region_LabelSet(region,"GeneralisedLaplaceEquation",err)
      !Finish the creation of the region
      CALL cmfe_Region_CreateFinish(region,err)

      !-----------------------------------------------------------------------------------------------------------
      ! BASIS
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of a basis (default is trilinear lagrange)
      CALL cmfe_Basis_Initialise(basis,err)
      CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
      SELECT CASE(interpolationType)
      CASE(1,2,3,4)
        CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
      CASE(7,8,9)
        CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_SIMPLEX_TYPE,err)
      CASE DEFAULT
        CALL HandleError("Invalid interpolation type.")
      END SELECT
      SELECT CASE(interpolationType)
      CASE(1)
        numberOfGaussXi=2
      CASE(2)
        numberOfGaussXi=3
      CASE(3,4)
        numberOfGaussXi=4
      CASE DEFAULT
        numberOfGaussXi=0 !Don't set number of Gauss points for tri/tet
      END SELECT
      !Set the basis to be a bi-interpolation basis
      CALL cmfe_Basis_NumberOfXiSet(basis,2,err)
      CALL cmfe_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType],err)
      IF(numberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
      ENDIF
      !Finish the creation of the basis
      CALL cmfe_Basis_CreateFinish(basis,err)

      !-----------------------------------------------------------------------------------------------------------
      ! MESH
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of a generated mesh in the region
      CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
      CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
      !Set up a regular x*y*z mesh
      CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
      !Set the default basis
      CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)
      !Define the mesh on the region
      CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
      !Finish the creation of a generated mesh in the region
      CALL cmfe_Mesh_Initialise(mesh,err)
      CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

      !Create a decomposition
      CALL cmfe_Decomposition_Initialise(decomposition,err)
      CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
      !Finish the decomposition
      CALL cmfe_Decomposition_CreateFinish(decomposition,err)

      !Decompose
      CALL cmfe_Decomposer_Initialise(decomposer,err)
      CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
      !Add in the decomposition
      CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
      !Finish the decomposer
      CALL cmfe_Decomposer_CreateFinish(decomposer,err)

      !Destory the mesh now that we have decomposed it
      !CALL cmfe_Mesh_Destroy(mesh,err)

      !-----------------------------------------------------------------------------------------------------------
      ! GEOMETRIC FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Start to create a default (geometric) field on the region
      CALL cmfe_Field_Initialise(geometricField,err)
      CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
      !Set the decomposition to use
      CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
      !Set the domain to be used by the field components.
      CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
      CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
      !Finish creating the field
      CALL cmfe_Field_CreateFinish(geometricField,err)

      !Update the geometric field parameters
      CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

      !-----------------------------------------------------------------------------------------------------------
      ! FIBRE FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Start to create a fibre field on the region
      CALL cmfe_Field_Initialise(fibreField,err)
      CALL cmfe_Field_CreateStart(FIBRE_FIELD_USER_NUMBER,region,fibreField,err)
      !Set the field type
      CALL cmfe_Field_TypeSet(fibreField,CMFE_FIELD_FIBRE_TYPE,err)
      !Set the decomposition to use
      CALL cmfe_Field_DecompositionSet(fibreField,decomposition,err)
      !Set the geometric field
      CALL cmfe_Field_GeometricFieldSet(fibreField,geometricField,err)
      !Set the domain to be used by the field components.
      CALL cmfe_Field_ComponentMeshComponentSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
      !Set the interpolation type to be constant
      CALL cmfe_Field_ComponentInterpolationSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,err)
      !Finish creating the field
      CALL cmfe_Field_CreateFinish(fibreField,err)

      !Initialise the field parameters to the fibre angle (in radians)
      CALL cmfe_Field_ComponentValuesInitialise(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
        & fibreAngleRadians,err)

      !-----------------------------------------------------------------------------------------------------------
      ! EQUATIONS SETS
      !-----------------------------------------------------------------------------------------------------------

      !Create the Generalised Laplace equations set
      CALL cmfe_EquationsSet_Initialise(equationsSet,err)
      CALL cmfe_Field_Initialise(equationsSetField,err)
      CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
        & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE], &
        & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
      !Finish creating the equations set
      CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

      !-----------------------------------------------------------------------------------------------------------
      ! DEPENDENT FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Create the equations set dependent field variables
      CALL cmfe_Field_Initialise(dependentField,err)
      CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
      !Finish the equations set dependent field variables
      CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

      !Initialise the field with an initial guess
      CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
        & 0.0_CMISSRP, err)

      !Get the number of DOFs
      CALL cmfe_Field_NumberOfGlobalDOFsGet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfDOFs,err)
      convergenceData(2,refinementIdx,interpolationType)=REAL(numberOfDOFs,CMISSRP)

      !-----------------------------------------------------------------------------------------------------------
      ! MATERIALS FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Create the equations set materials field variables
      CALL cmfe_Field_Initialise(materialsField,err)
      CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
      !Set the interpolation type to be constant
      CALL cmfe_Field_ComponentInterpolationSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,err)
      CALL cmfe_Field_ComponentInterpolationSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,err)
      CALL cmfe_Field_ComponentInterpolationSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,err)
      !Finish the equations set materials field variables
      CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)

      !Initialise the materials conductivity tensor
      CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,SIGMA_11,err)
      CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,SIGMA_22,err)
      CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,SIGMA_12,err)

      !-----------------------------------------------------------------------------------------------------------
      ! ANALYTIC FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Create the equations set analytic field variables
      CALL cmfe_Field_Initialise(analyticField,err)
      CALL cmfe_EquationsSet_AnalyticCreateStart(equationsSet,CMFE_EQUATIONS_SET_GENERALISED_LAPLACE_EQUATION_TWO_DIM_1, &
        & ANALYTIC_FIELD_USER_NUMBER,analyticField,Err)
      !Finish the equations set analytic field variables
      CALL cmfe_EquationsSet_AnalyticCreateFinish(equationsSet,err)

      !Don't set the analytic field values, they will default to the materials and fibre field values

      !-----------------------------------------------------------------------------------------------------------
      ! EQUATIONS
      !-----------------------------------------------------------------------------------------------------------

      !Create the equations set equations
      CALL cmfe_Equations_Initialise(equations,err)
      CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
      !Set the equations matrices sparsity type
      CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
      !CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_FULL_MATRICES,err)
      !Set the equations set output
      CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
      !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_TIMING_OUTPUT,err)
      !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
      !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
      !Finish the equations set equations
      CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

      !-----------------------------------------------------------------------------------------------------------
      ! PROBLEM
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of a problem.
      CALL cmfe_Problem_Initialise(problem,err)
      CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
        & CMFE_PROBLEM_LAPLACE_EQUATION_TYPE,CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE],problem,err)
      !Finish the creation of a problem.
      CALL cmfe_Problem_CreateFinish(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! CONTROL LOOP
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the problem control loop
      CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
      !Finish creating the problem control loop
      CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! SOLVER
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the problem solvers
      CALL cmfe_Solver_Initialise(solver,err)
      CALL cmfe_Problem_SolversCreateStart(problem,err)
      CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
      CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT,err)
      !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
      !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
      !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
      !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)

      CALL cmfe_Solver_LinearTypeSet(solver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
      CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(solver,1.0E-12_CMISSRP,err)
      CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(solver,1.0E-12_CMISSRP,err)

      !CALL cmfe_Solver_LinearTypeSet(solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
      !CALL cmfe_Solver_LinearTypeSet(solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
      !CALL cmfe_Solver_LibraryTypeSet(solver,CMFE_SOLVER_MUMPS_LIBRARY,err)
      !CALL cmfe_Solver_LibraryTypeSet(solver,CMFE_SOLVER_LAPACK_LIBRARY,err)
      !CALL cmfe_Solver_LibraryTypeSet(solver,CMFE_SOLVER_SUPERLU_LIBRARY,err)
      !CALL cmfe_Solver_LibraryTypeSet(solver,CMFE_SOLVER_PASTIX_LIBRARY,err)
      !Finish the creation of the problem solver
      CALL cmfe_Problem_SolversCreateFinish(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! SOLVER EQUATIONS
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the problem solver equations
      CALL cmfe_Solver_Initialise(solver,err)
      CALL cmfe_SolverEquations_Initialise(solverEquations,err)
      CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
      !Get the solve equations
      CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
      CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
      !Set the solver equations sparsity
      CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
      !CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_FULL_MATRICES,err)
      !Add in the equations set
      CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
      !Finish the creation of the problem solver equations
      CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! BOUNDARY CONDITIONS
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the equations set boundary conditions
      CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
      CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
      !Set the analytic boundary conditions
      CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
      !Finish the creation of the equations set boundary conditions
      CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

      !-----------------------------------------------------------------------------------------------------------
      ! SOLVE
      !-----------------------------------------------------------------------------------------------------------

      !Solve the problem
      CALL cmfe_Problem_Solve(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! OUTPUT
      !-----------------------------------------------------------------------------------------------------------

      !Perform Analytic analysis
      CALL cmfe_AnalyticAnalysis_Output(dependentField,filename,err)
      CALL cmfe_AnalyticAnalysis_RMSErrorGetNode(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_ANALYTIC_ABSOLUTE_ERROR_TYPE, &
        & localValues,localGhostValues,globalValues,err)
      convergenceData(3,refinementIdx,interpolationType)=globalValues(1)
      CALL cmfe_AnalyticAnalysis_IntegralNumericalValueGet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,integralValues,&
        & ghostIntegralValues,err)
      convergenceData(3,refinementIdx,interpolationType)=integralValues(1)

      !Export results
      CALL cmfe_Fields_Initialise(fields,err)
      CALL cmfe_Fields_Create(region,fields,err)
      CALL cmfe_Fields_NodesExport(fields,filename,"FORTRAN",err)
      CALL cmfe_Fields_ElementsExport(fields,filename,"FORTRAN",err)
      CALL cmfe_Fields_Finalise(fields,err)

      !Turn ouput off
      CALL cmfe_OutputSetOff(err)
      
      !Destroy the context
      CALL cmfe_Context_Destroy(context,err)
      
      DO convergenceDataIdx=1,4
        logConvergenceData(convergenceDataIdx,refinementIdx,interpolationType)= &
          & LOG10(convergenceData(convergenceDataIdx,refinementIdx,interpolationType))
        meanLogConvergenceData(convergenceDataIdx,refinementIdx)=meanLogConvergenceData(convergenceDataIdx,refinementIdx)+ &
          & logConvergenceData(convergenceDataIdx,refinementIdx,interpolationType)
      ENDDO !convergenceDataIdx
      
    ENDDO !refinementIdx

    !Calculate the slope of the convergence line
    meanLogConvergenceData(:,interpolationType)=meanLogConvergenceData(:,interpolationType)/REAL(numberOfRefinements,CMISSRP)
    DO refinementIdx=1,numberOfRefinements
      sumLogConvergenceData(1,interpolationType)=sumLogConvergenceData(1,interpolationType)+ &
        & (logConvergenceData(1,refinementIdx,interpolationType)-meanLogConvergenceData(1,interpolationType))**2
      sumLogConvergenceData(2,interpolationType)=sumLogConvergenceData(2,interpolationType)+ &
        & (logConvergenceData(2,refinementIdx,interpolationType)-meanLogConvergenceData(2,interpolationType))**2
      sumLogConvergenceData(3,interpolationType)=sumLogConvergenceData(3,interpolationType)+ &
        & (logConvergenceData(1,refinementIdx,interpolationType)-meanLogConvergenceData(1,interpolationType))* &
        & (logConvergenceData(3,refinementIdx,interpolationType)-meanLogConvergenceData(3,interpolationType))
      sumLogConvergenceData(4,interpolationType)=sumLogConvergenceData(4,interpolationType)+ &
        & (logConvergenceData(1,refinementIdx,interpolationType)-meanLogConvergenceData(1,interpolationType))* &
        & (logConvergenceData(4,refinementIdx,interpolationType)-meanLogConvergenceData(4,interpolationType))
      sumLogConvergenceData(5,interpolationType)=sumLogConvergenceData(5,interpolationType)+ &
        & (logConvergenceData(2,refinementIdx,interpolationType)-meanLogConvergenceData(2,interpolationType))* &
        & (logConvergenceData(3,refinementIdx,interpolationType)-meanLogConvergenceData(3,interpolationType))
      sumLogConvergenceData(6,interpolationType)=sumLogConvergenceData(6,interpolationType)+ &
        & (logConvergenceData(2,refinementIdx,interpolationType)-meanLogConvergenceData(2,interpolationType))* &
        & (logConvergenceData(4,refinementIdx,interpolationType)-meanLogConvergenceData(4,interpolationType))
    ENDDO ! refinementIdx
    slopeConvergenceData(1,interpolationType)=sumLogConvergenceData(3,interpolationType)/sumLogConvergenceData(1,interpolationType)
    slopeConvergenceData(2,interpolationType)=sumLogConvergenceData(4,interpolationType)/sumLogConvergenceData(1,interpolationType)
    slopeConvergenceData(3,interpolationType)=sumLogConvergenceData(5,interpolationType)/sumLogConvergenceData(2,interpolationType)
    slopeConvergenceData(4,interpolationType)=sumLogConvergenceData(6,interpolationType)/sumLogConvergenceData(2,interpolationType)

    WRITE(*,'(X)')
    WRITE(*,'("Interpolation type : ",I0)') interpolationType
    WRITE(*,'("  RMS error vs h slope            = ",E12.4)') slopeConvergenceData(1,interpolationType)
    WRITE(*,'("  Integral^2 error vs h slope     = ",E12.4)') slopeConvergenceData(2,interpolationType)
    WRITE(*,'("  RMS error vs #DOFs slope        = ",E12.4)') slopeConvergenceData(3,interpolationType)
    WRITE(*,'("  Integral^2 error vs #DOFs slope = ",E12.4)') slopeConvergenceData(4,interpolationType)

  ENDDO !interpolationType
  
  !Finialise OpenCMISS
  CALL cmfe_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."
  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError

END PROGRAM GeneralisedLaplaceAnalytic
