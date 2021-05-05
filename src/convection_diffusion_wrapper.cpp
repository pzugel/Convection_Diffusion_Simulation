#include "run_log.h"
#include "read_vtk.h"
#include "../convection_diffusion_base.h"
#include "../fe/convection_diffusion_fe.h"
#include "../fv/convection_diffusion_fv.h"
#include "../fv1/convection_diffusion_fv1.h"

#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "registry/class_name_provider.h"
#include "bridge/bridge.h"
#include "bridge/domain_bridges/refinement_bridge.cpp"
#include "lib_algebra/lib_algebra.h"


#include "lib_algebra/operator/linear_solver/lu.h"
#include "lib_algebra/operator/linear_solver/cg.h"
#include "lib_algebra/operator/linear_solver/bicgstab.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"

#include "lib_algebra/operator/preconditioner/jacobi.h"
#include "lib_algebra/operator/preconditioner/ilu.h"
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"
#include "lib_disc/operator/linear_operator/multi_grid_solver/mg_solver.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

#include "lib_algebra/operator/convergence_check.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/io/vtkoutput.h"

#include "lib_disc/domain_util.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_grid/subset_handler.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"
#include "common/util/smart_pointer.h"

#include <cstring>
#include <string>
#include <time.h>
#include <fstream>

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

//------------------------------------------------------------------------------
//******************************************************************************
//	Helper Functions (not exported)
//******************************************************************************
//------------------------------------------------------------------------------

namespace ug{

static RunLog* newRunLog = new RunLog();
static VTKReader* newVTKReader = new VTKReader();
static std::string lastVTKOutput = "";

//Produces 2D Vector
template<int dim>
ConstUserVector<dim>* createConstVector(double arr[2][1]){
	ConstUserVector<dim>* constVector = new ConstUserVector<dim>;
	constVector->set_entry(0,arr[0][0]);
	constVector->set_entry(1,arr[1][0]);
	return constVector;
}

//Produces 2D Matrix
template<int dim>
ConstUserMatrix<dim>* createConstMatrix(double arr[2][2]){
	ConstUserMatrix<dim>* constMatrix = new ConstUserMatrix<dim>;
	constMatrix->set_entry(0,0,arr[0][0]);
	constMatrix->set_entry(0,1,arr[0][1]);
	constMatrix->set_entry(1,0,arr[1][0]);
	constMatrix->set_entry(1,1,arr[1][1]);
	return constMatrix;
}

//------------------------------------------------------------------------------
//******************************************************************************
//	Main Setup Funtion (not exported)
//******************************************************************************
//------------------------------------------------------------------------------

template<typename TDomain>
int setupConvectionDiffusion(
		const char* fileName,
		const char* discType, 
		const char* upwindType, 
		const char* subsets,
		bool setDiffusionValue,
		number diffusionValue,
		double diffusionMatrix[2][2],	
		number reactionRateValue,
		number sourceValue,
		number massScaleValue,
		double velocityVector[2][1],
		const char* precondName,
		const char* solverName,
		int numRefs,
		int maxNumIter,
		number absTol)
{
	std::stringbuf UG_LOG_buffer(std::ios::out);
    auto UG_LOG_out = std::cout.rdbuf(std::addressof(UG_LOG_buffer));

	UG_LOG("\tNew Run - ConvectionDiffusion\n");
	const int dim = TDomain::dim;
	UG_LOG("** Domain ");
	UG_LOG(std::to_string(dim));
	UG_LOG("D\n");
	const char* functions = "c";

	//*********************************
	//**  Physical Parameters
	//*********************************
	
	//---------------------------------
	//-- Upwind
	//---------------------------------

	UG_LOG("** Creating Upwind ... ");
	ug::IConvectionShapes<dim>* upwinding;	
	std::string upwindTypeString = upwindType;
	if(upwindTypeString.compare("Full")==0)
	{
		upwinding = new ConvectionShapesFullUpwind<dim>();
	}
	else if(upwindTypeString.compare("Partial")==0)
	{
		upwinding = new ConvectionShapesPartialUpwind<dim>();
	}
	else if(upwindTypeString.compare("None")==0)
	{
		upwinding = new ConvectionShapesNoUpwind<dim>();
	}

	UG_LOG("(");
	UG_LOG(upwindType);
	UG_LOG(")");
	UG_LOG(" ... Done!\n");

	//---------------------------------
	//-- DiscType
	//---------------------------------

	ConvectionDiffusionPlugin::ConvectionDiffusionBase<TDomain>* elemDisc;
	std::string discTypeString = discType;
	UG_LOG("** Creating Element Discretization ... ");
	if(discTypeString.compare("Finite Elements")==0)
	{
		elemDisc = new ConvectionDiffusionPlugin::ConvectionDiffusionFE<TDomain>(functions, subsets);
	}
	else if(discTypeString.compare("Finite Volumes")==0)
	{
		elemDisc = new ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain>(functions, subsets);
		((ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain>*) elemDisc)->set_upwind(make_sp(upwinding));
	}
	
	UG_LOG("(");
	UG_LOG(discType);
	UG_LOG(")");
	UG_LOG(" ... Done!\n");

	//---------------------------------
	//-- Set Values
	//---------------------------------

	UG_LOG("** Setting Values ... ");
	//Reaction Rate / Source / Mass
	elemDisc->set_reaction_rate(reactionRateValue);
	elemDisc->set_source(sourceValue);
	elemDisc->set_mass_scale(massScaleValue);

	//Velocity
	ConstUserVector<dim>* velVec = createConstVector<dim>(velocityVector);
	elemDisc->set_velocity(make_sp(velVec));

	//Diffusion
	if(setDiffusionValue)
	{
		UG_LOG("(Diffusion: Value entry)");
		elemDisc->set_diffusion(diffusionValue);
	}
	else
	{
		UG_LOG("(Diffusion: Matrix entry)");
		ConstUserMatrix<dim>* difMat = createConstMatrix<dim>(diffusionMatrix);
		elemDisc->set_diffusion(make_sp(difMat));	
	}
	UG_LOG(" ... Done!\n");

	//*********************************
	//**  Assemble
	//*********************************

	UG_LOG("** Assemble");
	bridge::InitUG(dim, AlgebraType("CPU", 1));
	
	auto  dom = new Domain<dim, MultiGrid, MGSubsetHandler>();
	LoadDomain<TDomain>(*dom, fileName);
	
	auto dom_sptr = SmartPtr<TDomain>(dom);
	ApproximationSpace<TDomain>* approxSpace = new ApproximationSpace<TDomain>(dom_sptr);
	approxSpace->add(functions, "Lagrange", 1);
	
	auto approxSpaceSmartPtr = SmartPtr<ApproximationSpace<TDomain>>(approxSpace);
	DomainDiscretization<TDomain, CPUAlgebra>* domainDisc = new DomainDiscretization<TDomain, CPUAlgebra>(approxSpaceSmartPtr);

	//MGSubsetHandler sh = dom.subset_handler(); Not needed right now
	
	auto elemDiscSmartPtr = SmartPtr<IElemDisc<TDomain>>(elemDisc);
		
	domainDisc->add(elemDiscSmartPtr);

	//---------------------------------
	//-- Boundary
	//---------------------------------
	
	auto dirichletBND = new DirichletBoundary<TDomain, CPUAlgebra>(); 
	dirichletBND->add(1.0, functions, "Boundary");
	//dirichletBND->add(1.0, functions, "Inner"); //DO NOT ADD BOTH TO THE SAME DOMAIN DISC!
	domainDisc->add((SmartPtr<IDomainConstraint<TDomain, CPUAlgebra>>)make_sp(dirichletBND));

	UG_LOG(" ... Done!\n");

	//*********************************
	//**  Solver Setup
	//*********************************

	//---------------------------------
	//-- Refiner
	//---------------------------------
	
	UG_LOG("** Setting up Refiner");
	SmartPtr<IRefiner> refiner = GlobalDomainRefiner<TDomain>(dom);

	for(int i = 0; i < numRefs; ++i){
		refiner.get()->refine();
	}
	UG_LOG(" ... Done!\n");
	
	//---------------------------------
	//-- Init
	//---------------------------------
	
	UG_LOG("** Init");
	std::string precondNameString = precondName;
	if(precondNameString.compare("Geometric MultiGrid")==0){
		approxSpace->init_levels();
	}
	approxSpace->init_top_surface();
	//approxSpace->print_statistic();
	UG_LOG(" ... Done!\n");
	
	//OrderCuthillMcKee.invoke(approxSpace, true);

	auto A = new MatrixOperator<CPUAlgebra::matrix_type, CPUAlgebra::vector_type>();
	auto u = new GridFunction<TDomain, CPUAlgebra>(approxSpaceSmartPtr, true); 
	auto b = new GridFunction<TDomain, CPUAlgebra>(approxSpaceSmartPtr, true);

	//---------------------------------
	//-- Convergence Check
	//---------------------------------	
	
	auto convCheck = StdConvCheck<CPUAlgebra::vector_type>();
	convCheck.set_maximum_steps(maxNumIter);
	convCheck.set_minimum_defect(absTol);
	convCheck.set_reduction(1e-20);

	//---------------------------------
	//-- Solver
	//---------------------------------	
	
	std::string solverNameString = solverName;
	UG_LOG("** Setting up Solver ...");
	IPreconditionedLinearOperatorInverse<CPUAlgebra::vector_type>* solver;
	if(solverNameString.compare("BiCGStab")==0) 
	{
		solver = new BiCGStab<CPUAlgebra::vector_type>();
	}
	else if(solverNameString.compare("Conjugate Gradient")==0) 
	{
		solver = new CG<CPUAlgebra::vector_type>();
	}    
	else if(solverNameString.compare("Linear Solver")==0) 
	{
		solver = new LinearSolver<CPUAlgebra::vector_type>();
	}    
	else if(solverNameString.compare("LU")==0) 
	{
		//This case is handeled below
	}

	if(!(solverNameString.compare("LU")==0))
	{
		solver->set_convergence_check(make_sp(&convCheck));
	}
	UG_LOG("(");
	UG_LOG(solverNameString);
	UG_LOG(")");
	UG_LOG(" ... Done!\n");

	//---------------------------------
	//-- Preconditioner
	//---------------------------------

	UG_LOG("** Setting up Preconditioner ... ");
	if(!(precondNameString.compare("None")==0) && !(solverNameString.compare("LU")==0))
	{
		ILinearIterator<CPUAlgebra::vector_type>* precond;

		if(precondNameString.compare("Jacobi")==0)
		{	
			precond = new Jacobi<CPUAlgebra>();
		}        
		else if(precondNameString.compare("Gauss-Seidel")==0)
		{		
			precond = new GaussSeidel<CPUAlgebra>();
		} 
		else if(precondNameString.compare("Incomplete LU")==0)
		{		
			precond = new ILU<CPUAlgebra>();
		} 
		else if(precondNameString.compare("Geometric MultiGrid")==0)
		{				
			precond = new AssembledMultiGridCycle<TDomain, CPUAlgebra>(approxSpaceSmartPtr);
			((AssembledMultiGridCycle<TDomain, CPUAlgebra>*) precond)->set_discretization(make_sp(domainDisc));
			((AssembledMultiGridCycle<TDomain, CPUAlgebra>*) precond)->set_base_level(0);
			((AssembledMultiGridCycle<TDomain, CPUAlgebra>*) precond)->set_base_solver(make_sp(new LU<CPUAlgebra>()));
			((AssembledMultiGridCycle<TDomain, CPUAlgebra>*) precond)->set_smoother(make_sp(new ILU<CPUAlgebra>()));
			((AssembledMultiGridCycle<TDomain, CPUAlgebra>*) precond)->set_cycle_type(1);
			((AssembledMultiGridCycle<TDomain, CPUAlgebra>*) precond)->set_num_presmooth(3);
			((AssembledMultiGridCycle<TDomain, CPUAlgebra>*) precond)->set_num_postsmooth(3);
		}

		solver->set_preconditioner(make_sp(precond));
	}
	UG_LOG("(");
	UG_LOG(precondNameString);
	UG_LOG(")");
	UG_LOG(" ... Done!\n");

	//*********************************
	//**  Apply Solver
	//*********************************
	
	if(solverNameString.compare("LU")==0) 
	{
		try{
			ILinearOperatorInverse<CPUAlgebra::vector_type>* solver_lucase = new LU<CPUAlgebra>();
			solver_lucase->set_convergence_check(make_sp(&convCheck));
			domainDisc->assemble_linear(*A, *b);
			domainDisc->adjust_solution(*u);
			UG_LOG("** Starting Linear solver ... ");

			if(solver_lucase->init(make_sp(A)) == false)
			{
				UG_LOG("Initialization of Linear Solver failed.\n\n");
				return 1;
			}
			
			if(solver_lucase->apply_return_defect(*u,*b) == false)
			{
				UG_LOG("Linear Solver did not converge.\n\n");
				return 2;
			}
		} catch (ug::UGError e) {
			UG_LOG("Instance of UGError!\n\n");
			return 3;
		}
	}
	else
	{
		try{
			domainDisc->assemble_linear(*A, *b);
			domainDisc->adjust_solution(*u);
			UG_LOG("** Starting Linear solver ... ");

			if(solver->init(make_sp(A)) == false)
			{
				UG_LOG("Initialization of Linear Solver failed.\n\n");
				return 1;
			}
			
			if(solver->apply_return_defect(*u,*b) == false)
			{
				UG_LOG("Linear Solver did not converge.\n\n");
				return 2;
			}
		} catch (ug::UGError e) {
			UG_LOG("Instance of UGError!\n\n");
			return 3;
		}
	}
	UG_LOG("Done!\n");

	//*********************************
	//**  VTK Output
	//*********************************

	UG_LOG("** Writing VTK ... \n");

	VTKOutput<dim> out;
	out.set_binary(false);
	time_t rawtime;
  	struct tm * timeinfo;
  	char buffer[80];

  	time (&rawtime);
  	timeinfo = localtime(&rawtime);

  	strftime(buffer,sizeof(buffer),"%d%m%Y%H%M%S",timeinfo);
  	std::string str(buffer);

	const char* fileOut = "/tmp/tmpVTU_";
	std::string fileOut_timestamp = ((std::string) fileOut) + str + ".vtu";
	out.print(fileOut_timestamp.c_str(), *u);
	lastVTKOutput = ((std::string) fileOut) + str + ".vtu";

	if(newVTKReader->init(fileOut_timestamp.c_str()))
	{
		//Everything worked
		/*
		newRunLog->writeLog("minXValue : " + std::to_string(newVTKReader->minXValue) + "\n");
		newRunLog->writeLog("minYValue : " + std::to_string(newVTKReader->minYValue) + "\n");
		newRunLog->writeLog("maxXValue : " + std::to_string(newVTKReader->maxXValue) + "\n");
		newRunLog->writeLog("maxYValue : " + std::to_string(newVTKReader->maxYValue) + "\n");
		newRunLog->writeLog("numXCoords : " + std::to_string(newVTKReader->numXCoords) + "\n");
		newRunLog->writeLog("numYCoords : " + std::to_string(newVTKReader->numYCoords) + "\n");
		newRunLog->writeLog("numXCoordsTotal : " + std::to_string(newVTKReader->numXCoordsTotal) + "\n");
		newRunLog->writeLog("numYCoordsTotal : " + std::to_string(newVTKReader->numYCoordsTotal) + "\n");
		newRunLog->writeLog("numValues : " + std::to_string(newVTKReader->numValues) + "\n");
		newRunLog->writeLog("numPoints : " + std::to_string(newVTKReader->numPoints) + "\n");
		newRunLog->writeLog("numCells : " + std::to_string(newVTKReader->numCells) + "\n");
		newRunLog->writeLog("xOffset : " + std::to_string(newVTKReader->xOffset) + "\n");
		newRunLog->writeLog("yOffset : " + std::to_string(newVTKReader->yOffset) + "\n");
		*/
		//remove(fileOut_timestamp.c_str());
	}
	else
	{
		UG_LOG("File does not exist!\n");
		return 4;
	}


	UG_LOG("EOF\n\n");
	std::cout.rdbuf( UG_LOG_out );
    std::string UG_LOG_output = UG_LOG_buffer.str(); 
    newRunLog->writeLog(UG_LOG_output);
	return 0;
}

} //end namespace ug

//------------------------------------------------------------------------------
//******************************************************************************
//	Exported Functions
//******************************************************************************
//------------------------------------------------------------------------------

namespace ug{
namespace ConvectionDiffusionPlugin{

extern "C"{

//*********************************
//	Init
//*********************************

int wrapperInit(
		int dim, 
		const char* fileName,
		const char* discType, 
		const char* upwindType, 
		const char* subsets,
		bool setDiffusionValue,
		number diffusionValue,
		double diffusionMatrix[2][2],	
		number reactionRateValue,
		number sourceValue,
		number massScaleValue,
		double velocityVector[2][1],
		const char* precondName,
		const char* solverName,
		int numRefs,
		int maxNumIter,
		number absTol)
{
		return ug::setupConvectionDiffusion<ug::Domain2d>(
			fileName, 
			discType, 
			upwindType, 
			subsets, 
			setDiffusionValue, 
			diffusionValue, 
			diffusionMatrix, 
			reactionRateValue, 
			sourceValue, 
			massScaleValue, 
			velocityVector,
			precondName,
			solverName,
			numRefs,
			maxNumIter,
			absTol);
};

//*********************************
//	Other
//*********************************

void saveVTKFile(const char* path){
	std::ifstream  src(lastVTKOutput.c_str(), std::ios::binary);
	std::ofstream  dst(path,   std::ios::binary);

	dst << src.rdbuf();
}

char* displayLog(){
	return newRunLog->displayLog();
}

void clearLog(){
	newRunLog->clearLog();
}

double getVTK_minXValue(){
	return newVTKReader->minXValue;
}

double getVTK_minYValue(){
	return newVTKReader->minYValue;
}

double getVTK_maxXValue(){
	return newVTKReader->maxXValue;
}

double getVTK_maxYValue(){
	return newVTKReader->maxYValue;
}

int getVTK_numXCoords(){
	return newVTKReader->numXCoords;
}

int getVTK_numYCoords(){
	return newVTKReader->numYCoords;
}

int getVTK_numPoints(){
	return newVTKReader->numPoints;
}

int getVTK_numCells(){
	return newVTKReader->numCells;
}

double getVTK_xOffset(){
	return newVTKReader->xOffset;
}

double getVTK_yOffset(){
	return newVTKReader->yOffset;
}

void getVTK_Data(double (&vtkdata)[100][100]){
	newVTKReader->get_VTKData(vtkdata);
}

}// end extern "C"

}//end namespace ConvectionDiffusionPlugin
}//end namespace ug



/*
gcc -DUG_ALGEBRA -DUG_BRIDGE -DUG_CPU_1 -DUG_DEBUG -DUG_DIM_1 -DUG_DIM_2 -DUG_DIM_3 -DUG_DISC -DUG_FOR_LUA -DUG_GRID -DUG_PLUGINS -DUG_POSIX -DUG_PROFILER_SHINY_CHECK_CONSISTENCY -DUG_TARGET=\"ugshell\" -I/home/paul/ug4/externals/BoostForUG4 -I/home/paul/ug4/ugcore/ugbase -L/home/paul/ug4/lib -L/home/paul/ug4/bin/plugins -Wl,-rpath=/home/paul/ug4/lib -Wl,-rpath=/home/paul/ug4/bin/plugins -shared -o libConvectionDiffusionWrapper.so -fPIC -w convection_diffusion_wrapper.cpp -lug4 -lConvectionDiffusion
*/


