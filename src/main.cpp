#include "/home/paul/ug4/plugins/ConvectionDiffusion/convection_diffusion_base.h"
#include "/home/paul/ug4/plugins/ConvectionDiffusion/fe/convection_diffusion_fe.h"
#include "/home/paul/ug4/plugins/ConvectionDiffusion/fv/convection_diffusion_fv.h"
#include "/home/paul/ug4/plugins/ConvectionDiffusion/fv1/convection_diffusion_fv1.h"

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

#include <iostream>
#include <cstring>
#include <string>
#include <time.h>	

namespace ug{

template<int dim>
ConstUserVector<dim>* createConstVector(double arr[3][1]){
	ConstUserVector<dim>* constVector = new ConstUserVector<dim>;
	constVector->set_entry(0,arr[0][0]);
	if(dim==2)
	{
		constVector->set_entry(1,arr[1][0]);
	}
	if(dim==3)
	{
		constVector->set_entry(2,arr[2][0]);
	}
	return constVector;
}

template<int dim>
ConstUserMatrix<dim>* createConstMatrix(double arr[3][3]){
	ConstUserMatrix<dim>* constMatrix = new ConstUserMatrix<dim>;
	constMatrix->set_entry(0,0,arr[0][0]);
	if(dim==2)
	{
		constMatrix->set_entry(0,1,arr[0][1]);
		constMatrix->set_entry(1,0,arr[1][0]);
		constMatrix->set_entry(1,1,arr[1][1]);
	}
	if(dim==3)
	{
		constMatrix->set_entry(0,2,arr[0][2]);
		constMatrix->set_entry(1,2,arr[1][2]);
		constMatrix->set_entry(2,0,arr[2][0]);
		constMatrix->set_entry(2,1,arr[2][1]);
		constMatrix->set_entry(2,2,arr[2][2]);
	}
	return constMatrix;
}
}//end namespace ug

namespace ug{
namespace ConvectionDiffusionPlugin{

extern "C" int main(){
	std::cout<<"New Run: ConvectionDiffusion\n";

	const int dim = 2;
	const char* fileName = "/home/paul/Schreibtisch/vorlesung-2018/laplace_sample_grid_2d.ugx";
	const char* discType = "Finite Elements";
	const char* upwindType = "Full";
	const char* subsets = "Inner";
	bool setDiffusionValue = true;
	number diffusionValue = 1.0;
	double diffusionMatrix[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};	
	number reactionRateValue = 0.0;
	number sourceValue = 1.0;
	number massScaleValue = 0.0;
	double velocityVector[3][1] = {{0.0},{0.0},{0.0}};	

	const char* precondName = "Jacobi";
	const char* solverName = "Linear Solver";
	int numRefs = 3;
	int maxNumIter = 100;
	double absTol = 0.00001;

	std::cout<<"Dimension: " + std::to_string(dim) + "\n\n";
	
	const char* functions = "c";

	//************************************************************************************************************************************
	//**  Physical Parameters
	//************************************************************************************************************************************
	
	//---------------------------------
	//-- Upwind
	//---------------------------------

	std::cout<<"Create Upwind ...";
	ug::IConvectionShapes<dim>* upwinding;	
	std::string upwindTypeString = upwindType;
	if(upwindTypeString.compare("Full")==0)
	{
		std::cout << "(Full)...";
		upwinding = new ConvectionShapesFullUpwind<dim>();
	}
	else if(upwindTypeString.compare("Partial")==0)
	{
		std::cout << "(Partial)...";
		upwinding = new ConvectionShapesPartialUpwind<dim>();
	}
	else if(upwindTypeString.compare("None")==0)
	{
		std::cout << "(None)...";
		upwinding = new ConvectionShapesNoUpwind<dim>();
	}
	std::cout<<"Done!\n";

	//---------------------------------
	//-- DiscType
	//---------------------------------

	std::cout << "Creating Element Disc ...";

	ConvectionDiffusionBase<Domain2d>* elemDisc;
	std::string discTypeString = discType;
	if(discTypeString.compare("Finite Elements")==0)
	{
		std::cout << "(Finite Elements)...";
		elemDisc = new ConvectionDiffusionFE<Domain2d>(functions, subsets);
	}
	else if(discTypeString.compare("Finite Volumes")==0)
	{
		std::cout << "(Finite Volumes, Upwind:";
		std::cout << upwindTypeString;
		std::cout << ")...";
		elemDisc = new ConvectionDiffusionFV1<Domain2d>(functions, subsets);
		((ConvectionDiffusionFV1<Domain2d>*) elemDisc)->set_upwind(make_sp(upwinding));
	}
	std::cout << "Done!\n";	

	//---------------------------------
	//-- Set Values
	//---------------------------------

	std::cout << "Setting Values ...";
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
		std::cout << "(Diffusion: Value entry)...";
		elemDisc->set_diffusion(diffusionValue);
	}
	else
	{
		std::cout << "(Diffusion: Matrix entry)...";
		ConstUserMatrix<dim>* difMat = createConstMatrix<dim>(diffusionMatrix);
		elemDisc->set_diffusion(make_sp(difMat));	
	}
	std::cout << "Done!\n";	

	//************************************************************************************************************************************
	//**  Assemble
	//************************************************************************************************************************************

	std::cout << "Assemble ...";
	bridge::InitUG(dim, AlgebraType("CPU", 1));
	
	auto  dom = new Domain<dim, MultiGrid, MGSubsetHandler>();
	LoadDomain<Domain2d>(*dom, fileName);
	
	auto dom_sptr = SmartPtr<Domain2d>(dom);
	ApproximationSpace<Domain2d>* approxSpace = new ApproximationSpace<Domain2d>(dom_sptr);
	approxSpace->add(functions, "Lagrange", 1);
	
	auto approxSpaceSmartPtr = SmartPtr<ApproximationSpace<Domain2d>>(approxSpace);
	DomainDiscretization<Domain2d, CPUAlgebra>* domainDisc = new DomainDiscretization<Domain2d, CPUAlgebra>(approxSpaceSmartPtr);

	//MGSubsetHandler sh = dom.subset_handler(); Not needed right now
	
	auto elemDiscSmartPtr = SmartPtr<IElemDisc<Domain2d>>(elemDisc);
		
	//domainDisc->add(elemDiscSmartPtr);

	std::cout << "Done!\n";	

	//************************************************************************************************************************************
	//**  Solver Setup
	//************************************************************************************************************************************

	//---------------------------------
	//-- Refiner
	//---------------------------------
	
	std::cout << "Setting up Refiner ...";

	SmartPtr<IRefiner> refiner = GlobalDomainRefiner<Domain2d>(dom);

	for(int i = 0; i < numRefs; ++i){
		refiner.get()->refine();
	}
	std::cout << "Done!\n";

	//---------------------------------
	//-- Init
	//---------------------------------
	
	std::cout << "Init ...\n";
	std::string precondNameString = precondName;
	if(precondNameString.compare("Geometric MultiGrid")==0){
		approxSpace->init_levels();
	}
	approxSpace->init_top_surface();
	approxSpace->print_statistic();
	std::cout << "Done!\n";
	
	
	//OrderCuthillMcKee.invoke(approxSpace, true);

	auto A = new MatrixOperator<CPUAlgebra::matrix_type, CPUAlgebra::vector_type>();
	auto u = new GridFunction<Domain2d, CPUAlgebra>(approxSpaceSmartPtr, true); 
	auto b = new GridFunction<Domain2d, CPUAlgebra>(approxSpaceSmartPtr, true); 

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
	std::cout << "Setting up Solver ...";
	IPreconditionedLinearOperatorInverse<CPUAlgebra::vector_type>* solver;

	if(solverNameString.compare("BiCGStab")==0) 
	{
		std::cout << "(BiCGStab)...";
		solver = new BiCGStab<CPUAlgebra::vector_type>();
	}
	else if(solverNameString.compare("Conjugate Gradient")==0) 
	{
		std::cout << "(Conjugate Gradient)...";
		solver = new CG<CPUAlgebra::vector_type>();
	}    
	else if(solverNameString.compare("Linear Solver")==0) 
	{
		std::cout << "(Linear Solver)...";
		solver = new LinearSolver<CPUAlgebra::vector_type>();
	}    
	else if(solverNameString.compare("LU")==0) 
	{
		std::cout << "(LU)...";
		//This case is handeled below
	}
	
	if(!(solverNameString.compare("LU")==0))
	{
		solver->set_convergence_check(make_sp(&convCheck));
	}
	std::cout << "Done!\n";

	//---------------------------------
	//-- Preconditioner
	//---------------------------------

	std::cout << "Setting up Preconditioner ...";
	if(!(precondNameString.compare("None")==0) && !(solverNameString.compare("LU")==0)) {
    	ILinearIterator<CPUAlgebra::vector_type>* precond;
 
        if(precondNameString.compare("Jacobi")==0)
		{	
			std::cout << "(Jacobi)...";
			precond = new Jacobi<CPUAlgebra>();
		}        
		else if(precondNameString.compare("Gauss-Seidel")==0)
		{		
			std::cout << "(Gauss-Seidel)...";
			precond = new GaussSeidel<CPUAlgebra>();
		} 
        else if(precondNameString.compare("Incomplete LU")==0)
		{		
			std::cout << "(Incomplete LU)...";
			precond = new ILU<CPUAlgebra>();
		} 
        else if(precondNameString.compare("Geometric MultiGrid")==0)
		{	
			std::cout << "(Geometric MultiGrid)...";
			
			precond = new AssembledMultiGridCycle<Domain2d, CPUAlgebra>(approxSpaceSmartPtr);
			((AssembledMultiGridCycle<Domain2d, CPUAlgebra>*) precond)->set_discretization(make_sp(domainDisc));
			((AssembledMultiGridCycle<Domain2d, CPUAlgebra>*) precond)->set_base_level(0);
			((AssembledMultiGridCycle<Domain2d, CPUAlgebra>*) precond)->set_base_solver(make_sp(new LU<CPUAlgebra>()));
			((AssembledMultiGridCycle<Domain2d, CPUAlgebra>*) precond)->set_smoother(make_sp(new ILU<CPUAlgebra>()));
			((AssembledMultiGridCycle<Domain2d, CPUAlgebra>*) precond)->set_cycle_type(1);
			((AssembledMultiGridCycle<Domain2d, CPUAlgebra>*) precond)->set_num_presmooth(3);
			((AssembledMultiGridCycle<Domain2d, CPUAlgebra>*) precond)->set_num_postsmooth(3);
        }

		solver->set_preconditioner(make_sp(precond));
	}
	std::cout << "Done!\n";

	//************************************************************************************************************************************
	//**  Apply Solver
	//************************************************************************************************************************************

	if(solverNameString.compare("LU")==0) 
	{
		try{
			ILinearOperatorInverse<CPUAlgebra::vector_type>* solver_lucase = new LU<CPUAlgebra>();
			solver_lucase->set_convergence_check(make_sp(&convCheck));
			domainDisc->assemble_linear(*A, *b);
			domainDisc->adjust_solution(*u);
			std::cout << "Starting Linear solver ...\n";

			if(solver_lucase->init(make_sp(A)) == false)
			{
				std::cout << "Initialization of Linear Solver failed.";
				return 1;
			}
			
			if(solver_lucase->apply_return_defect(*u,*b) == false)
			{
				std::cout << "Linear Solver did not converge.\n";
				return 2;
			}
		} catch (ug::UGError e) {
			std::cout << "Instance of UGError!";
			return 3;
		}
	}
	else
	{
		try{
			domainDisc->assemble_linear(*A, *b);
			domainDisc->adjust_solution(*u);
			std::cout << "Starting Linear solver ...\n";

			if(solver->init(make_sp(A)) == false)
			{
				std::cout << "Initialization of Linear Solver failed.";
				return 1;
			}
			
			if(solver->apply_return_defect(*u,*b) == false)
			{
				std::cout << "Linear Solver did not converge.\n";
				return 2;
			}
		} catch (ug::UGError e) {
			std::cout << "Instance of UGError!";
			return 3;
		}
	}
	//ug::GetLogAssistant().logger(); TODO 
	//TODO VTKOutput	
	/*
	const char* filename = "vtuASCII";

	VTKOutput<dim> out;
	out.set_binary(false);
	out.print(filename, *dom);
	std::cout << dom->domain_info().to_string();

	auto grid = dom->grid();
	VertexIterator verticesBegin = grid->vertices_begin();
	VertexIterator verticesEnd = grid->vertices_end();
	size_t vertices = grid->num_vertices();
	size_t edges = grid->num_edges();
	size_t faces = grid->num_faces();
	size_t volumes = grid->num_volumes();
	std::cout << "#Surface Elements: " << dom->domain_info().num_surface_elements() << "\n";
	std::cout << "#Vertices: " << vertices << "\n";
	std::cout << "#Edges: " << edges << "\n";
	std::cout << "#Faces: " << faces << "\n";
	std::cout << "#Volumes: " << volumes << "\n";
	*/
	std::cout << "... EOF\n";
	return 0;
}

} //end namespace ConvectionDiffusionPlugin
} //end namespace ug
