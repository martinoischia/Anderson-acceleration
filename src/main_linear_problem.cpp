//! @file main_linear_problem.cpp
//! @brief Applying the accelerators to a large linear problem 
//!
//! Different matrixes from Matrix Market can be chosen through the input
//! data file, but in many cases the solving will fail.
//! The preconditioner (template argument to the Eigen iterative solver)
//! can also be changed (you need to recompile the code though).

#include "FixedPointIterator.hpp"
#include "GetPot"
#include "MM_readers3.3.hpp"
#include "simpleRichardson.hpp"
#include "Timing.hpp"
#include <unsupported/Eigen/IterativeSolvers>

int main(int argc, char** argv)
{
	std::cout << "\n********************EXAMPLE 2****************************\n\n\n\n";
	using namespace FixedPoint;
	GetPot get_filename (argc, argv);
	std::string filename = get_filename ("filename", "data.input");
	
	GetPot get_problem_data (filename.c_str ());
	
	using Vector = Traits::Vector; // Traits is a macro defined in Accelerators.hpp, if not defined in other ways
	using Matrix = Traits::Matrix;
	
	// reading the matrix and the RHS
	std::string matrix_filename = get_problem_data ( "files/MMMfile", "print me an error please" );
	std::string vector_filename = get_problem_data ( "files/MMVfile", "print me an error please" );
	std::cout << "***READING MATRIX\n";
	Timing t;
	Matrix MMM = Eigen::read_MM_Matrix <Matrix> ( matrix_filename );
	Vector MMV = Eigen::read_MM_Vector <Vector> ( vector_filename );
	t.write( "reading");
	std::size_t dimension = MMV.size() ;
	std::cout << "\nProblem dimension:" << dimension << std::endl;
	std::cout << "\nMatrix Frobenius norm: " << MMM.norm() << std::endl;
	std::cout << "\n***COMPUTING PRECONDITIONER" <<std::endl;
	// Eigen::BiCGSTAB < Matrix , Eigen::IncompleteLUT <double >> eigenSolver(MMM);
	Eigen::BiCGSTAB < Matrix > eigenSolver(MMM);
	eigenSolver.compute(MMM);
	t.write( "computing");
	
	#ifdef EigenSolver
		//Solving
		//Eigen implementation
		std::cout<<"\n*** WITH EIGEN ITERATIVE SOLVER"<<std::endl;
		if(eigenSolver.info()!=Eigen::Success) {
			// decomposition failed
			return 1;
		}
		auto & x = eigenSolver.solve(MMV);
		if(eigenSolver.info()!=Eigen::Success) {
			std::cout << "failed..\n";// solving failed
		}
		else{
			t.write( "solving");
			auto residuo = MMM*x-MMV;
			auto iteraz = eigenSolver.iterations() ;
			std::cout <<"\niterations (seems bugged to me)\n" << iteraz ;
			std::cout <<"\nnorm of RHS\n" << MMV.norm() ;
			std::cout<< "norm of residuum\n" << residuo.norm() << std::endl;
		}
	#endif
	
	//With the class FixedPointIterator
	std::cout<<"\n*** WITH RICHARDSON METHOD, USING EIGEN PRECONDITIONER:\n"; 
	//The iterator object
	FixedPointIterator FPI_2;
	double relaxationParameter = get_problem_data ("FPI_2_param/relaxationParameter", 1.);
	FPI_2.setIterator( std::make_unique <Iterator>(SimpleRichardson < decltype(eigenSolver) >(MMM, eigenSolver.preconditioner(), MMV, relaxationParameter), dimension)); 
	FPI_2.getOptions().maxIter = get_problem_data ("FPI_2_param/maxIter", 10000);
	FPI_2.getOptions().memory = get_problem_data ("FPI_2_param/memory", 8);
	FPI_2.getOptions().tolerance = get_problem_data ("FPI_2_param/tolerance", 1e-7);
	Timing t2;
	FPI_2.compute();
	t2.write ("Richardson");
	FPI_2.printResidualHistory();
	FPI_2.printResult();
	
	std::cout<<"\n*** APPLYING SECANT ACCELERATION:\n\n";
	FPI_2.reset();
	FPI_2.setIterator( std::make_unique <ASecantAccelerator>(SimpleRichardson < decltype(eigenSolver) >(MMM, eigenSolver.preconditioner(), MMV, relaxationParameter), dimension)); 
	Timing t3;
	FPI_2.compute();
	t3.write("secant");
	FPI_2.printResidualHistory();
	FPI_2.printResult();
	
	std::cout<<"\n*** APPLYING ANDERSON ACCELERATION:\n\n";
	FPI_2.reset();
	double mixingParameter = get_problem_data ("FPI_2_param/mixingParameter", 1.);
	std::size_t memory = get_problem_data ("FPI_2_param/AndersonMemory", 5);
	FPI_2.setIterator( std::make_unique <AndersonAccelerator>(SimpleRichardson < decltype(eigenSolver) >(MMM, eigenSolver.preconditioner(), MMV, relaxationParameter), dimension, mixingParameter, memory)); 
	Timing t4;
	FPI_2.compute();
	t4.write( "Anderson");
	FPI_2.printResidualHistory();
	FPI_2.printResult();
	FPI_2.reset();


	// Tentative Eigen unsupported GMRES usage, but doesn't work exactly like it should

	std::cout<<"\n*** APPLYING GMRES (is it the same if Anderson memory is infinity?):\n\n";

	// Eigen::GMRES< Matrix, Eigen::IncompleteLUT <double > > GMRESSolver( relaxationParameter*MMM);
	Eigen::GMRES< Matrix > GMRESSolver( relaxationParameter*MMM);
	std::cout << "tolerance " <<GMRESSolver.tolerance() << std::endl;
	Timing t5;
	GMRESSolver.set_restart(750);
	std::size_t GMRESiter = get_problem_data ("FPI_2_param/GMRESiter", 5);	
	GMRESSolver.setMaxIterations(GMRESiter);
	Vector sol = GMRESSolver.solve(relaxationParameter*MMV);
	t4.write( "GMRES");
		std::cout << "info " <<GMRESSolver.info() << std::endl;

	Vector residuum = MMM*sol-MMV;
	// for (int i = 0; i < dimension; ++i) std::cout << residuum(i) << "\n";
	std::cout << "tentative error     " << GMRESSolver.error() << std::endl;
	std::cout << "\n#iterations:     " << GMRESSolver.iterations() << std::endl;
	std::cout<< "norm of residuum\n" << residuum.norm() << std::endl;
	
	
}				