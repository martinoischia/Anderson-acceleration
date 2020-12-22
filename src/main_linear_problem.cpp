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
int main(int argc, char** argv)
{
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
	std::cout << "\n***COMPUTING PRECONDITIONER" <<std::endl;
	std::size_t dimension = MMV.size() ;
	Eigen::BiCGSTAB < Matrix , Eigen::IncompleteLUT <double >> eigenSolver(MMM);
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
	FPI_2.setIterator( std::make_unique <Iterator>(SimpleRichardson < decltype(eigenSolver) >(MMM, eigenSolver.preconditioner(), MMV, 0.2), dimension)); 
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
	FPI_2.setIterator( std::make_unique <ASecantAccelerator>(SimpleRichardson < decltype(eigenSolver) >(MMM, eigenSolver.preconditioner(), MMV, 0.5), dimension)); 
	Timing t3;
	FPI_2.compute();
	t3.write("secant");
	FPI_2.printResidualHistory();
	FPI_2.printResult();
	
	std::cout<<"\n*** APPLYING ANDERSON ACCELERATION:\n\n";
	FPI_2.reset();
	double mixingParameter = get_problem_data ("FPI_2_param/mixingParameter", 1.);
	std::size_t memory = get_problem_data ("FPI_2_param/AndersonMemory", 5);
	FPI_2.setIterator( std::make_unique <AndersonAccelerator>(SimpleRichardson < decltype(eigenSolver) >(MMM, eigenSolver.preconditioner(), MMV), dimension)); 
	Timing t4;
	FPI_2.compute();
	t4.write( "Anderson");
	FPI_2.printResidualHistory();
	FPI_2.printResult();
}				