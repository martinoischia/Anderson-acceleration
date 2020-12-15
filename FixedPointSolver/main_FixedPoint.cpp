#include "FixedPointIterator.hpp"
#include "GetPot"
#include "MM_readers3.3.hpp"
#include "AcceleratorFactory.hpp"

int main(int argc, char** argv)
{
	using namespace FixedPoint;
	GetPot get_filename (argc, argv);
	std::string filename = get_filename ("filename", "data.txt");
	
	GetPot get_problem_data (filename.c_str ());
	
	double lambda = get_problem_data ("lambda",  18.);
	using Vector = Traits::Vector; // Traits is a macro defined in Accelerators.hpp, if not defined in other ways
	using Matrix = Traits::Matrix;
	// a function from R^2 to R^2... 
	auto phi_1 = [lambda] (Vector const & x){
		Vector vec(2);
		vec.coeffRef(0) = lambda*std::sin( x.coeff(0) );
		vec.coeffRef(1) = std::cos( x.coeff(1) );
		return vec; 
		};
	// ... and a starting point
	Vector startingPoint_1(2);
	startingPoint_1.coeffRef(0) = 5;
	startingPoint_1.coeffRef(1) = 7;
	
	//The iterator object
	FixedPointIterator FPI_1;
	FPI_1.setIterator( makeIterator (NoAccelerator, phi_1, 2) ); // care that phi_1 has been moved
	
	//Solving
	std::cout<<"*** WITH BASIC METHOD:\n";
	FPI_1.compute(startingPoint_1);
	FPI_1.printResult();
	
	// Now with the acceleration provided by the alternate secant method
	FPI_1.reset();
	FPI_1.setIterator( makeIterator (ASecantAccel, FPI_1.getIterator().getIterationFunction(), 2) );
	
	//Solving
	std::cout<<"*** WITH SECANT ACCELERATION:\n";
	FPI_1.compute(startingPoint_1);
	FPI_1.printResult();
	
	// Now we solve iteratively a sparse linear system.
	#include<Eigen/IterativeLinearSolvers>
	// reading the matrix and the RHS
	std::string matrix_filename = get_problem_data ( "MMMfile", "cavity15.mtx" );
	std::string vector_filename = get_problem_data ( "MMVfile", "cavity15_rhs1.mtx" );
	Matrix MMM = Eigen::read_MM_Matrix <Matrix> ( matrix_filename );
	Vector MMV = Eigen::read_MM_Vector <Vector> ( vector_filename );
	
	// Eigen::BiCGSTAB <SparseMatrix<double> > solver;
	
	// // a starting point
	// Vector startingPoint_2 = Vector.setRandom (MMV.size());
	
	
	// FixedPointIterator FPI_2 { makeIterator (NoAccelerator, ???)};
	
	// //Solving
	// std::cout<<"*** WITH BASIC METHOD:\n";
	// FPI_2.printResult(FPI_2.compute(startingPoint_2));
	
	// // Now with the acceleration provided by the alternate secant method
	
	// FPI_2.setAccelerator( makeIterator (ASecantAccel, FPI_2.getIterator().getIterationFunction()) );
	
	// //Solving
	// std::cout<<"*** WITH SECANT ACCELERATION:\n";
	// FPI_2.printResult(FPI_2.compute(startingPoint_2));
}



