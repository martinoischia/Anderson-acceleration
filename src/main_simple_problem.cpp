//! @file main_simple_problem.cpp
//! @brief A simple example for testing
//!
//! It contains a simple fixed-point problem that is solved with different acceleration methods.
//! Many parameters can be read from input through the GetPot utility.
#include "FixedPointIterator.hpp"
#include "GetPot"

int main(int argc, char** argv)
{
	std::cout << "\n********************EXAMPLE 1****************************\n\n";	
	using namespace FixedPoint;
	GetPot get_filename (argc, argv);
	std::string filename = get_filename ("filename", "data.input");
	
	GetPot get_problem_data (filename.c_str ());
	
	using Vector = Traits::Vector; // Traits is a macro defined in Accelerators.hpp, if not defined in other ways
	using Matrix = Traits::Matrix;
	
	
	// Two simple fixed point problems: I commented the one in dimension 2
	// and used the one in dimension 3
	
	// double lambda = get_problem_data ("FPI_1_param/lambda",  18.);
	// // a function from R^2 to R^2... 
	// std::size_t dimension = 2;
	// auto phi_1 = [lambda, dimension] (Vector const & x){
	// Vector vec(dimension);
	// vec.coeffRef(0) = lambda*std::sin( x.coeff(0) );
	// vec.coeffRef(1) = std::cos( x.coeff(1) );
	// return vec; 
	// };
	// // ... and a starting point
	// Vector startingPoint_1(dimension);
	
	// a function from R^3 to R^3...
	std::size_t dimension = 3;
	auto phi_1 = [dimension] (Vector const & x){
	Vector vec(dimension);
	vec.coeffRef(0) = -1./81.*std::cos( x.coeff(0) ) + 1./9.* x.coeff(1)*x.coeff(1) + 1./3. * std::sin(x.coeff(2));
	vec.coeffRef(1) = 1./3.*std::sin( x.coeff(0) ) +1./3.*std::cos(x.coeff(2));
	vec.coeffRef(2) = -1./9.*std::cos( x.coeff(0) ) + 1./3.* x.coeff(1) + 1./6. * std::sin(x.coeff(2));
	return vec; 
	};
	// ... and a starting point
	Vector startingPoint_1(dimension);
	startingPoint_1.coeffRef(0) = get_problem_data ("FPI_1_param/x00", 5.);
	startingPoint_1.coeffRef(1) = get_problem_data ("FPI_1_param/x01", 7.);
	startingPoint_1.coeffRef(2) = 1.;
	
	//The iterator object
	FixedPointIterator FPI_1;
	FPI_1.setIterator( std::make_unique <Iterator> (phi_1, dimension) ); // care that phi_1 has been moved
	FPI_1.getOptions().maxIter = get_problem_data ("FPI_1_param/maxIter", 100);
	FPI_1.getOptions().memory = get_problem_data ("FPI_1_param/memory", 100);
	FPI_1.getOptions().tolerance = get_problem_data ("FPI_1_param/tolerance", 0.01);
	
	//Solving
	std::cout<<"\n*** WITH BASIC METHOD:\n\n";
	FPI_1.compute(startingPoint_1);
	FPI_1.printResidualHistory();
	FPI_1.printResult();
	
	// Now with the acceleration provided by the alternate secant method
	FPI_1.reset();
	FPI_1.setIterator( std::make_unique <ASecantAccelerator> (std::move (FPI_1.getIterator().getIterationFunction()), dimension) );
	
	//Solving
	
	std::cout<<"\n*** WITH SECANT ACCELERATION:\n\n";
	FPI_1.compute(startingPoint_1);
	FPI_1.printResidualHistory();
	FPI_1.printResult();
	
	// Now with the Anderson Accelerator
	FPI_1.reset();
	double mixingParameter = get_problem_data ("FPI_1_param/mixingParameter", 1.);
	std::size_t memory = get_problem_data ("FPI_1_param/AndersonMemory", 5);
	FPI_1.setIterator( std::make_unique <AndersonAccelerator> (std::move (FPI_1.getIterator().getIterationFunction()), dimension, mixingParameter, memory) );
	
	//Solving
	std::cout<<"\n*** WITH ANDERSON ACCELERATION:\n\n";
	FPI_1.compute(startingPoint_1);
	FPI_1.printResidualHistory();
	FPI_1.printResult();	
	}