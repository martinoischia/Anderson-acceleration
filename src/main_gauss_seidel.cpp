//! @file main_gauss_seidel.cpp
//! @brief A simple discretization of a differential problem
//!
//! Solves 1D heat equation on a bar through Gauss-Seidel fixed-point iterations.
//! Implements also an alternating Anderson accelerator.
#include "gauss_seidel.cpp"
#include "FixedPointIterator.hpp"
#include "GetPot"
#include "Timing.hpp"
#include <iomanip>

int main()
{
	std::cout << "\n********************EXAMPLE 3****************************\n\n";
	std::cout << std::fixed;
	using namespace FixedPoint;
	using Vector = Traits::Vector;
	GetPot get_problem_data ("data.input");
	int dimension = get_problem_data ("FPI_3_param/M", 7);
	GaussSeidelIterator solver(dimension, Vector::Zero(dimension), 0.);
	
	// Initial vector: a linear variation of T
	Vector theta(dimension);
	for(int m=0; m<dimension; ++m)
	theta[m]=(1.-(m+1)*solver.h)*(solver.To-solver.Te)/solver.Te;
	// Analitic solution
	Vector true_sol(dimension);
	for(int m=0; m<dimension; m++)
	true_sol[m]=solver.Te+(solver.To-solver.Te)*cosh(sqrt(solver.act)*(1-(m+1)*solver.h))/cosh(sqrt(solver.act));
	
	//The iterator object
	FixedPointIterator FPI_3;
	FPI_3.setIterator( std::make_unique <Iterator> (std::move(solver), dimension) );
	FPI_3.getOptions().maxIter = get_problem_data ("FPI_3_param/maxIter", 100);
	FPI_3.getOptions().memory = get_problem_data ("FPI_3_param/memory", 100);
	FPI_3.getOptions().tolerance = get_problem_data ("FPI_3_param/tolerance", 0.01);
	FPI_3.getOptions().printMemory = get_problem_data ("FPI_3_param/printMemory", 100);
	//Solving
	std::cout<<"\n*** WITH BASIC METHOD:\n\n";
	Timing t;
	FPI_3.compute(theta);
	t.write("fixed-point iterations");
	FPI_3.printResidualHistory();
	std::cout<<"\n L^infinity norm of the error: " << ((FPI_3.getHistory().back()+Vector::Ones(dimension))*solver.Te-true_sol).maxCoeff()/26. <<std::endl;
	FPI_3.printResult();
	
	// Now with the Anderson Accelerator
	FPI_3.reset();
	double mixingParameter = get_problem_data ("FPI_3_param/mixingParameter", 1.);
	std::size_t memory = get_problem_data ("FPI_3_param/AndersonMemory", 5);
	FPI_3.setIterator( std::make_unique <AndersonAccelerator> (std::move (FPI_3.getIterator().getIterationFunction()), dimension, mixingParameter, memory) );
	//Solving
	std::cout<<"\n*** WITH ANDERSON ACCELERATION:\n\n";
	t.reset();
	bool good = FPI_3.compute(theta);
	t.write("Anderson method");
	FPI_3.printResidualHistory();
	std::cout<<"\n L^infinity norm of the error: " << ((FPI_3.getHistory().back()+Vector::Ones(dimension))*solver.Te-true_sol).maxCoeff()/26. <<std::endl;
	FPI_3.printResult();
	// Just for my postprocessing
	if ( good ) std::cout << "true" <<std::endl;
	
	
	// Now alternating Gauss-Seidel iterations to Anderson iterations
	// I swap unique pointers to Iterator objects
	FPI_3.reset();
	FPI_3.setIterator( std::make_unique <Iterator> (std::move (FPI_3.getIterator().getIterationFunction()), dimension) );
	GaussSeidelIterator solver_copy(dimension, Vector::Zero(dimension), 0.);
	std::unique_ptr<Iterator> IteratorToSwap = std::make_unique <AndersonAccelerator> (std::move (solver_copy), dimension, mixingParameter, memory);
	double alternatingParameter = get_problem_data ("FPI_3_param/alternatingParameter", 10);
	//Solving
	std::cout<<"\n*** WITH ALTERNATING ANDERSON METHOD:\n\n";
	auto & stop = FPI_3.getOptions().maxIter;
	auto maxIter = get_problem_data ("FPI_3_param/maxIter", 100);
	bool converged = false;
	unsigned int Iterations = 0;
	t.reset();
	while( !converged and Iterations < maxIter ){
		stop = alternatingParameter;
		FPI_3.getIteration() = 0;
		converged = FPI_3.compute( theta );
		Iterations += FPI_3.getIteration();
		if( !converged ) 
		{
			FPI_3.getIterator_p().swap( IteratorToSwap );
			stop = 1;
			FPI_3.getIteration() = 0;
			FPI_3.getIterator().SetUp( FPI_3.getHistory() );
			converged = FPI_3.compute();
			++Iterations;
			FPI_3.getIterator_p().swap( IteratorToSwap );
		}
	}
	t.write("Alternating Anderson method");
	FPI_3.printResidualHistory();
	std::cout<<"\n L^infinity norm of the error: " << ((FPI_3.getHistory().back()+Vector::Ones(dimension))*solver.Te-true_sol).maxCoeff()/26. <<std::endl;
	std::cout.precision(10);
	std::cout <<"\n The fixed-point iteration ";
	if (converged)
	{
		std::cout<<" converged ";
	}
	else
	{
		std::cout << "has not converged ";
		
	}
	std::cout <<" in "<< Iterations<<" Iterations."<<std::endl;
	std::cout << "________________________________________________________"<< std::endl;	
}	