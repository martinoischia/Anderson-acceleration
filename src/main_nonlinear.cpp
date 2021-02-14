//! @file main_nonlinear.cpp
//! @brief Discretization of a nonlinear differential problem
//!
//! Solves 1D heat equation on a bar with a thermal radiation term (~)
#include "FixedPointIterator.hpp"
#include "GetPot"
#include "Timing.hpp"
#include <stdexcept>
#include "gauss_seidel.cpp"
// #include <iomanip>

namespace FixedPoint
{	
	//! @brief It is basically the same alternating solver that was used in main_gauss_seidel.cpp
	//!
	//! Solves the linear system associated to the single iteration of the 
	//! nonlinear problem with Gauss-Seidel method
	class AlternatingSolver: public Traits
	{
		private: 
		//! Newton solver or not
		bool Newton;
		
		public:
		//! Constructor
		AlternatingSolver ( bool N = false ) : Newton (N) {};
		//! Call operator
		Vector operator ()(Vector const & x) const {
			GetPot get_problem_data ("data.input");
			int dimension = get_problem_data ("FPI_3_param/M", 7);
			GaussSeidelIterator solver(dimension, x, 1.);
			if (Newton) {
				solver.rhs = solver.f() ;
				solver.diagonal  = 2*Vector::Ones(dimension).array()+
				solver.h*solver.h*( solver.act*Vector::Ones(dimension).array() +
				4.*solver.sigma*x.array()*x.array()*x.array()) ;
			}
			// Initial vector: a linear variation of T
			Vector theta(dimension);
			for(int m=0; m<dimension; ++m)
			theta[m]=(1.-(m+1)*solver.h)*(solver.To-solver.Te)/solver.Te;
			
			//The iterator object
			FixedPointIterator FPI_3;
			FPI_3.setIterator( std::make_unique <Iterator> (std::move(solver), dimension) );
			FPI_3.getOptions().memory = get_problem_data ("FPI_3_param/memory", 100);
			FPI_3.getOptions().tolerance = get_problem_data ("FPI_3_param/tolerance", 0.01);
			double mixingParameter = get_problem_data ("FPI_3_param/mixingParameter", 1.);
			std::size_t memory = get_problem_data ("FPI_3_param/AndersonMemory", 5);
			
			// Alternating Gauss-Seidel iterations to Anderson iterations
			GaussSeidelIterator solver_copy(dimension, x, 1.);
			if (Newton) {
				solver_copy.rhs = solver_copy.f() ;
				solver_copy.diagonal  = 2*Vector::Ones(dimension).array()+
				solver_copy.h*solver_copy.h*( solver_copy.act*Vector::Ones(dimension).array() +
				4.*solver_copy.sigma*x.array()*x.array()*x.array()) ;
				
			}
			std::unique_ptr<Iterator> IteratorToSwap = std::make_unique <AndersonAccelerator> (std::move (solver_copy), dimension, mixingParameter, memory);
			double alternatingParameter = get_problem_data ("FPI_3_param/alternatingParameter", 10);
			//Solving
			auto & stop = FPI_3.getOptions().maxIter;
			auto maxIter = get_problem_data ("FPI_3_param/maxIter", 100);
			bool converged = false;
			unsigned int Iterations = 0;
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
			if (!converged) throw std::runtime_error( "the linear step did non converge" ) ;
			// std::cout <<"inner iterations " << Iterations << std::endl;
			return (Newton? (x - FPI_3.getHistory().back()) : FPI_3.getHistory().back());
		}
		
	};
}
int main()
{
	std::cout << "\n********************EXAMPLE 4****************************\n\n";
	// std::cout << std::fixed;
	using namespace FixedPoint;
	using Vector = Traits::Vector;
	GetPot get_problem_data ("data.input");
	int dimension = get_problem_data ("FPI_3_param/M", 7);
	
	// Initial vector: constant T
	Vector theta(dimension);
	for(int m=0; m<dimension; ++m)
	theta[m]=(GaussSeidelIterator::To-GaussSeidelIterator::Te)/GaussSeidelIterator::Te;
	//The iterator object
	FixedPointIterator FPI_3;
	FPI_3.setIterator( std::make_unique <Iterator> ( AlternatingSolver(), dimension) );
	FPI_3.getOptions().maxIter = get_problem_data ("newton_param/maxIter", 100);
	FPI_3.getOptions().memory = get_problem_data ("newton_param/memory", 100);
	FPI_3.getOptions().tolerance = get_problem_data ("newton_param/tolerance", 0.01);
	FPI_3.getOptions().printMemory = get_problem_data ("newton_param/printMemory", 100);
	
	//Solving
	std::cout<<"\n*** WITH BASIC METHOD:\n\n";
	Timing t;
	FPI_3.compute(theta);
	t.write("fixed-point iterations");
	// FPI_3.printHistory();
	FPI_3.printResidualHistory();
	FPI_3.printResult();
	Vector v1 = FPI_3.getHistory().back();
	
	// Now with Newton method
	FPI_3.reset();
	FPI_3.setIterator( std::make_unique <Iterator> ( AlternatingSolver (true), dimension ) );
	//Solving
	std::cout<<"\n*** WITH NEWTON METHOD:\n\n";
	t.reset();
	FPI_3.compute(theta);
	t.write("Newton method");
	FPI_3.printResidualHistory();
	FPI_3.printResult();
	Vector v2 = FPI_3.getHistory().back();
	
	// Now with the Anderson Accelerator
	FPI_3.reset();
	double mixingParameter = get_problem_data ("newton_param/mixingParameter", 1.);
	std::size_t memory = get_problem_data ("newton_param/AndersonMemory", 5);
	FPI_3.setIterator( std::make_unique <AndersonAccelerator> ( AlternatingSolver(), dimension, mixingParameter, memory) );
	//Solving
	std::cout<<"\n*** WITH ANDERSON ACCELERATION:\n\n";
	t.reset();
	bool good = FPI_3.compute(theta);
	t.write("Anderson method");
	FPI_3.printResidualHistory();
	FPI_3.printResult();
	Vector v3 = FPI_3.getHistory().back();
	
	// Now with alternating method
	FPI_3.reset();
	FPI_3.setIterator( std::make_unique <Iterator> ( AlternatingSolver(), dimension) );
	std::unique_ptr<Iterator> IteratorToSwap = std::make_unique <AndersonAccelerator> ( AlternatingSolver(), dimension, mixingParameter, memory);
	double alternatingParameter = get_problem_data ("newton_param/alternatingParameter", 10);
	//Solving
	std::cout<<"\n*** WITH ALTERNATING ANDERSON METHOD:\n\n";
	auto & stop = FPI_3.getOptions().maxIter;
	auto maxIter = get_problem_data ("newton_param/maxIter", 100);
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
	Vector v4 = FPI_3.getHistory().back();
	auto res = [dimension](const Vector & x){
		Vector y(dimension);
		y[0]= (2+1./dimension*1./dimension*( GaussSeidelIterator::act + x[0]*x[0]*x[0]))*x[0]-x[1]-(GaussSeidelIterator::To-GaussSeidelIterator::Te)/GaussSeidelIterator::Te;
		for(int m=1; m < dimension-1; ++m){
			y[m]= (2+1./dimension*1./dimension*( GaussSeidelIterator::act + x[m]*x[m]*x[m]))*x[m] - x[m-1] - x[m+1] ;
		}
		y[dimension-1] = x[dimension-1]-x[dimension-2] ;
		return y.norm();
	};
	
	
	// Now with Anderson applied to Newton method
	FPI_3.reset();
	stop = maxIter;
	FPI_3.setIterator( std::make_unique <AndersonAccelerator> ( AlternatingSolver (true), dimension, mixingParameter, memory ) );
	//Solving
	std::cout<<"\n*** WITH ANDERSON-NEWTON METHOD:\n\n";
	t.reset();
	FPI_3.compute(theta);
	t.write("Anderson-Newton method");
	FPI_3.printResidualHistory();
	FPI_3.printResult();
	Vector v5 = FPI_3.getHistory().back();
	
	std::cout << "proper residual:" << std::endl;
	std::cout << "base case " << res(v1) << std::endl;
	std::cout << "Newton " << res(v2) << std::endl;
	std::cout << "Anderson " << res(v3) << std::endl;
	std::cout << "Alternating Anderson " << res(v4) << std::endl;
	std::cout << "Anderson-Newton " << res(v5) << std::endl;
}			