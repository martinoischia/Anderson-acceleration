//! @file main_nonlinear.cpp
//! @brief Discretization of a nonlinear differential problem
//!
//! Solves 1D heat equation on a bar with a thermal radiation term (~)
#include "FixedPointIterator.hpp"
#include "GetPot"
#include "Timing.hpp"
#include <stdexcept>
// #include <iomanip>

//! @brief Gauss-Seidel solver
//!
//! It is taken from the HeatExchange directory in the repository PACS of Luca Formaggia,
//! have a look there to understand what is going on.
namespace FixedPoint
{
	
	class LinearSolver2
	{
		public:
		using Vector = Traits::Vector;
		//constructor
		LinearSolver2(const int size, const Vector & v): M(size), x(v){};
		//call operator
		Vector operator ()(const Vector & z) const
		{
			Vector y(M);
			y[0] = ((To-Te)/Te + z[1])/( 2.+h*h* ( act+ sigma *x[0]*x[0]*x[0] ));
			for(int m=1; m < M-1; ++m)
			{  
				y[m]  = (y[m-1]+z[m+1])/(2.+h*h* ( act + sigma *x[m]*x[m]*x[m] ));
			}      
			y[M-1] = 1./( 1.+ sigma * h*h*x[M-1]*x[M-1]*x[M-1] ) * y[M-2];
			return y;
		}
		static constexpr double L= 40.;  // Bar length
		static constexpr double a1=4.; // First longitudinal dimension
		static constexpr double a2=50; //  Second longitudinal dimension
		static constexpr double To=46.; // Dirichlet condition
		static constexpr double Te=20.; // External temperature (Centigrades)
		static constexpr double k=0.164;  // Thermal conductivity
		static constexpr double hc=200.0e-6; // Convection coefficient
		static constexpr double act=2.*(a1+a2)*hc*L*L/(k*a1*a2);
		static constexpr double sigma = 1.; //Radiation term
		const int M; // Number of grid elements
		//! Precomputed coefficient for adimensional form of equation
		// mesh size
		const double h=1./M;
		const Vector & x;
	};
	
	
	class FixedPointLinearSolver
	{
		public:
		using Vector = Traits::Vector;
		Vector operator ()(Vector const & x) const {
			using namespace FixedPoint;
			GetPot get_problem_data ("data.input");
			int dimension = get_problem_data ("FPI_3_param/M", 7);
			LinearSolver2 solver(dimension, x);
			
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
			LinearSolver2 solver_copy(dimension, x);
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
			return FPI_3.getHistory().back();
		}
		
	};
	

	
	
	
	class LinearSolver3
	{
		public:
		using Vector = Traits::Vector;
		//constructor
		LinearSolver3(const int size, const Vector & v): M(size), x(v){};
		//call operator
		Vector operator ()(const Vector & z) const
		{
			Vector rhs = f();
			Vector y(M);
			y[0] = ( rhs[0] + z[1])/( 2.+h*h* ( act+ 4*sigma *x[0]*x[0]*x[0] ));
			for(int m=1; m < M-1; ++m)
			{  
				y[m]  = (y[m-1]+z[m+1]+rhs[m])/(2.+h*h* ( act + 4*sigma *x[m]*x[m]*x[m] ));
			}      
			y[M-1] = 1./( 1.+ sigma * x[M-1]*x[M-1]*x[M-1] ) * (y[M-2]+rhs[M-1]);
			return y;
		}
		static constexpr double L= 40.;  // Bar length
		static constexpr double a1=4.; // First longitudinal dimension
		static constexpr double a2=50; //  Second longitudinal dimension
		static constexpr double To=46.; // Dirichlet condition
		static constexpr double Te=20.; // External temperature (Centigrades)
		static constexpr double k=0.164;  // Thermal conductivity
		static constexpr double hc=200.0e-6; // Convection coefficient
		static constexpr double act=2.*(a1+a2)*hc*L*L/(k*a1*a2);
		static constexpr double sigma = 1.; //Radiation term
		const int M; // Number of grid elements
		//! Precomputed coefficient for adimensional form of equation
		// mesh size
		const double h=1./M;
		const Vector & x;
		
		private:
		Vector f () const{
			Vector y(M);
			y[0]= (2.+h*h* (act + sigma* x[0]*x[0]*x[0]))*x[0]-x[1]-(To-Te)/Te;
			for(int m=1; m < M-1; ++m){
				y[m]= (2.+h*h*(act + sigma* x[m]*x[m]*x[m]))*x[m]-x[m-1] - x[m+1];
			}
			y[M-1] = (1.+ sigma*h*h* x[M-1]*x[M-1]*x[M-1])*x[M-1]-x[M-2];
			return y;
		}
	};	
	
		class NewtonFixedPointLinearSolver
	{
		public:
		using Vector = Traits::Vector;
		Vector operator ()(Vector const & x) const {
			using namespace FixedPoint;
			GetPot get_problem_data ("data.input");
			int dimension = get_problem_data ("FPI_3_param/M", 7);
			LinearSolver3 solver(dimension, x);
			
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
			LinearSolver3 solver_copy(dimension, x);
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
			return (x - FPI_3.getHistory().back());
		}
		
	};
}
int main()
{
	std::cout << "\n********************EXAMPLE 4****************************\n\n";
	// std::cout << std::fixed;
	using namespace FixedPoint;
	using Vector = Traits::Vector;
	using Matrix = Traits::Matrix;
	GetPot get_problem_data ("data.input");
	int dimension = get_problem_data ("FPI_3_param/M", 7);
	
	// Initial vector: a linear variation of T
	Vector theta(dimension);
	{
		LinearSolver2 solver (dimension, Vector{} );
		for(int m=0; m<dimension; ++m)
		theta[m]=(solver.To-solver.Te)/solver.Te;
	}
	//The iterator object
	FixedPointIterator FPI_3;
	FPI_3.setIterator( std::make_unique <Iterator> ( FixedPointLinearSolver{}, dimension) );
	FPI_3.getOptions().maxIter = get_problem_data ("FPI_3_param/maxIter", 100);
	FPI_3.getOptions().memory = get_problem_data ("FPI_3_param/memory", 100);
	FPI_3.getOptions().tolerance = get_problem_data ("FPI_3_param/tolerance", 0.01);
	FPI_3.getOptions().printMemory = get_problem_data ("FPI_3_param/printMemory", 100);
	
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
	FPI_3.setIterator( std::make_unique <Iterator> (NewtonFixedPointLinearSolver{}, dimension ));
	
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
	double mixingParameter = get_problem_data ("FPI_3_param/mixingParameter", 1.);
	std::size_t memory = get_problem_data ("FPI_3_param/AndersonMemory", 5);
	FPI_3.setIterator( std::make_unique <AndersonAccelerator> ( FixedPointLinearSolver{}, dimension, mixingParameter, memory) );
	
	//Solving
	std::cout<<"\n*** WITH ANDERSON ACCELERATION:\n\n";
	t.reset();
	bool good = FPI_3.compute(theta);
	t.write("Anderson method");
	FPI_3.printResidualHistory();
	FPI_3.printResult();
	Vector v3 = FPI_3.getHistory().back();

	// Now alternating Gauss-Seidel iterations to Anderson iterations
	// I swap unique pointers to Iterator objects
	FPI_3.reset();
	FPI_3.setIterator( std::make_unique <Iterator> (std::move (FPI_3.getIterator().getIterationFunction()), dimension) );
	std::unique_ptr<Iterator> IteratorToSwap = std::make_unique <AndersonAccelerator> (FixedPointLinearSolver{}, dimension, mixingParameter, memory);
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
	std::cout << "norm of sol1-sol2 " << (v1-v2).norm() << std::endl;
	std::cout << "norm 1-3 " << (v1-v3).norm() << std::endl;
	std::cout << "norm 1-4 " << (v1-v4).norm() << std::endl;
}			