#include "FixedPointIterator.hpp"
#include "GetPot"

class LinearSolver{
	public:
	using Vector = FixedPoint::Traits::Vector;
	//constructor
	LinearSolver(const int size): M(size){};
	//call operator
	Vector operator ()(const Vector x){
		Vector y(M);
		y[0] = ((To-Te)/Te + x[1])/(2.+h*h*act);
		for(int m=1; m < M-1; ++m)
		{  
			y[m]  = (y[m-1]+x[m+1])/(2.+h*h*act);
		}      
      y[M-1] = y[M-2];
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
	const int M; // Number of grid elements
	//! Precomputed coefficient for adimensional form of equation
	// mesh size
	const double h=1./M;
	
};

int main()
{
	std::cout << "\n********************Gauss Seidel****************************\n\n";
	using namespace FixedPoint;
	using Vector = Traits::Vector;
	using Matrix = Traits::Matrix;
	GetPot get_problem_data ("data.input");
	int dimension = get_problem_data ("FPI_3_param/M", 7);
	LinearSolver solver(dimension);
	
	// Gauss Siedel is initialised with a linear variation
	// of T
	Vector theta(dimension);
	for(int m=0; m<dimension; ++m)
	theta[m]=(1.-(m+1)*solver.h)*(solver.To-solver.Te)/solver.Te;
	// Analitic solution
	Vector thetaa(dimension);
	for(int m=0; m<dimension; m++)
	thetaa[m]=solver.Te+(solver.To-solver.Te)*cosh(sqrt(solver.act)*(1-(m+1)*solver.h))/cosh(sqrt(solver.act));

	//The iterator object
	FixedPointIterator FPI_3;
	FPI_3.setIterator( std::make_unique <Iterator> (std::move(solver), dimension) );
	FPI_3.getOptions().maxIter = get_problem_data ("FPI_3_param/maxIter", 100);
	FPI_3.getOptions().memory = get_problem_data ("FPI_3_param/memory", 100);
	FPI_3.getOptions().tolerance = get_problem_data ("FPI_3_param/tolerance", 0.01);
	
	//Solving
	std::cout<<"\n*** WITH BASIC METHOD:\n\n";
	FPI_3.compute(theta);
	FPI_3.printResidualHistory();
	std::cout<<"\n norm of true error: " << ((FPI_3.getHistory().back()+Vector::Ones(dimension))*solver.Te-thetaa).norm() <<std::endl;
	FPI_3.printResult();
	// Now with the Anderson Accelerator
	FPI_3.reset();
	double mixingParameter = get_problem_data ("FPI_3_param/mixingParameter", 1.);
	std::size_t memory = get_problem_data ("FPI_3_param/AndersonMemory", 5);
	FPI_3.setIterator( std::make_unique <AndersonAccelerator> (std::move (FPI_3.getIterator().getIterationFunction()), dimension, mixingParameter, memory) );
	
	//Solving
	std::cout<<"\n*** WITH ANDERSON ACCELERATION:\n\n";
	FPI_3.compute(theta);
	FPI_3.printResidualHistory();
	std::cout<<"\n norm of true error: " << ((FPI_3.getHistory().back()+Vector::Ones(dimension))*solver.Te-thetaa).norm() <<std::endl;
	FPI_3.printResult();	

}