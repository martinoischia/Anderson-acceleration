//! @file simpleRichardson.hpp
//! @brief Makes an iterative solver for linear systems compatible for being accelerated
//!
//! Implements stationary Richardson iterative algorithm.


#include<Eigen/IterativeLinearSolvers>

namespace FixedPoint {
	//! Variable members are saved as references because Eigen has deleted
	//! copy and move constructor of preconditioners,
	//! making it not usable for my FixedPoint classes
	//! @tparam Solver_t an Eigen sparse iterative solver
	template <typename Solver_t>
	class SimpleRichardson: public Traits
	{
		public:
		
		using Precond_t = decltype(std::declval < Solver_t >().preconditioner());
		
		SimpleRichardson( Matrix & m, Precond_t & p, Vector & v, double par = 1. ):
		A (m),
		precond (p), 
		b (v),
		param (par)
		{}
		
		Vector operator()( Vector const & x){
			auto residual = A*x - b;
			return x - param*precond.solve(residual);			
		}
		
		private:
		
		//! relaxation parameter
		double param;
		Precond_t & precond;
		Matrix & A;
		Vector & b ;
	};
}					