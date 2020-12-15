#ifndef SRC_NONLINSYSSOLVER_ACCELERATORS_HPP_
	#define SRC_NONLINSYSSOLVER_ACCELERATORS_HPP_
	
	#include <deque>
	#include <cassert>
	#include "FixedPointTraits.hpp"
	
	#ifndef Traits
		#define Traits SparseTraits
	#endif
	namespace FixedPoint
	{		
		class Iterator : public Traits
		{
			public:
			
			template <class IterationFun>
			Iterator(IterationFun&& IF , std::size_t dim): phi(std::forward<IterationFun>(IF)), dimension (dim) {}
			
			virtual Vector operator()(const std::deque < Vector > & past){
				assert (!past.empty());
				return phi( past.back() );
			}
			
			virtual ~Iterator(){};
			
			IterationFunction & getIterationFunction(){ return phi;}
			IterationFunction getIterationFunction() const { return phi;}
			
			std::size_t getDimension() const { return dimension ;}
			
			
			protected:
			IterationFunction phi;
			std::size_t dimension;
		};
		
		//! Alternate secant method
		/*!
			* This is a two-level acceleration method equivalent to a two level Anderson method
			* \f[
			* x_{n+1}=\phi(x_n)-\frac{(\Delta x_n -\Delta x_{n-1})\cdot\Delta x_n}{||\Delta x_n -\Delta x_{n-1}||^2}(\phi(x_n)-\phi(x_{n-1})
			* \f]
			*
			* V. Eyert, A Comparative Study on Methods for Convergence Acceleration of Iterative Vector Sequences, Journal of Computational Physics 124 (1996) 271–285.
			* H. Fang, Y. Saad, Two classes of multisecant methods for nonlinear accel- eration, Numerical Linear Algebra with Applications 16 (2009) 197–221.
			*
			* \tparam IteratorFunction. The iterator function. Note that at this level the iterator function
			* is used only as a trait to extract the type used for the arguments. No other requirements are made
			*
		*/
		
		
		class ASecantAccelerator : public Iterator
		{
			public:
			
			template <class IterationFun>
			ASecantAccelerator(IterationFun&& IF, std::size_t dim):
				Iterator(std::forward<IterationFun>(IF), dim), deltaXNew(dim), deltaXOld(dim), phiOld(dim) {}
			
			Vector operator()(const std::deque < Vector > &) override;
			
			void reset()
			{
				firstTime=true;
			}
			
			private:
			
			Vector deltaXNew;
			Vector deltaXOld;
			Vector phiOld;
			bool firstTime=true;
		};
		
		class GeneralAndersonAccelerator : public Iterator
		{
			public:
			template <class IterationFun>
			GeneralAndersonAccelerator(IterationFun&& IF, std::size_t dim, int memory = 5): 
			Iterator(std::forward<IterationFun>(IF), dim ), m(memory) {}
			
			Vector operator()(const std::deque < Vector > &) override;
			
			private:
			
			int m;
			
			
		};
		
	}
#endif /* SRC_NONLINSYSSOLVER_ACCELERATORS_HPP_ */					