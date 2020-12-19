#ifndef SRC_NONLINSYSSOLVER_ACCELERATORS_HPP_
	#define SRC_NONLINSYSSOLVER_ACCELERATORS_HPP_
	
	#include <deque>
	#include <cassert>
	#include "FixedPointTraits.hpp"
	#include <algorithm>
	
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
			
			virtual inline Vector operator()(const std::deque < Vector > & past);
			
			virtual void reset(){}
			
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
			Iterator(std::forward<IterationFun>(IF), dim), deltaXOld(dim), phiOld(dim), firstTime(true) {}
			
			Vector operator()(const std::deque < Vector > &) override;
			
			void reset() override {firstTime=true;}
			
			private:
			
			Vector deltaXOld;
			Vector phiOld;
			
			bool firstTime;
		};
		
		class AndersonAccelerator : public Iterator
		{
			public:
			
			template <class IterationFun>
			AndersonAccelerator(IterationFun&& IF, std::size_t dim, double p =1, std::size_t m = 10): 
			Iterator(std::forward<IterationFun>(IF), dim ), mixingParameter(p), memory(std::min (m, dim)) {}
			
			Vector operator()(const std::deque < Vector > &) override;
			
			void reset() {
				evaluated = false ; 
				evaluation_history.clear();
				evaluation_history.shrink_to_fit();
			}
			
			private:
			
			double mixingParameter;
			
			std::size_t memory;	
			
			std::deque < Vector > evaluation_history;
			
			bool evaluated = false ;
			
			void evaluate (const std::deque < Vector > & past) {
				evaluation_history.resize( std::min (past.size(), memory) );
				std::transform( past.cbegin() , past.cend() , evaluation_history.begin() , phi ) ;
			}
			
		};
		
		
		Traits::Vector Iterator::operator()(const std::deque < Vector > & past)
		{
			assert (!past.empty());
			return phi( past.back() );
		}
		
	}
#endif /* SRC_NONLINSYSSOLVER_ACCELERATORS_HPP_ */					