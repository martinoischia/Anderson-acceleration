//! @file Accelerators.hpp
//! @brief Header file for the polymorphic family of iterators

#ifndef SRC_ACCELERATORS_HPP_
	#define SRC_ACCELERATORS_HPP_
	
	#include <deque>
	#include <cassert>
	#include "FixedPointTraits.hpp"
	#include <algorithm>
	
	#ifndef Traits
		//! All the classes in this project derive from Traits to have uniform types.
		#define Traits SparseTraits
	#endif
	namespace FixedPoint
	{
		//! @brief A function object producing iterations.
		
		//! This class is the parent class for the subsequent classes: what it does
		//! is providing, through the call operator, the next value of a (hopefully converging) sequence, based on an
		//! iteration function *phi* and previous iterations data.
		class Iterator : public Traits
		{
			public:
			
			//! The constructor
			template <class IterationFun>
			Iterator(IterationFun && IF , std::size_t dim): phi(std::forward<IterationFun>(IF)), dimension (dim) {}
			
			//! Function call operator, is defined at the end of this file
			
			//! @param past the sequence of past iterates used to compute a new value
			//! @return the computed value
			virtual inline Vector operator()(const std::deque < Vector > & past);
			
			virtual void reset(){}
			
			virtual ~Iterator(){};
			
			IterationFunction & getIterationFunction(){ return phi;}
			IterationFunction getIterationFunction() const { return phi;}
			
			std::size_t getDimension() const { return dimension ;}
			
			
			protected:
			IterationFunction phi;
			//! Fixed-point problems go from \f$\mathbb{R}^{dim}\f$ to \f$\mathbb{R}^{dim}\f$
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
		*/
		
		
		class ASecantAccelerator : public Iterator
		{
			public:
			
			//! The constructor just forwards arguments to the parent class constructor and
			//! allocates memory for the vectors.
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
		
		//! Anderson Accelerator
		
		//! Anderson, Donald, Iterative Procedures for Nonlinear Integral Equations J. ACM 12 (1965) 547-560, doi 10.1145/321296.321305
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
			//! A sort of relaxation parameter.
			double mixingParameter;
			
			//! How many past iterates are considered at each step.
			std::size_t memory;	
			
			//! Contains the values of the iterator function evaluated in the past iterations.
			std::deque < Vector > evaluation_history;
			
			bool evaluated = false ;
			
			//! Fills evaluation_history
			//! @param past the past iteration values
			void evaluate (const std::deque < Vector > & past) {
				evaluation_history.resize( std::min (past.size(), memory) );
				std::transform( past.cbegin() , past.cend() , evaluation_history.begin() , phi ) ;
			}
			
		};
		
		//! Here only the last value of the past container is used
		Traits::Vector Iterator::operator()(const std::deque < Vector > & past)
		{
			assert (!past.empty());
			return phi( past.back() );
		}
		
	}
#endif /* SRC_ACCELERATORS_HPP_ */					