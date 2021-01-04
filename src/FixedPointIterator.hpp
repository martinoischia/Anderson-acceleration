//! @file FixedPointIterator.hpp
//! @brief The class for solving a fixed-point problem
//!
//! Contains several options and the history of the iteration process.

#ifndef SRC_FIXEDPOINTITERATOR_HPP_
	#define SRC_FIXEDPOINTITERATOR_HPP_
	
	#include <iostream>
	#include "Accelerators.hpp"
	#include <limits>
	
	namespace FixedPoint
	{
		struct FixedPointOptions
		{
			double tolerance=1.e-8;
			unsigned int maxIter=100;
			//! number of data to keep track
			unsigned int memory = maxIter;
			//! number of data to print with printHistory function
			unsigned int print_memory = 30;
		};
		
		
		class FixedPointIterator: public Traits
		{
			public:
			
			//TOFIX
			// This code has lead me to an Internal Compiler Error. Making the class default constructable seems good.
			// FixedPointIterator(std::unique_ptr <Iterator> && iter, const FixedPointOptions& opt=FixedPointOptions{}):
			// iterator(std::move(iter)), options{opt}, iteration{0}, history{} {}
			
			//! Default constructor
			FixedPointIterator():iterator(nullptr), options{FixedPointOptions{}}, iteration{0}, history{} {}
			
			
			FixedPointOptions getOptions() const {
				return options;
			}
			
			FixedPointOptions & getOptions () {
				return options;
			}
			
			void setIterator( std::unique_ptr<Iterator> && iter){iterator=std::move(iter);}
			
			const Iterator & getIterator() const {return *iterator;}
			Iterator & getIterator(){return *iterator;}
			
			//! Calculates the sequence of iterates until convergence or maximum iteration is reached
			
			//! If it starts from the first iteration, a default initial vector of zeros is used.
			//! @return The last computed vector
			Vector compute ();
			//! @param x0 suggested initial value: it is used only if there are no past iterations.
			//! @return The last computed vector
			Vector compute ( Vector const & x0 ) {if ( history.empty() ) history.emplace_back ( x0 ); return compute(); }
			
			void reset() { 
				iterator->reset();
				history.clear();
				history.shrink_to_fit();
				iteration = 0;				
			}
			
			void printResult ( std::ostream & OS = std::cout) const;
			
			void printHistory ( std::ostream & OS = std::cout) const {
				unsigned int m = std::min (history.size(), static_cast<long unsigned int>(options.print_memory));
				std::cout<< "Last " << m << " values:\n";
				for (int i = 0 ; i < m ; ++i) Traits::print ( history [i] , OS) ;
			}
			
			//! probably to be redesigned but useful
			void printResidualHistory ( std::ostream & OS = std::cout) const {
				unsigned int m = std::min (history.size(), static_cast<long unsigned int>(options.print_memory));
				std::cout<< "Last " << m << " residual norms:\n";
				for (int i = history.size()-m ; i < history.size() ; ++i) OS << ((iterator->getIterationFunction())( history [i]) - history [i]).norm() << std::endl ;
			}
			
			const std::deque < Vector > & getHistory () const { return history ;}
			
			private:
			
			//! A pointer to a function object that produces the iterates
			std::unique_ptr<Iterator> iterator;
			FixedPointOptions options;
			//! History of iterations
			
			//! The maximum number it can contain is established by the memory option.
			std::deque < Vector > history;
			//! Number of computed iterations
			unsigned int iteration;
		};
		
		
	}
#endif /* SRC_FIXEDPOINTITERATOR_HPP_ */
