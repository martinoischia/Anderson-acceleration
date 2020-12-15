#ifndef SRC_NONLINSYSSOLVER_FIXEDPOINTITERATION_HPP_
	#define SRC_NONLINSYSSOLVER_FIXEDPOINTITERATION_HPP_

	#include <iostream>
	#include "Accelerators.hpp"
	#include <limits>
	
	namespace FixedPoint
	{
		struct FixedPointOptions
		{
			double tolerance=1.e-8;
			unsigned int maxIter=100;
			unsigned int memory = maxIter;
		};
		
		
		class FixedPointIterator: public Traits
		{
			public:
			
			//TOFIX
			// This code has lead me to an Internal Compiler Error. Making the class default constructable seems good.
			// FixedPointIterator(std::unique_ptr <Iterator> && iter, const FixedPointOptions& opt=FixedPointOptions{}):
				// iterator(std::move(iter)), options{opt}, iteration{0}, history{} {}
			
			
			FixedPointIterator():iterator(nullptr), options{FixedPointOptions{}}, iteration{0}, history{} {}
			
			
			FixedPointOptions getOptions() const {
				return options;
			}
			
			void setOptions(const FixedPointOptions& options) {
				this->options = options;
			}
			
			void setIterator( std::unique_ptr<Iterator> && iter){iterator=std::move(iter);}
			
			const Iterator & getIterator() const {return *iterator;}
			Iterator & getIterator(){return *iterator;}
			
			Vector compute ();
			Vector compute ( Vector const & x0 ) {if ( history.empty() ) history.emplace_back ( x0 ); return compute(); }
			
			void reset() { 
				history.clear();
				history.shrink_to_fit();
				iteration = 0;				
				}
			
			void printResult ( std::ostream & OS = std::cout) const;
			
			void printHistory ( std::ostream & OS = std::cout) const {
				for (int i = 0 ; i < history.size() ; ++i) Traits::print ( history [i] , OS) ;
			}
			
			const std::deque < Vector > & getHistory () const { return history ;}
			
			private:
			
			std::unique_ptr<Iterator> iterator;
			FixedPointOptions options;
			std::deque < Vector > history;
			unsigned int iteration;
		};
		
		
	}
#endif /* SRC_NONLINSYSSOLVER_FIXEDPOINTITERATION_HPP_ */
