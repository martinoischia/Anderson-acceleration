//! @file FixedPointIterator.cpp

#include "FixedPointIterator.hpp"

namespace FixedPoint {
	
	Traits::Vector FixedPointIterator::compute()
	{
		if (history.empty()) history.emplace_back ( Vector::Zero ( iterator-> getDimension()));
		
		double currentDistance = std::numeric_limits<double>::max();
		
		while(iteration < options.maxIter && currentDistance > options.tolerance)
		{
			auto & Old = history[history.size()-1] ;
			auto & New = history.emplace_back ( (*iterator)(history) );
			currentDistance = distance( New, Old );
			++iteration;
			if (history.size() > options.memory) history.pop_front();
		}		
		
		return history.back();
	}
	
	void FixedPointIterator::printResult (std::ostream & OS) const
	{
		auto & current = history.back();
		auto & previous = history[history.size()-2];
		auto distance = Traits::distance ( current , previous ) ;
		auto residual = Traits::distance ((iterator->getIterationFunction())(current), current );
		bool converged = (distance < options.tolerance) ;
		
		OS.precision(10);
		OS <<"\n The fixed-point iteration ";
		if (converged)
		{
			OS<<" converged ";
		}
		else
		{
			OS << "has not converged ";
			
		}
		OS <<" in "<< iteration<<" Iterations."<<std::endl;
		OS << "Residual has value " << residual << std::endl;
		OS << "___________________________________________________"<< std::endl;
	}
	
}							