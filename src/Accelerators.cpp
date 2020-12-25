//! @file Accelerators.cpp
//! @brief Implementation of the acceleration algorithms
#include "Accelerators.hpp"
#include <Eigen/QR>

namespace FixedPoint
{	
	Traits::Vector  ASecantAccelerator::operator()(const std::deque < Vector > & past)
	{
		assert (!past.empty());
		Vector xNew (dimension);
		if (firstTime)
		{
			xNew = phi (past.back());
			phiOld = xNew;
			deltaXOld = xNew - past.back();
			firstTime = false;
		}
		else
		{
			Vector tmp (dimension), deltaXNew (dimension);
			xNew = phi (past.back());
			// compute \Delta x_n = phi(x_n)-x_n
			deltaXNew = xNew - past.back();
			
			// \Delta x_n - \Delta x_{n-1}
			tmp = deltaXNew - deltaXOld;
			// Get ||\Delta x_n - \Delta x_{n-1}||^2
			double norm2=tmp.squaredNorm();
			//to do: check norm2
			double factor = tmp.dot(deltaXNew);
			factor /= norm2;
			// Save phi(X_n) for later use
			tmp=xNew;
			// scale (phi(x_n)-\phi(x_{n-1}) factor and subtract
			// 
			for (std::size_t i=0;i<dimension;++i)
			xNew.coeffRef(i) -= factor * (xNew.coeffRef(i) - phiOld.coeffRef(i));
			phiOld=tmp;
			deltaXOld=deltaXNew;
		}
		
		return xNew ;
		
	}
	
	Traits::Vector  AndersonAccelerator::operator()(const std::deque < Vector > & past){
		
		assert (!past.empty());
		if ( !evaluated ) 
		{
			evaluate ( past );
			evaluated = true;
		}
		
		if ( evaluation_history.size() == 1 ) // evaluation_history size is kept controlled by the (Anderson) memory parameter, unlike past.size()
		{
			Vector solution = (1-mixingParameter) * past [0] +  mixingParameter * evaluation_history.back(); // if memory == 1 it's called simple mixing
			evaluation_history.emplace_back( phi (solution) );
			if ( memory== 1) evaluation_history.pop_front();			
			return ( solution ) ;
		}
		
		else
		{	
			int m_k = evaluation_history.size()-1;
			
			//build "delta x" matrix
			
			AndersonMatrix X ( dimension, m_k );
			for (int j = 0; j < m_k; ++j )
			{
				for (int i = 0; i < dimension; ++i )
				{
					X (i,j) = past [ past.size()-m_k+j ](i)- past [ past.size()-m_k+j-1 ](i) ;
				}		
			}
			
			//Build "delta residual" matrix
			AndersonMatrix F ( dimension, m_k );
			for (int j = 0; j < m_k; ++j )
			{
				for (int i = 0; i < dimension; ++i )
				{
					F (i,j) = 
					evaluation_history [ evaluation_history.size()-m_k+j ](i)- past [ past.size()-m_k+j ](i) -
					evaluation_history [ evaluation_history.size()-m_k+j-1 ](i) + past [ past.size()-m_k+j-1 ](i);
				}		
			}
			
			//Calculating the solution using QR decomposition for solving least square problem
			
			Vector f = evaluation_history.back() - past.back() ;
			Vector solution ;
			solution = past.back() + mixingParameter*f - ( X + mixingParameter*F)*F.colPivHouseholderQr().solve(f) ;
			
			evaluation_history.emplace_back( phi ( solution ));
			if (evaluation_history.size() > memory) 
			evaluation_history.pop_front();
			
			return solution;
		}
	}
	
}	




