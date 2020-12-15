#include "Accelerators.hpp"
namespace FixedPoint
{
	Traits::Vector  ASecantAccelerator::operator()(const std::deque < Vector > & past)
	{
		Vector xNew(dimension);
		
		if (firstTime)
		{
			xNew = phi (past.back());
			phiOld = xNew;
			deltaXOld = xNew - past.back();
			firstTime = false;
		}
		else
		{
			xNew = phi (past.back());
			// compute \Delta x_n = phi(x_n)-x_n
			deltaXNew = xNew - past.back();
			
			Vector tmp (dimension);
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
	
	// Traits::Vector  GeneralAndersonAccelerator::operator()(Traits::Vector const & x0){
	// auto dimension = x0.size() ;
	
	// if (history.empty())		{ history.emplace_back( x0 ) ;}
	// if (history.size == 1 ; ) { history.emplace_back ( phi (x0)) ;}
	
	// int k = history.size() - 1;
	// int m_k = std::min ( m, k )
	// for ()
	
	// using ResidualMatrix = EigenTraits::Matrix;
	// ResidualMatrix F ( size, m );
	
	// Eigen::Map < Vector > v (F.data(), size);
	// v = history[1] - history[0] ;
	
	
	// int iter = history.size() ;
	
	// do{
	
	// new (&v) Eigen::Map < Vector> (F.data() + size, size);
	// // v = 
	// ++iter;
	// }
	// while (iter <= m) ;
	
	// do
	
	// std::size_t index ;
	
	// return x1;
	// }
	
}	




