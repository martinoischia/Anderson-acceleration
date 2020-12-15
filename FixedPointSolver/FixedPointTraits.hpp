#ifndef NONLINSYSSOLVER_FIXEDPOINTTRAITS_HPP_
	#define NONLINSYSSOLVER_FIXEDPOINTTRAITS_HPP_
	#include <memory>
	#include <vector>
	#include <Eigen/Core>
	#include <Eigen/SparseCore>
	namespace FixedPoint
	{
		class VectorTraits // it should be extended, for now I have just made the main
		//compliant with the Eigen related traits
		{
			public:
			using Vector = std::vector <double>;
			using IterationFunction = std::function <Vector (Vector const &)>;
			double a;
			static std::ostream & print (const Vector& v, std::ostream & OS)
			{
				for (const auto & i: v) OS<< i << " ";
				OS << std::endl ;
				return OS;
			}
			
			static double distance(Vector const & current, Vector const & previous)
			{
				double res{0.0};
				for (std::size_t i=0; i<current.size();++i)
				res+=(previous[i]-current[i])*(previous[i]-current[i]);
				return std::sqrt(res);
			}
		};
		
		struct EigenTraits
		{
			using Vector = Eigen::Matrix<double,Eigen::Dynamic,1>;
			using IterationFunction = std::function <Vector (Vector const &)>;
			using Matrix = Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic>;
			
			static double distance (Vector const & current, Vector const & previous)
			{
				return (current-previous).norm() ;
			}
			
			static std::ostream & print (const Vector& v, std::ostream & OS)
			{
				for (std::size_t i = 0; i < v.size() ; ++i) OS << *(v.data()+i) << " ";
				OS << std::endl ;
				return OS;
			}
			
		};
		
		struct SparseTraits
		{
			using Vector = Eigen::Matrix<double,Eigen::Dynamic,1>;
			using IterationFunction= std::function <Vector (Vector const &)>;			
			using Matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
			
			static double distance (Vector const & current, Vector const & previous)
			{
				return (current-previous).norm() ;
			}
			
			static std::ostream & print (const Vector& v, std::ostream & OS)
			{
				for (std::size_t i = 0; i < v.size() ; ++i) OS << *(v.data()+i) << " ";
				OS << std::endl ;
				return OS;
			}
			
			//This one in case I want to use a sparse vector
			
			// static std::ostream & print (const Vector& v, std::ostream & OS)
			// {
				// for (Vector::InnerIterator it(v); it; ++it)
				// {
					// OS << it.index() <<" ";
					// OS << it.value() <<std::endl;
					// return OS;
				// }
			// }
			
		};
		
	}
	
	
	
	
#endif /* NONLINSYSSOLVER_FIXEDPOINTTRAITS_HPP_ */
