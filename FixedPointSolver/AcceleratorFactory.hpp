#ifndef ACCELERATOR_FACTORY_HPP_
	#define ACCELERATOR_FACTORY_HPP_
	
	#include "Accelerators.hpp"
	
	namespace FixedPoint
	{
		enum IteratorType { NoAccelerator = 0, ASecantAccel = 1};
		
		//! A simple factory that returns an Iterator polymorphic object wrapped in a unique_prt
		//!
		//! /param kind The kind of Iterator class you want
		//! /param args Optional arguments to be forwarded to the constructor
		//!
		template <typename ...Args>
		std::unique_ptr<FixedPoint::Iterator> makeIterator (IteratorType ID, Args&&... args)
		{
			switch(ID)
			{
				case NoAccelerator:
				return std::make_unique<FixedPoint::Iterator>(std::forward<Args>(args)...);
				break;
				
				case ASecantAccel:
				return std::make_unique<FixedPoint::ASecantAccelerator>(std::forward<Args>(args)...);
				break;
				
				default:
				throw std::runtime_error("Error in IteratorType: You must specify a valid type");
			}
		}
	}
	
	#endif	