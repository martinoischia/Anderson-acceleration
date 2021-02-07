#ifndef __TIMING__HPP_____
	#define __TIMING__HPP_____
	#include <iostream>
	#include <sys/time.h>
	#include <string>
	
	class Timing
	{
		private:
		
		double st;
		double gettime(){return clock();}
		
		public:
		Timing(){st = gettime();}
		void reset(){st = gettime();}
		double time(){return (gettime()- st)/CLOCKS_PER_SEC;}
		void write(const std::string & mes){
			std::cout << "TIME: " + mes + " is done in " << time() <<" s" << std::endl;
			st = gettime();
		}
	};
	
#endif 			