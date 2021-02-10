//! @file main_test2_param.cpp
//! @brief Produces a table of execution times of main_gauss_seidel
//!
//! Compiling is trivial. As you can read, here we are varying the alternatingParameter
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <vector>
#include <sstream>

// Getting as string the output of a call to a shell command
std::string exec ( const char* cmd ) 
{
	std::array<char, 128> buffer;
	std::string result;
	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
	if (!pipe) {
		throw std::runtime_error("popen() failed!");
	}
	while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
		result += buffer.data();
	}
	return result;
}


int main(){
	std::string num_as_string;
	constexpr int stat_size = 10;
	std::vector <double > alternating {2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 50, 60, 80, 100, 150, 200};
	
	std::ostringstream tmp2;
	
	
	for (int i = 0; i < alternating.size(); ++i){
		std::cout << alternating[i];
		tmp2 << "alternatingParameter = " << alternating[i] << std::flush;
		std::system ( ("sed -i '45s/.*/" + tmp2.str() + "/' data.input").c_str() );
		tmp2.str("");
		tmp2.clear();	
		
		double sum = 0;
		
		for (int k = 0 ; k < stat_size ; ++k) {
			num_as_string = exec("./main_gauss_seidel | grep 'TIME: Alternating Anderson method' | tr -dc '.0-9'");
			sum += stod (num_as_string) ;
		}
		
		sum/=stat_size;
		std::cout <<" " << sum ;	
		std::cout << std::endl;
	}		
}		