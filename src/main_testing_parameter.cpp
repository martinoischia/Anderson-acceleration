#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <vector>
#include <sstream>

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
	std::vector <double > beta { 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1., 2., 5., 1e1, 50, 1e2, 1e3 };
	std::vector <double > memory {2, 4, 6, 8, 10, 15, 20, 25};
	
	std::ostringstream tmp;
	std::ostringstream tmp2;
	
	std::cout << "empty";
	for (int j = 0; j < beta.size(); ++j) std::cout <<" " << beta[j];
	std::cout << std::endl;
	
	for (int i = 0; i < memory.size(); ++i){
		std::cout << memory[i];
		tmp2 << "AndersonMemory = " << memory[i] << std::flush;
		std::system ( ("sed -i '43s/.*/" + tmp2.str() + "/' data.input").c_str() );
		tmp2.str("");
		tmp2.clear();	
		
		for (int j = 0; j < beta.size(); ++j){
			
			tmp << "mixingParameter = " << beta[j] << std::flush;
			std::system ( ("sed -i '42s/.*/" + tmp.str() + "/' data.input").c_str() );
			tmp.str("");
			tmp.clear();
			
			double sum = 0;
			if ( ! stoi (exec("./main_gauss_seidel | grep true >/dev/null; echo $?")) ){
				
				for (int k = 0 ; k < stat_size ; ++k) {
					num_as_string = exec("./main_gauss_seidel | grep 'TIME: Anderson method' | tr -dc '.0-9'");
					sum += stod (num_as_string) ;
				}
				
				sum/=stat_size;
				std::cout <<" " << sum ;	
			}
			else std::cout <<" nan" ;
		}		
		std::cout << std::endl;
	}
	
	std::system ( "sed -i '42s/.*/mixingParameter = 1/' data.input" );
	std::system ( "sed -i '43s/.*/AndersonMemory = 10/' data.input" );
	
	
}		