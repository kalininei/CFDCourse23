#include "common.hpp"
#include "prog_common.hpp"

int main(){
	try{
		//

		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << "  " << e.what() << std::endl;
	}
}

