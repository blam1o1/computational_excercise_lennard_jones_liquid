#include <iostream>
#include <fstream>
#include "../include/params.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <parameter file>" << std::endl;
        return 1;
    }
    Params params = read_params(argv[1]);
    check_params(params);
    print_params(params, std::cout);
    return 0;   
}