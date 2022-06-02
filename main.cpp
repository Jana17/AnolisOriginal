//
//  main.cpp
//  Jana
//
//  Created by thijsjanzen on 01/03/2022.
//

#include "Anolis.h"

int main() {

    Param parameters;
  //  parameters.set_parameters("ParameterFile.txt");

    std::cout << parameters.number_of_timesteps << "\n";
    
    auto t1 = std::chrono::system_clock::now();

    for (int seed = parameters.start_seed; seed < parameters.end_seed; ++seed) {
        Simulation sim(parameters, seed);
        sim.run();
        auto t2 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = t2 - t1;
        t1 = t2;
        std::cout << "This took: " << elapsed_seconds.count() << " seconds\n";
    }
    return 0;
}
