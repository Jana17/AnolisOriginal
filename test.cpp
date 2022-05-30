#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <iostream>
#include "Anolis.h"
#include "catch.hpp"

#include <array>

TEST_CASE ("random number methods") {
    rnd_j rndgen(1);
    rnd_j rndgen2(1);
    
    auto v1 = rndgen.uniform();
    auto v2 = rndgen2.uniform();
    
    REQUIRE(v1 == v2);
    
    std::array<size_t, 2> count = {0, 0};
    for (size_t r = 0; r < 1000; ++r) {
        if (rndgen.flip_coin()) {
            count[1]++;
        } else {
            count[0]++;
        }
    }
    
    CHECK(count[0] - count[1] < 100); // give some leeway
    
    
    
    
    
    
}

