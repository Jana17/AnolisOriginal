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
    
    // mutate trait
    rndgen.set_mutate_trait_dist(0.01);
    rndgen.set_mutate_prob(1.0);
    
    for (size_t r = 0; r < 100; ++r) {
        double old_trait = rndgen.uniform();
        double new_trait = rndgen.mutate_trait(old_trait);
        CHECK(new_trait >= 0.0);
        CHECK(new_trait <= 1.0);
        CHECK(new_trait != old_trait);
    }
    
    // draw random niche
    
    size_t num_niches = 6;
    size_t num_repl = 100;
    
    std::array< int, 600 > found_niches;
    for (size_t current_niche = 0; current_niche < num_niches; ++current_niche) {
        for (size_t r = 0; r < num_repl; ++r) {
            auto new_niche = rndgen.draw_random_niche(current_niche, num_niches);
            CHECK(new_niche != current_niche);
            CHECK(new_niche < num_niches);
            found_niches[current_niche * num_repl + r] = new_niche;
        }
    }
    
    auto max_niche = *std::max_element(found_niches.begin(), found_niches.end());
    auto min_niche = *std::min_element(found_niches.begin(), found_niches.end());
    CHECK(max_niche == num_niches - 1);
    CHECK(min_niche == 0);
    
    std::array<int, 6> hist = {0, 0, 0, 0, 0, 0};
    for(const auto& i : found_niches) {
        hist[i]++;
    }
    
    for (size_t i = 0; i < 6; ++i) {
        CHECK(std::abs(hist[i] - 100) < 50); // this is very generous! Should never break.
    }
}

TEST_CASE ("migration") {
    Param P;
    std::vector<double> trait_goals = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double sigma = 1.0;
    rnd_j rndgen(1);
    
    Individual Offspring(trait_goals, sigma, P, female, rndgen);
    
    for (double proportion_dispersal_density = 0.0; proportion_dispersal_density <= 1.0;
         proportion_dispersal_density += 0.01) {
    
        double proportion_available = 0.5;
        Offspring.fit_to_niche = 0.9;
        double OldNicheFit = Offspring.fit_to_niche;
        double basal_dispersal = 0.01;
        
        double mismatch_dispersal = std::exp(-5 * pow(OldNicheFit, 2));
        //old function:
         double density_dispersal = std::exp(0.5 * pow(proportion_available, 2)) - 1;
         double prob_disperse = std::max(basal_dispersal,
                                         ((1 - proportion_dispersal_density) * mismatch_dispersal
                                        + proportion_dispersal_density * density_dispersal));
        
        double new_prob = Offspring.calc_migration_prob(proportion_available, proportion_dispersal_density, basal_dispersal);
        
        REQUIRE(std::abs(new_prob - prob_disperse) < 1e-5);
    }
}
