//
//  main.cpp
//  Jana
//
//  Created by thijsjanzen on 01/03/2022.
//

#include <stdio.h>
#include <iostream>

#include <fstream>

#include <vector>
#include "rndutils.hpp"

using reng_type = rndutils::default_engine;
enum sex {male, female};

enum output_type {only_average, only_females};

struct Param {
    int number_of_timesteps = 20000;
    int size_per_niche = 1000;
    
    int num_niches = 6;
    int num_traits = 6;
    
    int pop_size_max = 1000;
    double lambda = 0.5;
    
    double basal_death_rate = 0.05;
    double basal_birth_rate = 0.8;
    double basal_migration_rate = 0.01;
    
    double sigma = 0.01;
    double d_sigma = 0.001;
    
    double mu = 0.01;
    double recom_rate = 0.01;
    
    size_t seed = 42;
    
    bool use_random_niches = false;
    
    int init_males = 100;
    int init_females = 100;
    
    output_type chosen_output_type = only_average;
    std::string only_average_file_name = "averages.txt";
    
};

class rnd_j {
public:
    reng_type rndgen;
    
    rnd_j() {
        std::random_device rd;
        reng_type rndgen_t(rd());
        rndgen = rndgen_t;
    }
    
    rnd_j(size_t seed) {
        rndgen = reng_type(seed);
    }
    
    rnd_j(const Param& P) {
        rndgen = reng_type(P.seed);
        set_mutate_trait_dist(P.sigma);
        set_dewlap_dist(P.d_sigma);
        set_mutate_prob(P.mu);
    }
    
    // true or false with 50/50 probability:
    bool flip_coin() {
        return coin_flip(rndgen);
    }
    
    bool bernouilli(double p) {
        return std::bernoulli_distribution(p)(rndgen);
    }
    
    sex get_random_sex() {
        if (flip_coin()) {
            return female;
        }
        return male;
    }
    
    double mutate_trait(double old_trait_value) {
        auto new_trait_value = old_trait_value;
        if (bernouilli(mutate_prob)) {
            new_trait_value += cauchy_dist(rndgen);
            new_trait_value = std::max(new_trait_value, 0.0);
            new_trait_value = std::min(new_trait_value, 1.0);
        }
        return new_trait_value;
    }
    
    // picks a random number in [0, n-1], useful when picking randomly from a vector.
    size_t random_number(int n)    {
        return std::uniform_int_distribution<> (0, n - 1)(rndgen);
    }

    double dewlap_noise() {
        return dewlap_dist(rndgen);
    }

    size_t draw_random_niche(size_t current_niche, size_t num_niches) {
        std::uniform_int_distribution<> niche_number(0, num_niches - 1);
        size_t new_niche = niche_number(rndgen);
        while(new_niche == current_niche) {
            new_niche = niche_number(rndgen);
        }
        return new_niche;
    }
    
    
    //// setters:

    void set_mutate_trait_dist(double s) {
        cauchy_dist = std::cauchy_distribution<double>(0.0 ,s);
    }
    
    void set_dewlap_dist(double s) {
        dewlap_dist = std::normal_distribution<double>(0.0, s);
    }
    
    void set_mutate_prob(double m) {
        mutate_prob = m;
    }
    
private:
    std::cauchy_distribution<double> cauchy_dist;
    std::bernoulli_distribution coin_flip = std::bernoulli_distribution(0.5); // we can keep this fixed.
    std::normal_distribution<double> dewlap_dist;
    double mutate_prob;
};






struct Trait {
    double a;
    double b;
    double c;
    double phenotype;
    
    Trait(double A, double B, double C) : a(A), b(B), c(C) {
    }

    void set_phenotype(sex S, double sigma) {
        if (S == female) {
            phenotype = sigma * a + (1 - sigma) * b;
        }
        if (S == male) {
            phenotype = sigma * a + (1 - sigma) * c;
        }
    }

    void mutate(rnd_j& rnd) {
        a = rnd.mutate_trait(a);
        b = rnd.mutate_trait(b);
        c = rnd.mutate_trait(c);
    }
};



struct Individual {
    std::vector< Trait > traits;
    const sex S;
    double resource_level;
    double mismatch;
    double fit_to_niche;
    double carotenoid_investment;
    double dewlap;
    
    Individual(const std::vector< int >& initial_trait_values,
               double sigma,
               sex initial_sex) : S(initial_sex) {
        for (int i = 0; i < initial_trait_values.size(); ++i) {
            traits.push_back(Trait(initial_trait_values[i],
                                   initial_trait_values[i],
                                   initial_trait_values[i]));
        }
    }
    
    Individual(const Individual& parent1, const Individual& parent2,
               double sigma, double recom_rate, rnd_j& rnd) : S(rnd.get_random_sex()) {

        // recombination
        for (size_t i = 0; i < parent1.traits.size(); ++i) {
            if (!rnd.bernouilli(recom_rate)) {
                traits.push_back(parent1.traits[i]);
            } else {
                traits.push_back(parent2.traits[i]);
            }
        }
        for (auto& i : traits) {
            i.mutate(rnd);
            i.set_phenotype(S, sigma);
        }
    }
    
    double calculate_match_to_niche(const std::vector<int>& selection_goals) {
        double fit = 0.0;
        for (size_t i = 0; i < selection_goals.size(); ++i) {
            auto d = selection_goals[i] - traits[i].phenotype;
            fit += std::abs(d);
        }
        return fit;
    }
    
    double allocate_resources(double niche_fit, rnd_j& rnd) {
        dewlap = carotenoid_investment * niche_fit + rnd.dewlap_noise();
        return niche_fit * (1 - carotenoid_investment);
    }
    
    void calculate_resources(const std::vector<int>& selection_goals,
                             int max_mismatch,
                             rnd_j rnd) {
        mismatch = calculate_match_to_niche(selection_goals);
        fit_to_niche = (max_mismatch - mismatch) / max_mismatch;
        if (S == female) {
            resource_level = fit_to_niche;
        }
        if (S == male) {
            resource_level = allocate_resources(fit_to_niche,
                                                rnd);
        }
    }
    
    bool will_migrate(double p, double lambda, double min_rate, rnd_j& rnd) {
        double disp_fit = std::exp(-5 * fit_to_niche * fit_to_niche);
        double disp_dens = std::exp(0.5 * p * p) - 1;
        double prob_disp = lambda * disp_fit + (1 - lambda) * disp_dens;
        prob_disp = std::max(prob_disp, min_rate);
        return rnd.bernouilli(prob_disp);
    }
    
    std::vector<double> collect_information() const {
        std::vector<double> info = {static_cast<double>(S),
                                    resource_level,
                                    fit_to_niche,
                                    carotenoid_investment,
                                    dewlap};
        for (auto i : traits) {
            info.push_back(i.a);
            info.push_back(i.b);
            info.push_back(i.c);
            info.push_back(i.phenotype);
        }
        return info;
    }
}; 



struct Niche {
    
    
    Niche(const std::vector< int >& goals_from_data,
          int num_males,
          int num_females,
          Param P) : selection_goals(goals_from_data),
                          num_traits(goals_from_data.size()),
    death_rate(P.basal_death_rate) {
       
        for (int i = 0; i < num_males; ++i) {
            males.push_back(Individual(goals_from_data, P.sigma, male));
        }
        for (int i = 0; i < num_females; ++i) {
            females.push_back(Individual(goals_from_data, P.sigma, male));
        }
    }
    
    std::vector< Individual > males;
    std::vector< Individual > females;
    
    std::vector< Individual > migrants; // both males and females, dispersing AWAY.
    
    
    void survival(std::vector< Individual>& v, rnd_j& rnd) {
        std::vector< Individual > alive;
        for (const auto& i : v) {
            double indiv_deathrate = 0.75 * death_rate +
                                     0.25 * std::exp(-2 * pow(i.resource_level, 8));
            if (!rnd.bernouilli(indiv_deathrate)) {
                alive.push_back(i);
            }
        }
        std::swap(v, alive);
    }
    
    
    void viability_selection(rnd_j& rnd) {
        survival(males, rnd);
        survival(females, rnd);
    }
    
    void reproduction(rnd_j& rnd, const Param& P) {
        std::vector< Individual > kids;
        
        if (males.empty()) return; // no reproduction
        
        double p = 1.0 * (P.pop_size_max - males.size() + females.size()) / P.pop_size_max;
        if (p < 0.0) p = 0.0;
        for (auto& mother : females) {
            double prob_repro = p * (0.8 * P.basal_birth_rate + 0.2 * mother.resource_level);
            if (rnd.bernouilli(prob_repro)) {
                // I first code here random mating. Non-random mating we need to look at a bit more closely!
                auto father = males[ rnd.random_number(males.size()) ];
                
                auto offspring = Individual(mother, father, P.sigma, P.recom_rate, rnd);
                
                if (offspring.will_migrate(p, P.lambda, P.basal_migration_rate, rnd)) {
                    migrants.push_back(offspring);
                } else {
                    kids.push_back(offspring);
                }
            }
        }
        
        for (auto& i : kids) {
            if (i.S == female) females.push_back(i);
            if (i.S == male)   males.push_back(i);
        }
        kids.clear(); // just to be sure.
    }
    
    void add_individual(const Individual& new_individual) {
        if (new_individual.S == female) {
            females.push_back(new_individual);
        } else {
            males.push_back(new_individual);
        }
    }
    
    const std::vector< int > selection_goals;
    const int num_traits;
    const double death_rate; // in case niches differ in basal death rate
};



struct Output {
    Output(const Param& P) {
        o = P.chosen_output_type;
        file_name = get_file_name(P, o);
    }
    
    Output() {
        
    }
    
    std::string file_name;
    output_type o;
    
    void update(const std::vector< Niche >& world, int t) {
        switch(o) {
            case only_average:
                output_averages(world, t);
                break;
            case only_females:
                // output_females(world, t);
                break;
        }
    }
    
    std::string make_file_name(std::string base, const Param& P) {
        base += "_" + std::to_string(P.recom_rate);
        base += "_" + std::to_string(P.basal_birth_rate);
        return base;
    }
    
    std::string get_file_name(const Param& P, output_type o) {
        switch (o) {
            case only_average:
                return make_file_name(P.only_average_file_name, P);
                break;
            case default:
                return "test.txt";
                break;
        }
    }
    
    void output_averages(const std::vector< Niche >& world, int t) {
        std::ofstream out_file(file_name.c_str(), std::ios::app);
        out_file << t << "\t";
        for (size_t i = 0; i < world.size(); ++i) {
            std::vector< std::vector< double > > individual_info;
            for (const auto& j : world[i].males) {
                individual_info.push_back(j.collect_information());
            }
            for (const auto& j : world[i].females) {
                individual_info.push_back(j.collect_information());
            }
            std::vector< double > mean_values = get_mean_values(individual_info);
            std::vector< double > sd_values   = get_sd_values(individual_info, mean_values);
            out_file << i << "\t"; // instead of i, trait name could also work
            // readr::read_tsv();
            for (size_t j = 0; j < mean_values.size(); ++j) {
                out_file << mean_values[i] << "\t" << sd_values[i] << "\t";
            }
            out_file << "\n";
        }
        out_file.close();
    }
    
    std::vector<double> get_mean_values(const std::vector< std::vector< double >>& v) {
        std::vector<double> m(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
          /*  double s = std::accumulate(v.begin(), v.end(), 0.0,
                                       [&](const std::vector<double>& a,
                                           const std::vector<double>& b) {
                                           return a[i] + b[i];
                                        });*/
            
            ////              a, b, c  = i
            ///indiv 1         b
            ///indiv 2         b
            ///indiv 3         b
            ///indiv 4         b
            /// = j
            
            double s = 0.0;
            /*
            for (size_t j = 0; j < v.size(); ++j) {
                s += v[j][i];
            }
            
            for (std::vector<std::vector<double>>::const_iterator j = v.begin();
                j != v.end(); ++j) {
                s += (*j)[i];
            }*/
            
            for (auto j = v.begin(); j != v.end(); ++j) {
                s += (*j)[i];
            }
            /*
            for (const auto& j : v) {
                s += j[i];
            }*/
            m[i] = s * 1.0 / v[i].size();
        }
        return m;
    }
    std::vector<double> get_sd_values(const std::vector< std::vector< double >>& v,
                                      const std::vector<double> m) {
        std::vector<double> sd(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            double s = 0.0;
            for (auto j = v.begin(); j != v.end(); ++j) {
                s += ((*j)[i] - m[i]) * ((*j)[i] - m[i]);
            }
            sd[i] = sqrt( s * 1.0 / v.size());
        }
        return sd;
    }
    
};

struct Simulation {
    Simulation() {
      // later on, we can have a look at code to read parameters from file.
      //  parameters = read_parameters();
        
        master_random_generator = rnd_j(parameters);
        record = Output(parameters);
       
        std::vector< std::vector< int > > goals;
        if (parameters.use_random_niches) {
            goals =  create_random_goals();
        } else {
            // code below not rwitten yet.
           // goals  = read_niches_from_file("niche_goals.txt");
        }
        
        for (int i = 0; i < parameters.num_niches; ++i) {
            int m = 0;
            int f = 0;
            if (i == 0) {
                m = parameters.init_males;
                f = parameters.init_females;
            }
            
            world.push_back( Niche(goals[i], m, f, parameters) );
        }
    }

    void run() {
        for (size_t t = 0; t < parameters.number_of_timesteps; ++t) {
            for (size_t i = 0; i < world.size(); ++i) {
                world[i].viability_selection(master_random_generator);
                world[i].reproduction(master_random_generator, parameters);
            }
            distribute_migrants();
            record.update(world, t);
        }
    }

    void distribute_migrants() {
        for (size_t current_niche = 0; current_niche < world.size(); ++current_niche) {
            
            for (const auto& i : world[current_niche].migrants) {
                size_t other_niche = master_random_generator.draw_random_niche(current_niche, world.size());
                world[other_niche].add_individual(i);
            }
        }
    }
    
    std::vector< std::vector< int >> create_random_goals() {
        std::vector< std::vector< int >> new_goals;
        for (int i = 0; i < parameters.num_niches; ++i) {
            std::vector< int > niche_goals(parameters.num_traits);
            for (int j = 0; j < parameters.num_traits; ++j) {
                niche_goals[j] = master_random_generator.random_number(100);
            }
            new_goals.push_back(niche_goals);
        }
        return new_goals;
    }
    
    Param parameters;
    rnd_j master_random_generator;
    Output record;
    
    std::vector< Niche > world;
};


int main() {
    Simulation sim;
    sim.run();
    return 0;
}
