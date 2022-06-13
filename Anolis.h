//
//  Anolis.h
//  Anolis_original
//
//  Created by thijsjanzen on 30/05/2022.
//

#ifndef Anolis_h
#define Anolis_h

#include <stdio.h>
#include <iostream>

#include <fstream>

#include <vector>
#include <string>
#include <array>

#include <mutex>

#include <chrono>

#include "rndutils.hpp"

enum sex { male, female };
enum output_type { only_average, only_females, indiv_data_end, only_dispersers, extinction_metrics, cInvest };

using reng_type = rndutils::default_engine; // alternative is mersenne twister

struct Param {
    size_t number_of_timesteps = 1001;
    int save_interval = 100; //e.g. save output only every 10 timesteps

    size_t num_niches = 6;
    size_t num_traits = 6;
    int initial_niche = 5;//which niche do you start in?

    int pop_size_max = 1000;
    double lambda = 0.0; // propensity dispersal

    double basal_death_rate = 0.05; // deathrate
    double basal_birth_rate = 0.8; // birthrate
    double basal_migration_rate = 0.01; // basal_dispersal

    double sigma = 0.0; // should be  0, 0.5 or 1. // PropSexDep

    double gamma = 3.0; // standard deviation of cauchy distribution of mutation.
                        // also known as muteffect_size

    double d_sigma = 0.0001;

    double mu = 0.01; //mutrate
    double recom_rate = 0.0;

    int start_seed = 0;
    int end_seed = 100;
    int used_seed;

    bool use_random_niches = false;

    int init_males = 500;
    int init_females = 500;
    double init_investment = 0.1;

    int SexSel = 1; //Between how many males can the female choose? The higher this variable is, the stronger sexual selection

    output_type chosen_output_type = cInvest;
    std::string only_average_file_name = "averages.txt";
    std::string extinction_metrics_file_name = "extinction_metrics";
    std::string cInvest_file_name = "cInvest";
    std::string final_indiv_file_name = "final_indiv_data.txt";
    std::string only_dispersers_file_name = "dispersers.txt";

    std::string niche_file_name = "niche_goals_largescale.txt";

    std::vector<double> MaxMism{ 410, 330, 420, 460, 460, 500 };//this is of course a very inelegant way of implementing the mismatch, but we are not keeping the mismatch anyways, just to check if we replicate the results

    void set_parameters(std::string parFileName) {

        std::ifstream ifs(parFileName.c_str());
        if (!ifs.is_open()) { std::cerr << "Unable to open parfile " << parFileName << '\n'; exit(EXIT_FAILURE); }
        for (;;) {
            std::string parId;
            ifs >> parId;
            if (ifs.good()) {
                if (parId == "maxtime") { ifs >> number_of_timesteps; }
                else if (parId == "saveinterval") { ifs >> save_interval; }
                else if (parId == "num_niches") { ifs >> num_niches; }
                else if (parId == "num_traits") { ifs >> num_traits; }
                else if (parId == "start_niche_ID") { ifs >> initial_niche; }
                else if (parId == "popsize_per_niche") { ifs >> pop_size_max; }
                else if (parId == "lambda") { ifs >> lambda; }
                else if (parId == "basal_d") { ifs >> basal_death_rate; }
                else if (parId == "basal_b") { ifs >> basal_birth_rate; }
                else if (parId == "basal_mig") { ifs >> basal_migration_rate; }
                else if (parId == "sigma") { ifs >> sigma; }
                else if (parId == "d_sigma") { ifs >> d_sigma; }
                else if (parId == "mutrate") { ifs >> mu; }
                else if (parId == "recombirate") { ifs >> recom_rate; }
                else if (parId == "use_random_niches") { ifs >> use_random_niches; }
                else if (parId == "init_males") { ifs >> init_males; }
                else if (parId == "init_females") { ifs >> init_females; }
                else if (parId == "init_investment") { ifs >> init_investment; }
                else if (parId == "SexSel") { ifs >> SexSel; }


                else { std::cerr << "unknown parname in file"; exit(EXIT_FAILURE); }

            }
            else break;
        }
        ifs.close();
    }
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
        rndgen = reng_type(P.used_seed);
        set_mutate_trait_dist(P.gamma);
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

    double mutate_trait(double old_trait_value, double max_trait) {
        auto new_trait_value = old_trait_value;
        if (bernouilli(mutate_prob)) {
            double trait_change =  cauchy_dist(rndgen);
            if (max_trait == 1.0) trait_change = trait_change / 100.0;
            
            new_trait_value = std::clamp(old_trait_value + trait_change, 0.0, max_trait);
        }
        return new_trait_value;
    }

    // picks a random number in [0, n-1], useful when picking randomly from a vector.
    size_t random_number(size_t n) {
        if (n == 1) return 0;
        return std::uniform_int_distribution<>(0, static_cast<int>(n) - 1)(rndgen);
    }

    double dewlap_noise() {
        return dewlap_dist(rndgen);
    }

    size_t draw_random_niche(size_t current_niche, size_t num_niches) {
        std::uniform_int_distribution<> niche_number(0, static_cast<int>(num_niches) - 1);
        size_t new_niche = niche_number(rndgen);
        while (new_niche == current_niche) {
            new_niche = niche_number(rndgen);
        }
        return new_niche;
    }

    double uniform() {
        return unif_dist(rndgen);
    }

    //// setters:
    void set_mutate_trait_dist(double s) {
        cauchy_dist = std::cauchy_distribution<double>(0.0, s);
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
    std::uniform_real_distribution<> unif_dist =
        std::uniform_real_distribution<>(0, 1.0);

    double mutate_prob;
};


struct Trait {
    double a;
    double b;
    double c;
    double phenotype;

    Trait() {};
    
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
        a = rnd.mutate_trait(a, 100.0);
        b = rnd.mutate_trait(b, 100.0);
        c = rnd.mutate_trait(c, 100.0);
    }
};

struct store_info { // used for sexual selection
    size_t index;
    double dew_lap;
};

struct Individual {
    std::vector< Trait > traits;
    sex S;
    double resource_level;
    double mismatch;
    double fit_to_niche;
    double carotenoid_investment;
    double dewlap;
    bool dispersed = false; // I added this in order to be able to output e.g. specifically only dispersers

    int prev_niche;
    int niche;
    double prev_mismatch; // mismatch before dispersal, i.e. had they not dispersed

    int age;
    int LRS;

    Individual(const std::vector< double >& trait_goals,
               double sigma,
               const Param& P,
               sex initial_sex,
               rnd_j& rnd) : S(initial_sex), niche(P.initial_niche) {
        
        std::uniform_real_distribution<double> initializeWithNoise(-10, 10);
        
        traits.resize(trait_goals.size());
        
        for (size_t i = 0; i < traits.size(); ++i) {
            
            std::vector<double> new_traits {trait_goals[i], trait_goals[i], trait_goals[i]};
            for (size_t j = 0; j < new_traits.size(); ++j) {
                new_traits[j] += initializeWithNoise(rnd.rndgen);
                new_traits[j] = std::clamp(new_traits[j], 0.0, 100.0);
            }
            
            traits[i] = Trait(new_traits[0], new_traits[1], new_traits[2]);
            
            traits[i].set_phenotype(S, sigma);
        }
        niche = P.initial_niche;
        carotenoid_investment = P.init_investment;
        calculate_resources(trait_goals, P.MaxMism[P.initial_niche], rnd);
        age = 0;
        LRS = 0;
    }

    Individual(const Individual& parent1,
               const Individual& parent2,
               const std::vector<double>& trait_goals,
               Param P,
               rnd_j& rnd) : S(rnd.get_random_sex()), niche(parent1.niche) {

        // recombination
        traits = parent1.traits;
        for (size_t i = 0; i < parent1.traits.size(); ++i) {
            if (rnd.bernouilli(P.recom_rate)) {
                traits[i] = parent2.traits[i];
            }
        }

        carotenoid_investment = parent1.carotenoid_investment;
        if (rnd.bernouilli(P.recom_rate)) {
            carotenoid_investment = parent2.carotenoid_investment;
        }

        //mutation&setting phenotype
        for (auto& i : traits) {
            i.mutate(rnd);
            i.set_phenotype(S, P.sigma);
        }
        
        carotenoid_investment = rnd.mutate_trait(carotenoid_investment, 1.0); // is this correct? TJ: yes

        niche = parent1.niche;
        
        calculate_resources(trait_goals, P.MaxMism[niche], rnd);
        age = 0;
        LRS = 0;
    }

    double calculate_match_to_niche(const std::vector<double>& selection_goals) {
            double current_mismatch = 0.0;
           
            // legacy code!!!!!! only for initial testing.
            switch( niche ) {
                case 0: {
                    current_mismatch = std::abs(traits[0].phenotype - selection_goals[0]) +
                                       0 +
                                       std::abs(traits[2].phenotype - selection_goals[2]) +
                                       std::abs(traits[3].phenotype - selection_goals[3]) +
                                       std::abs(traits[4].phenotype - selection_goals[4]) +
                                       std::abs(traits[5].phenotype - selection_goals[5]);
                    break;
                }
                case 1: {
                    current_mismatch = std::abs(traits[0].phenotype - selection_goals[0]) +
                                       0 +
                                       std::abs(traits[2].phenotype - selection_goals[2]) +
                                       std::abs(traits[3].phenotype - selection_goals[3]) +
                                       std::abs(traits[4].phenotype - selection_goals[4]) +
                                       std::abs(traits[5].phenotype - selection_goals[5]);
                    break;
                }
                case 2: {
                    current_mismatch = std::abs(traits[0].phenotype - selection_goals[0]) +
                                       std::abs(traits[1].phenotype - selection_goals[1])+
                                       std::abs(traits[2].phenotype - selection_goals[2]) +
                                       std::abs(traits[3].phenotype - selection_goals[3]) +
                                       std::abs(traits[4].phenotype - selection_goals[4]) +
                                       std::abs(traits[5].phenotype - selection_goals[5]);
                    break;
                }
                case 3: {
                    double size_mismatch;
                    if (traits[0].phenotype < 10) { size_mismatch = 10 - traits[0].phenotype; }
                    if (traits[0].phenotype > 50) { size_mismatch = traits[0].phenotype - 50; }
                    else { size_mismatch = 0; }
                    
                    current_mismatch =  size_mismatch +
                                        std::abs(traits[1].phenotype - selection_goals[1])+
                                        std::abs(traits[2].phenotype - selection_goals[2]) +
                                        std::abs(traits[3].phenotype - selection_goals[3]) +
                                        std::abs(traits[4].phenotype - selection_goals[4]) +
                                        std::abs(traits[5].phenotype - selection_goals[5]);
                    break;
                }
                case 4: {
                    current_mismatch = std::abs(traits[0].phenotype - selection_goals[0]) +
                                       std::abs(traits[1].phenotype - selection_goals[1])+
                                       std::abs(traits[2].phenotype - selection_goals[2]) +
                                       std::abs(traits[3].phenotype - selection_goals[3]) +
                                       std::abs(traits[4].phenotype - selection_goals[4]) +
                                       std::abs(traits[5].phenotype - selection_goals[5]);
                    break;
                }
                case 5: {
                    double size_mismatch;
                    if (traits[0].phenotype < 10) { size_mismatch = 10 - traits[0].phenotype; }
                    if (traits[0].phenotype > 50) { size_mismatch = traits[0].phenotype - 50; }
                    else { size_mismatch = 0; }
                    
                    current_mismatch =  size_mismatch +
                                        std::abs(traits[1].phenotype - selection_goals[1])+
                                        std::abs(traits[2].phenotype - selection_goals[2]) +
                                        std::abs(traits[3].phenotype - selection_goals[3]) +
                                        std::abs(traits[4].phenotype - selection_goals[4]) +
                                        std::abs(traits[5].phenotype - selection_goals[5]);
                    break;
                }
                default: {
                    throw std::runtime_error("Error! Individual in non-existing niche!");
                }
            }
            
            return current_mismatch;
        }

    double allocate_resources(double niche_fit, rnd_j& rnd) {
        dewlap = carotenoid_investment * niche_fit + rnd.dewlap_noise();
        return niche_fit * (1 - carotenoid_investment);
    }

    void calculate_resources(const std::vector<double>& selection_goals,
                                 int max_mismatch,
                                 rnd_j& rnd) {
            mismatch = calculate_match_to_niche(selection_goals);
            fit_to_niche = (max_mismatch - mismatch) / max_mismatch;
            if (S == female) {
                resource_level = fit_to_niche;
                dewlap = 0.0; // not really used, but just to be sure.
            }
            if (S == male) {
                dewlap = carotenoid_investment * fit_to_niche + rnd.dewlap_noise();
                resource_level = resource_level - carotenoid_investment * resource_level;
            }
            resource_level = std::clamp(resource_level, 0.0, 1.0);
        }

    double calc_migration_prob(double p, double lambda, double min_rate) {
            double disp_fit = std::exp(-5 * fit_to_niche * fit_to_niche);
            double disp_dens = std::exp(0.5 * p * p) - 1;
            double prob_disp = (1 - lambda) * disp_fit + lambda * disp_dens;
            prob_disp = std::max(prob_disp, min_rate);
            return prob_disp;
    }
    
    bool will_migrate(double p, double lambda, double min_rate, rnd_j& rnd) {
        auto prob_disp = calc_migration_prob(p, lambda, min_rate);
        return rnd.bernouilli(prob_disp);
    }

    std::vector<double> collect_information() const {
        std::vector<double> info = { static_cast<double>(S),
                                    resource_level,
                                    fit_to_niche,
                                    carotenoid_investment,
                                    dewlap };
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

    Niche(const std::vector< double >& goals_from_data,
        int num_males,
        int num_females,
        Param P,
        rnd_j& rnd) : selection_goals(goals_from_data),
        num_traits(goals_from_data.size()),
        death_rate(P.basal_death_rate) {

        for (int i = 0; i < num_males; ++i) {
            males.push_back(Individual(goals_from_data, P.sigma, P, male, rnd));
        }
        for (int i = 0; i < num_females; ++i) {
            females.push_back(Individual(goals_from_data, P.sigma, P, female, rnd));
        }
    }

    std::vector< Individual > males;
    std::vector< Individual > females;

    std::vector< Individual > migrants; // both males and females, dispersing AWAY.


    void survival(std::vector< Individual>& v, rnd_j& rnd) {
        if (v.empty()) return;


        // if order of individuals is important:

      /*
        std::vector< Individual > alive;
        for (const auto& i : v) {
            double indiv_deathrate = 0.75 * death_rate +
                                     0.25 * std::exp(-2 * pow(i.resource_level, 8));
            if (!rnd.bernouilli(indiv_deathrate)) {
                alive.push_back(i);
            }
        }
        std::swap(v, alive);
      */
      // else (~2x faster!) :

        for (size_t i = 0; i < v.size(); ) {
            double indiv_deathrate = 0.75 * death_rate +
                0.25 * std::exp(-2 * pow(v[i].fit_to_niche, 8)); // THIS IS WRONG FOR MALES
            if (rnd.bernouilli(indiv_deathrate)) {
                v[i] = v.back();
                v.pop_back();
            }
            else {
                v[i].age = v[i].age + 1;
                ++i;
            }
        }
    }


    void viability_selection(rnd_j& rnd,
                             int& num_dead_females,
                             int& num_dead_males) {
        num_dead_males = 0;
        num_dead_females = 0;
        
        std::vector<Individual> tmp_pop;
        for (const auto& ind : females) {
            double indiv_deathrate = 0.75 * death_rate +
                                     0.25 * std::exp(-2 * std::pow(ind.fit_to_niche, 8));

            std::bernoulli_distribution DeathEvent(indiv_deathrate);
            if (!DeathEvent(rnd.rndgen)) {
                tmp_pop.emplace_back(ind);
            } else {
                num_dead_females++;
            }
        }
        std::swap(females, tmp_pop);

        tmp_pop.clear();
        for (const auto& ind : males) {
            double indiv_deathrate = 0.75 * death_rate +
            0.25 * std::exp(-2 * pow(ind.resource_level, 8));
            std::bernoulli_distribution DeathEvent(indiv_deathrate);
            if (!DeathEvent(rnd.rndgen)) {
                tmp_pop.emplace_back(ind);
            } else {
                num_dead_males++;
            }
        }
        std::swap(males, tmp_pop);
    }


    bool is_in_candidates(const std::vector<store_info>& v, size_t index) {
        for (const auto& i : v) {
            if (i.index == index) return true;
        }
        return false;
    }

    size_t get_father_index(const std::vector< Individual>& males,
        int sexsel,
        rnd_j& rnd) {
        if (sexsel == 1) {
            return static_cast<int>(rnd.random_number(males.size()));
        }

        // if the number of males to pick from is larger than the
        // number of available males, limit choice to all available males.
        if (sexsel > static_cast<int>(males.size())) sexsel = static_cast<int>(males.size());

        std::vector< store_info > candidates(sexsel);

        if (sexsel >= static_cast<int>(males.size())) { // border case where there are few males
            for (size_t i = 0; i < males.size(); ++i) {
                candidates[i].index = i;
                candidates[i].dew_lap = males[i].dewlap;
            }
        }
        else {
            for (int i = 0; i < sexsel; ++i) {
                size_t index = rnd.random_number(males.size());
                while (is_in_candidates(candidates, index)) {
                    index = rnd.random_number(males.size());
                }

                candidates[i].index = static_cast<int>(index);
                candidates[i].dew_lap = males[index].dewlap;
            }
        }

        double s = 0.0;
        // For later (if Franjo really really really really wants this):
        // if we want to implement a soft max function
        // we need to iterate twice over the dewlaps to correct them:
        // for (auto& i : candidates) {
        //   i.dewlap = exp(i.dewlap * b); s += i.dewlap;
        // }
        // for (auto& i : candidates) {
        //   i.dewlap *= 1.0 / s;
        // }
        // s = 1.0;

        for (const auto& i : candidates) {
            s += i.dew_lap;
        }

        double r = s * rnd.uniform();
        size_t picked_indiv = 0;
        for (; picked_indiv < candidates.size(); ++picked_indiv) {
            r -= candidates[picked_indiv].dew_lap;
            if (r <= 0.0) break;
        }

        return picked_indiv;
    }

    void reproduction(rnd_j& rnd, const Param& P,
                      int& num_kids,
                      int& num_migrants_sent) {
        
        num_kids = 0;
        num_migrants_sent = 0;
        
        std::vector< Individual > kids;

        if (males.empty()) return; // no reproduction

        int current_pop_size = static_cast<int>(males.size()) + static_cast<int>(females.size());

        double p = 1.0 * (P.pop_size_max - current_pop_size) / P.pop_size_max;
        if (p < 0.0) p = 0.0;
        
        for (auto& mother : females) {
            double prob_repro = p * (0.8 * P.basal_birth_rate + 0.2 * mother.resource_level);
            if (rnd.bernouilli(prob_repro)) {

                auto father_index = get_father_index(males, P.SexSel, rnd);

                auto offspring = Individual(mother, males[father_index],                               selection_goals, P, rnd);

                mother.LRS = mother.LRS + 1;
                males[father_index].LRS = males[father_index].LRS + 1;

                if (offspring.will_migrate(p, P.lambda, P.basal_migration_rate, rnd)) {
                    offspring.dispersed = true;
                    offspring.prev_niche = mother.niche;
                    offspring.prev_mismatch = offspring.mismatch;
                    migrants.push_back(std::move(offspring));
                    num_migrants_sent++;
                } else {
                    kids.push_back(std::move(offspring));
                    num_kids++;
                }
            }
        }

        for (auto& i : kids) {
            if (i.S == female) females.push_back(std::move(i));
            if (i.S == male)   males.push_back(std::move(i));
        }
    }

    void add_individual(const Individual& new_individual) {
        if (new_individual.S == female) {
            females.push_back((new_individual));
        }
        else {
            males.push_back((new_individual));
        }
    }

    const std::vector< double > selection_goals;
    const size_t num_traits;
    const double death_rate; // in case niches differ in basal death rate
};


struct Output {
    Output(const Param& P) {
        o = P.chosen_output_type;
        file_name = get_file_name(P, o);
    }

    Output() {

    }

    std::vector< std::array< std::string, 4 > > track_record;
    
    std::string file_name;
    output_type o;
    
    

    void update(const std::vector< Niche >& world, size_t t, const Param& P) {
        if (t % P.save_interval != 0) return;

        switch (o) {
        case only_average:
            output_averages(world, t, P);
            break;
        case only_females:
            // output_females(world, t);
            break;
        case indiv_data_end:
            output_indivData_end(world, t, P);
        case only_dispersers:
            output_dispersers(world, t, P);
            break;
        case extinction_metrics:
            output_extinction_metrics(world, t, P);
            break;
        case cInvest:
            output_cInvest(world, t, P);
            break;
        }
    }



    std::string make_file_name(std::string base, const Param& P) {
        base += "_Sigma_" + std::to_string(P.sigma);
        base += "_RecRate_" + std::to_string(P.recom_rate);
        base += "_DewNoise_" + std::to_string(P.d_sigma);
        base += "_SexSel_" + std::to_string(P.SexSel);
        base += "_Seed_" + std::to_string(P.used_seed);
        base += ".csv";

        return base;
    }

    std::string get_file_name(const Param& P, output_type o) {
        switch (o) {
        case only_average:
            return make_file_name(P.only_average_file_name, P);
            break;
        case indiv_data_end:
            return make_file_name(P.final_indiv_file_name, P);
            break;
        case only_dispersers:
            return make_file_name(P.only_dispersers_file_name, P);
        case extinction_metrics:
            return make_file_name(P.extinction_metrics_file_name, P);
            break;
        case cInvest:
            return make_file_name(P.cInvest_file_name, P);
            break;
        default:
            return "test.txt";
            break;
        }
    }

    void output_averages(const std::vector< Niche >& world, size_t t,
        const Param& P) {
        static std::once_flag header_written;
       
        //I tried to add column names here but not sure if it worked
        std::call_once(header_written, [&]() {
            std::ofstream out_file(file_name.c_str());
            
            out_file << "Time" << "\t" << "numMales" << "\t" << "numFemales" << "\t" << "Niche" << "\t";
            for (size_t TraitNr = 0; TraitNr < P.num_traits; TraitNr++) {
                out_file << "Trait_" << TraitNr << "\t" << "T" << TraitNr << "_Avg" << "\t" << "T" << TraitNr << "_Stdev" << "\t";
            }
            out_file << "\n";
            out_file.close();
        });
        
        std::ofstream out_file(file_name.c_str(), std::ios::app);
    

        for (size_t i = 0; i < world.size(); ++i) {
            std::vector< std::vector< double > > individual_info;
            individual_info.reserve(world[i].males.size() + world[i].females.size());
            for (const auto& j : world[i].males) {
                individual_info.push_back(j.collect_information());
            }

            for (const auto& j : world[i].females) {
                individual_info.push_back(j.collect_information());
            }

            out_file << t << "\t";
            out_file << world[i].males.size() << "\t" << world[i].females.size() << "\t";
            out_file << i << "\t"; // instead of i, trait name could also work

            if (individual_info.empty()) {
                size_t num_columns = (5 + P.num_traits * 4) * 2;
                for (size_t i = 0; i < num_columns; ++i) {
                    out_file << "NA" << "\t";
                }
                out_file << "\n";

            }
            else {
                std::vector< double > mean_values = get_mean_values(individual_info);
                std::vector< double > sd_values = get_sd_values(individual_info, mean_values);
                assert(mean_values.size() == sd_values.size());
                for (size_t j = 0; j < mean_values.size(); ++j) {
                    out_file << mean_values[i] << "\t" << sd_values[i] << "\t";
                }
                out_file << "\n";
            }
        }
        out_file.close();
    }


    void output_indivData_end(const std::vector< Niche >& world, size_t t,
        const Param& P) {
        static std::once_flag header_written;
       
        //I tried to add column names here but not sure if it worked
        std::call_once(header_written, [&]() {
            std::ofstream out_file(file_name.c_str());
        
            out_file << "Time" << "\t" << "Sex" << "\t" << "Resources" << "\t" << "Mismatch" << "\t" << "C_investment" << "\t" << "Dewlap" << "\t";
            for (size_t TraitNr = 0; TraitNr < P.num_traits; TraitNr++) {//Do I need to pass Param to the function if I use it here?
                out_file << TraitNr << "_A" << "\t";
                out_file << TraitNr << "_B" << "\t";
                out_file << TraitNr << "_C" << "\t";
                out_file << TraitNr << "_Phen" << "\t";
            }
            out_file << "\n";
            out_file.close();
        });
        
        std::ofstream out_file(file_name.c_str(), std::ios::app);
        
        if (static_cast<int>(t) == (P.number_of_timesteps - 1)) {

            for (size_t i = 0; i < world.size(); ++i) {
                //for (const auto& j : world[i].males) {
                for (size_t j = 0; j < world[i].males.size(); ++j) {
                    out_file << t << "\t" << world[i].males[j].S << "\t" << world[i].males[j].resource_level << "\t" << world[i].males[j].mismatch << "\t"
                        << world[i].males[j].carotenoid_investment << "\t" << world[i].males[j].dewlap << "\t";
                    for (size_t TraitNr = 0; TraitNr < P.num_traits; TraitNr++) {
                        out_file << world[i].males[j].traits[TraitNr].a << "\t";
                        out_file << world[i].males[j].traits[TraitNr].b << "\t";
                        out_file << world[i].males[j].traits[TraitNr].c << "\t";
                        out_file << world[i].males[j].traits[TraitNr].phenotype << "\t";
                    }
                    out_file << "\n";
                }
                for (size_t j = 0; j < world[i].females.size(); ++j) {
                    // for (const auto& j : world[i].females) {
                    out_file << t << "\t" << world[i].females[j].S << "\t" << world[i].females[j].resource_level << "\t" << world[i].females[j].mismatch << "\t"
                        << world[i].females[j].carotenoid_investment << "\t" << "NA" << "\t";
                    for (size_t TraitNr = 0; TraitNr < P.num_traits; TraitNr++) {
                        out_file << world[i].females[j].traits[TraitNr].a << "\t";
                        out_file << world[i].females[j].traits[TraitNr].b << "\t";
                        out_file << world[i].females[j].traits[TraitNr].c << "\t";
                        out_file << world[i].females[j].traits[TraitNr].phenotype << "\t";
                    }
                    out_file << "\n";
                }

            }
        }
        out_file.close();
    }

    void output_dispersers(const std::vector< Niche >& world, size_t t, const Param& P) {
        //I tried to add column names here but not sure if it worked
        static std::once_flag header_written;
        
        std::call_once(header_written, [&]() {
            std::ofstream out_file(file_name.c_str());
            
            out_file << "Time" << "\t" << "Sex" << "\t" << "Resources" << "\t" << "C_investment" << "\t" << "Dewlap" << "\t";
            for (size_t TraitNr = 0; TraitNr < P.num_traits; TraitNr++) {//Do I need to pass Param to the function if I use it here?
                out_file << TraitNr << "_A" << "\t";
                out_file << TraitNr << "_B" << "\t";
                out_file << TraitNr << "_C" << "\t";
                out_file << TraitNr << "_Phen" << "\t";
            }
            out_file << "NicheBefore" << "\t" << "NicheAfter" << "\t" << "MismatchBefore" << "\t" << "MismatchAfter" << "\t";
            out_file << "\n";
            out_file.close();
        });
                       
       std::ofstream out_file(file_name.c_str(), std::ios::app);
                       
        for (size_t i = 0; i < world.size(); ++i) {
            // for (const auto& j : world[i].males) {
            for (size_t j = 0; j < world[i].males.size(); ++j) {
                out_file << t << "\t" << world[i].males[j].S << "\t" << world[i].males[j].resource_level << "\t"
                    << world[i].males[j].carotenoid_investment << "\t" << world[i].males[j].dewlap << "\t";
                for (size_t TraitNr = 0; TraitNr < P.num_traits; TraitNr++) {
                    out_file << world[i].males[j].traits[TraitNr].a << "\t";
                    out_file << world[i].males[j].traits[TraitNr].b << "\t";
                    out_file << world[i].males[j].traits[TraitNr].c << "\t";
                    out_file << world[i].males[j].traits[TraitNr].phenotype << "\t";
                }
                out_file << world[i].males[j].prev_niche << "\t" << world[i].males[j].niche << "\t" << world[i].males[j].prev_mismatch << "\t" << world[i].males[j].mismatch << "\t";
                out_file << "\n";
            }
            // for (const auto& j : world[i].females) {
            for (size_t j = 0; j < world[i].females.size(); ++j) {
                out_file << t << "\t" << world[i].females[j].S << "\t" << world[i].females[j].resource_level << "\t"
                    << world[i].females[j].carotenoid_investment << "\t" << world[i].females[j].dewlap << "\t";
                for (size_t TraitNr = 0; TraitNr < P.num_traits; TraitNr++) {
                    out_file << world[i].females[j].traits[TraitNr].a << "\t";
                    out_file << world[i].females[j].traits[TraitNr].b << "\t";
                    out_file << world[i].females[j].traits[TraitNr].c << "\t";
                    out_file << world[i].females[j].traits[TraitNr].phenotype << "\t";
                }
                out_file << world[i].females[j].prev_niche << "\t" << world[i].females[j].niche << "\t" << world[i].females[j].prev_mismatch << "\t" << world[i].females[j].mismatch << "\t";
                out_file << "\n";
            }
        }
        out_file.close();
    }



    void output_extinction_metrics(const std::vector< Niche >& world,
                                   size_t t,
                                   const Param& P) {
        
        static std::once_flag header_written;
       
        //I tried to add column names here but not sure if it worked
        std::call_once(header_written, [&]() {
            std::ofstream out_file(file_name.c_str());
         
            out_file << "Time" << "\t";
            for (size_t NicheNr = 0; NicheNr < P.num_niches; NicheNr++) {
                out_file
                    << "Niche_ID" << "\t"
                    << "Niche_" << NicheNr << "_NrM" << "\t"
                    << "Niche_" << NicheNr << "_NrF" << "\t";
            }
            out_file << "\n";
            out_file.close();
        });
        
        std::ofstream out_file(file_name.c_str(), std::ios::app);

        out_file << t << "\t";
        for (size_t i = 0; i < world.size(); ++i) {
            out_file << i << "\t";
            out_file << world[i].males.size() << "\t" << world[i].females.size() << "\t";
        }
        out_file << "\n";

    }


    void output_cInvest(const std::vector< Niche >& world, size_t t,
        const Param& P) {
        static std::once_flag header_written;
       
        //I tried to add column names here but not sure if it worked
        std::call_once(header_written, [&]() {
            std::ofstream out_file(file_name.c_str());
         
            out_file << "Time" << "\t" << "Niche" << "\t" << "Sex" << "\t" << "cInvest" << "\t" << "Dispersed" << "\t" << "Age" << "\t" << "LRS" << "\t";
            out_file << "\n";
            out_file.close();
        });
        
        std::ofstream out_file(file_name.c_str(), std::ios::app);

        if (t % P.save_interval == 0) {
            for (size_t i = 0; i < world.size(); ++i) {
                for (size_t j = 0; j < world[i].males.size(); ++j) {
                    out_file << t << "\t" << i << "\t" << world[i].males[j].S << "\t" << world[i].males[j].carotenoid_investment << "\t"
                        << world[i].males[j].dispersed << "\t" << world[i].males[j].age << "\t" << world[i].males[j].LRS << "\t";
                    out_file << "\n";
                }
                for (size_t j = 0; j < world[i].females.size(); ++j) {
                    out_file << t << "\t" << i << "\t" << world[i].females[j].S << "\t" << world[i].females[j].carotenoid_investment << "\t"
                        << world[i].females[j].dispersed << "\t" << world[i].females[j].age << "\t" << world[i].females[j].LRS << "\t";
                    out_file << "\n";
                }
            }
        }
        out_file.close();
    }


    std::vector<double> get_mean_values(const std::vector< std::vector< double >>& v) {
        std::vector<double> m(v[0].size());
        for (size_t i = 0; i < v[0].size(); ++i) {
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
        std::vector<double> sd(v[0].size());
        for (size_t i = 0; i < v[0].size(); ++i) {
            double s = 0.0;
            for (auto j = v.begin(); j != v.end(); ++j) {
                s += ((*j)[i] - m[i]) * ((*j)[i] - m[i]);
            }
            sd[i] = sqrt(s * 1.0 / v.size());
        }
        return sd;
    }
    
    void update_dead_indivs(size_t t, int niche_no, int num_dead_females, int num_dead_males) {
        
        std::array< std::string, 4> add1 = {std::to_string(t), std::to_string(niche_no), "dead_females", std::to_string(num_dead_females)};
        std::array< std::string, 4> add2 = {std::to_string(t), std::to_string(niche_no), "dead_males", std::to_string(num_dead_males)};
        track_record.push_back(add1);
        track_record.push_back(add2);
    }
    
    void update_offspring(size_t t, int niche_no, int num_local_kids, int num_migrants_sent) {
        std::array< std::string, 4> add1 = {std::to_string(t), std::to_string(niche_no), "num_offspring", std::to_string(num_local_kids)};
        std::array< std::string, 4> add2 = {std::to_string(t), std::to_string(niche_no), "num_migrants_sent", std::to_string(num_migrants_sent)};
        track_record.push_back(add1);
        track_record.push_back(add2);
    }
    
    void write_track_record(std::string file_name) {
        std::ofstream out_file(file_name.c_str());
        out_file << "t" << "\t" << "niche" << "\t" << "statistic" << "\t" << "value" << "\n";
        
        for (const auto& i : track_record) {
            for (const auto& j : i ) {
                out_file << j << "\t";
            }
            out_file << "\n";
        }
        out_file.close();
    }
    
    std::vector< std::array< std::string, 4>> record_fit;
    
    void update_niche_fit(size_t t, const std::vector<Niche>& world) {
        for (size_t niche = 0; niche < world.size(); ++niche) {
            for (const auto& i : world[niche].females) {
                std::array<std::string, 4 > to_add = {std::to_string(t),
                                                      std::to_string(niche),
                                                      "female",
                                                     std::to_string(i.fit_to_niche)};
                record_fit.push_back(to_add);
            }
            
            for (const auto& i : world[niche].males) {
                std::array<std::string, 4 > to_add = {std::to_string(t),
                                                      std::to_string(niche),
                                                      "male",
                                                     std::to_string(i.fit_to_niche)};
                record_fit.push_back(to_add);
            }
        }
        return;
    }
    
    void write_niche_fit(const std::string file_name) {
        std::ofstream out_file(file_name.c_str());
        out_file << "t" << "\t" << "niche" << "\t" << "sex" << "\t" << "fit_to_niche" << "\n";
        
        for (const auto& i : record_fit) {
            for (const auto& j : i ) {
                out_file << j << "\t";
            }
            out_file << "\n";
        }
        out_file.close();
    }

};


struct Simulation {
    Simulation(const Param& p, int seed) : parameters(p) {

        parameters.used_seed = seed;
        std::cout << "Seed: " << seed << std::endl;
        master_random_generator = rnd_j(parameters);
        record = Output(parameters);

        std::vector< std::vector< double > > goals;
        if (parameters.use_random_niches) {
            goals = create_random_goals();
        }
        else {
            goals = read_niches_from_file();
        }

        for (int i = 0; i < parameters.num_niches; ++i) {
            int m = 0;
            int f = 0;
            if (i == parameters.initial_niche) {
                m = parameters.init_males;
                f = parameters.init_females;
            }

            world.push_back(Niche(goals[i], m, f, parameters, master_random_generator));
        }
    }

    void run() {
        auto clock_start = std::chrono::system_clock::now();
        
        int num_dead_females, num_dead_males, num_local_kids, num_migrants_sent;
        
        
        for (size_t t = 0; t < parameters.number_of_timesteps; ++t) {
            
            record.update_niche_fit(t, world);
            
            for (size_t i = 0; i < world.size(); ++i) {
                world[i].viability_selection(master_random_generator, num_dead_females, num_dead_males);
          //      record.update_dead_indivs(t, i, num_dead_females, num_dead_males);
                
                world[i].reproduction(master_random_generator, parameters, num_local_kids, num_migrants_sent);
         //       record.update_offspring(t, i, num_local_kids, num_migrants_sent);
            }
            distribute_migrants();
            record.update(world, t, parameters);


            if (t % parameters.save_interval == 0) {
                std::cout << "Time: " << t << " ";
                for (size_t i = 0; i < world.size(); ++i) {
                    std::cout << world[i].males.size() + world[i].females.size() << " ";
                }
                
                auto clock_now = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = clock_now - clock_start;
                std::cout << "this took: " << elapsed_seconds.count() << "seconds\n";
                clock_start = clock_now;
            }
        }
      //  record.write_track_record("track_record.txt");
      //  record.write_niche_fit("niche_fit.txt");
    }

    void distribute_migrants() {
        for (size_t current_niche = 0; current_niche < world.size(); ++current_niche) {

            for (auto& i : world[current_niche].migrants) {
                size_t other_niche = master_random_generator.draw_random_niche(current_niche, world.size());
                i.niche = static_cast<int>(other_niche);
                i.calculate_resources(world[other_niche].selection_goals,
                    parameters.MaxMism[other_niche],
                    master_random_generator); // update resources because of change of niche.

                world[other_niche].add_individual(i);
            }
            world[current_niche].migrants.clear();
        }
    }

    std::vector< std::vector< double >> create_random_goals() {
        std::vector< std::vector< double >> new_goals;
        for (size_t i = 0; i < parameters.num_niches; ++i) {
            std::vector< double > niche_goals(parameters.num_traits);
            for (size_t j = 0; j < parameters.num_traits; ++j) {
                niche_goals[j] = static_cast<double>(master_random_generator.random_number(100)); //why 100 shouldn't this be a parameter??
            }
            new_goals.push_back(niche_goals);
        }
        return new_goals;
    }

    std::vector<std::vector<double>> read_niches_from_file() {

        std::ifstream ifs(parameters.niche_file_name);
        if (!ifs.is_open()) { std::cerr << "Unable to open selection goals file " << parameters.niche_file_name << '\n'; exit(EXIT_FAILURE); }
        std::vector< std::vector< double >> new_goals;
        for (size_t niche_nr = 0; niche_nr < parameters.num_niches; niche_nr++) {
            std::vector<double> new_niche;
            for (size_t trait_nr = 0; trait_nr < parameters.num_traits; trait_nr++) {
                std::string parId_IS;
                ifs >> parId_IS;
                std::string parID_SHOULD = "Niche_" + std::to_string(niche_nr) + "_Trait_" + std::to_string(trait_nr);
                //if (parId_IS == parID_SHOULD) {
                double next_trait;
                ifs >> next_trait;
                new_niche.push_back(next_trait);
                //}
                /*else {
                    std::cerr << "unknown selection goal in file";
                    std::cerr << "found: " << parId_IS << std::endl;
                    std::cerr << "expected: " << parID_SHOULD << std::endl;
                    exit(EXIT_FAILURE);
                }*/
            }
            new_goals.push_back(new_niche);
        }

        ifs.close();
        return new_goals;
    }

    Param parameters;
    rnd_j master_random_generator;
    Output record;

    std::vector< Niche > world;
};


#endif /* Anolis_h */
