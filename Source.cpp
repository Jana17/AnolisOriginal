#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip>
#include <utility>
#include <algorithm>
#include "rndutils.hpp"   // some performance gains


/// Modified from: Anolis_Feb_14_2022_NewNoise; change: now the order in which the niches reproduce is randomized
/// Investment does not affect noise, as it is also done in Anolis_Feb_14_2022_NewNoise
/// Noise only on sex sel (as it also is in Anolis_Jan_26_2022)
/// Mut steps smaller for dewlap investment (as it is also in Anolis_Feb_03_2022_NewMut)
/// Addendum 18/02: I added Individual.get_size() etc. functions, because I wanted to plot individual phenotypes and like a complete moron put those into the "private" section
/// Addendum 21/02: I added output bits from extinction metrics bc i want to look at them
/// Addendum 22/02: I added output bits to get niche fit of dispersers (File name: "DisperserFit....")

// we can easily switch random number generators
// using reng_type = std::mt19937_64;
using reng_type = rndutils::default_engine;


constexpr int maxseed = 30;
//int maxtime = 2000;
constexpr int saveinterval = 1000;

//int maxseed = 50;
constexpr int maxtime = 20001;
//int saveinterval = 100;


constexpr int popsize = 1000; //per niche
constexpr double mutrate = 0.01;
constexpr double muteffect_size = 3;
constexpr double deathrate = 0.05;
constexpr double birthrate = 0.8;
constexpr double basal_dispersal = 0.01;


//Ideal Phenotypes: size, front limb, hind limb, lamellae, tail, movement. 
//if length does not matter, put zero - but later, it is not treated as "ideal length 0" but as "doesn't influence fitness"
// std::max(100 - trait, trait - 0);
constexpr std::array<double, 6> grass_bush{ 10,0,90,50,90,10 };
constexpr std::array<double, 6> trunk_ground{ 50,0,90,50,50,10 };
constexpr std::array<double, 6> trunk{ 10,50,50,50,10,90 };
constexpr std::array<double, 6> trunk_crown{ 50,10,10,90,50,90 };//size is given as intermediate but can also be small
constexpr std::array<double, 6> crown_giant{ 90,10,10,50,50,10 };
constexpr std::array<double, 6> twig{ 50,10,10,10,10,90 }; // size is given as intermediate but can also be small
//Maximal mismatches (given that all traits are bounded between 0 and 100):
constexpr double maxmis_BG = 90.0 + 0.0 + 90.0 + 50.0 + 90.0 + 90.0;
constexpr double maxmis_TG = 50.0 + 0.0 + 90.0 + 50.0 + 50.0 + 90.0;
constexpr double maxmis_T = 90.0 + 50.0 + 50.0 + 50.0 + 90.0 + 90.0;
constexpr double maxmis_TC = 50.0 + 90.0 + 90.0 + 90.0 + 50.0 + 90.0;
constexpr double maxmis_CG = 90.0 + 90.0 + 90.0 + 50.0 + 50.0 + 90.0;
constexpr double maxmis_Tw = 50.0 + 90.0 + 90.0 + 90.0 + 90.0 + 90.0;


class Individual {
public:
    double size;
    double front_limb;
    double hind_limb;
    double lamellae;
    double tail;
    double movement;

    double size_male;
    double front_limb_male;
    double hind_limb_male;
    double lamellae_male;
    double tail_male;
    double movement_male;

    double size_female;
    double front_limb_female;
    double hind_limb_female;
    double lamellae_female;
    double tail_female;
    double movement_female;

    int sex;

    double carotenoid_reserve; //depends on fit to niche
    double carotenoid_investment; //evolvable
    double dewlap_colour; //depends on reserve and investment

    void set_phen(double proportion_sexindep);
    void set_niche_fit();
    void set_niche(int NicheNr);
    double get_niche_fit() const;
    double get_phen_size() const;
    double get_phen_front_limb() const;
    double get_phen_hind_limb() const;
    double get_phen_lamellae() const;
    double get_phen_tail() const;
    double get_phen_movement() const;
    int get_niche() const;

private:
    double phen_size;
    double phen_front_limb;
    double phen_hind_limb;
    double phen_lamellae;
    double phen_tail;
    double phen_movement;
    int niche;
    double niche_fit;
};


void Individual::set_phen(double d_proportion_sexindep) {
    if (sex == 0) {
        phen_size = d_proportion_sexindep * size + (1 - d_proportion_sexindep) * size_female;
        phen_front_limb = d_proportion_sexindep * front_limb + (1 - d_proportion_sexindep) * front_limb_female;
        phen_hind_limb = d_proportion_sexindep * hind_limb + (1 - d_proportion_sexindep) * hind_limb_female;
        phen_lamellae = d_proportion_sexindep * lamellae + (1 - d_proportion_sexindep) * lamellae_female;
        phen_tail = d_proportion_sexindep * tail + (1 - d_proportion_sexindep) * tail_female;
        phen_movement = d_proportion_sexindep * movement + (1 - d_proportion_sexindep) * movement_female;
    }
    else {
        phen_size = d_proportion_sexindep * size + (1 - d_proportion_sexindep) * size_male;
        phen_front_limb = d_proportion_sexindep * front_limb + (1 - d_proportion_sexindep) * front_limb_male;
        phen_hind_limb = d_proportion_sexindep * hind_limb + (1 - d_proportion_sexindep) * hind_limb_male;
        phen_lamellae = d_proportion_sexindep * lamellae + (1 - d_proportion_sexindep) * lamellae_male;
        phen_tail = d_proportion_sexindep * tail + (1 - d_proportion_sexindep) * tail_male;
        phen_movement = d_proportion_sexindep * movement + (1 - d_proportion_sexindep) * movement_male;
    }
}

void Individual::set_niche(int nicheNr) {
    niche = nicheNr;
}

void Individual::set_niche_fit() {

    double max_mismatch;
    double current_mismatch;

    switch (niche) {
    case 0: {
        //Here, the front limb length does not matter. Therefore, the mismatch is always zero.
        current_mismatch = std::abs(phen_size - grass_bush[0]) + 0 + std::abs(phen_hind_limb - grass_bush[2]) + std::abs(phen_lamellae - grass_bush[3]) + std::abs(phen_tail - grass_bush[4]) + std::abs(phen_movement - grass_bush[5]);
        max_mismatch = maxmis_BG;
        break;
    }
    case 1: {
        current_mismatch = std::abs(phen_size - trunk_ground[0]) + 0 + std::abs(phen_hind_limb - trunk_ground[2]) + std::abs(phen_lamellae - trunk_ground[3]) + std::abs(phen_tail - trunk_ground[4]) + std::abs(phen_movement - trunk_ground[5]);
        max_mismatch = maxmis_TG;
        break;
    }
    case 2: {
        current_mismatch = std::abs(phen_size - trunk[0]) + std::abs(phen_front_limb - trunk[1]) + std::abs(phen_hind_limb - trunk[2]) + std::abs(phen_lamellae - trunk[3]) + std::abs(phen_tail - trunk[4]) + std::abs(phen_movement - trunk[5]);
        max_mismatch = maxmis_T;
        break;
    }
    case 3: {
        //now for trunk-crown we need to be careful: they can be small (10) or intermediate (50) in size, it does not matter.
        double size_mismatch;
        if (phen_size < 10) { size_mismatch = 10 - phen_size; }
        if (phen_size > 50) { size_mismatch = phen_size - 50; }
        else { size_mismatch = 0; }
        current_mismatch = size_mismatch + std::abs(phen_front_limb - trunk_crown[1]) + std::abs(phen_hind_limb - trunk_crown[2]) + std::abs(phen_lamellae - trunk_crown[3]) + std::abs(phen_tail - trunk_crown[4]) + std::abs(phen_movement - trunk_crown[5]);
        max_mismatch = maxmis_TC;
        break;
    }
    case 4: {
        current_mismatch = std::abs(phen_size - crown_giant[0]) + std::abs(phen_front_limb - crown_giant[1]) + std::abs(phen_hind_limb - crown_giant[2]) + std::abs(phen_lamellae - crown_giant[3]) + std::abs(phen_tail - crown_giant[4]) + std::abs(phen_movement - crown_giant[5]);
        max_mismatch = maxmis_CG;
        break;
    }
    case 5: {
        double size_mismatch_b;
        if (phen_size < 10) { size_mismatch_b = 10 - phen_size; }
        if (phen_size > 50) { size_mismatch_b = phen_size - 50; }
        else { size_mismatch_b = 0; }
        current_mismatch = size_mismatch_b + std::abs(phen_front_limb - twig[1]) + std::abs(phen_hind_limb - twig[2]) + std::abs(phen_lamellae - twig[3]) + std::abs(phen_tail - twig[4]) + std::abs(phen_movement - twig[5]);
        max_mismatch = maxmis_Tw;
        break;
    }
    default: {
        throw std::runtime_error("Error! Individual in non-existing niche!");
    }
    }

    double niche_fit_total = max_mismatch - current_mismatch;
    niche_fit = niche_fit_total / max_mismatch;

    //std::cout << "Niche is " << niche << std::endl;
    //std::cout << "niche mismatch: " << niche_mismatch << std::endl;
    //std::cout << "niche fit total: " << niche_fit_total << std::endl;
    //std::cout << "niche fit: " << niche_fit << std::endl;

}

double Individual::get_niche_fit() const { return niche_fit; }
double Individual::get_phen_size() const { return phen_size; }
double Individual::get_phen_front_limb() const { return phen_front_limb; }
double Individual::get_phen_hind_limb() const { return phen_hind_limb; }
double Individual::get_phen_lamellae() const { return phen_lamellae; }
double Individual::get_phen_tail() const { return phen_tail; }
double Individual::get_phen_movement() const { return phen_movement; }

int Individual::get_niche() const { return niche; }

bool randomBernoulli(double p, reng_type& rng) {
    std::bernoulli_distribution randomBernoulli(p);
    return randomBernoulli(rng);
}

std::vector<Individual> InitializePop(const Param& P) {
    Individual NewIndividual(Param);
    
}

//Initializing the population - only in one niche (initialize at "twig")
std::vector<Individual> InitializePop(int initialPopsize,
                                      double init_size,
                                      double init_front,
                                      double init_hind,
                                      double init_lam,
                                      double init_tail,
                                      double init_move,
                                      double prop_sexindep,
                                      int init_sex,
                                      std::uniform_real_distribution<double>& Noise,
                                      reng_type& rng) {

    std::vector<Individual> Population;
    for (int i = 0; i < initialPopsize; i++) {
        Individual NewIndividual;
        // this should go into Individual::Individual(...)
        NewIndividual.size = init_size + Noise(rng);
        NewIndividual.size_female = init_size + Noise(rng);
        NewIndividual.size_male = init_size + Noise(rng);
        NewIndividual.front_limb = init_front + Noise(rng);
        NewIndividual.front_limb_female = init_front + Noise(rng);
        NewIndividual.front_limb_male = init_front + Noise(rng);
        NewIndividual.hind_limb = init_hind + Noise(rng);
        NewIndividual.hind_limb_female = init_hind + Noise(rng);
        NewIndividual.hind_limb_male = init_hind + Noise(rng);
        NewIndividual.lamellae = init_lam + Noise(rng);
        NewIndividual.lamellae_female = init_lam + Noise(rng);
        NewIndividual.lamellae_male = init_lam + Noise(rng);
        NewIndividual.tail = init_tail + Noise(rng);
        NewIndividual.tail_female = init_tail + Noise(rng);
        NewIndividual.tail_male = init_tail + Noise(rng);
        NewIndividual.movement = init_move + Noise(rng);
        NewIndividual.movement_female = init_move + Noise(rng);
        NewIndividual.movement_male = init_move + Noise(rng);
        NewIndividual.sex = init_sex;
        NewIndividual.carotenoid_investment = 0.1;
        NewIndividual.set_phen(prop_sexindep);
        NewIndividual.set_niche(5);
        NewIndividual.set_niche_fit();
        Population.push_back(NewIndividual);
    }

    return Population;
}

//void mutate_trait(double& trait, double change) {
//    trait = std::clamp(trait + change, 0, 100);
//}

void do_mutations(Individual& CurrentIndiv,
                  const std::bernoulli_distribution& MutDistrib,
                  const std::cauchy_distribution<double>& muteffectDistrib,
                  reng_type& rng) {
    // would help if traits are iterable.
    auto mutate_trait = [&](double& trait, double max_trait = 100) {
        if (MutDistrib(rng)) {
            double trait_change = muteffectDistrib(rng);
            if (max_trait == 1.0) { trait_change = trait_change / 100.0; }
            trait = std::clamp(trait + trait_change, 0.0, max_trait);
        }
    };

    mutate_trait(CurrentIndiv.size);
    mutate_trait(CurrentIndiv.size_female);
    mutate_trait(CurrentIndiv.size_male);
    mutate_trait(CurrentIndiv.front_limb);
    mutate_trait(CurrentIndiv.front_limb_female);
    mutate_trait(CurrentIndiv.front_limb_male);
    mutate_trait(CurrentIndiv.hind_limb);
    mutate_trait(CurrentIndiv.hind_limb_female);
    mutate_trait(CurrentIndiv.hind_limb_male);
    mutate_trait(CurrentIndiv.lamellae);
    mutate_trait(CurrentIndiv.lamellae_female);
    mutate_trait(CurrentIndiv.lamellae_male);
    mutate_trait(CurrentIndiv.tail);
    mutate_trait(CurrentIndiv.tail_female);
    mutate_trait(CurrentIndiv.tail_male);
    mutate_trait(CurrentIndiv.movement);
    mutate_trait(CurrentIndiv.movement_female);
    mutate_trait(CurrentIndiv.movement_male);

    //add mutation of investment - this one is bounded between 1 and 0
    mutate_trait(CurrentIndiv.carotenoid_investment, 1.0);
}


// just to give things names...
struct pop_t {
    std::vector<Individual> females;
    std::vector<Individual> males;
};


std::vector<double> summarize_pop(const std::vector<pop_t>& Populations) {
    std::vector<double> pop_summary;
    std::vector<double> dewlaps;
    for (int k = 0; k < Populations.size(); k++) {

        double average_fit_female = 0;
        double average_fit_male = 0;
        double average_dewlap_investment = 0;
        for (int l = 0; l < Populations[k].females.size(); l++) {
            average_fit_female = average_fit_female + Populations[k].females[l].get_niche_fit();
            average_dewlap_investment = average_dewlap_investment + Populations[k].females[l].carotenoid_investment;
        }
        for (int l = 0; l < Populations[k].males.size(); l++) {
            average_fit_male = average_fit_male + Populations[k].males[l].get_niche_fit();
            average_dewlap_investment = average_dewlap_investment + Populations[k].males[l].carotenoid_investment;
        }

        double popsize_female = static_cast<double>(Populations[k].females.size());
        double popsize_male = static_cast<double>(Populations[k].males.size());
        average_fit_female = average_fit_female / popsize_female;
        average_fit_male = average_fit_male / popsize_male;
        average_dewlap_investment = average_dewlap_investment / (popsize_female + popsize_male);

        pop_summary.push_back(average_fit_female);
        pop_summary.push_back(average_fit_male);
        pop_summary.push_back(popsize_female);
        pop_summary.push_back(popsize_male);
        dewlaps.push_back(average_dewlap_investment);


    }
    for (int k = 0; k < dewlaps.size(); k++) {
        pop_summary.push_back(dewlaps[k]);
    }

    return pop_summary;

}

Individual make_Offspring(const Individual& Mother,
    const Individual& Father,
    const std::bernoulli_distribution& Recombi,
    const std::bernoulli_distribution& Rando,
    reng_type& rng) {

    const bool mp = Rando(rng);
    const Individual& MainParent = mp ? Mother : Father;
    const Individual& RecombiOnlyParent = mp ? Father : Mother;
    Individual Offspring(MainParent);   // inherit from MainParent if no recombination

    //Size
    if (Recombi(rng)) {
        Offspring.size = RecombiOnlyParent.size;
        Offspring.size_female = RecombiOnlyParent.size_female;
        Offspring.size_male = RecombiOnlyParent.size_male;
    }
    //Front limb
    if (Recombi(rng)) {
        Offspring.front_limb = RecombiOnlyParent.front_limb;
        Offspring.front_limb_female = RecombiOnlyParent.front_limb_female;
        Offspring.front_limb_male = RecombiOnlyParent.front_limb_male;
    }
    //Hind limb
    if (Recombi(rng)) {
        Offspring.hind_limb = RecombiOnlyParent.hind_limb;
        Offspring.hind_limb_female = RecombiOnlyParent.hind_limb_female;
        Offspring.hind_limb_male = RecombiOnlyParent.hind_limb_male;
    }
    //Number of lamellae
    if (Recombi(rng)) {
        Offspring.lamellae = RecombiOnlyParent.lamellae;
        Offspring.lamellae_female = RecombiOnlyParent.lamellae_female;
        Offspring.lamellae_male = RecombiOnlyParent.lamellae_male;
    }
    //Tail
    if (Recombi(rng)) {
        Offspring.tail = RecombiOnlyParent.tail;
        Offspring.tail_female = RecombiOnlyParent.tail_female;
        Offspring.tail_male = RecombiOnlyParent.tail_male;
    }
    //Movement
    if (Recombi(rng)) {
        Offspring.movement = RecombiOnlyParent.movement;
        Offspring.movement_female = RecombiOnlyParent.movement_female;
        Offspring.movement_male = RecombiOnlyParent.movement_male;
    }
    //carotenoid investment
    if (Recombi(rng)) {
        Offspring.carotenoid_investment = RecombiOnlyParent.carotenoid_investment;
    }

    Offspring.sex = Rando(rng) ? 0 : 1;

    // are we really shure we didn't forgot something?
    return Offspring;
}

void run() {

    std::cout << maxmis_BG << std::endl;
    std::cout << maxmis_TG << std::endl;
    std::cout << maxmis_T << std::endl;
    std::cout << maxmis_TC << std::endl;
    std::cout << maxmis_CG << std::endl;
    std::cout << maxmis_Tw << std::endl;



    reng_type rng; //create rng

    double proportion_indep_of_sex;
    double recombi_rate;
    double proportion_dispersal_density;
    double carotenoid_noise;


    // pull temporaries out of inner loop:
    std::vector<Individual> PossibleMates;
    std::vector<double> Attractiveness;
    std::vector<Individual> tmp_pop;    // s. vibility selection
    rndutils::mutable_discrete_distribution<int, rndutils::all_zero_policy_uni> weightedLottery;   // s. reproduction

    for (int PropSexDep = 1; PropSexDep < 2; PropSexDep++) {
        if (PropSexDep == 0) { proportion_indep_of_sex = 0.0; }
        //else if (PropSexDep == 1) { proportion_indep_of_sex = 0.5; }
        else { proportion_indep_of_sex = 1.0; }

        for (int PropDensityDispersal = 0; PropDensityDispersal < 1; PropDensityDispersal++) {
            if (PropDensityDispersal == 0) { proportion_dispersal_density = 0.0; }
            //else if (PropDensityDispersal == 1) { proportion_dispersal_density = 0.5; }
            else { proportion_dispersal_density = 1.0; }

            for (int RecRate = 0; RecRate < 2; RecRate++) {
                if (RecRate == 0) { recombi_rate = 0.0; }
                else if (RecRate == 1) { recombi_rate = 0.1; }
                else { recombi_rate = 0.5; }

                for (int Noise_Carotenoids = 0; Noise_Carotenoids < 3; Noise_Carotenoids++) {
                    if (Noise_Carotenoids == 0) { carotenoid_noise = 0.0001; }
                    else if (Noise_Carotenoids == 1) { carotenoid_noise = 0.1; }
                    else if (Noise_Carotenoids == 2) { carotenoid_noise = 0.3; }
                    else { carotenoid_noise = 1.0; }

                    for (int SexSel = 1; SexSel < 4; SexSel++) {
                        for (int seed = 0; seed < maxseed; seed++) {
                            std::cout << "PropSexDep: " << PropSexDep << ", PropDensDisp: " << PropDensityDispersal << ", recombi: " << RecRate <<
                                ", Noise: " << Noise_Carotenoids << ", SexSel: " << SexSel << ", Seed: " << seed << std::endl;
                            rng.seed(seed);//set seed
                            

                            //Output file stream
                            std::ofstream ofs1("TraitValues_Anolis_PropSexDep_" + std::to_string(proportion_indep_of_sex) +
                                "_PropDensityDisp_" + std::to_string(proportion_dispersal_density) + "_RecRate_" + std::to_string(recombi_rate) +
                                "_GGNoise_" + std::to_string(carotenoid_noise) + "_SexSel_" + std::to_string(SexSel) +
                                "_Seed_" + std::to_string(seed) + ".csv");
                            //For normal data:
                            /*ofs1 << "Time" << ","
                                << "Average_niche_fit_GB_female" << "," << "Average_niche_fit_GB_male" << "," << "PopSize_GB_female" << "," << "PopSize_GB_male" << ","
                                << "Average_niche_fit_TG_female" << "," << "Average_niche_fit_TG_male" << "," << "PopSize_TG_female" << "," << "PopSize_TG_male" << ","
                                << "Average_niche_fit_T_female" << "," << "Average_niche_fit_T_male" << "," << "PopSize_T_female" << "," << "PopSize_T_male" << ","
                                << "Average_niche_fit_TC_female" << "," << "Average_niche_fit_TC_male" << "," << "PopSize_TC_female" << "," << "PopSize_TC_male" << ","
                                << "Average_niche_fit_CG_female" << "," << "Average_niche_fit_CG_male" << "," << "PopSize_CG_female" << "," << "PopSize_CG_male" << ","
                                << "Average_niche_fit_Tw_female" << "," << "Average_niche_fit_Tw_male" << "," << "PopSize_Tw_female" << "," << "PopSize_Tw_male" << ","
                                << "Average_dewlap_investment_GB" << "," << "Average_dewlap_investment_TG" << "," << "Average_dewlap_investment_T" << ","
                                << "Average_dewlap_investment_TC" << "," << "Average_dewlap_investment_CG" << "," << "Average_dewlap_investment_Tw" << ","
                                << "Seed" << ","
                                << "newbirth_count" << "," << "disperse_count" << "," << "stay_count" << "," << "prop_off_disperse" << ","
                                << "extinction_count" << "," << "recolonize_count";
                            ofs1 << '\n';*/

                            //For extinction metrics
                            /*ofs1 << "Time" << ","
                                << "Time_Since_Extinction_GB" << "," << "ExtinctionCount_GB" << "," << "Extinction_GB" << "," << "Entered_GB" << ","
                                << "Time_Since_Extinction_TG" << "," << "ExtinctionCount_TG" << "," << "Extinction_TG" << "," << "Entered_TG" << ","
                                << "Time_Since_Extinction_T" << "," << "ExtinctionCount_T" << "," << "Extinction_T" << "," << "Entered_T" << ","
                                << "Time_Since_Extinction_TC" << "," << "ExtinctionCount_TC" << "," << "Extinction_TC" << "," << "Entered_TC" << ","
                                << "Time_Since_Extinction_CG" << "," << "ExtinctionCount_CG" << "," << "Extinction_CG" << "," << "Entered_CG" << ","
                                << "Time_Since_Extinction_Tw" << "," << "ExtinctionCount_Tw" << "," << "Extinction_Tw" << "," << "Entered_Tw" << ","
                                << "Seed";
                            ofs1 << '\n';*/

                            //For disperser fit
                            /*ofs1 << "Time" << ","
                                << "NicheOld" << "," << "NicheNew" << "," << "NicheFitOld" << "," << "NicheFitNew" << "," << "Founder" << "," << "Sex" << ","
                                << "Seed";
                            ofs1 << '\n';*/

                            //For trait values fit
                            //TraitIDs: same order as in class Individual: size, front_limb, hind_limb, lamellae, tail,movement
                            ofs1 << "Time" << ","
                                << "TraitID" << "," << "TraitAverage" << "," << "TraitStDev" << "," << "NicheNr" << "," << "PopSizeM" << "," << "PopSizeF" << ","
                                << "Seed";
                            ofs1 << '\n';


                            //For individual data:
                            /*ofs1 << "Time" << ","
                                << "Seed" << ","
                                << "Niche" << ","
                                << "Niche_Fit" << ","
                                << "Dewlap_investment" << ","
                                << "Dewlap_colour" << ","
                                << "MPopSize" << ","
                                << "FPopSize";
                            ofs1 << '\n';*/

                            //For full individual data:
                            /*ofs1 << "Time" << ","
                                << "Seed" << ","
                                << "Niche" << ","
                                << "Niche_Fit" << ","
                                << "Size" << ","
                                << "Front_limb" << ","
                                << "Hind_limb" << ","
                                << "Lamellae" << ","
                                << "Tail" << ","
                                << "Movement" << ","
                                << "Dewlap_investment" << ","
                                << "Dewlap_colour" << ","
                                << "Sex" << ","
                                << "MPopSize" << ","
                                << "FPopSize";
                            ofs1 << '\n';*/


                            std::bernoulli_distribution Mutation(mutrate);
                            std::bernoulli_distribution Rando(0.5);
                            std::cauchy_distribution<double> muteffect(0.0, muteffect_size);
                            std::uniform_real_distribution<double> initializeWithNoise(-10, 10);
                            std::normal_distribution<double> NoiseInCarotenoids(0, carotenoid_noise);
                            //std::uniform_real<double> NoiseInCarotenoids(-carotenoid_noise, carotenoid_noise); //tried out: 20th Feb, for full noise
                            std::uniform_int_distribution<int> RandomHabitat(0, 5);
                            std::bernoulli_distribution Recombi(recombi_rate);


                            std::vector<Individual> Pop_twig_female = InitializePop(popsize / 2, twig[0], twig[1], twig[2], twig[3], twig[4], twig[5], proportion_indep_of_sex, 0, initializeWithNoise, rng);
                            std::vector<Individual> Pop_twig_male = InitializePop(popsize / 2, twig[0], twig[1], twig[2], twig[3], twig[4], twig[5], proportion_indep_of_sex, 1, initializeWithNoise, rng);
                            //std::vector<Individual> Pop_grass_bush_female = InitializePop(popsize / 2, grass_bush[0], grass_bush[1], grass_bush[2], grass_bush[3], grass_bush[4], grass_bush[5], proportion_indep_of_sex, 0, initializeWithNoise, rng);
                            //std::vector<Individual> Pop_grass_bush_male = InitializePop(popsize / 2, grass_bush[0], grass_bush[1], grass_bush[2], grass_bush[3], grass_bush[4], grass_bush[5], proportion_indep_of_sex, 1, initializeWithNoise, rng);

                            std::vector<pop_t> WholePop(6);
                            WholePop[5].females = Pop_twig_female;
                            WholePop[5].males = Pop_twig_female;

                            for (int i = 0; i < WholePop[5].males.size(); i++) {

                                WholePop[5].males[i].carotenoid_reserve = WholePop[5].males[i].get_niche_fit();
                                if (WholePop[5].males[i].carotenoid_reserve < 0) { WholePop[5].males[i].carotenoid_reserve = 0; }
                                if (WholePop[5].males[i].carotenoid_reserve > 1) { WholePop[5].males[i].carotenoid_reserve = 1; }

                                WholePop[5].males[i].dewlap_colour = (WholePop[5].males[i].carotenoid_investment * WholePop[5].males[i].carotenoid_reserve) + NoiseInCarotenoids(rng);
                                if (WholePop[5].males[i].dewlap_colour < 0) {
                                    WholePop[5].males[i].dewlap_colour = 0;
                                }
                                if (WholePop[5].males[i].dewlap_colour < 0) {
                                    std::cout << "Dewlap colour below zero! Investment: " << WholePop[5].males[i].carotenoid_investment << ", Reserve: " << WholePop[5].males[i].carotenoid_reserve << std::endl;
                                }
                                if (WholePop[5].males[i].dewlap_colour > 1) {
                                    WholePop[5].males[i].dewlap_colour = 1;
                                }


                                WholePop[5].males[i].carotenoid_reserve = WholePop[5].males[i].carotenoid_reserve - (WholePop[5].males[i].carotenoid_investment * WholePop[5].males[i].carotenoid_reserve);
                                if (WholePop[5].males[i].carotenoid_reserve < 0) {
                                    std::cout << "Reserve below zero! Investment: " << WholePop[5].males[i].carotenoid_investment << ", Reserve: " << WholePop[5].males[i].carotenoid_reserve << std::endl;
                                }
                            }

                            int extinct_count = 0;
                            int recolonize_count = 0;

                            /// Following variables added (21st Feb) to output time to extinction
                            std::vector<int> TimeSinceExtinction(6, 0);
                            std::vector<int> ExtinctionCount(6, 0);
                            std::vector<bool> Extinction(6, false);
                            std::vector<bool> Entered{ false,false,false,false,false,true };


                            for (int time = 0; time < maxtime; time++) {
                                //if (time % saveinterval == 0) { std::cout << " time is " << time << std::endl; }

                                //vibility selection
                                for (int n = 0; n < WholePop.size(); n++) {
                                    bool niche_occupied_before = false;
                                    double entire_niche_size_before = WholePop[n].females.size() + WholePop[n].males.size();
                                    if (entire_niche_size_before != 0) {
                                        niche_occupied_before = true;
                                    }
                                    tmp_pop.clear();
                                    for (const auto& ind : WholePop[n].females) {
                                        double indiv_deathrate = 0.75 * deathrate + 0.25 * std::exp(-2 * std::pow(ind.get_niche_fit(), 8));
                                        std::bernoulli_distribution DeathEvent(indiv_deathrate);
                                        if (!DeathEvent(rng)) {
                                            tmp_pop.emplace_back(ind);
                                        }
                                    }
                                    std::swap(WholePop[n].females, tmp_pop);

                                    tmp_pop.clear();
                                    for (const auto& ind : WholePop[n].males) {
                                        double indiv_deathrate = 0.75 * deathrate + 0.25 * std::exp(-2 * pow(ind.carotenoid_reserve, 8));
                                        std::bernoulli_distribution DeathEvent(indiv_deathrate);
                                        if (!DeathEvent(rng)) {
                                            tmp_pop.emplace_back(ind);
                                        }
                                    }
                                    std::swap(WholePop[n].males, tmp_pop);

                                    bool niche_occupied_after = false;
                                    double entire_niche_size_after = WholePop[n].females.size() + WholePop[n].males.size();
                                    if (entire_niche_size_after != 0) {
                                        niche_occupied_after = true;
                                    }


                                    /// Following if statement added (21st Feb) to output time to extinction
                                    if (Entered[n] == true) {
                                        //std::cout << "we're in 'if entered'" << std::endl;
                                        if (niche_occupied_before != niche_occupied_after) {
                                            //std::cout << "we're in 'if extinction'" << std::endl;
                                            //Niche was occupied before but now isn't
                                            //therefore: extinction!
                                            //sanity check:
                                            if (niche_occupied_before != true) { std::cout << "Niche wasn't occupied before viability sel, erraneous extinction!" << std::endl; }
                                            if (niche_occupied_after == true) { std::cout << "Niche is still occupied, erraneous extinction!" << std::endl; }
                                            //update the vectors accordingly
                                            TimeSinceExtinction[n] = 0;
                                            ExtinctionCount[n] = ExtinctionCount[n] + 1;
                                            Extinction[n] = true;
                                            Entered[n] = false;
                                        }
                                        else {
                                            //std::cout << "we're in 'else'" << std::endl;
                                            //No extinction, niche still occupied!
                                            //sanity check:
                                            if (niche_occupied_before != true) { std::cout << "Niche wasn't occupied before viability sel, error!" << std::endl; }
                                            if (niche_occupied_after != true) { std::cout << "Niche is not occupied anymore, overlooked extinction!" << std::endl; }
                                            //update the vectors accordingly
                                            //std::cout << "Time since extinction: " <<TimeSinceExtinction[n] << std::endl;
                                            TimeSinceExtinction[n] = TimeSinceExtinction[n] + 1;
                                            Extinction[n] = false;
                                        }
                                    }
                                    else { Extinction[n] = false; }




                                    if (niche_occupied_before != niche_occupied_after) {
                                        extinct_count = extinct_count + 1;
                                    }
                                    /*
                                    * Don't want this:
                                    * std::vector::erase is expensive
                                    * mutating a container while iterating over it is ugly, slow and dangerous
                                                                        for (int m = 0; m < WholePop[n].females.size(); m++) {
                                                                            //old function:
                                                                            //double indiv_deathrate = 0.75 * deathrate + 0.25 * std::exp(-2 * pow(WholePop[n].females[m].get_niche_fit(), 10));
                                                                            //new function:
                                                                            double indiv_deathrate = 0.75 * deathrate + 0.25 * std::exp(-2 * pow(WholePop[n].females[m].get_niche_fit(), 10));
                                                                            std::bernoulli_distribution DeathEvent(indiv_deathrate);
                                                                            if (DeathEvent(rng)) {
                                                                                WholePop[n].females.erase(WholePop[n].females.begin() + m);
                                                                                m = m - 1;
                                                                            }
                                                                        }
                                                                        for (int m = 0; m < WholePop[n].malessize(); m++) {
                                                                            //old function:
                                                                            //double indiv_deathrate = 0.75 * deathrate + 0.25 * std::exp(-2 * pow(WholePop[n].males[m].carotenoid_reserve, 10));
                                                                            //new function:
                                                                            double indiv_deathrate = 0.75 * deathrate + 0.25 * std::exp(-2 * pow(WholePop[n].males[m].carotenoid_reserve, 10));
                                                                            std::bernoulli_distribution DeathEvent(indiv_deathrate);
                                                                            if (DeathEvent(rng)) {
                                                                                WholePop[n].maleserase(WholePop[n].males.begin() + m);
                                                                                m = m - 1;
                                                                            }
                                                                        }
                                    */
                                }//end viability sel


                                //SAVE POPULATION!
                                //save averages
                                /*if (time % saveinterval == 0) {
                                    std::vector<double> PopSummary = summarize_pop(WholePop);
                                    ofs1 << time << ", "
                                        << PopSummary[0] << ", " << PopSummary[1] << ", " << PopSummary[2] << ", " << PopSummary[3] << ", "
                                        << PopSummary[4] << ", " << PopSummary[5] << ", " << PopSummary[6] << ", " << PopSummary[7] << ", "
                                        << PopSummary[8] << ", " << PopSummary[9] << ", " << PopSummary[10] << ", " << PopSummary[11] << ", "
                                        << PopSummary[12] << ", " << PopSummary[13] << ", " << PopSummary[14] << ", " << PopSummary[15] << ", "
                                        << PopSummary[16] << ", " << PopSummary[17] << ", " << PopSummary[18] << ", " << PopSummary[19] << ", "
                                        << PopSummary[20] << ", " << PopSummary[21] << ", " << PopSummary[22] << ", " << PopSummary[23] << ", "
                                        << PopSummary[24] << ", " << PopSummary[25] << ", " << PopSummary[26] << ", " << PopSummary[27] << ", " << PopSummary[28] << ", " << PopSummary[29] << ", "
                                        << seed << ", ";
                                }*/
                                //save individual data
                                 /*if (time % saveinterval == 0) {
                                    for (int whichniche = 0; whichniche < WholePop.size(); whichniche++) {
                                        bool is_established = false;
                                        if (WholePop[whichniche].males.size() > 50) { is_established = true; }
                                        for (int indivnr = 0; indivnr < WholePop[whichniche].males.size(); indivnr++) {
                                            ofs1 << time << ", "
                                                << seed << ", "
                                                << whichniche << ","
                                                << WholePop[whichniche].males[indivnr].get_niche_fit() << ", "
                                                << WholePop[whichniche].males[indivnr].carotenoid_investment << ", "
                                                << WholePop[whichniche].males[indivnr].dewlap_colour << ", "
                                                << is_established << '\n';
                                        }
                                    }
                                }*/
                                //save full indivi data
                                /*if (time % saveinterval == 0) {
                                    for (int whichniche = 0; whichniche < WholePop.size(); whichniche++) {
                                        for (int indivnr = 0; indivnr < WholePop[whichniche].males.size(); indivnr++) {
                                            ofs1 << time << ", "
                                                << seed << ", "
                                                << WholePop[whichniche].males[indivnr].get_niche() << ", "
                                                << WholePop[whichniche].males[indivnr].get_niche_fit() << ", "
                                                << WholePop[whichniche].males[indivnr].get_phen_size() << ", "
                                                << WholePop[whichniche].males[indivnr].get_phen_front_limb() << ", "
                                                << WholePop[whichniche].males[indivnr].get_phen_hind_limb() << ", "
                                                << WholePop[whichniche].males[indivnr].get_phen_lamellae() << ", "
                                                << WholePop[whichniche].males[indivnr].get_phen_tail() << ", "
                                                << WholePop[whichniche].males[indivnr].get_phen_movement() << ", "
                                                << WholePop[whichniche].males[indivnr].carotenoid_investment << ", "
                                                << WholePop[whichniche].males[indivnr].dewlap_colour << ", "
                                                << WholePop[whichniche].males[indivnr].sex << ", "
                                                << WholePop[whichniche].males.size() << ", "
                                                << WholePop[whichniche].females.size() << '\n';
                                        }
                                        for (int indivnr = 0; indivnr < WholePop[whichniche].females.size(); indivnr++) {
                                            ofs1 << time << ", "
                                                << seed << ", "
                                                << WholePop[whichniche].females[indivnr].get_niche() << ", "
                                                << WholePop[whichniche].females[indivnr].get_niche_fit() << ", "
                                                << WholePop[whichniche].females[indivnr].get_phen_size() << ", "
                                                << WholePop[whichniche].females[indivnr].get_phen_front_limb() << ", "
                                                << WholePop[whichniche].females[indivnr].get_phen_hind_limb() << ", "
                                                << WholePop[whichniche].females[indivnr].get_phen_lamellae() << ", "
                                                << WholePop[whichniche].females[indivnr].get_phen_tail() << ", "
                                                << WholePop[whichniche].females[indivnr].get_phen_movement() << ", "
                                                << WholePop[whichniche].females[indivnr].carotenoid_investment << ", "
                                                << WholePop[whichniche].females[indivnr].dewlap_colour << ", "
                                                << WholePop[whichniche].females[indivnr].sex << ", "
                                                << WholePop[whichniche].males.size() << ", "
                                                << WholePop[whichniche].females.size() << '\n';
                                        }
                                    }
                                }*/

                                //save trait averages &stdevs - here only for sex-indep expressed genes
                                if (time % saveinterval == 0) {
                                    for (int whichniche = 0; whichniche < WholePop.size(); whichniche++) {
                                        //averages
                                        double average_size = 0;
                                        double average_front_limb = 0;
                                        double average_hind_limb = 0;
                                        double average_lamellae = 0;
                                        double average_tail = 0;
                                        double average_movement = 0;
                                        for (int l = 0; l < WholePop[whichniche].females.size(); l++) {
                                            average_size = average_size + WholePop[whichniche].females[l].size;
                                            average_front_limb = average_front_limb + WholePop[whichniche].females[l].front_limb;
                                            average_hind_limb = average_hind_limb + WholePop[whichniche].females[l].hind_limb;
                                            average_lamellae = average_lamellae + WholePop[whichniche].females[l].lamellae;
                                            average_tail = average_tail + WholePop[whichniche].females[l].tail;
                                            average_movement = average_movement + WholePop[whichniche].females[l].movement;
                                        }
                                        for (int l = 0; l < WholePop[whichniche].males.size(); l++) {
                                            average_size = average_size + WholePop[whichniche].males[l].size;
                                            average_front_limb = average_front_limb + WholePop[whichniche].males[l].front_limb;
                                            average_hind_limb = average_hind_limb + WholePop[whichniche].males[l].hind_limb;
                                            average_lamellae = average_lamellae + WholePop[whichniche].males[l].lamellae;
                                            average_tail = average_tail + WholePop[whichniche].males[l].tail;
                                            average_movement = average_movement + WholePop[whichniche].males[l].movement;
                                        }
                                        double popsize_female = static_cast<double>(WholePop[whichniche].females.size());
                                        double popsize_male = static_cast<double>(WholePop[whichniche].males.size());
                                        double whole_popsize = popsize_female + popsize_male;
                                        average_size = average_size / whole_popsize;
                                        average_front_limb = average_front_limb / whole_popsize;
                                        average_hind_limb = average_hind_limb / whole_popsize;
                                        average_lamellae = average_lamellae / whole_popsize;
                                        average_tail = average_tail / whole_popsize;
                                        average_movement = average_movement / whole_popsize;
                                        //stdevs
                                        double dev_sum_size = 0;
                                        double dev_sum_front_limb = 0;
                                        double dev_sum_hind_limb = 0;
                                        double dev_sum_lamellae = 0;
                                        double dev_sum_tail = 0;
                                        double dev_sum_movement = 0;
                                        for (int l = 0; l < WholePop[whichniche].females.size(); l++) {
                                            dev_sum_size += pow((WholePop[whichniche].females[l].size - average_size), 2);
                                            dev_sum_front_limb += pow((WholePop[whichniche].females[l].front_limb - average_front_limb), 2);
                                            dev_sum_hind_limb += pow((WholePop[whichniche].females[l].hind_limb - average_hind_limb), 2);
                                            dev_sum_lamellae += pow((WholePop[whichniche].females[l].lamellae - average_lamellae), 2);
                                            dev_sum_tail += pow((WholePop[whichniche].females[l].tail - average_tail), 2);
                                            dev_sum_movement += pow((WholePop[whichniche].females[l].movement - average_movement), 2);
                                        }
                                        for (int l = 0; l < WholePop[whichniche].males.size(); l++) {
                                            dev_sum_size += pow((WholePop[whichniche].males[l].size - average_size), 2);
                                            dev_sum_front_limb += pow((WholePop[whichniche].males[l].front_limb - average_front_limb), 2);
                                            dev_sum_hind_limb += pow((WholePop[whichniche].males[l].hind_limb - average_hind_limb), 2);
                                            dev_sum_lamellae += pow((WholePop[whichniche].males[l].lamellae - average_lamellae), 2);
                                            dev_sum_tail += pow((WholePop[whichniche].males[l].tail - average_tail), 2);
                                            dev_sum_movement += pow((WholePop[whichniche].males[l].movement - average_movement), 2);
                                        }
                                        double std_dev_size = sqrt(dev_sum_size / whole_popsize);
                                        double std_dev_front_limb = sqrt(dev_sum_front_limb / whole_popsize);
                                        double std_dev_hind_limb = sqrt(dev_sum_hind_limb / whole_popsize);
                                        double std_dev_lamellae = sqrt(dev_sum_lamellae / whole_popsize);
                                        double std_dev_tail = sqrt(dev_sum_tail / whole_popsize);
                                        double std_dev_movement = sqrt(dev_sum_movement / whole_popsize);

                                        ofs1 << time << ", " << 1 << ", "<< average_size << ", " << std_dev_size << ", " << whichniche << ", "<< popsize_male << ", "<< popsize_female << ", " << seed << '\n';
                                        ofs1 << time << ", " << 2 << ", " << average_front_limb << ", " << std_dev_front_limb << ", " << whichniche << ", " << popsize_male << ", " << popsize_female << ", " << seed << '\n';
                                        ofs1 << time << ", " << 3 << ", " << average_hind_limb << ", " << std_dev_hind_limb << ", " << whichniche << ", " << popsize_male << ", " << popsize_female << ", " << seed << '\n';
                                        ofs1 << time << ", " << 4 << ", " << average_lamellae << ", " << std_dev_lamellae << ", " << whichniche << ", " << popsize_male << ", " << popsize_female << ", " << seed << '\n';
                                        ofs1 << time << ", " << 5 << ", " << average_tail << ", " << std_dev_tail << ", " << whichniche << ", " << popsize_male << ", " << popsize_female << ", " << seed << '\n';
                                        ofs1 << time << ", " << 6 << ", " << average_movement << ", " << std_dev_movement << ", " << whichniche << ", " << popsize_male << ", " << popsize_female << ", " << seed << '\n';

                                    }//end which niche loop
                                }//end saving loop

                                //reproduction
                                int disperse_count = 0;
                                int stay_count = 0;
                                int newbirth_count = 0;
                                double prop_offs_disperse = 0.0;
                                std::vector<int> ReproNumbers{ 0,1,2,3,4,5 };
                                std::shuffle(ReproNumbers.begin(), ReproNumbers.end(), rng);
                                int n;
                                
                                for (size_t n_old = 0; n_old < WholePop.size(); n_old++) {
                                    n = ReproNumbers[n_old];
                                    size_t CurrentSizeFemale = WholePop[n].females.size();
                                    size_t CurrentSizeMale = WholePop[n].males.size();

                                    if (CurrentSizeMale != 0) {
                                        double available_spaces = static_cast<double>(popsize) - static_cast<double>(WholePop[n].females.size()) - static_cast<double>(WholePop[n].males.size());
                                        if (available_spaces < 0) { available_spaces = 0; }
                                        double proportion_available = static_cast<double>(available_spaces) / static_cast<double>(popsize);
                                        auto mate_choice = std::uniform_int_distribution<size_t>(0, CurrentSizeMale - 1);

                                        for (size_t m = 0; m < CurrentSizeFemale; m++) {
                                            double indiv_birthrate = (0.8 * birthrate + 0.2 * WholePop[n].females[m].get_niche_fit()) * proportion_available;
                                            std::bernoulli_distribution BirthEvent(indiv_birthrate);
                                            if (BirthEvent(rng)) {
                                                if (n == 5) { newbirth_count = newbirth_count + 1; }
                                                const Individual& Mother = WholePop[n].females[m];
                                                Individual Father = WholePop[n].males[mate_choice(rng)];  // default: random mate choice

                                                //sexual selection: assume female choice based on dewlap. Best of n: female meets n males.
                                                //think of cartenoids: better males can make prettier dewlaps.
                                                //specifically, if you have more food, you can invest more into dewlaps bc you have enough carotenoid resources.
                                                if (SexSel != 1) {
                                                    PossibleMates.clear();
                                                    Attractiveness.clear();
                                                    bool AllAttractivenessZero = true;
                                                    Individual PossibleFather;
                                                    for (int mates = 0; mates < SexSel; mates++) {
                                                        PossibleFather = WholePop[n].males[mate_choice(rng)];
                                                        if (PossibleFather.dewlap_colour < 0) {
                                                            std::cout << "Negative colour!" << std::endl;
                                                            std::cout << "Colour: " << PossibleFather.dewlap_colour << ", Resources: " << PossibleFather.carotenoid_reserve << std::endl;
                                                        }
                                                        PossibleMates.push_back(PossibleFather);
                                                        Attractiveness.push_back(PossibleFather.dewlap_colour);
                                                        if (PossibleFather.dewlap_colour != 0) { AllAttractivenessZero = false; }
                                                    }
                                                    if (AllAttractivenessZero == true) {
                                                        Father = PossibleMates[0];
                                                    }
                                                    else {
                                                        weightedLottery.mutate(Attractiveness.begin(), Attractiveness.end());
                                                        Father = PossibleMates[weightedLottery(rng)];
                                                    }
                                                }

                                                Individual Offspring = make_Offspring(Mother, Father, Recombi, Rando, rng);
                                                do_mutations(Offspring, Mutation, muteffect, rng);


                                                Offspring.set_phen(proportion_indep_of_sex);
                                                int parent_niche = Mother.get_niche(); //Both parents have the same niche anyways
                                                Offspring.set_niche(parent_niche);

                                                //dispersal!
                                                Offspring.set_niche_fit();
                                                double OldNicheFit = Offspring.get_niche_fit();
                                                double mismatch_dispersal = std::exp(-5 * pow(Offspring.get_niche_fit(), 2));
                                                //old function: 
                                                //double density_dispersal = std::exp(0.2 * pow(proportion_available, 4)) - 1;
                                                double density_dispersal = std::exp(0.5 * pow(proportion_available, 2)) - 1;
                                                //double density_dispersal = std::exp(0.2 * pow(proportion_available, 2)) - 1;
                                                double prob_disperse = std::max(basal_dispersal, ((1 - proportion_dispersal_density) * mismatch_dispersal
                                                    + proportion_dispersal_density * density_dispersal));

                                                bool Founder = false;
                                                if (std::bernoulli_distribution(prob_disperse)(rng)) {
                                                    if (n == 5) { disperse_count = disperse_count + 1; }
                                                    int NewNiche = RandomHabitat(rng);
                                                    double entire_new_niche_size = WholePop[NewNiche].females.size() + WholePop[NewNiche].males.size();
                                                    if (entire_new_niche_size == 0) {
                                                        recolonize_count = recolonize_count + 1;
                                                        Founder = true;
                                                    }
                                                    // This statement was added (21st Feb) to output time to extinction. 
                                                    // Unlike the version called "Extinction metrics, I decided to not have it in the if bracket, because it should not matter whether the niche previously was empty or not, entered is entered.
                                                    Entered[NewNiche] = true;
                                                    Offspring.set_niche(NewNiche);
                                                    Offspring.set_niche_fit();

                                                    /*ofs1 << time << ", "
                                                        << parent_niche << ", " << Offspring.get_niche() << ", " << OldNicheFit << ", " << Offspring.get_niche_fit() << ", " << Founder << ", " << Offspring.sex << ", "
                                                        << seed;
                                                    ofs1 << '\n';*/

                                                }
                                                else { if (n == 5) { stay_count = stay_count + 1; } }

                                                //make dewlap
                                                Offspring.carotenoid_reserve = Offspring.get_niche_fit();
                                                if (Offspring.carotenoid_reserve < 0) { Offspring.carotenoid_reserve = 0; }
                                                if (Offspring.carotenoid_reserve > 1) { Offspring.carotenoid_reserve = 1; }
                                                Offspring.dewlap_colour = (Offspring.carotenoid_investment * Offspring.carotenoid_reserve) + NoiseInCarotenoids(rng);
                                                if (Offspring.dewlap_colour < 0) {
                                                    Offspring.dewlap_colour = 0;
                                                }
                                                if (Offspring.dewlap_colour < 0) {
                                                    std::cout << "Dewlap colour below zero! Investment: " << Offspring.carotenoid_investment << ", Reserve: " << Offspring.carotenoid_reserve << std::endl;
                                                }
                                                if (Offspring.dewlap_colour > 1) {
                                                    Offspring.dewlap_colour = 1;
                                                }

                                                Offspring.carotenoid_reserve = Offspring.carotenoid_reserve - (Offspring.carotenoid_investment * Offspring.carotenoid_reserve);
                                                if (Offspring.carotenoid_reserve < 0) {
                                                    std::cout << "Reserve colour below zero! Investment: " << Offspring.carotenoid_investment << ", Reserve: " << Offspring.carotenoid_reserve << std::endl;
                                                }

                                                //saving new offspring
                                                /*if (time % saveinterval == 0) {
                                                    ofs1 << time << ", "
                                                        << seed << ", "
                                                        << Offspring.get_niche() << ","
                                                        << Offspring.get_niche_fit() << ", "
                                                        << Offspring.carotenoid_investment << ", "
                                                        << Offspring.dewlap_colour << ", "
                                                        << Offspring.sex << '\n';
                                                }*/
                                                
                                                //Use this individual data to check dewlap investement (18th Feb)
                                                /*if (Offspring.get_niche() == parent_niche) {
                                                    ofs1 << time << ", "
                                                        << seed << ", "
                                                        << Offspring.get_niche() << ","
                                                        << Offspring.get_niche_fit() << ", "
                                                        << Offspring.carotenoid_investment << ", "
                                                        << Offspring.dewlap_colour << ", "
                                                        << WholePop[Offspring.get_niche()].males.size() << ", " 
                                                        << WholePop[Offspring.get_niche()].females.size() << '\n';
                                                }*/
                                                

                                               

                                                if (Offspring.sex == 0) {
                                                    WholePop[Offspring.get_niche()].females.push_back(Offspring);
                                                }
                                                else {
                                                    WholePop[Offspring.get_niche()].males.push_back(Offspring);
                                                }

                                            }//end of "If there is a birth event" statement

                                        }//going through all females that are currently in this niche

                                    }//end of if statement for "if there is still a male in this niche"

                                }//end repro//end of reproduction loop

                                //save some more output!
                               /* ofs1 << time << ", "
                                    << TimeSinceExtinction[0] << ", " << ExtinctionCount[0] << ", " << Extinction[0] << ", " << Entered[0] << ", "
                                    << TimeSinceExtinction[1] << ", " << ExtinctionCount[1] << ", " << Extinction[1] << ", " << Entered[1] << ", "
                                    << TimeSinceExtinction[2] << ", " << ExtinctionCount[2] << ", " << Extinction[2] << ", " << Entered[2] << ", "
                                    << TimeSinceExtinction[3] << ", " << ExtinctionCount[3] << ", " << Extinction[3] << ", " << Entered[3] << ", "
                                    << TimeSinceExtinction[4] << ", " << ExtinctionCount[4] << ", " << Extinction[4] << ", " << Entered[4] << ", "
                                    << TimeSinceExtinction[5] << ", " << ExtinctionCount[5] << ", " << Extinction[5] << ", " << Entered[5] << ", "
                                    << seed << '\n';*/

                                /*prop_offs_disperse = static_cast<double>(disperse_count) / static_cast<double>(newbirth_count);
                                if (time % saveinterval == 0) {
                                    ofs1 << newbirth_count << ", " << disperse_count << ", " << stay_count << ", " << prop_offs_disperse
                                        << ", " << extinct_count << ", " << recolonize_count << '\n';
                                }*/

                            }//end of time loop
                        }//end of seed loop
                    }//end Sex Sel loop

                }//End Noise in Carotenoids loop

            }//End Recombi rate loop

        }//End dispersal type loop

    }//End sex dep expression loop
}


int main() {
    try {
        run();
        return 0;
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << '\n';
    }
    catch (...) {
        std::cerr << "unknown exception\n";
    }
    return -1;
}
