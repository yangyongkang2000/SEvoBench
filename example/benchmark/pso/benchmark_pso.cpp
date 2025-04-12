//-----------------------------------------------------------------------------
// RealPSO.cpp
//-----------------------------------------------------------------------------
//*
// An instance of a VERY simple Real-coded Particle Swarm Optimization Algorithm
//
//-----------------------------------------------------------------------------
#include "SEvoBench/sevobench.hpp"
#include "eo"
#include "eoModifyStandardVelocity.h"
#include <chrono>
#include <iostream>
#include <sstream>
#include <stdexcept>
// Use functions from namespace std
// using namespace std; // Do not do this, this shadows EO's `apply` function.

//-----------------------------------------------------------------------------
typedef eoMinimizingFitness FitT;
typedef eoRealParticle<FitT> Particle;
//-----------------------------------------------------------------------------

// EVALFUNC
//-----------------------------------------------------------------------------
// a simple fitness function that computes the euclidian norm of a real vector
FitT real_value(const Particle &t) {
  return std::inner_product(t.begin(), t.end(), t.begin(), double(0));
}

template <typename T> inline auto sphere(std::span<const T> t) noexcept {
  ;
  return std::inner_product(t.begin(), t.end(), t.begin(), T(0));
}

template <int Dim, int Pop_Size, int Max> inline auto sevobench_pso() noexcept {
  using namespace sevobench;
  using namespace sevobench::pso_module;
  constexpr double lb = -100;
  constexpr double ub = 100;
  constexpr double v_min = -(ub - lb) * 0.2;
  constexpr double v_max = (ub - lb) * 0.2;
  tool::rng RNG;
  population<double> pop(Pop_Size, Dim, lb, ub);
  pso_velocity<double> vec(Pop_Size, std::vector<double>(Dim));
  for (auto &v : vec)
    for (auto &e : v)
      e = RNG.rand_float(v_min, v_max);
  evolutionary_algorithm alg((Max + 1) * Pop_Size);
  alg.set_max_iterator(Max);
  auto pso = pso_algorithm_builder<double>()
                 .topology(std::make_unique<gbest_topology<double>>())
                 .update(std::make_unique<inertia_weight_update<double>>())
                 .constraint_handler(pso_constraint<double>{
                     .vc = std::make_unique<spso_velocity_constraint<double>>(
                         v_min, v_max)})
                 .build();
  pso.run(pop, vec, sphere<double>, lb, ub, alg);
  return pso.topology()->best_value();
}

template <int Dim, int Pop_Size, int Max> inline auto paradiseo_pso() noexcept {
  // PARAMETRES
  // all parameters are hard-coded!

  constexpr auto MAX_GEN = Max;
  constexpr auto VEC_SIZE = Dim;
  constexpr auto POP_SIZE = Pop_Size;

  constexpr double POS_INIT_MIN = -100;
  constexpr double POS_INIT_MAX = 100;

  constexpr double VELOCITY_INIT_MIN = -0.2 * (POS_INIT_MAX - POS_INIT_MIN);
  constexpr double VELOCITY_INIT_MAX = 0.2 * (POS_INIT_MAX - POS_INIT_MIN);

  constexpr double VELOCITY_MIN = -0.2 * (POS_INIT_MAX - POS_INIT_MIN);
  constexpr double VELOCITY_MAX = 0.2 * (POS_INIT_MAX - POS_INIT_MIN);

  constexpr double INERTIA = 0.4;
  constexpr double LEARNING_FACTOR1 = 2;
  constexpr double LEARNING_FACTOR2 = 2;

  //////////////////////////
  //  RANDOM SEED
  //////////////////////////
  // reproducible random seed: if you don't change SEED above,
  //  you'll aways get the same result, NOT a random run
  std::random_device rd;
  rng.reseed(rd());

  /// SWARM
  // population <=> swarm
  eoPop<Particle> pop;

  /// EVALUATION
  // Evaluation: from a plain C++ fn to an EvalFunc Object
  eoEvalFuncPtr<Particle, FitT, const Particle &> eval(real_value);

  ///////////////
  /// TOPOLOGY
  //////////////
  // linear topology
  eoStarTopology<Particle> topology;

  /////////////////////
  // INITIALIZATION
  ////////////////////
  // position initialization
  eoUniformGenerator<double> uGen(POS_INIT_MIN, POS_INIT_MAX);
  eoInitFixedLength<Particle> random(VEC_SIZE, uGen);
  pop.append(POP_SIZE, random);

  // velocities initialization component
  eoUniformGenerator<double> sGen(VELOCITY_INIT_MIN, VELOCITY_INIT_MAX);
  eoVelocityInitFixedLength<Particle> veloRandom(VEC_SIZE, sGen);

  // first best position initialization component
  eoFirstIsBestInit<Particle> localInit;

  // Create an eoInitialier that:
  // 		- performs a first evaluation of the particles
  //  	- initializes the velocities
  //  	- the first best positions of each particle
  // 		- setups the topology
  eoInitializer<Particle> fullInit(eval, veloRandom, localInit, topology, pop);

  // Full initialization here to be able to print the initial population
  // Else: give the "init" component in the eoEasyPSO constructor
  fullInit();

  ///////////////
  /// VELOCITY
  //////////////
  // Create the bounds for the velocity not go to far away
  eoRealVectorBounds bnds(VEC_SIZE, VELOCITY_MIN, VELOCITY_MAX);

  // the velocity itself that needs the topology and a few constants
  eoModifyStandardVelocity<Particle> velocity(
      topology, INERTIA, LEARNING_FACTOR1, LEARNING_FACTOR2, bnds);

  ///////////////
  /// FLIGHT
  //////////////
  // flight
  eoStandardFlight<Particle> flight;

  ////////////////////////
  /// STOPPING CRITERIA
  ///////////////////////
  // the algo will run for MAX_GEN iterations
  eoGenContinue<Particle> genCont(MAX_GEN);

  // GENERATION
  /////////////////////////////////////////
  // the algorithm
  ////////////////////////////////////////
  // standard PSO requires
  // stopping criteria, evaluation,velocity, flight

  eoEasyPSO<Particle> pso(genCont, eval, velocity, flight);

  // Apply the algo to the swarm - that's it!
  rng.reseed(rd());
  pso(pop);

  // OUTPUT
  // Print (sorted) intial population
  return topology.globalBest().fitness();
}

// A main that catches the exceptions

int main(int argc, char **argv) {
  constexpr int Dim = 30;
  constexpr int Pop_Size = 30;
  constexpr int Max = 1000;
  constexpr int Runs = 10;
  double fit[Runs];
  double time[Runs];
  double fit1[Runs];
  double time1[Runs];
  for (int i = 0; i < Runs; i++) {
    {
      auto t0 = std::chrono::steady_clock::now();
      fit[i] = sevobench_pso<Dim, Pop_Size, Max>();
      auto t1 = std::chrono::steady_clock::now();
      time[i] = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                    .count();
    }
    {
      auto t0 = std::chrono::steady_clock::now();
      fit1[i] = paradiseo_pso<Dim, Pop_Size, Max>();
      auto t1 = std::chrono::steady_clock::now();
      time1[i] = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                     .count();
    }
  }
  auto ms = sevobench::tool::mean_std(fit, fit + Runs);
  auto ms1 = sevobench::tool::mean_std(fit1, fit1 + Runs);
  auto ms2 = sevobench::tool::mean_std(time, time + Runs);
  auto ms3 = sevobench::tool::mean_std(time1, time1 + Runs);
  std::printf("sevobench pso:mean:%f,std:%f,mean_time:%f\n", ms[0], ms[1],
              ms2[0]);
  std::printf("paradiseo pso:mean:%f,std:%f,mean_time:%f\n", ms1[0], ms1[1],
              ms3[0]);
  return 0;
}
//-----------------------------------------------------------------------------
