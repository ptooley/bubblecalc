#ifndef TYPEDEFS_HPP
#define TYEPDEFS_HPP

#include <vector>
#include <memory>

#include "distribution_function.hpp"

// Forwards decls

class Bubble;
class Electron;

// Typedefs

typedef double num;
typedef DistributionFunction<num> dist_fn;
typedef std::unique_ptr<dist_fn> dfn_uptr;

typedef NormalDistribution<num> norm_dist;
typedef ConstantDistribution<num> const_dist;
typedef UniformDistribution<num> uniform_dist;
typedef Sin2Distribution<num> sin2_dist;
typedef LinearSpacedDistribution<num> linspace_dist;

#endif

