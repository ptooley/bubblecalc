#ifndef BUBBLE_HPP
#define BUBBLE_HPP

#include <cmath>
#include <string>
#include <vector>
#include <iostream>

#include "hdf5.h"

#include "typedefs.hpp"
#include "electron.hpp"

class Bubble{
  friend class Electron;
  public:
    Bubble(num lambda_0, num eta2);

    void preallocate_hint(size_t n_electron);
    void add_electron(num t0, num x0, num y0, num z0, num px0, num py0, num pz0, num q0);
    void calculate_tracks(double dt, double t_end, double t_start);
    herr_t save_tracks(const std::string &filename, bool radt_format);
    const void add_electron_distribution(size_t n_electron, dfn_uptr& t_fn,
                                         dfn_uptr& xi_fn, 
                                         dfn_uptr& y_fn, 
                                         dfn_uptr& z_fn,
                                         dfn_uptr& px_fn, 
                                         dfn_uptr& py_fn, 
                                         dfn_uptr& pz_fn,
                                         dfn_uptr& q_fn);

  const num get_gamma(){ return(gamma);}
  protected:
    num lambda_0, omega_0;
    num eta2, eta;
    num omega_p2, omega_p, k_p;
    num beta, gamma;
    std::vector<Electron> electrons;

};

#endif
