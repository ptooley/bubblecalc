#ifndef ELECTRON_HPP
#define ELECTRON_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "hdf5.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "utils.hpp"
#include "consts.hpp"
#include "typedefs.hpp"
#include "bubble.hpp"

class Electron{
  public:
		Electron(Bubble* _bubble, size_t _pid) : bubble(_bubble), pid(_pid){};

    void initial_conditions(num t0, num x0, num y0, num z0, num px0, num py0, num pz0, num q0);
		void calculate_track(size_t calc_steps);
    void radt_convert();
    void radt_revert();
    herr_t save_track(hid_t file_handle, bool radt_format, herr_t dset_props);

  private:
    Bubble* bubble;
    size_t pid;
    bool radt_format = false;

    num q = 0;
    std::vector<num> t = {0.0};
    std::vector<num> xi = {0.0};
    std::vector<num> x = {0.0};
    std::vector<num> y = {0.0};
    std::vector<num> z = {0.0};
    std::vector<num> px = {0.0};
    std::vector<num> py = {0.0};
    std::vector<num> pz = {0.0};
    std::vector<num> gamma = {0.0};

		void enlarge_as_needed(size_t);
		void destructive_resize(size_t);
    herr_t write_dset(hid_t group_h, const char* name, hsize_t* dsize, num* data,
                      herr_t dset_props);
    herr_t fill_dset(hid_t group_h, const char* name, hsize_t* dsize, num* data,
                      herr_t dset_props);
};

#endif
