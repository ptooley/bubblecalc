#include "bubble.hpp"

Bubble::Bubble(num _lambda, num _eta2){
  eta2 = _eta2;
  lambda_0 = _lambda;
  omega_0 = 2 * k_pi * k_c / lambda_0;
  eta = std::sqrt(eta2);
  omega_p = eta * omega_0;
  omega_p2 = eta2 * omega_0 * omega_0;
  k_p = omega_p / k_c;
  beta = std::sqrt(1 - eta2);
  gamma = 1 / eta;
};

void Bubble::preallocate_hint(size_t n_electron){
  electrons.reserve(n_electron);
}

void Bubble::add_electron(num t0, num x0, num y0, num z0, num px0, num py0, num pz0, num q0){
  size_t curr_count = electrons.size();
  electrons.push_back(Electron(this, curr_count));  
  electrons.back().initial_conditions(t0, x0, y0, z0, px0, py0, pz0, q0);
}

void Bubble::calculate_tracks(double dt, double t_end, double t_start){

  size_t write_count = 0;
  size_t total_size = electrons.size();
  size_t update_freq = 10;

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (size_t i = 0; i < electrons.size(); i++){
    if(++write_count % update_freq == 0){
      std::cout << "Progress " << write_count << "/" << total_size << "\r" << std::flush;
    }
    electrons[i].calculate_track(dt, t_end, t_start);
  }
}

herr_t Bubble::save_tracks(std::string filename, bool radt_format){
  hid_t file_h, fapl;
  herr_t status, dset_props;
  dset_props = H5P_DEFAULT;

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_libver_bounds(fapl, H5F_LIBVER_18, H5F_LIBVER_LATEST);
  file_h = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

  std::cout << "Writing \"" << filename << "\":" << std::endl;

  for(auto it = electrons.begin(); it != electrons.end(); it++){
    status = it->save_track(file_h, radt_format, dset_props);
  }
  status = H5Fclose(file_h);
  return(status);
}
    
const void Bubble::add_electron_distribution(size_t n_electron, dfn_uptr& t_fn, dfn_uptr& xi_fn, 
                                             dfn_uptr& y_fn, dfn_uptr& z_fn, dfn_uptr& px_fn, 
                                             dfn_uptr& py_fn, dfn_uptr& pz_fn, dfn_uptr& q_fn){
  
  electrons.reserve(electrons.size() + n_electron);

  std::vector<num> t0 = (*t_fn)(n_electron);
  std::vector<num> xi0 = (*xi_fn)(n_electron);
  std::vector<num> y0 = (*y_fn)(n_electron);
  std::vector<num> z0 = (*z_fn)(n_electron);
  std::vector<num> px0 = (*px_fn)(n_electron);
  std::vector<num> py0 = (*py_fn)(n_electron);
  std::vector<num> pz0 = (*pz_fn)(n_electron);
  std::vector<num> q = (*q_fn)(n_electron);

  for (size_t i = 0; i < n_electron; i++){
    add_electron(t0[i], xi0[i], y0[i], z0[i], px0[i], py0[i], pz0[i], q[i]);
  }
}
