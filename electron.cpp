#include "electron.hpp"

struct int_params{
  num c_vphi;
  num c_05momega2;
  num c_05momega2vphi;
  num c_025m2omega2vphi;
  num c_05m2omega2vphi;
};

int integrator (double t, const double y[], double f[], void* params){
  int_params *p = (int_params*) params;
  double gamma = std::sqrt(1 + (std::pow(y[2],2) + std::pow(y[3],3))/k_me2c2);
  f[0] = y[2] / (gamma * k_me);
  f[1] = y[3] / (gamma * k_me);
  f[2] = - (p->c_05momega2 * (y[0] - p->c_vphi*t))  - (p->c_025m2omega2vphi * y[3] * y[1] /gamma);
  f[3] = - p->c_05momega2 * y[1] - p->c_05m2omega2vphi * y[2] * y[1] / gamma;
  return(GSL_SUCCESS);
}

int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params){
  int_params *p = (int_params*)params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
  gsl_matrix *m = &dfdy_mat.matrix;
  double gamma = std::sqrt(1 + std::pow(y[2],2) + std::pow(y[3],3));

  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 0.0);
  gsl_matrix_set (m, 0, 2, k_mec2/gamma);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 1, 0, 0.0);
  gsl_matrix_set (m, 1, 1, 0.0);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, k_mec2/gamma);
  gsl_matrix_set (m, 2, 0, -p->c_05momega2);
  gsl_matrix_set (m, 2, 1, -p->c_025m2omega2vphi * y[3]/gamma);
  gsl_matrix_set (m, 2, 2, 0.0);
  gsl_matrix_set (m, 2, 3, -p->c_025m2omega2vphi * y[1]/gamma);
  gsl_matrix_set (m, 3, 0, 0.0);
  gsl_matrix_set (m, 3, 1, -p->c_05momega2 - p->c_05m2omega2vphi * y[2]/gamma);
  gsl_matrix_set (m, 3, 2, -p->c_05m2omega2vphi * y[1]/gamma);
  gsl_matrix_set (m, 3, 3, 0.0);

  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = p->c_05momega2vphi;
  dfdt[3] = 0.0;

  return(GSL_SUCCESS);
}


void Electron::initial_conditions(num t0, num x0, num y0, num z0, num px0, num py0, num pz0,
                                  num q0){
  num gamma0 = std::sqrt( 1 + (px0*px0 + py0*py0 + pz0*pz0)/k_me2c2);
  t[0] = t0;
  x[0] = x0 + k_c*bubble->beta*t0;
  xi[0] = x0;
  y[0] = y0;
  z[0] = z0;
  px[0] = px0;
  py[0] = py0;
  pz[0] = pz0;
  gamma[0] = gamma0;
  q = q0;
}


void Electron::calculate_track(size_t calc_steps){
  const size_t ndims = 4;
  const double relerr = 0.0;
  const double abserr = 1e-6;
	enlarge_as_needed(calc_steps);

  num omega_p2 = bubble->omega_p2;
  num beta_b = bubble->beta;
  num gamma_b = bubble->gamma;

  num t_d = (-2.*std::pow(gamma_b,2.0) * xi[0])/k_c;
  num dt = t_d / calc_steps;

  int_params params;
  params.c_vphi = k_c * beta_b;
  params.c_05momega2 = 0.5*k_me*omega_p2;
  params.c_05momega2vphi = 0.5*k_me*omega_p2* params.c_vphi;
  params.c_05m2omega2vphi = k_me*params.c_05momega2vphi;
  params.c_025m2omega2vphi = 0.5*k_me*params.c_05momega2vphi; 

  gsl_odeiv2_system sys = {integrator, jacobian, ndims, &params};
  gsl_odeiv2_driver* driver = 
    gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, dt, abserr, relerr);

  double ti = t[0];
  double t_int = ti;
  double u[4] = {x[0], y[0], px[0], py[0]};
	for (size_t s = 1; s <= calc_steps; s++){
      double tt = ti + s * (t_d - ti) / calc_steps;
      int status = gsl_odeiv2_driver_apply(driver, &t_int, tt, u);
      if (status != GSL_SUCCESS){
        std::cout << "Fell over! " << status <<std::endl;
      }
      t[s] = t_int;
			xi[s] = u[0] - t_int*k_c*beta_b;
			x[s] = u[0];
			y[s] = u[1];
			z[s] = 0;
			px[s] = u[2];
			py[s] = u[3];
			pz[s] = 0;
      gamma[s] = std::sqrt(1.  + std::pow(px[s]/k_mec,2.) + std::pow(py[s]/k_mec,2.) 
                 + std::pow(pz[s]/k_mec,2.));
  }
}


herr_t Electron::save_track(hid_t file_handle, bool radt_format = false,
                            herr_t dset_props = H5P_DEFAULT){
  if (radt_format) {
    radt_convert(); // these are very quick if format is already correct
  } else {
    radt_revert();
  }
  hsize_t arrsize = t.size();
  hsize_t* dsize = &arrsize;
  herr_t status = 0;
  hid_t group_h = H5Gcreate(file_handle, std::to_string(pid).c_str(),
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = write_dset(group_h, "t", dsize, t.data(), dset_props);
  status = write_dset(group_h, "x1", dsize, x.data(), dset_props);
  status = write_dset(group_h, "x2", dsize, y.data(), dset_props);
  status = write_dset(group_h, "x3", dsize, z.data(), dset_props);
  status = write_dset(group_h, "p1", dsize, px.data(), dset_props);
  status = write_dset(group_h, "p2", dsize, py.data(), dset_props);
  status = write_dset(group_h, "p3", dsize, pz.data(), dset_props);
  status = write_dset(group_h, "ene", dsize, gamma.data(), dset_props);
  status = fill_dset(group_h, "q", dsize, &q, dset_props);

  status = H5Gclose(group_h);

  return(status);
}


herr_t Electron::write_dset(hid_t group_h, const char* name, hsize_t* dsize, num* data,
                            herr_t dset_props){
    herr_t status;
    hid_t dspace_h, dset_h;

    dspace_h = H5Screate_simple(1,dsize, NULL);
    dset_h = H5Dcreate(group_h, name, H5T_NATIVE_DOUBLE, dspace_h, H5P_DEFAULT, dset_props, H5P_DEFAULT); 
    status = H5Dwrite(dset_h, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    status = H5Dclose(dset_h);
    status = H5Sclose(dspace_h);
    
    return(status);
}

herr_t Electron::fill_dset(hid_t group_h, const char* name, hsize_t* dsize, num* data,
                            herr_t dset_props){
    herr_t status;
    hid_t dspace_h, dset_h;
    std::vector<num> dbuf;
    dbuf.resize(*dsize);
    std::fill(dbuf.begin(), dbuf.end(), *data);
    dspace_h = H5Screate_simple(1,dsize, NULL);
    dset_h = H5Dcreate(group_h, name, H5T_NATIVE_DOUBLE, dspace_h, H5P_DEFAULT, dset_props, H5P_DEFAULT); 
    status = H5Dwrite(dset_h, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dbuf.data());
    status = H5Dclose(dset_h);
    status = H5Sclose(dspace_h);
    return(status);
}

void Electron::radt_convert(){
  if (!radt_format){
    radt_format = true;
    q /= k_e;
    iterable_multiply(&t, k_c*1e6);
    iterable_multiply(&x, 1e6);
    iterable_multiply(&y, 1e6);
    iterable_multiply(&z, 1e6);
    iterable_multiply(&px, 1/(k_mec));
    iterable_multiply(&py, 1/(k_mec));
    iterable_multiply(&pz, 1/(k_mec));
    iterable_add(&gamma, -1.);
  }
}

void Electron::radt_revert(){
  if (radt_format){
    radt_format = false;
    q *= k_e;
    iterable_multiply(&t, 1e-6/k_c);
    iterable_multiply(&x, 1e-6);
    iterable_multiply(&y, 1e-6);
    iterable_multiply(&z, 1e-6);
    iterable_multiply(&px, k_mec);
    iterable_multiply(&py, k_mec);
    iterable_multiply(&pz, k_mec);
    iterable_add(&gamma, 1.);
  }
}

void Electron::enlarge_as_needed(size_t new_steps){
	if (new_steps < t.size()) return;
	t.resize(new_steps);
	x.resize(new_steps);
	xi.resize(new_steps);
	y.resize(new_steps);
	z.resize(new_steps);
	px.resize(new_steps);
	py.resize(new_steps);
	pz.resize(new_steps);
	gamma.resize(new_steps);
}


void Electron::destructive_resize(size_t new_steps){
	t.resize(new_steps);
	x.resize(new_steps);
	xi.resize(new_steps);
	y.resize(new_steps);
	z.resize(new_steps);
	px.resize(new_steps);
	py.resize(new_steps);
	pz.resize(new_steps);
	gamma.resize(new_steps);
}


