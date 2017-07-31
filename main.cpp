#include "main.hpp"

dfn_uptr dist_fn_helper(double x, double dx, std::string dist){
  dfn_uptr dfn;
  if(dist == "uniform"){
    dfn = std::make_unique<uniform_dist>(x, dx);
  } else if(dist == "normal"){
    dfn = std::make_unique<norm_dist>(x, dx);
  } else if(dist == "constant"){
    dfn = std::make_unique<const_dist>(x);
  } else if(dist == ""){
    std::cout << "WARNING: No distribution specified, falling back to ";
    if(dx <= 0){
      std::cout << "constant distribution.\n";
      dfn = std::make_unique<const_dist>(x);
    } else {
      std::cout << "normal distribution.\n";
      dfn = std::make_unique<norm_dist>(x, dx);
    }
  } else {
    std::ostringstream err;
    err << "Unknown distribution type " << dist << "check your config file.\n";
    throw std::runtime_error(err.str());
  } 
  return(dfn);
}
    

int main(int argc, char* argv[]){
  
  std::vector<std::string> args(argv, argv + argc); //nifty pointer arithmetic
  
  if(argc < 2){
    std::cout << "Usage: " << args[0] << " configfile\n";
    return(0);
  }

  std::string config_file(args[1]);

  ConfigParser* config;

  try{
    config = new ConfigParser(config_file);
  }
  catch (const std::exception &e){
    std::cout << "Fatal: " << e.what() << std::endl;
    return(-1);
  }

  Bubble bubble = Bubble(config->lambda0, config->eta2);

  for(const auto it : config->bunches){
    dfn_uptr t_dist = dist_fn_helper(it.t, it.dt, it.tdist);
    dfn_uptr xi_dist = dist_fn_helper(it.x, it.dx, it.xdist);
    dfn_uptr y_dist = dist_fn_helper(it.y, it.dy, it.ydist);
    dfn_uptr z_dist = dist_fn_helper(it.z, it.dz, it.zdist);
    dfn_uptr px_dist = dist_fn_helper(it.px, it.dpx, it.pxdist);
    dfn_uptr py_dist = dist_fn_helper(it.py, it.dpy, it.pydist);
    dfn_uptr pz_dist = dist_fn_helper(it.pz, it.dpz, it.pzdist);
    dfn_uptr q_dist = dist_fn_helper(it.q/it.npart, 0., "uniform");
  
    bubble.add_electron_distribution(it.npart, t_dist, xi_dist, y_dist, z_dist, px_dist,
                                     py_dist, pz_dist, q_dist);
  }

  bubble.calculate_tracks(config->nsteps);
  bubble.save_tracks(config->outfile, config->radt_mode);

  std::cout << "Calculation complete." << std::endl;

  return(0);
}
