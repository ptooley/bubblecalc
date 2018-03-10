#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include "consts.hpp"

namespace pt = boost::property_tree;

struct bunchdata{
  double x, y, z, px, py, pz, t, q;
  double dx, dy, dz, dpx, dpy, dpz, dt;
  std::string xdist, ydist, zdist, pxdist, pydist, pzdist, tdist;
  size_t npart;
};

class ConfigParser{
public:
  std::vector<bunchdata> bunches;
  double eta2;
  double lambda0;
  double int_time;
  std::string outfile;
  bool radt_mode;
  size_t nsteps;

  ConfigParser(std::string filename): filename(filename){
    parse_config();
  }
    
  void parse_config(){
    try{
      pt::read_json(filename, config);
    }
    catch (const std::exception& e){
      std::ostringstream err;
      err << "Failed to parse " << filename << " check your syntax!\n";
      throw std::runtime_error(err.str());
    }

    for(const auto &it : config){
      if( it.first == "Bunch"){ handle_bunch_block(it.first, it.second);}
      else if( it.first == "Output"){ handle_output_block(it.first, it.second);}
      else if( it.first == "Bubble"){ handle_bubble_block(it.first, it.second);}
      else { handle_unknown_block(it.first, it.second);}
    }
  }

  void handle_bunch_block(std::string name, pt::ptree block){
    bunchdata currbunch;
    std::cout << "Parsing Bunch Spec:\n\n";
    try{
      currbunch.npart = block.get<size_t>("npart");
    } catch (const std::exception& e){
      std::string err = "Error: Number of particles (npart) not specified in bunch.\n";
      throw std::runtime_error(err);
    }
    try{
      currbunch.q = block.get<double>("q");
    } catch (const std::exception& e){
      currbunch.q = (double)currbunch.npart * k_e;
      std::cout << "Warning: No bunch charge specified, assuming npart * e = "
        << currbunch.q << "C\n";
    }
    currbunch.xdist = block.get("xdist", "");
    currbunch.ydist = block.get("ydist", "");
    currbunch.zdist = block.get("zdist", "");
    currbunch.pxdist = block.get("pxdist", "");
    currbunch.pydist = block.get("pydist", "");
    currbunch.pzdist = block.get("pzdist", "");
    currbunch.tdist = block.get("tdist", "");

    currbunch.x = block.get("x", 0.);
    currbunch.y = block.get("y", 0.);
    currbunch.z = block.get("z", 0.);
    currbunch.px = block.get("px", 0.);
    currbunch.py = block.get("py", 0.);
    currbunch.pz = block.get("pz", 0.);
    currbunch.t = block.get("t", 0.);

    currbunch.dpx = block.get("dpx", 0.);
    currbunch.dpy = block.get("dpy", 0.);
    currbunch.dpz = block.get("dpz", 0.);
    currbunch.dx = block.get("dx", 0.);
    currbunch.dy = block.get("dy", 0.);
    currbunch.dz = block.get("dz", 0.);
    currbunch.dt = block.get("dt", 0.);
  
    std::cout << "Read bunch parameters:"
      << "\n\tnpart = " << currbunch.npart 
      << "\n\tq = " << currbunch.q
      << "\n\tt = " << currbunch.t << "\tdt = " << currbunch.dt << "\tdist: " << currbunch.tdist
      << "\n\tx = " << currbunch.x << "\tdx = " << currbunch.dx << "\tdist: " << currbunch.xdist
      << "\n\ty = " << currbunch.y << "\tdy = " << currbunch.dy << "\tdist: " << currbunch.ydist
      << "\n\tz = " << currbunch.z << "\tdz = " << currbunch.dz << "\tdist: " << currbunch.zdist
      << "\n\tpx = " << currbunch.px << "\tdpx = " << currbunch.dpx << "\tdist: " << currbunch.pxdist
      << "\n\tpy = " << currbunch.py << "\tdpy = " << currbunch.dpy << "\tdist: " << currbunch.pydist
      << "\n\tpz = " << currbunch.pz << "\tdpz = " << currbunch.dpz << "\tdist: " << currbunch.pzdist
      << std::endl << std::endl;

    this->bunches.push_back(currbunch);
  }

  void handle_output_block(std::string name, pt::ptree block){
    try{
      this->nsteps = block.get<size_t>("nsteps");
    } catch (const std::exception& e){
      std::string err = "Error: Number of timesteps (nsteps) not specified.\n";
      throw std::runtime_error(err);
    }
    try{
      this->outfile = block.get<std::string>("filename");
    } catch (const std::exception& e){
      std::string err = "Error: filename not specified.\n";
      throw std::runtime_error(err);
    }
    try{
      this->radt_mode = block.get<bool>("radt_mode");
    } catch (const std::exception& e){
      this->radt_mode = false;
    }

    this->int_time = block.get("int_time", 1.0);

    if(this->radt_mode){
      std::cout << "Writing to " << this->outfile 
        << " using radt normalized units.";
    } else {
      std::cout << "Writing to " << filename
        << " using SI units.";
    }
    std::cout << "\n\n";
  }

  void handle_bubble_block(std::string name, pt::ptree block){
    try{
      this->eta2 = block.get<double>("eta2");
    } catch (const std::exception& e){
      std::string err = "Error: eta2 not specified.\n";
      throw std::runtime_error(err);
    }
    try{
      this->lambda0 = block.get<double>("lambda0");
    } catch (const std::exception& e){
      std::string err = "Error: lambda0 not specified.\n";
      throw std::runtime_error(err);
    }
    std::cout << "Read Bubble Parameters:"
      << "\n\teta2 = " << this->eta2
      << "\n\tlambda0 = " << this->lambda0
      << std::endl << std::endl;
  }

  void handle_unknown_block(std::string name, pt::ptree block){
    std::cout << "Warning: Unknown block" << name << "ignoring...\n";
  }

protected:
  pt::ptree config;
  std::string filename;

};

