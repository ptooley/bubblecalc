#ifndef DIST_FN_HPP
#define DIST_FN_HPP

#include <vector>
#include <random>
#include <iostream>

template <class T>
class DistributionFunction{
  public:

  DistributionFunction(){
    gen = std::mt19937(rd());
  }

  //private:
  virtual std::vector<T> operator()(size_t n_values) = 0;

  protected:
  static constexpr double k_pi = 3.141592653589793238462643383279502884;
  std::random_device rd;
  std::mt19937 gen;

};


template <class T>
class NormalDistribution : public DistributionFunction<T>{
  public:
  NormalDistribution(T center, T width): c(center), w(width){
    rng = std::normal_distribution<T>(center, width);
  }

  std::vector<T> operator()(size_t n_values){
    std::vector<T> dist_data;
    dist_data.resize(n_values);
    for (size_t i = 0; i < n_values; i++){
      dist_data[i] = rng(this->gen);
    }
    return(dist_data);
  }
   
  private:
  T c, w, denom;
  std::normal_distribution<T> rng;
};


template <class T>
class ConstantDistribution : public DistributionFunction<T>{
  public:
  ConstantDistribution(T constant): c(constant){};

  std::vector<T> operator()(size_t n_values){
    std::vector<T> dist_data;
    dist_data.resize(n_values);
    std::fill(dist_data.begin(), dist_data.end(), c);
    return(dist_data);
  }
   
  private:
  T c;
};

template <class T>
class UniformDistribution: public DistributionFunction<T>{
  public:
  UniformDistribution(T center, T width): c(center), w(width){
    a = c - w/2.;
    b = c + w/2.;
    rng = std::uniform_real_distribution<T>(a,b);
  }

  std::vector<T> operator()(size_t n_values){
    std::vector<T> dist_data;
    dist_data.resize(n_values);
    for (size_t i = 0; i < n_values; i++){
      dist_data[i] = rng(this->gen);
    }
    return(dist_data);
  }
   
  private:
  T c, w, a, b;
  std::uniform_real_distribution<T> rng;
};

template <class T>
class Sin2Distribution: public DistributionFunction<T>{
  public:
  Sin2Distribution(T center, T width): c(center), w(width){
    a = c - w/2.;
    b = c + w/2.;
    rng_x = std::uniform_real_distribution<T>(a,b);
    rng_y = std::uniform_real_distribution<T>(0.0,1.0);
  }

  std::vector<T> operator()(size_t n_values){
    std::vector<T> dist_data;
    dist_data.resize(n_values);
    T xpoint, ypoint, ythres;
    for (size_t i = 0; i < n_values; i++){
      do{
        xpoint = rng_x(this->gen);
        ypoint = rng_y(this->gen);
        ythres = 0.5 + 0.5*cos(2*M_PI*(xpoint - this->c)/this->w);
      }while(ypoint > ythres);
      dist_data[i] = xpoint;
    }
    return(dist_data);
  }
   
  private:
  T c, w, a, b;
  std::uniform_real_distribution<T> rng_x;
  std::uniform_real_distribution<T> rng_y;
};
#endif
