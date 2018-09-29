#include "distribution_function.hpp"

int main(int argc, char* argv[]){
  NormalDistribution<double> foo(0,1);

  std::vector<double> bar = foo((size_t)100);

  for (auto it = bar.begin(); it != bar.end(); it++){
    std::cout << *it << std::endl;
  }


  return(0);
}
