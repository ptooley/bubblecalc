#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>

template<class V, class N>
void iterable_multiply(V iterable, N operand){
  for(auto it = iterable->begin(); it != iterable->end(); it++){
//    std::cout << *it << " -> ";
    *it = *it * operand;
//    std::cout << *it << std::endl;
  }
}

template<class V, class N>
void iterable_add(V iterable, N operand){
  for(auto it = iterable->begin(); it != iterable->end(); it++){
    *it = *it + operand;
  }
}

#endif
