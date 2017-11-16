#ifndef UTILITIES_H
#define UTILITIES_H
#include<vector>
#include<iostream>


template<typename T>
void print_vec(std::vector<T> & vec){
  for(typename std::vector<T>::iterator iter = vec.begin(); iter!=vec.end();++iter)
    std::cout<<*iter<<"\t";
  std::cout<<std::endl;
}


#endif
