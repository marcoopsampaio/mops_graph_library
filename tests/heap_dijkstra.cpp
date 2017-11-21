#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../src/graph.hpp"
#include "../src/utilities.hpp"

using namespace std;

int main(){

  cout << "\n***************************************************"
       << "\n**                                               **"
       << "\n** Testing adapted heap for Dijkstra's algorithm **"
       << "\n**                                               **"
       << "\n***************************************************\n" << endl;
  
  // Test of myHeap

  unsigned int dmax{5};
  cout << "Creating heap with max size " << dmax << endl;
  
  myHeap h1(dmax);
  
  cout << "\nInsert 1..." << endl;
  h1.insert(std::tuple<unsigned int, unsigned int, unsigned int>(7,0,3));
  h1.print();
  
  cout << "\nInsert 2..." << endl;
  h1.insert(std::tuple<unsigned int, unsigned int, unsigned int>(10,1,4));
  h1.print();

  cout << "\nInsert 3..." << endl;
  h1.insert(std::tuple<unsigned int, unsigned int, unsigned int>(1,2,5));
  h1.print();

  cout << "\nInsert 4..." << endl;
  h1.insert(std::tuple<unsigned int, unsigned int, unsigned int>(5,3,6));
  h1.print();

  cout << "\nInsert 5..." << endl;
  h1.insert(std::tuple<unsigned int, unsigned int, unsigned int>(6,4,6));
  h1.print();

  cout << "\nRemove 1..." << endl;
  h1.remove(0);
  h1.print();

  cout << "\nRemove 2..." << endl;
  h1.remove(3);
  h1.print();

  cout << "\nPost remove insert..." << endl;
  h1.insert(std::tuple<unsigned int, unsigned int, unsigned int>(3,1,7));
  h1.print();

  cout << "\nExtract min..." << endl;
  std::tuple<unsigned int, unsigned int, unsigned int>
    mintuple = h1.extract_min();
  cout << "mintuple" << endl;
  cout << std::get<0>(mintuple)<<endl;
  cout << std::get<1>(mintuple)<<endl;
  cout << std::get<2>(mintuple)<<endl << endl;
  
  h1.print();


}
