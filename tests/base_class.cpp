#include <iostream>
#include <fstream>
#include <string>
#include "../src/graph.hpp"

using namespace std;

int main(){

  cout << "\n **********************************************"
       << "\n ***                                        ***"
       << "\n *** This is a test of the graph base class ***"
       << "\n ***                                        ***"
       << "\n **********************************************\n\n"
       << flush;
  
  //string filename1("tests/graphs/graph1_dir_edges.txt");
  string filename1("tests/graphs/graph3_dir_edges_small.txt");
  ifstream fgraph1(filename1);
  
  cout << "1.1 - Testing read of directed version of graph from file "
       << filename1 << endl;
  Graph g1(fgraph1,0,1);
  g1.print_graph(cout);

  cout << "1.2 - Testing read of undirected version of graph from file: "
       << filename1 << endl;
  fgraph1.clear();
  fgraph1.seekg(0);
  Graph g2(fgraph1,0,0);
  g2.print_graph(cout);
  
}
