#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../src/graph.hpp"
#include "../src/utilities.hpp"

using namespace std;

void run_all_tests(string filename, unsigned int node1, unsigned int node2,
		   bool directed,unsigned int nprint, bool printgraph)
{
  ifstream graphin(filename);
  
  cout << ">>> Testing Dijkstra's like algorithms -----------------------\n\n"
       << "Graph filename is " << filename << "\n";
  if(directed)
    cout << "Loading graph as DIRECTED";
  else
    cout << "Loading graph as UNDIRECTED";
  cout << "\n\n" << flush;
  
  cout << "> Loading graph into memory ... " << endl;
  Graph g0(graphin, 1, directed, 1); // undirected weighted read from nodes list
  
  std::clog<<" ... done with loading graph! \n "<<endl;
  cout << "number of nodes is " << g0.n() << endl;
  cout << "number of edges is " << g0.m() << "\n" << endl;
  
  if(printgraph)
    g0.print_graph(cout, true);
 
  cout << "> Applying DFS to find all connected components (regardless of "
    "directionality) ... \n" << endl;
  
  vector<bool> g0explored;
  vector<unsigned int> cc0_nodes;
  std::vector<std::vector<unsigned int> > all0_ccs;
  DFS_all_connected_components(g0, all0_ccs);
  
  set<unsigned int> cc0_sizes;
  multiset<unsigned int> cc0_sizes_multi;
  for(std::vector<std::vector<unsigned int> >::iterator it =
	all0_ccs.begin(); it != all0_ccs.end(); ++it)
    {
      cc0_sizes.insert(it->size());
      cc0_sizes_multi.insert(it->size());
    }
  
  cout << "Found " << all0_ccs.size() << " connected components with "
       << cc0_sizes.size() << " different sizes. \n \n"
       << "Summary of top " << nprint << " (size, multiplicity) pairs: \n";
  
  set<unsigned int>::iterator it0 = cc0_sizes.end();
  for(unsigned int i = 0; i != nprint && i < cc0_sizes.size(); ++i)
    {
      --it0;
      cout << "(" << *it0 << ", " << cc0_sizes_multi.count(*it0) <<  ")"
	   << endl;
    }
  cout << endl;
  
  cout << "> Applying DFS to find all nodes reachable from "
       << node1 << "... \n"
       << endl;
  std::vector<unsigned int> cfrom_nodes;
  DFS_reachable_from(g0, node1, cfrom_nodes);

  cout << " Printing "
       << ((nprint < cfrom_nodes.size()) ? nprint : cfrom_nodes.size())
       << " nodes reachable: \n";
  
  for(unsigned int i = 0; i < nprint && i < cfrom_nodes.size(); ++i)
    cout << cfrom_nodes[i] << ", \t";
  cout << "\n" << endl;
  /*
  std::multiset<unsigned int> all_SCC_sizes;
  std::set<unsigned int> all_SCC_setsizes;

  SCC_Kosaraju(g0, all_SCC_sizes, all_SCC_setsizes);

  std::cout << " all_SCC_sizes.size() is "
	    << all_SCC_sizes.size() << std::endl;
  std::cout << " all_SCC_setsizes.size() is "
	    << all_SCC_setsizes.size() << std::endl;

  std::cout << "Printing out multiplicity, SCC size pairs for the top 5:"
	    << std::endl;
  short int i{0};
  for(std::set<unsigned int>::iterator it = --all_SCC_setsizes.end();
      it != --all_SCC_setsizes.begin() && i != 5; --it)
    {
      std::cout << all_SCC_sizes.count(*it) << " of size "<< *it<< std::endl;
      ++i;
    }*/
}

/*bool comp(std::pair<unsigned int , unsigned int> p1,
	  std::pair<unsigned int , unsigned int> p2)
{ // WARNING: Here I invert the std heap behaviour to keep the smallest value
  // at the top pf the heap!
  return p1.second > p2.second;
}//*/

int main(){

    
  /*cout << "WARNING: you may have to set $ulimit -s 80000, for this to work!"
    "type Y and enter to continue!"
       << endl;
  string dummy;
  cin >> dummy;*/

  cout << "\n*************************************"
       << "\n***                               ***"
       << "\n*** Tests of dijkstra's algorithm ***"
       << "\n***                               ***"
       << "\n***                               ***"
       << "\n*************************************\n\n"
       << flush;

  // Small connected graph
  run_all_tests("tests/graphs/dijkstra1.txt", 1, 4, false, 6, true);
  
  //run_all_tests("tests/graphs/dijkstraData.txt", 1, 4, false, 200, true);
  
  /*/ Small disconnected graph
  run_all_tests("tests/graphs/graph1_cut1.txt", 1, 4, false, 6, true);

  // Slightly larger connected graph
  run_all_tests("tests/graphs/graph2_dir_edges.txt", 7, 2, false, 13, true);

  // Slightly larger connected graph considering directed edges
  run_all_tests("tests/graphs/graph2_dir_edges.txt", 5, 11, true, 13, true);
  
  // Slightly larger connected graph -- disconnected version
  run_all_tests("tests/graphs/graph2_cut1.txt", 6, 13, false, 13, true);

  // Large graph
  run_all_tests("tests/graphs/SCC_graph_dir_edges.txt",
		1, 501,false, 50, false);*/
  
}
