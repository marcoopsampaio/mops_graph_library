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
  
  cout << ">>> Testing DFS algorithms -------------------------------------\n\n"
       << "Graph filename is " << filename << "\n";
  if(directed)
    cout << "Loading graph as DIRECTED";
  else
    cout << "Loading graph as UNDIRECTED";
  cout << "\n\n" << flush;
  
  cout << "> Loading graph into memory ... " << endl;
  Graph g0(graphin, 0, directed); // undirected read from edges list
  std::clog<<" ... done with loading graph! \n "<<endl;
  if(printgraph)
    g0.print_graph(cout);

  /*////////
  vector<bool> g0explored;
  vector<unsigned int> cc0_nodes;
  std::vector<std::vector<unsigned int> > all0_ccs;
  //cout << "----->" << std::numeric_limits<std::deque<unsigned int>::size_type>::max() << endl;
  //exit(0);
  DFS_connected(g0, 1, g0explored, cc0_nodes);
  //DFS_reachable_from(g0, 3, cc0_nodes);
  // BFS_all_connected_components(g0, all0_ccs);*/


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
  cout << "number of nodes is "<<g0.n()<<endl;
  cout << "number of edges is "<<g0.m()<<endl;
  /*
  cout << "> Applying BFS to find all distances from node " << node1 << "... \n"
       << endl;
  std::vector<unsigned int> dists;
  std::vector<unsigned int> cc_nodes;
  unsigned int nlayers = BFS_distances_from_source(g0, node1, dists, cc_nodes);

  cout << "Number of nodes in first "
       << ((nprint < nlayers) ? nprint : nlayers)
       << " layers: \n";

  for(unsigned int ilayer = 0, iccnodes = 0;
      ilayer < nprint && iccnodes != cc_nodes.size(); ++ilayer)
    {
      unsigned int icount = 0;
      while(iccnodes != cc_nodes.size() &&
	    dists[g0.index_of_id(cc_nodes[iccnodes])] == ilayer)
	{
	  ++icount;
	  ++iccnodes;
	}
      cout << icount << " nodes in layer " << ilayer << "\n";
    }
  cout << endl;

  cout << "Distance between nodes " << node1 << " and " << node2
       << " is " << BFS_shortest_path_length(g0, node1, node2)
       << "\n" << endl;*/
}

int main(){

  cout << "\n*************************************"
       << "\n***                               ***"
       << "\n*** Tests of DFS based algorithms ***"
       << "\n***     (undirected graphs)       ***"
       << "\n***                               ***"
       << "\n*************************************\n\n"
       << flush;

  // Small connected graph
  run_all_tests("tests/graphs/graph1_dir_edges.txt", 1, 4, false, 10, true);
  
  // Small disconnected graph
  run_all_tests("tests/graphs/graph1_cut1.txt", 1, 4, false, 10, true);

  // Slightly larger connected graph
  run_all_tests("tests/graphs/graph2_dir_edges.txt", 7, 2, false, 10, true);

  // Slightly larger connected graph considering directed edges
  run_all_tests("tests/graphs/graph2_dir_edges.txt", 5, 11, true, 10, true);
  
  // Slightly larger connected graph -- disconnected version
  run_all_tests("tests/graphs/graph2_cut1.txt", 6, 13, false, 10, true);

  // Large graph
  run_all_tests("tests/graphs/SCC_graph_dir_edges.txt",
		1, 501,false, 50, false);
  
}
