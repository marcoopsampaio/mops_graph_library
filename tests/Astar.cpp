#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../src/graph.hpp"
#include "../src/utilities.hpp"

using namespace std;

void run_all_tests(string filename, unsigned int node1,
		   vector< pair<unsigned int, unsigned int > > pairs_nodes,
		   bool directed,unsigned int nprint, bool printgraph)
{
  ifstream graphin(filename);
  
  cout << ">>> Testing A* like algorithms -----------------------\n\n"
       << "Graph filename is " << filename << "\n";
  if(directed)
    cout << "Loading graph as DIRECTED";
  else
    cout << "Loading graph as UNDIRECTED";
  cout << "\n\n" << flush;
  
  cout << "> Loading graph into memory ... " << endl;
  Graph g0(graphin, 2, directed, 1); // undirected weighted read from nodes list
  
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



  cout << "> Testing now several pairs of points for A* algorithm" << std::endl;

  for(unsigned int i = 0; i != pairs_nodes.size(); ++i)
    {
      std::vector<unsigned int> path;
      double dist{Astar_shortest_path(g0, pairs_nodes[i].first,
				      pairs_nodes[i].second, path)};
      
    /*Dijkstra_shortest_paths(g0, 0, dists, paths);*/
      cout << "> Printing out A* algorithm results\n" << endl
	   << ">> Shortest path with length " << dist << " between node "
	   << pairs_nodes[i].first << " and "
	   << "target node  "<< pairs_nodes[i].second << ": ";
      
      for(unsigned int j = 0; j != path.size(); ++j)
	if(j == 0)
	  cout << path[j];
	else
	  cout << ", " << path[j] ;
      
      cout << endl;
    }

}

int main(){

    
  /*cout << "WARNING: you may have to set $ulimit -s 80000, for this to work!"
    "type Y and enter to continue!"
       << endl;
  string dummy;
  cin >> dummy;*/

  cout << "\n*************************************"
       << "\n***                               ***"
       << "\n***  Tests of A* like algorithms  ***"
       << "\n***                               ***"
       << "\n***                               ***"
       << "\n*************************************\n\n"
       << flush;

  vector<pair<unsigned int, unsigned int> > pairs_nodes;
  pairs_nodes.push_back(make_pair(0,10));
  pairs_nodes.push_back(make_pair(8,2));
  pairs_nodes.push_back(make_pair(2,8));
  // Small connected graph
  run_all_tests("tests/graphs/graph_geo1.txt", 0, pairs_nodes, false, 6, true);

  
}
