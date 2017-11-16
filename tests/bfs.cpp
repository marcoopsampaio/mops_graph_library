#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../src/graph.hpp"
#include "../src/utilities.hpp"

using namespace std;

void run_all_tests(string filename, unsigned int node1, unsigned int node2,
		   unsigned int nprint, bool printgraph)
{
  ifstream graphin(filename);
  
  cout << ">>> Testing BFS algorithms -------------------------------------\n\n"
       << "Graph filename is " << filename << "\n\n"
       << flush;
  cout << "> Loading graph into memory ... " << endl;
  Graph g0(graphin); // undirected read from edges list
  std::clog<<" ... done with loading graph! \n "<<endl;
  if(printgraph)
    g0.print_graph(cout);
  
  cout << "> Applying BFS to find all connected components ... \n" << endl;
  
  vector<bool> g0explored;
  vector<unsigned int> cc0_nodes;
  std::vector<std::vector<unsigned int> > all0_ccs;
  BFS_all_connected_components(g0, all0_ccs);

  set<unsigned int> cc0_sizes;
  multiset<unsigned int> cc0_sizes_multi;
  for(std::vector<std::vector<unsigned int> >::iterator it = all0_ccs.begin();
      it != all0_ccs.end(); ++it){
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
  
  cout << "> Applying BFS to find all distances from a source node ... \n"
       << endl;
  std::vector<unsigned int> dists;
  std::vector<unsigned int> cc_nodes;
  BFS_distances_from_source(g0, g0.id_of_index(0), dists, cc_nodes);

  cout << "Number of nodes in first " << nprint << " layers: \n";

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
       << " is " << BFS_shortest_path_length(g0, node1, node2) << endl;
}

int main(){

  cout << "\n*************************************"
       << "\n***                               ***"
       << "\n*** Tests of BFS based algorithms ***"
       << "\n***     (undirected graphs)       ***"
       << "\n***                               ***"
       << "\n*************************************\n\n"
       << flush;

  // Small connected graph
  run_all_tests("tests/graphs/graph1_dir_edges.txt", 1, 4, 5, true);
  
  // Small disconnected graph
  run_all_tests("tests/graphs/graph1_cut1.txt", 1, 4, 5, true);

  // Slightly larger connected graph
  run_all_tests("tests/graphs/graph2_dir_edges.txt", 7, 2, 5, true);
  
  // Slightly larger connected graph -- disconnected version
  run_all_tests("tests/graphs/graph2_cut1.txt", 6, 13, 5, true);

  // Large graph
  run_all_tests("tests/graphs/SCC_graph_dir_edges.txt", 1, 501, 100, false);


  

  /*string filename; // to store filename for each graph
  ifstream graphin;

  filename = "tests/graphs/graph1_dir_edges.txt";
  graphin.close();
  graphin.open(filename);
  cout << ">>> Testing BFS algorithms on a small connected graph \n\n"
       << "Graph filename is " << filename << "\n\n"
       << flush;
  cout << "> Loading graph into memory ... " << endl;
  Graph g0(graphin); // undirected read from edges list
  std::clog<<" ... done with loading graph! \n "<<endl;
  g0.print_graph(cout);
  
  cout << "\n> Applying BFS to find all connected components ... \n" << endl;
  
  vector<bool> g0explored;
  vector<unsigned int> cc0_nodes;
  std::vector<std::vector<unsigned int> > all0_ccs;
  BFS_all_connected_components(g0, all0_ccs);

  set<unsigned int> cc0_sizes;
  multiset<unsigned int> cc0_sizes_multi;
  for(std::vector<std::vector<unsigned int> >::iterator it = all0_ccs.begin();
      it != all0_ccs.end(); ++it){
    cc0_sizes.insert(it->size());
    cc0_sizes_multi.insert(it->size());
  }
  
  cout << "Found " << all0_ccs.size() << " connected components with "
       << cc0_sizes.size() << " different sizes. "
       << endl <<  "(Size, number of cc graphs) " << endl;

  set<unsigned int>::iterator it0 = cc0_sizes.end();
  for(unsigned int i = 0; i != 50 && i < cc0_sizes.size(); ++i)
    {
      --it0;
      cout << "(" << *it0 << ", " << cc0_sizes_multi.count(*it0) <<  ")" << endl;
    }*/
  
  //vector<bool> g1explored;
  //vector<unsigned int> cc_nodes;

  
  /*exit(0);

  cout << ">>> Testing BFS algorithms on a small disconnected graph \n\n"
       << "> Graph filename is " << "\n\n"
       << flush;


  //string filename1("tests/graphs/graph1_cut1.txt");
  
  string filename1("tests/graphs/SCC_graph_dir_edges.txt");  
  //string filename1("tests/graphs/graph_test5_scc.txt");
  ifstream fgraph1(filename1);
  
  cout << "1.1 - Testing bfs on undirected version of graph from following file"
       << filename1 << endl;
  std::clog << "\nLoading graph into memory ... \n";
  Graph g1(fgraph1,0,0);
  vector<bool> g1explored;
  vector<unsigned int> cc_nodes;
  //g1.print_graph(cout);
  std::clog<<" Done with loading graph! \n "<<endl;*/


  /*BFS_connected(g1, 1, g1explored, cc_nodes);
  print_vec(cc_nodes);
  BFS_connected(g1, 6, g1explored, cc_nodes);
  print_vec(cc_nodes);
  BFS_connected(g1, 9, g1explored, cc_nodes);
  print_vec(cc_nodes);*/
  //////////gout.print_graph(cout);

  /*std::vector<std::vector<unsigned int> > all_ccs;
  BFS_all_connected_components(g1, all_ccs);

  
  set<unsigned int> cc_sizes;
  multiset<unsigned int> cc_sizes_multi;
  for(std::vector<std::vector<unsigned int> >::iterator it = all_ccs.begin();
      it != all_ccs.end(); ++it){
    cc_sizes.insert(it->size());
    cc_sizes_multi.insert(it->size());
  }
  
  cout << "Found " << all_ccs.size() << " connected components with "
       << cc_sizes.size() << " different sizes. "
       << endl <<  "(Size, number of cc graphs) " << endl;

  set<unsigned int>::iterator it = cc_sizes.end();
  for(unsigned int i = 0; i != 50 && i < cc_sizes.size(); ++i)
    {
      --it;
      cout << "(" << *it << ", " << cc_sizes_multi.count(*it) <<  ")" << endl;
    }
  
  exit(0);
    
  filename1 = "tests/graphs/graph2_dir_edges.txt";
  fgraph1.close();
  fgraph1.open(filename1);
  
  cout << "1.2 - Testing bfs on undirected version of graph from following file"
       << filename1 << endl;
  std::clog << "\nLoading graph into memory ... \n";
  Graph g2(fgraph1,0,0);
  vector<bool> g2explored; // Initialise 
  std::clog<<" Done with loading graph! \n "<<endl;

  BFS_connected(g2, 1, g2explored, cc_nodes);*/
  
}
