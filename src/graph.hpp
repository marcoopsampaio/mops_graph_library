#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <list>
#include <limits>
#include <sstream>
#include <cmath>
#include <set>
#include <queue>
#include <stack>
#include <numeric>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Classes & Types /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class Edge{
public:
  unsigned int first,second;
  const bool directed;
  Edge(unsigned int first_in,unsigned int second_in):first{first_in},
						     second{second_in},
						     directed{false}{};
  Edge(unsigned int first_in,
       unsigned int second_in, bool directed_in):first{first_in},
						     second{second_in},
						     directed{directed_in}{};
};

typedef std::vector<std::list<unsigned int> > NodesVec;
typedef std::vector<Edge> EdgesVec;

class Graph
{
  NodesVec nodes_; // Nodes in adjacency list
  EdgesVec edges_; // Edges in adjacency list
  // Map to translate internal node representation to external
  std::vector<unsigned int> id_nodes_;
  std::map<unsigned int,unsigned int> index_of_ids_;

  // Utilities for constructor
  void stream_read0_graph_by_edges(std::istream & in_stream, bool directed = 0);
  void stream_read1_graph_by_nodes(std::istream & in_stream, bool directed = 0);
  
public:
  // Constructors
  Graph(){}; // start empty (to add nodes and edges in algorithms)
  Graph(std::istream & in_stream, bool input_type = 0,
	bool directed = 0); // from stream
  ////////Graph(EdgesVec & edges_in, bool directed = 0); // from vector of edges
  
  // Reset the graph to an empty state
  void clear();
  
  // retrieve node id
  unsigned int id_of_index(unsigned int index){return id_nodes_[index];};  
  // retrieve index of given node id
  unsigned int index_of_id(unsigned int id){return index_of_ids_.at(id);};
  // retrieve number of nodes
  unsigned int n(){return nodes_.size();};
  // retrieve number of nodes
  unsigned int m(){return edges_.size();};
  
  // Print out utilities
  void print_graph(std::ostream & os, bool internal = 0);

  
  friend void BFS_connected(Graph & gin, unsigned int start_node,
			    std::vector<bool> & explored,
			    std::vector<unsigned int> & cc_nodes);
  friend void BFS_reachable_from(Graph & gin,  unsigned int start_node,
				 std::vector<unsigned int> & cfrom_nodes);
  friend unsigned int BFS_distances_from_source(Graph & gin,
						unsigned int start_node,
						std::vector<unsigned int>
						& dists,
						std::vector<unsigned int>
						& cc_nodes);
  
  friend unsigned int BFS_shortest_path_length(Graph & gin,
					       unsigned int start_node,
					       unsigned int end_node);

  friend void DFS_connected_internal_ids_(Graph & gin, unsigned int start_node,
			    std::vector<bool> & explored,
			    std::vector<unsigned int> & cc_nodes);
  
  friend void DFS_reachable_from_internal_ids_(Graph & gin,
					       unsigned int start_node,
					       std::vector<bool> & explored,
					       std::vector<unsigned int>
					       & cfrom_nodes);

  friend void DFS_Kosaraju_reversed_(Graph & gin, unsigned int start_node);

  friend void DFS_Kosaraju_direct_(Graph & gin, unsigned int start_node);
};

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Algorithms ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------------
  BFS based algorithms
  ----------------------------------------------------------------------------*/

// Explore a connected component of a graph and return its number of nodes. 
// explored is either populated (if empty) or updated (if passed from
// previous calls). This ignores edge direction.
void BFS_connected(Graph & gin,unsigned int start_node,
		   std::vector<bool> & explored,
		   std::vector<unsigned int> & cc_nodes);

// Find all connected components of a graph and return a vector of vectors  of
// nodes in each subgraph
void BFS_all_connected_components(Graph & gin,
				  std::vector<std::vector<unsigned int> >
				  & all_ccs);
// Explore nodes reachable from a given node (connected component if undirected)
void BFS_reachable_from(Graph & gin,  unsigned int start_node,
			std::vector<unsigned int> & cfrom_nodes);

// Find all shortest distances from a source node
unsigned int BFS_distances_from_source(Graph & gin, unsigned int start_node,
			       std::vector<unsigned int> & dists,
			       std::vector<unsigned int> & cc_nodes);

// Find shortest distance between two nodes
unsigned int BFS_shortest_path_length(Graph & gin, unsigned int start_node,
				      unsigned int end_node);


/*------------------------------------------------------------------------------
  DFS based algorithms
  ----------------------------------------------------------------------------*/


// Explore a connected component of a graph (ignores direction of edges)
void DFS_connected(Graph & gin, unsigned int start_node,
		   std::vector<bool> & explored,
		   std::vector<unsigned int> & cc_nodes);

// Find all connected components of a graph and return a vector of vectors  of
// nodes in each subgraph
void DFS_all_connected_components(Graph & gin,
				  std::vector<std::vector<unsigned int> >
				  & all_ccs);

// Explore nodes reachable from a given node (connected component if undirected)
void DFS_reachable_from(Graph & gin,  unsigned int start_node,
			std::vector<unsigned int> & cfrom_nodes);

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// TO DO: make function that returns connected graph and then make the function
// below
// Apply DFS, find topological sorting and return true is graph is acyclic
//////bool DFS_topological_sort(Graph & gin, unsigned int start_node,
//////			  std::vector<unsigned int> & idx_topo);

/////////


/*/ Find the sizes of the top 5 strongly connected components of a graph
void SCC_Kosaraju_sizes_top5(Graph & gin, std::set<unsigned int> & dims_top5);*/

// Find all the strongly connected components of a graph
void SCC_Kosaraju(Graph & gin, std::multiset<unsigned int> & all_SCC_sizes,
		  std::set<unsigned int> & all_SCC_setsizes);

static std::vector<bool> explored1;
static std::vector<bool> explored2;
static std::vector<unsigned int> finishing_vec;
static unsigned int leader;
static unsigned int leader_count;

#endif
