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
  friend void BFS_distances_from_source(Graph & gin, unsigned int start_node,
					std::vector<unsigned int> & dists,
					std::vector<unsigned int> & cc_nodes);
  friend unsigned int BFS_shortest_path_length(Graph & gin,
					       unsigned int start_node,
					       unsigned int end_node);
};

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Algorithms ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------------
  BFS based algorithms
  ----------------------------------------------------------------------------*/

// Explore a connected component of a graph and return its number of nodes. 
// explored is either populated (if empty) or updated (if passed from
// previous calls)
void BFS_connected(Graph & gin,unsigned int start_node,
		   std::vector<bool> & explored,
		   std::vector<unsigned int> & cc_nodes);

// Find all connected components of a graph and return a vector of vectors  of
// nodes in each subgraph
void BFS_all_connected_components(Graph & gin,
				  std::vector<std::vector<unsigned int> >
				  & all_ccs);

// Find all shortest distances from a source node
void BFS_distances_from_source(Graph & gin, unsigned int start_node,
			       std::vector<unsigned int> & dists,
			       std::vector<unsigned int> & cc_nodes);

// Find shortest distance between two nodes
unsigned int BFS_shortest_path_length(Graph & gin, unsigned int start_node,
				      unsigned int end_node);


/*------------------------------------------------------------------------------
  DFS based algorithms
  ----------------------------------------------------------------------------*/


// Explore a connected component of a graph and return it as new graph
void DFS_connected(Graph & gin, Graph & gout, unsigned int start_node);

// Find all connected components of a graph and return as a vector of graphs
void DFS_all_connected_components(Graph & gin, std::vector<Graph> & all_gout);

// Apply DFS, find topological sorting and check graphi is acyclic 
void DFS_topological_sort(Graph & gin,
			  std::vector<unsigned int> & topological_index);

// Find the sizes of the top 5 strongly connected components of a graph
void SCC_Kosaraju_sizes_top5(Graph & gin, std::set<unsigned int> & dims_top5);

// Find all the strongly connected components of a graph
void SCC_Kosaraju(Graph & gin, std::vector<Graph> & all_SCCs);

#endif
