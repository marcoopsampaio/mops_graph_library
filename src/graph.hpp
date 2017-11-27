#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <list>
#include <limits>
#include <sstream>
#include <cmath>
#include <set>
#include <queue>
#include <stack>
#include <numeric>
#include <algorithm>
#include <tuple>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Classes & Types /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//typedef std::tuple<unsigned int, unsigned int, unsigned int> Triplet;

class Triplet
{
public:
  double t0;
  unsigned int t1, t2;
  Triplet(){};
  Triplet(double t0_in, unsigned int t1_in,
	  unsigned int t2_in): t0{t0_in}, t1{t1_in}, t2{t2_in}{}; 
};

class myHeap
{
public:
  // Contains score, node number and current edge associated with score
  std::vector<Triplet > heap_data;
  // Vector to hold the location of each node in the heap vector;
  std::vector<unsigned int> pos_node_heap;
  unsigned int heapsize;

  void insert(Triplet triplet);
  Triplet extract_min();
  void remove(unsigned int pos_del);

  // Constructor
  myHeap(unsigned int n_nodes): heap_data(n_nodes),
				pos_node_heap(n_nodes,
					      std::numeric_limits<unsigned
					      int>::max()),
				heapsize{0}{};

  void print();
  
};

class Edge
{
public:
  unsigned int first, second;
  double weight;
  const bool directed;
  Edge(unsigned int first_in,unsigned int second_in):
    first{first_in}, second{second_in}, weight{1.}, directed{false}{};
  
  Edge(unsigned int first_in,unsigned int second_in, bool directed_in):
    first{first_in}, second{second_in}, weight{1.}, directed{directed_in}{};
  
  Edge(unsigned int first_in,unsigned int second_in, double weight_in):
    first{first_in},second{second_in}, weight{weight_in}, directed{false}{};
  
  Edge(unsigned int first_in, unsigned int second_in, double weight_in,
       bool directed_in):first{first_in}, second{second_in}, weight{weight_in},
			 directed{directed_in}{};

};

typedef std::vector<std::list<unsigned int> > NodesVec;
typedef std::vector<std::pair<double, double> > NodesCoords;
typedef std::vector<Edge> EdgesVec;

class Graph
{
  NodesVec nodes_; // Nodes in adjacency list
  EdgesVec edges_; // Edges in adjacency list
  NodesCoords nodes_coords_; // 2D coordinates of nodes
  // Map to translate internal node representation to external
  std::vector<unsigned int> id_nodes_;
  std::unordered_map<unsigned int,unsigned int> index_of_ids_;

  // Utilities for constructor
  void stream_read0_graph_by_edges(std::istream & in_stream, bool directed = 0);
  void stream_read1_graph_by_nodes(std::istream & in_stream, bool directed = 0,
				   bool weighted = 0);
  void stream_read2_graph_by_nodes_and_edges(std::istream & in_stream,
					     bool directed = 0,
					     bool weighted = 0);
  
public:
  // Constructors
  Graph(){}; // start empty (to add nodes and edges in algorithms)
  Graph(std::istream & in_stream, char input_type = 0,
	bool directed = 0, bool weighted = 0); // from stream
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
  void print_graph(std::ostream & os, bool weighted = 0, bool internal = 0);

  
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

  friend void Dijkstra_shortest_paths(Graph & gin,  unsigned int start_node,
				      std::vector<double> & dists,
				      std::vector< std::vector<unsigned int> >
				      & paths);
  friend double Astar_shortest_path(Graph & gin,  unsigned int start_node,
				    unsigned int target_node,
				    std::vector<unsigned int> & path);
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

/*------------------------------------------------------------------------------
  Dijkstra's based algorithms
  ----------------------------------------------------------------------------*/

// Find all shortest paths from a source node and corresponding distances
void Dijkstra_shortest_paths(Graph & gin, unsigned int start_node,
			std::vector<unsigned int> & dists,
			std::vector< std::vector<unsigned int> > & paths);

double Astar_shortest_path(Graph & gin, unsigned int start_node,
			   unsigned int target_node,
			   std::vector<unsigned int> & path);

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
