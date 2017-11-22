#include "graph.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Classes & Types /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/*------------------------------------------------------------------------------
  Graph class
  ----------------------------------------------------------------------------*/

//////////////////
// Constructors //
//////////////////


// Constructor from stream
Graph::Graph(std::istream & in_stream, bool input_type, bool directed,
	     bool weighted)
{
  if(input_type)
    this->stream_read1_graph_by_nodes(in_stream, directed, weighted);
  else
    this->stream_read0_graph_by_edges(in_stream, directed);
}


/////////////////////
// Graph Utilities //
/////////////////////

// Read graph from list of edges
void Graph::stream_read0_graph_by_edges(std::istream & in_stream,
					bool directed)
{ 
  if(!in_stream.good())
    {
      std::cerr << "error: from void Graph::stream_read0_graph_by_edges("
	"std::istream & in_stream, bool directed): bad in_stream, please"
	"check your filename. Crashing program in error." << std::endl;
      exit(EXIT_FAILURE);
    }
  
  
  std::string line; // string to read lines of file one by one
  std::set<unsigned int> nodes_set; // auxiliary set to re-map node numbers
  
  // Read in all edges
  while(std::getline(in_stream,line))
    {
      std::stringstream ss(line); // current edge info
      unsigned int i,j; // integers to read node indices
      if(ss >> i && ss >> j) // read pair of nodes
	{ 
	  this->edges_.push_back(Edge(i,j,directed)); // Add edge
	  // If not in map add to map and then increment
	  if(nodes_set.find(i) == nodes_set.end())
	    nodes_set.insert(i);
	  if(nodes_set.find(j) == nodes_set.end())
	    nodes_set.insert(j);
	}
    }
  
  // Populate indexing variables
  for(std::set<unsigned int>::iterator it = nodes_set.begin();
      it != nodes_set.end(); ++it)
    {
      index_of_ids_[*it] = id_nodes_.size();
      id_nodes_.push_back(*it);
    }
  
  // Set internal nodes vector to correct size
  nodes_.resize(nodes_set.size());
  
  // Populate nodes_  and re-label edges
  for(EdgesVec::size_type i = 0; i != edges_.size(); ++i)
    {
      // Re-numbering nodes in edges to range from 0 to nodes_set.size()-1
      edges_[i].first = index_of_ids_[edges_[i].first];
      edges_[i].second = index_of_ids_[edges_[i].second];
      // Add edge to the two nodes in the edge
      nodes_[edges_[i].first].push_back(i);
      nodes_[edges_[i].second].push_back(i);
    }
  
}


// Read graph from list of nodes
void Graph::stream_read1_graph_by_nodes(std::istream & in_stream,
					bool directed, bool weighted)
{
  /*******************************************
    Variables to read graph from input stream
  ********************************************/
  std::map<unsigned int, std::list< std::pair<unsigned int, unsigned int> > >
    graph_in_by_nodes; //input graph
  std::string line; // string to read lines of file one by one
  unsigned int i; // vertex raw index

  /*******************************************
    Read the graph
  ********************************************/
  while(std::getline(in_stream,line))
    {
      std::stringstream ss(line); // current vertex info
      unsigned int j, weight; // integer to read incident vertex in edge
      char dummy;
      ss >> i; // read current vertex index
      // list of vertices reached from i0
      std::list<std::pair<unsigned int, unsigned int> > jlist; 
      while(ss >> j)
	{
	  ss >> dummy;
	  ss >> weight;
	  // read vertices into list
	  jlist.push_back(std::pair<unsigned int, unsigned int>(j,weight)); 
	}
      graph_in_by_nodes[i] = jlist; //add line to graph
    }
  
  // new index to run from zero to graph_in_by_nodes.size() -1
  unsigned int i0{0}; 
  std::map<unsigned int, unsigned int> indx, indx0; // map to save standardized
  // index
  
  /*******************************************
    Re-map indices
  ********************************************/
  for(std::map<unsigned int,
	std::list<std::pair<unsigned int, unsigned int> > >::iterator
	it = graph_in_by_nodes.begin(); it != graph_in_by_nodes.end(); ++it)
    {
      id_nodes_.push_back(it->first);
      index_of_ids_[it->first] = i0;
      ++i0;
    }
  this->nodes_.clear(); // Make sure nodes_ is empty.
  this->nodes_.resize(i0); // prepare nodes data structure to hold edges
  
  /*******************************************
    Create maps representation of the graph
    and check that the graph is consistent using 
    the adjacency matrix
  ********************************************/

  // Adjacency matrix initialised to zero
  std::vector< std::vector<unsigned int> > 
    Adj(graph_in_by_nodes.size(),
	std::vector<unsigned int>(graph_in_by_nodes.size(),0));

  unsigned int iedge = 0;
  for(unsigned int i0 = 0; i0 != graph_in_by_nodes.size(); ++i0)
    {
      for(std::list< std::pair<unsigned int, unsigned int> >::iterator
	    it = graph_in_by_nodes.at(id_nodes_[i0]).begin();
	  it != graph_in_by_nodes.at(id_nodes_[i0]).end(); ++it)
	{
	    
	  if(index_of_ids_[it->first] > i0) // avoid double counting of edges
	    {
	      this->edges_.push_back(Edge(i0, index_of_ids_[it->first], it->second,
					  directed)); // Add edge
	      nodes_[i0].push_back(iedge); // add edge to first node
	      // add edge to second node
	      nodes_[index_of_ids_[it->first]].push_back(iedge); 
	      ++iedge; // next edge
	    }
	  Adj[i0][index_of_ids_[it->first]] += 1; // add to adjacency matrix
	}
    }
  
  // Check for inconsistencies
  for(unsigned int i = 0; i != Adj.size();++i)
    for(unsigned int j = i+1; j != Adj.size(); ++j)
      if(Adj[i][j] != Adj[j][i]) // Exit program in error 
	{
	  std::cerr << std::endl
		    << "Error from void Graph::stream_read1_graph_by_nodes"
	    "(std::istream & in_stream, bool directed, bool weighted):"
		    << " Inconsistent adjacency matrix."<<std::endl;
	  exit(EXIT_FAILURE);
	}	
  
}

void Graph::clear()
{
  nodes_.clear();
  edges_.clear();
  id_nodes_.clear();
  index_of_ids_.clear();
}


////////////
// Output //
////////////

void Graph::print_graph(std::ostream & os, bool weighted, bool internal)
{
  os << " *** Printing out adjacency list data for this graph *** \n\n";
  // Print out nodes list
  os << " --- Nodes with list of edges"<<std::endl;
  for(NodesVec::size_type inode = 0; inode != nodes_.size(); ++inode)
    {
      os << " node ";
      internal ? os << inode : os << id_nodes_[inode];
      os << ": ";
      for(std::list<unsigned int>::iterator it = nodes_[inode].begin();
	  it != nodes_[inode].end(); ++it)
	os << *it <<"\t";
      os << std::endl;
    }
  
  // Print out edges list
  os << "\n --- Edges"<<std::endl;
  for(EdgesVec::size_type i = 0; i != edges_.size(); ++i)
    {
      os << " edge " << i << ": ";
      unsigned int i0, i1;
      if(internal){
	i0 = edges_[i].first;
	i1 = edges_[i].second;
      }
      else{
	i0 = id_nodes_[edges_[i].first];
	i1 = id_nodes_[edges_[i].second];
      }
      os << i0;
      edges_[i].directed ? os <<" ---> " : os <<" --- ";
      os << i1;
      weighted ? os << " | weight: " << edges_[i].weight : os ; 
      os << std::endl;
    }
  os << std::endl;
}

// Heap stuff

void myHeap::insert(Triplet triplet)
{
  heap_data[heapsize] = triplet; // add new element to heap (breaks heap)
  // store iterator pointing to element just added
  pos_node_heap[triplet.t1] = heapsize;
  //bubble up element to the correct position in heap
  unsigned int pos_child = heapsize++; // get heapsize and increase it
  int pos_parent = (pos_child+1)/2-1;
  
  while(pos_parent >= 0 && heap_data[pos_parent].t0 > heap_data[pos_child].t0)
    {
      // Swap positions
      unsigned int pos_temp = pos_node_heap[heap_data[pos_child].t1];
      
      pos_node_heap[heap_data[pos_child].t1] =
	pos_node_heap[heap_data[pos_parent].t1];
      pos_node_heap[heap_data[pos_parent].t1] = pos_temp;

      // Swap triplet data
      Triplet temp{heap_data[pos_child]};
      
      heap_data[pos_child] = heap_data[pos_parent];
      heap_data[pos_parent] = temp;
      
      // Recompute parent
      pos_child = pos_parent;
      pos_parent = (pos_child+1)/2-1;
    }
  
}

void myHeap::remove(unsigned int pos_del)
{
  
  // Swap with last element and reduce heapsize
  if(pos_del == (heapsize - 1))
    {
      // Destroy position of node in heap
      pos_node_heap[heap_data[pos_del].t1]
	= std::numeric_limits<unsigned int>::max();
      --heapsize;
    }
  else if(pos_del< heapsize)
    {
      // Destroy position of node in heap
      pos_node_heap[heap_data[pos_del].t1]
	= std::numeric_limits<unsigned int>::max();
      // Put last element in heap position of deleted element
      heap_data[pos_del] = heap_data[--heapsize];
      // correct stored position
      pos_node_heap[heap_data[heapsize].t1] = pos_del;
      unsigned int score{heap_data[pos_del].t0};
      // Bubble down if necessary to restore heap property
      unsigned int pos_child1{2*pos_del+1}, pos_child2{2*(pos_del+1)};
      unsigned int pos{pos_del};
      while(pos_child1 < heapsize)
	{
	  if(pos_child2 < heapsize)
	    {
	      unsigned int score1{heap_data[pos_child1].t0};
	      unsigned int score2{heap_data[pos_child2].t0};
	      if(score > score1 && score1 <= score2)
		{
		  //Swap with child1
		  pos_node_heap[heap_data[pos_child1].t1] = pos;
		  pos_node_heap[heap_data[pos].t1] = pos_child1;

		  // Swap triplet data
		  Triplet temp{heap_data[pos_child1]};
		  heap_data[pos_child1] = heap_data[pos];
		  heap_data[pos] = temp;
		  pos = pos_child1;
		  pos_child1 = 2*pos+1;
		  pos_child2 = 2*pos+2;
		}
	      else if (score > score2 && score2 <= score1)
		{
		  //Swap with child2
		  pos_node_heap[heap_data[pos_child2].t1] = pos;
		  pos_node_heap[heap_data[pos].t1] = pos_child2;

		  // Swap triplet data
		  Triplet temp{heap_data[pos_child2]};
		  heap_data[pos_child2] = heap_data[pos];
		  heap_data[pos] = temp;
		  pos = pos_child2;
		  pos_child1 = 2*pos+1;
		  pos_child2 = 2*pos+2;
		}
	      else
		break; // Nothing left to be done. Heap property restored!
	    }
	  else
	    {// Case with only one child
	      unsigned int score1{heap_data[pos_child1].t0};
	      if(score > score1)
		{
		  //Swap with child1
		  pos_node_heap[heap_data[pos_child1].t1] = pos;
		  pos_node_heap[heap_data[pos].t1] = pos_child1;

		  // Swap triplet data
		  Triplet temp{heap_data[pos_child1]};
		  heap_data[pos_child1] = heap_data[pos];
		  heap_data[pos] = temp;
		  pos = pos_child1;
		}
	      else
		break; // Nothing left to be done. Heap property restored!
	    }
	}
    }
  else
    {
      std::cerr << "error: from void myHeap::remove(unsigned int pos_del):"
	" attempting to delete a node not in the heap" << std::endl;
      exit(EXIT_FAILURE);
    }
  
}


Triplet myHeap::extract_min()
{
  Triplet out{heap_data[0]};
  this->remove(0);
  return out;
}

void myHeap::print()
{
  std::cout << "> Printing out data stored in this myHeap... \n" << std::endl;
  std::cout << ">> (score, node, edge) \n";
  for(unsigned int i = 0; i != heap_data.size(); ++i)
    std::cout << heap_data[i].t0 << "\t"
	      << heap_data[i].t1 << "\t"
	      << heap_data[i].t2 << std::endl;
  std::cout << ">> (node, position) \n";
  for(unsigned int i = 0; i != heap_data.size(); ++i)
    std::cout << i << "\t"
	 << pos_node_heap[i] <<  std::endl;
    
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Algorithms ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------------
  BFS based algorithms
  ----------------------------------------------------------------------------*/

// Explore a connected component of a graph and return its number of nodes. 
// explored is either populated (if empty) or updated (if passed from
// previous calls)
void BFS_connected(Graph & gin,  unsigned int start_node,
		   std::vector<bool> & explored,
		   std::vector<unsigned int> & cc_nodes)
{
  std::queue<unsigned int> q; // auxiliary queue for algorithm
  q.push(gin.index_of_id(start_node));
  
  cc_nodes.clear(); // clear list of nodes in the connected component
  cc_nodes.push_back(start_node); // add starting node to connected component
  
  if(explored.size() == 0) // prepared the "explored" vector
    explored.resize(gin.n(),false);
  explored[q.front()]=true; // mark start node explored
  
  while(q.size() != 0)
    {
      unsigned int inow = q.front();
      q.pop();
      // Take first in queue and loop over all its edges
      for(std::list<unsigned int>::iterator it = gin.nodes_[inow].begin();
	  it != gin.nodes_[inow].end(); ++it)
	{
	  if(gin.edges_[*it].first == inow && !explored[gin.edges_[*it].second])
	    {
	      unsigned int iadd = gin.edges_[*it].second;
	      explored[iadd] = true;
	      q.push(iadd);
	      cc_nodes.push_back(gin.id_of_index(iadd));
	    }
	  else if(gin.edges_[*it].second == inow
		  && !explored[gin.edges_[*it].first])
	    {
	      unsigned int iadd = gin.edges_[*it].first;
	      explored[iadd] = true;
	      q.push(iadd);
	      cc_nodes.push_back(gin.id_of_index(iadd));
	    }
	    
	}  
    }
  
#ifdef VERBOSE
  std::cout << "* INFO * \nFrom  void BFS_connected(Graph & gin, "
    "std::vector<bool> & explored, std::vector<unsigned int> & cc_nodes,"
    "unsigned int start_node): BFS found "
	    << cc_nodes.size()
	    << " nodes out of " << gin.n()
	    << " nodes, starting from node id "
	    << start_node << ".\n"
	    <<std::endl;
#endif
}

// Explore nodes reachable from a given node (connected component if undirected)
void BFS_reachable_from(Graph & gin,  unsigned int start_node,
		   std::vector<unsigned int> & cfrom_nodes)
{
  std::queue<unsigned int> q; // auxiliary queue for algorithm
  q.push(gin.index_of_id(start_node));
  
  cfrom_nodes.clear(); // clear list of nodes in the connected component
  cfrom_nodes.push_back(start_node); // add starting node to connected component
  
  std::vector<bool> explored(gin.n(),false);
  explored[q.front()]=true; // mark start node explored
  
  while(q.size() != 0)
    {
      unsigned int inow = q.front();
      q.pop();
      // Take first in queue and loop over all its edges
      for(std::list<unsigned int>::iterator it = gin.nodes_[inow].begin();
	  it != gin.nodes_[inow].end(); ++it)
	{
	  if(gin.edges_[*it].first == inow && !explored[gin.edges_[*it].second])
	    {
	      unsigned int iadd = gin.edges_[*it].second;
	      explored[iadd] = true;
	      q.push(iadd);
	      cfrom_nodes.push_back(gin.id_of_index(iadd));
	    }
	  else if(!gin.edges_[*it].directed && gin.edges_[*it].second == inow
		  && !explored[gin.edges_[*it].first])
	    {
	      unsigned int iadd = gin.edges_[*it].first;
	      explored[iadd] = true;
	      q.push(iadd);
	      cfrom_nodes.push_back(gin.id_of_index(iadd));
	    }
	    
	}  
    }
  
#ifdef VERBOSE
  std::cout << "* INFO * \nFrom  void BFS_reachable_from(Graph & gin,  unsigned"
    " int start_node, std::vector<unsigned int> & cfrom_nodes): BFS found "
	    << cfrom_nodes.size()
	    << " nodes out of " << gin.n()
	    << " nodes, starting from node id "
	    << start_node << ".\n"
	    <<std::endl;
#endif
}

// Find all connected components of a graph and return a vector of vectors  of
// nodes in each subgraph
void BFS_all_connected_components(Graph & gin,
				  std::vector<std::vector<unsigned int> >
				  & all_ccs)
{
  std::vector<bool> explored; // To keep track of explored nodes
  std::vector<unsigned int> cc_nodes; // To save nodes in 1 connected component
  all_ccs.clear(); // Clear output listing of connected components

  // Find first connected component (also populates "explored" appropriately)
  BFS_connected(gin,gin.id_of_index(0),explored,cc_nodes);
  all_ccs.push_back(cc_nodes); // save connected component
  
  // Loop over all other nodes...
  for(unsigned int start_node = 1; start_node != gin.n(); ++start_node)
    {
      if(!explored[start_node]) // ... and start a new search if not visited yet
	{
	  BFS_connected(gin,gin.id_of_index(start_node),explored,cc_nodes);
	  all_ccs.push_back(cc_nodes); // save connected component
	}
    }
}

// Find all shortest distances from a source node
unsigned int BFS_distances_from_source(Graph & gin, unsigned int start_node,
				       std::vector<unsigned int> & dists,
				       std::vector<unsigned int> & cc_nodes)
{ 
  std::queue<unsigned int> q; // auxiliary queue for algorithm
  q.push(gin.index_of_id(start_node));
  
  cc_nodes.clear(); // clear list of nodes in the connected component
  cc_nodes.push_back(start_node); // add starting node to connected component
  
  std::vector<bool> explored(gin.n(),false);
  explored[q.front()]=true; // mark start node explored
    
  dists.clear(); // reset dists
  dists.resize(gin.n(),std::numeric_limits<unsigned int>::max()); 
  dists[q.front()] = 0; // set first node to have zero distance
  
  unsigned int inow;
  while(q.size() != 0)
    {
      inow = q.front();
      q.pop();
      // Take first in queue and loop over all its edges
      for(std::list<unsigned int>::iterator it = gin.nodes_[inow].begin();
	  it != gin.nodes_[inow].end(); ++it)
	{
	  if(gin.edges_[*it].first == inow && !explored[gin.edges_[*it].second])
	    {
	      unsigned int iadd = gin.edges_[*it].second;
	      explored[iadd] = true;
	      dists[iadd] = dists[inow] + 1;
	      q.push(iadd);
	      cc_nodes.push_back(gin.id_of_index(iadd));
	    }
	  else if(!gin.edges_[*it].directed && gin.edges_[*it].second == inow
		  && !explored[gin.edges_[*it].first])
	    {
	      unsigned int iadd = gin.edges_[*it].first;
	      explored[iadd] = true;
	      dists[iadd] = dists[inow] + 1;
	      q.push(iadd);
	      cc_nodes.push_back(gin.id_of_index(iadd));
	    }
	    
	}  
    }
  
#ifdef VERBOSE
  std::cout << "* INFO * \nFrom  unsigned int BFS_distances_from_source(Graph &"
    " gin, unsigned int start_node, std::vector<unsigned int> & dists): "
    "BFS found "
	    << cc_nodes.size()
	    << " nodes out of " << gin.n()
	    << " nodes, starting from node id "
	    << start_node << ".\n"
	    << std::endl
	    << "--- Distances from source node also found!" << std::endl;
#endif
  return dists[inow] + 1; // return number of layers
}

// Find shortest distance between two nodes
// (returns std::numeric_limits<unsigned int>::max() if no path exists)
unsigned int BFS_shortest_path_length(Graph & gin, unsigned int start_node,
				      unsigned int end_node)
{
  std::queue<unsigned int> q; // auxiliary queue for algorithm
  q.push(gin.index_of_id(start_node));
  unsigned int target{gin.index_of_id(end_node)};
  
  std::vector<bool> explored(gin.n(),false);
  explored[q.front()]=true; // mark start node explored

  std::vector<unsigned int> dists;
  dists.resize(gin.n(),0); 
  dists[q.front()] = 0; // set first node to have zero distance
  
  while(q.size() != 0)
    {
      unsigned int inow = q.front();
      q.pop();
      // Take first in queue and loop over all its edges
      for(std::list<unsigned int>::iterator it = gin.nodes_[inow].begin();
	  it != gin.nodes_[inow].end(); ++it)
	{
	  if(gin.edges_[*it].first == target) 
	    return dists[inow] + 1; // undirected or directed cases
	  else if(!gin.edges_[*it].directed && gin.edges_[*it].second == target)
	    return dists[inow] + 1; // directed only
	  if(gin.edges_[*it].first == inow && !explored[gin.edges_[*it].second])
	    {
	      unsigned int iadd = gin.edges_[*it].second;
	      explored[iadd] = true;
	      dists[iadd] = dists[inow] + 1;
	      q.push(iadd);
	    }
	  else if(!gin.edges_[*it].directed && gin.edges_[*it].second == inow
		  && !explored[gin.edges_[*it].first])
	    {
	      unsigned int iadd = gin.edges_[*it].first;
	      explored[iadd] = true;
	      dists[iadd] = dists[inow] + 1;
	      q.push(iadd);
	    }
	  
	}  
    }
  return std::numeric_limits<unsigned int>::max();
}


//------------------------------------------------------------------------------
//--- DFS based algorithms 
//------------------------------------------------------------------------------

/*/ WARNING: This passes a copy of the graph because edges need to be destroyed.
// Explore a connected component of a graph and return its number of nodes. 
// explored is updated (if passed from previous calls)
void DFS_connected_internal_ids_(Graph & gin,  unsigned int start_node,
				 std::vector<bool> & explored,
				 std::vector<unsigned int> & cc_nodes)
{
  std::stack<unsigned int> sk; // auxiliary queue for algorithm
  sk.push(start_node);
  
  cc_nodes.clear(); // clear list of nodes in the connected component
  cc_nodes.push_back(start_node); // add starting node to connected component
  
  if(explored.size() == 0) // prepared the "explored" vector
    explored.resize(gin.m(),false);
  //explored[sk.top()]=true; // mark start node explored
  unsigned int i{0};
  
  std::cout <<  " sktop!! "<< sk.top()  << std::endl;
  
  while(sk.size() != 0)
    {
      unsigned int inow = sk.top();
      std::cout << " still in while " << std::endl;
      // Explore one edge
      if(gin.nodes_[inow].size() != 0)
	{
	  std::cout << " point 1 " << std::endl;
	  unsigned int iedge = gin.nodes_[inow].front();
	  if(!explored[iedge] && gin.edges_[iedge].first == inow)
	    //gin.edges_[iedge].first == inow
	    //&& !explored[gin.edges_[iedge].second])
	    //&& !explored[gin.edges_[iedge].second])
	    {
	      std::cout << " point 1.1 " << std::endl;
	      unsigned int iadd = gin.edges_[iedge].second;
	      std::cout << " point 1.1.1 " << std::endl;
	      explored[iedge] = true;
	      std::cout << " point 1.1.2 " << std::endl;
	      //explored[iadd] = true;  // Mark new node explored
	      sk.push(iadd); // Add node to stack
	      std::cout << " point 1.1.3 " << std::endl;
	      
	      // Delete edge from the node to avoid re-checks in future
	      gin.nodes_[inow].erase(gin.nodes_[inow].begin());
	      //cc_nodes.push_back(gin.id_of_index(iadd));
	    }
	  else if(!explored[iedge] && gin.edges_[iedge].second == inow)
	    //else if(gin.edges_[iedge].second == inow
	    // && !explored[gin.edges_[iedge].first])
	    {
	      std::cout << " point 1.2 " << std::endl;
	      unsigned int iadd = gin.edges_[iedge].first;
	      std::cout << " point 1.2.1 " << std::endl;
	      //explored[iadd] = true;  // Mark new node explored
	      explored[iedge] = true;  // Mark new node explored
	      std::cout << " point 1.2.2 " << std::endl;
	      sk.push(iadd); // Add node to stack
	      // Delete edge from this node
	      std::cout << " point 1.2.3 " << std::endl;
	      gin.nodes_[inow].erase(gin.nodes_[inow].begin());
	      //cc_nodes.push_back(gin.id_of_index(iadd));
	    }
	  else{
	    std::cout << " point 1.3 " << std::endl;
	    gin.nodes_[inow].erase(gin.nodes_[inow].begin());
	  }
	}
      else
	sk.pop(); // There is nothing else to explore so remove from stack
      if(sk.size())
	std::cout << " stack size "<< sk.size()<< "top is " << sk.top() << std::endl;
      //std::cout << " ??? "<< sk.size() << " sktop? "<< sk.top() << " gin.nodes_[inow].front() " << gin.nodes_[inow].front() << std::endl;
      ///++i;
      //if(i>50)
	exit(0);
    }

  
#ifdef VERBOSE
  std::cout << "* INFO * \nFrom  void DFS_connected_internal_ids_(Graph & gin,"
    " unsigned int start_node, std::vector<bool> & explored,"
    " std::vector<unsigned int> & cc_nodes) DFS found "
	    << cc_nodes.size()
	    << " nodes out of " << gin.n()
	    << " nodes, starting from node with internal index "
	    << start_node << ".\n"
	    <<std::endl;
#endif
}*/

// THIS IS A RECURSIVE VERSION
// Explore a connected component of a graph and return its number of nodes. 
// explored is updated (if passed from previous calls)
void DFS_connected_internal_ids_(Graph & gin,  unsigned int start_node,
				 std::vector<bool> & explored,
				 std::vector<unsigned int> & cc_nodes)
{
  explored[start_node]=true; // mark start node explored
  cc_nodes.push_back(gin.id_of_index(start_node));

  for(std::list<unsigned int>::iterator it = gin.nodes_[start_node].begin();
      it != gin.nodes_[start_node].end(); ++it)
    {
      if(gin.edges_[*it].first == start_node
	 && !explored[gin.edges_[*it].second])
	DFS_connected_internal_ids_(gin,
				    gin.edges_[*it].second, explored, cc_nodes);
      else if(gin.edges_[*it].second == start_node
	      && !explored[gin.edges_[*it].first])
	DFS_connected_internal_ids_(gin,
				    gin.edges_[*it].first, explored, cc_nodes);
    }
  
#ifdef VERBOSE
  std::cout << "* INFO * \nFrom  void DFS_connected_internal_ids_(Graph & gin,"
    " unsigned int start_node, std::vector<bool> & explored,"
    " std::vector<unsigned int> & cc_nodes) DFS found "
	    << cc_nodes.size()
	    << " nodes out of " << gin.n()
	    << " nodes, starting from node with internal index "
	    << start_node << ".\n"
	    <<std::endl;
#endif
}

void DFS_connected(Graph & gin,  unsigned int start_node,
		   std::vector<bool> & explored,
		   std::vector<unsigned int> & cc_nodes)
{
  // Prepare "explored"
  if(explored.size() != gin.n())
    {
      explored.clear();
      explored.resize(gin.n(),false);
    }
  cc_nodes.clear();
  //Run using internal ids
  DFS_connected_internal_ids_(gin, gin.index_of_id(start_node),
			      explored, cc_nodes);
}

// THIS IS A RECURSIVE VERSION
// Explore a nodes reachable from a start node of a graph
void DFS_reachable_from_internal_ids_(Graph & gin,  unsigned int start_node,
				 std::vector<bool> & explored,
				 std::vector<unsigned int> & cfrom_nodes)
{

  explored[start_node]=true; // mark start node explored
  cfrom_nodes.push_back(gin.id_of_index(start_node));

  for(std::list<unsigned int>::iterator it = gin.nodes_[start_node].begin();
      it != gin.nodes_[start_node].end(); ++it)
    { 
      if(gin.edges_[*it].first == start_node
	 && !explored[gin.edges_[*it].second])
	DFS_reachable_from_internal_ids_(gin,
					 gin.edges_[*it].second, explored,
					 cfrom_nodes);
      else if(!gin.edges_[*it].directed && gin.edges_[*it].second == start_node
	      && !explored[gin.edges_[*it].first])
	DFS_reachable_from_internal_ids_(gin,
					 gin.edges_[*it].first, explored,
					 cfrom_nodes);
    }
  
#ifdef VERBOSE
  std::cout << "* INFO * \nFrom  DFS_reachable_from_internal_ids_(Graph & gin,"
    " unsigned int start_node, std::vector<bool> & explored,"
    " std::vector<unsigned int> & cfrom_nodes) DFS found "
	    << cfrom_nodes.size()
	    << " nodes out of " << gin.n()
	    << " nodes, starting from node with internal index "
	    << start_node << ".\n"
	    <<std::endl;
#endif
}

// Find all connected components of a graph and return a vector of vectors  of
// nodes in each subgraph
void DFS_all_connected_components(Graph & gin,
				  std::vector<std::vector<unsigned int> >
				  & all_ccs)
{
  std::vector<bool> explored; // To keep track of explored nodes
  std::vector<unsigned int> cc_nodes; // To save nodes in 1 connected component
  all_ccs.clear(); // Clear output listing of connected components

  // Find first connected component (also populates "explored" appropriately)
  DFS_connected(gin,gin.id_of_index(0),explored,cc_nodes);
  all_ccs.push_back(cc_nodes); // save connected component
  
  // Loop over all other nodes...
  for(unsigned int start_node = 1; start_node != gin.n(); ++start_node)
    {
      if(!explored[start_node]) // ... and start a new search if not visited yet
	{
	  DFS_connected(gin,gin.id_of_index(start_node),explored,cc_nodes);
	  all_ccs.push_back(cc_nodes); // save connected component
	}
    }
}

// Explore nodes reachable from a given node (connected component if undirected)
void DFS_reachable_from(Graph & gin,  unsigned int start_node,
			std::vector<unsigned int> & cfrom_nodes)
{
  // Prepare "explored"
  std::vector<bool> explored(gin.n(),false);
  //Run using internal ids
  DFS_reachable_from_internal_ids_(gin, gin.index_of_id(start_node),
				   explored, cfrom_nodes);
}

void DFS_Kosaraju_reversed_(Graph & gin, unsigned int start_node)
{
  explored1[start_node] = true;
  for(std::list<unsigned int>::iterator it = gin.nodes_[start_node].begin();
      it != gin.nodes_[start_node].end(); ++it)
    if(gin.edges_[*it].second == start_node
       && !explored1[gin.edges_[*it].first])
      DFS_Kosaraju_reversed_(gin, gin.edges_[*it].first);
  finishing_vec.push_back(start_node); // add finished node to vector
}

void DFS_Kosaraju_direct_(Graph & gin, unsigned int start_node)
{
  explored2[start_node] = true;
  ++leader_count; // increment count for this leader
  for(std::list<unsigned int>::iterator it = gin.nodes_[start_node].begin();
      it != gin.nodes_[start_node].end(); ++it)
    if(gin.edges_[*it].first == start_node
       && !explored2[gin.edges_[*it].second])
      DFS_Kosaraju_direct_(gin, gin.edges_[*it].second);
}


void SCC_Kosaraju(Graph & gin, std::multiset<unsigned int> & all_SCC_sizes,
		  std::set<unsigned int> & all_SCC_setsizes)
{
  // Perform first pass on the whole reversed graph
  explored1.resize(gin.n(),false);
  // loop over all nodes
  for(unsigned int i = 0; i != gin.n(); ++i)
      if(!explored1[i])
	DFS_Kosaraju_reversed_(gin, i);
  
  // Peform second pass and store multiplicities of leaders in the multiset
  explored2.resize(gin.n(),false);
  // loop over all nodes
  for(std::vector<unsigned int>::iterator it = --finishing_vec.end();
      it != --finishing_vec.begin(); --it)
    if(!explored2[*it])
      {
	leader = *it;
	leader_count = 0;
	DFS_Kosaraju_direct_(gin, *it);
	all_SCC_sizes.insert(leader_count);
	all_SCC_setsizes.insert(leader_count);
      }
  
  finishing_vec.clear();
  explored1.clear();
  explored2.clear();
  
}

/*------------------------------------------------------------------------------
  Dijkstra's based algorithms
  ----------------------------------------------------------------------------*/

// WARNING: Implemented as undirected for now!
// Find all shortest paths from a source node and corresponding distances
void Dijkstra_shortest_paths(Graph & gin,  unsigned int start_node,
			std::vector<unsigned int> & dists,
			std::vector< std::vector<unsigned int> > & paths)
{
  // Prepare data structures

  // Auxiliary integer to test if node is in heap
  unsigned int not_in_heap{std::numeric_limits<unsigned int>::max()};
  
  // Prepare vector to hold distances of each node to start_node
  if(dists.size() != 0) 
    dists.clear();
  dists.resize(gin.n(),1000000);
  // Set start_node to zero distance
  dists[start_node] = 0; 
  
  // Prepare vector to hold shortest paths between each node and start_node
  if(paths.size() != 0) 
    paths.clear();
  paths.resize(gin.n());
  paths[start_node].push_back(start_node);

  // Prepare boolean to check if node is in explored set
  std::vector<bool> explored(gin.n(),false);
  explored[start_node] = true;
  
  // Prepare myHeap object
  myHeap nodesHeap(gin.n());
  // Deal separately with start node
  // loop over start_node edges to populate nodesHeap
  for(std::list<unsigned int>::iterator
	it_edge = gin.nodes_[start_node].begin();
      it_edge != gin.nodes_[start_node].end(); ++it_edge)
    {// Assuming only one edge between any pair of nodes
      // Get other node number
      unsigned int node2{gin.edges_[*it_edge].second};
      if(gin.edges_[*it_edge].first != start_node)
	node2 = gin.edges_[*it_edge].first;
      // Insert cost, node and edge in Heap
      nodesHeap.insert(Triplet(gin.edges_[*it_edge].weight, node2, *it_edge));
      ////std::cout << "Initial heap building..." << std::endl;
      ///nodesHeap.print();
    }

  for(unsigned int i = 1; i != gin.n(); ++i)
    {   
      ////nodesHeap.print();
      // Extract min and add to explored nodes
      Triplet node_min{nodesHeap.extract_min()};
      // Mark node explored
      explored[node_min.t1] = true;

      // Save distance
      dists[node_min.t1] = node_min.t0;

      // Save path
      //// First work out origin node
      unsigned int node1{gin.edges_[node_min.t2].first};
      if(gin.edges_[node_min.t2].first == node_min.t1)
	node1 = gin.edges_[node_min.t2].second;
      //// Now copy path
      paths[node_min.t1] = paths[node1];
      //// and add extra edge
      paths[node_min.t1].push_back(node_min.t1);
      
      // loop over edges of added node to enlarge heap
      for(std::list<unsigned int>::iterator
	    it_edge = gin.nodes_[node_min.t1].begin();
	  it_edge != gin.nodes_[node_min.t1].end(); ++it_edge)
	{
	  // Get other node number
	  unsigned int node2{gin.edges_[*it_edge].second};
	  
	  if(gin.edges_[*it_edge].second == node_min.t1)
	    node2 = gin.edges_[*it_edge].first;
	  if(!explored[node2])
	    {
	      // score of new path
	      unsigned int score_try{dists[node_min.t1]
		  + gin.edges_[*it_edge].weight};
	      // position in heap (if valid)
	      unsigned int pos2 = nodesHeap.pos_node_heap[node2];
	      
	      // If in heap and if we should replace
	      if( pos2 != not_in_heap
		  && score_try < nodesHeap.heap_data[pos2].t0)
		{
		  nodesHeap.remove(pos2);
		  nodesHeap.insert(Triplet(score_try,node2, *it_edge));
		}
	      else if(pos2 == not_in_heap) // If not in heap add it
		nodesHeap.insert(Triplet(score_try,node2, *it_edge));
	    }
	}	
    }

      

}
