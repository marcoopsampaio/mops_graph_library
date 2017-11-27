Graph utilities library
=======================
### by Marco O. P. Sampaio

This is a simple library with a Graph class and various algorithms. To compile the library you will need a ```C++11``` compliant compiler. To compile all the sources run

```bash
$ make
```
Test program binaries are then found in ```bin/``` and must be run from the ```Makefile``` directory like this (replace ```<binaryname>``` with relevant test binary)
```bash
$ bin/<binaryname>
```
the compiled static library is found in ```lib/libgraph.o```.
To clean the object files generated during compilation run
```bash
$ make clean
```
or if you need to clean completely the compilation run
```bash
$ make cleanall
```

Algorithms implemented so far include:
* **Breadth First Search (BFS)**: To compute paths and distances, and connected components of a graph.
* **Depth First Search:** Also so compute paths and connectedness, but also all strongly connected components for directed graphs.
* **Heap class:** For Dijkstra's algorithm
* **Dijkstra's algorithm:** Shortest paths on weighted graphs
* __A* algorithm:__ Shortests paths on a weighted graph with a euclidean heuristic.

