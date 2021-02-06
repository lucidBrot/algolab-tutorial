## What?

* Advanced Intro to boost graph library
* Intro to CGAL
* Common Pitfalls
* Further Things You Should Know

## Why?

The Algolab Lecture at ETH aims to teach you how to put some algorithmic building blocks together and to deal with intentionally confusing exercise statements.
The exercises they provide are great, but the guidance less so.

I have especially found the tutorial sessions on what the libraries we are supposed to use offer and how to use them quite lacking. Since I haven't found many good comprehensive resources for boost beginners out there, I have kept track of my own progress.

> "You don't need to understand that, you can copy paste it from the provided example codes during the exam."

## How?

The intended audience are master students taking the Algorithms Lab course. I am assuming generic coding skills similar to the ones I had at the start of the semester.

Check the chapters out by clicking the links in the [What?-section](#what?)

# How To Boost

The [boost](https://boost.org) library offers various algorithms that are already implemented and you only need to pick and combine. It is very useful to be aware that there is a [BGL table of contents](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/table_of_contents.html) which is available during the exam.

It's important to know what they already have, so you don't end up implementing dijstra's shortest path algorithm yourself... again.

## Algorithms

Some algorithms you should know exist follow. I'll get to how to use them in a sec.

* Maximum Flow and Matching Algorithms
  1. [`edmonds_karp_max_flow`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/edmonds_karp_max_flow.html)
  2. [`push_relabel_max_flow`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/push_relabel_max_flow.html)
  3. [`edmonds_maximum_cardinality_matching`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/maximum_matching.html)
  4. [`maximum_weighted_matching`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/maximum_weighted_matching.html)
* Minimum Spanning Tree Algorithms
  1. [`kruskal_minimum_spanning_tree`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/kruskal_min_spanning_tree.html)
  2. [`prim_minimum_spanning_tree`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/prim_minimum_spanning_tree.html)
* Minimum Cost Maximum Flow Algorithms
  1. [`cycle_canceling`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/cycle_canceling.html)
  2. [`successive_shortest_path_nonnegative_weights`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/successive_shortest_path_nonnegative_weights.html)
* Connected Components Algorithms
  1. [`connected_components`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/connected_components.html)
  2. [`strong_components`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/strong_components.html)
  3. [`biconnected_components`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/biconnected_components.html)
  4. [`articulation_points`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/biconnected_components.html#sec:articulation_points)
* [`dijkstra_shortest_paths`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/dijkstra_shortest_paths.html)

## Building A Graph

There are different graph types, but for Algolab you will use the same one almost always: [adjacency_list](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/using_adjacency_list.html) with underlying storage `boost::vecS` for the Vertex List and Edge List.

There is also the option of using an `adjacency_matrix` instead of an `adjacency_list` but I have actually never done that during the whole course.

The type of the Graph looks somewhat scary when you're not used to it:

```c++
typedef boost::adjacency_list<
               boost::vecS,
               boost::vecS,
			   boost::undirectedS> graph_t;
```

The options for the underlying types of the Edge List and Vertex List are the first two template arguments. They have an effect on speed, memory consumption, and on which operations break the iterators over vertices or edges. On [this page](https://www.boost.org/doc/libs/1_58_0/libs/graph/doc/adjacency_list.html) you can find a table that tells you more on that. Some facts that I found somewhere:

* If you don't want parallel edges (two edges with the same source and target) you could use a `setS` or `hash_setS` type for the Edge List. If you don't care, use a sequence based type like `listS`, `vecS`, or `slistS` because `add_edge()` is faster there.
  However, finding an edge is faster in AssociativeContainers and hence clearing a vertex also is.

* The iterator over edges is fastest with `vecS` and slowest with `hash_setS`.

* My notes say

  > Use `listS`for VertexList  (the second type argument). It uses more space than `vecS` but is constant time instead of `O(V+E)` for `remove_vertex()` and it does not cause invalidation of vertex descriptors on deletion.

  although looking back, I think I never heeded my own advice and just used `vecS`for everything.

### Directedness

A graph is either `directedS`, `undirectedS`, or `bidirectionalS`. The last one stores all the information twice. Only in that last case can you get the *incoming* edges of a vertex. But you don't need that often.

### Weights / Internal Properties

Some algorithms work with special properties like "Vertex Color" or "Edge Weight" etc. The simplest way to make them work is by using a predefined type in the graph typedef to tell boost that we need a graph with that property. For example, let's try run [kruskal's algorithm](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/kruskal_min_spanning_tree.html) for building a minimum spanning tree.

For this we make use of some (optional) template arguments that come after the ones I've already mentioned.

```c++
typedef boost::adjacency_list<
               boost::vecS,
               boost::vecS, boost::undirectedS,
               boost::no_property,
               boost::property<boost::edge_weight_t, int> 
            	> weighted_graph_t;
```

We define the `Vertex Property` to be unset by specifying `boost::no_property` because we don't care about that. And for the edges, we want a weight.

Note that this snippet uses `int` as a type for the edge weight. But sometimes we need a larger type. There's no problem in changing this to a `long`!

The `boost::property` syntax here seems unneccessary bloat at first, but once we would like to have more than just one property associated with the Edge, it makes sense: The `boost::property` actually takes three template arguments if you give it three. The first is the kind of property (`edge_weight_t` here), the second is the type of that property's value (`int` here) and the third, which I have not yet shown, is... *another property*. That's right, template recursion!

We could then start nesting them like this:

```c++
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property,
  boost::property<boost::edge_capacity_t, long,
    boost::property<boost::edge_residual_capacity_t, long,
      boost::property<boost::edge_reverse_t,
                      traits::edge_descriptor>>>> graph_t;
```

But you shouldn't have to worry too much about writing these yourself. The spirit of Algolab is to copy-paste such examples from the slides or the code examples.

### Inserting Vertices

When you have your typedef the way you want it, it's quite easy to instantiate a graph variable of that type with a given number of vertices.

```c++
int main(){
    graph_t G(100);
    return 0;
}
```

Every vertex is identified by some vertex handle. The type for that can be typedef'd e.g. with

```c++
typedef boost::graph_traits<graph_t>::vertex_descriptor          vertex_desc;
```

but in practise as long as your underlying types are always `vecS`, these are actually just integer values. Makes sense, because the vertices are stored in a vector, right?

When you delete a vertex from the graph, all the vertices after that one are shifted around, so deleting a vertex can invalidate your vertex descriptor. The solution for algolab is to just not do that.

### Inserting Edges

```c++
int n = 100;
graph_t G(n);
// add an edge from vertex 0 to vertex 4
boost::add_edge(0, 4, G);
// adding an edge automatically adds the vertex if necessary. 
// But why would you not set the correct graph size from the start?
boost::add_edge(0, 150, G);
```

> ```
> std::pair<edge_descriptor, bool>
> add_edge(vertex_descriptor u, vertex_descriptor v,
>          adjacency_list& g)
> ```
>
> Adds edge *(u,v)* to the graph and returns the edge descriptor for the new edge. For graphs that do not allow parallel edges, if the edge is already in the graph then a duplicate will not be added and the `bool` flag will be `false`. When the flag is `false`, the returned edge descriptor points to the already existing edge.

The `add_edge` function returns an edge descriptor. That allows us to later access that edge. For example if we want to modify its weight.

If we have our weighted graph again, with only one edge property, then we can set that property right when adding the edge.

```c++
boost::add_edge(u, v, weight, G);
```

Other times, when we have many properties per edge, I find it easier to do it in a different way. We can kindly ask boost to give us a "property map" which we can then use to look up a specific property for an edge descriptor.

```c++
// type of edge descriptor
typedef boost::graph_traits<graph_t>::edge_descriptor          edge_desc;
// ...
auto weight_map = boost::get(boost::edge_weight, G);
// Note how we get the property map for edge_weight,
// not for edge_weight_t. Because the _t was only in the typedefs

auto edge_one_pair = boost::add_edge(0, 1, G);
auto edge_two_pair = boost::add_edge(2, 2, G);
// these are pairs of an edge descriptor and a boolean
// but "auto" can save you a lot of typing and confusion

// set weight to 12 for the first edge
edge_desc e1 = edge_one_pair.first;
weight_map[e1] = 12;
```

### Custom Properties

Sometimes we want to store additional information with an edge or a vertex. For example some `id` or whether we have already visited that vertex/edge.

The simplest way to do that is to specify a struct or class as the last of the recursive properties:

```c++
typedef struct {
    uint id;
} EdgeData;

class VertexData {
	public:
    	bool visited;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS,
			boost::undirectedS,
	        VertexData,
			boost::property<
                boost::edge_weight_t, uint, EdgeData>
                >
  weighted_graph;
```

Then, we can access these values using the very simple syntax

```c++
G[my_vertex_descriptor].visited = true;
uint value = G[my_edge_descriptor].id;
```

### Iterating

There are functions to get all incoming edges, all outgoing edges, all neighbouring vertices, all vertices, and so on. All of these give you the descriptors you need for accessing such things.

When you have an edge, you can get the vertices it connects with `boost::source(edge, G)` and `boost::target(edge, G)`.

One way to iterate over all outgoing edges from a vertex (you can also use this on an undirected graph):

```c++
auto it_pair = boost::out_edges(my_vertex_desc, G);
for ( auto it = it_pair.first; it != it_pair.second; it++){
    // do stuff
}
```

I've come to prefer this alternative way:

```c++
#include <boost/range/iterator_range.hpp>
for ( auto it : boost::make_iterator_range(boost::out_edges(my_vertex_desc, G))){
    // do stuff
}
```

I can't tell you why I prefer that, though.