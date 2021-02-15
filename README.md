## What?

* Advanced Intro to boost graph library
* Intro to CGAL
* Common Pitfalls and Further Things You Should Know

## Why?

The Algolab Lecture at ETH aims to teach you how to put some algorithmic building blocks together and to deal with intentionally confusing exercise statements.
The exercises they provide are great, but the guidance less so.

I have especially found the tutorial sessions on what the libraries we are supposed to use offer and how to use them quite lacking. Since I haven't found many good comprehensive resources for boost beginners out there, I have kept track of my own progress.

> "You don't need to understand that, you can copy paste it from the provided example codes during the exam."

## How?

The intended audience are master students taking the Algorithms Lab course. I am assuming generic coding skills similar to the ones I had at the start of the semester.

I originally planned to go through all my notes again and write them up in better words and more detail, but I think I won't be doing that. See the file `codingtricks.md` for my raw notes that cover more than just the basics in the sections below.

# How To Boost

The [boost](https://boost.org) library offers various algorithms that are already implemented and you only need to pick and combine. It is very useful to be aware that there is a [BGL table of contents](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/table_of_contents.html) which is available during the exam.

It's important to know what they already have, so you don't end up implementing dijstra's shortest path algorithm yourself... again.

If you like videos, there's a [good one here](https://www.youtube.com/watch?v=GSp2531Wti4) that covers the basics.

## Algorithms

Some algorithms you should know exist follow. I'll get to how to use them later.

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

Once you have built a graph with such internal properties as required by the function you want to call, you can call it quite simply and it will automatically use these properties. For example after a quick look at the [kruskal](https://stackoverflow.com/a/56124666/2550406) docs, I add the include from "Where Defined" and call it like this:

```c++
#include <boost/graph/kruskal_min_spanning_tree.hpp>
// [Graph Type Definition Here]

// now inside some function:
boost::kruskal_minimum_spanning_tree ( 
    G /* our graph */,
    out_it /* where it should write the output to */ )
```

But we don't know yet what that `out_it` is supposed to be. More on that in the section on [calling kruskal](#calling-kruskal) below.

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

I can't tell you why I prefer that, though. Probably because there are fewer ways to mess this one up.

## Calling

You know now how to build a graph using internal properties. They will automatically be used when you call a BGL function like [`kruskal_minimum_spanning_tree`](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/kruskal_min_spanning_tree.html). But what about the other arguments?

The first argument is the graph, the second is an Output Iterator according to the docs, so ignoring all the optional parameters we could call it like this:

```c++
std::vector<edge_desc> mst;
boost::kruskal_minimum_spanning_tree(G, std::back_inserter(mst));
```

### Named Parameters

Many boost functions have in their declaration an optional parameter of type`bgl_named_params`:

```c++
template <  class Graph, class OutputIterator,
			class P, class T, class R >  OutputIterator

kruskal_minimum_spanning_tree(
    Graph& g, OutputIterator tree_edges,
    const bgl_named_params<P, T, R>& params = all defaults);
```

You can ignore that option if you want to.

What they do is outlined further down in the docs and they either specify `IN` (allow you to give additional inputs) or `OUT` (allow you to read additional outputs) or `UTIL`(used internally).

For example you might want to get a distance map from [Prim's algorithm](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/prim_minimum_spanning_tree.html) so that you know the weight of the edges in the MST. That is supposed to be a Property Map that hands you the distance from a given vertex to its direct parent in the MST. We can create such a map in multiple ways - see below. Or perhaps you want to tell Prim's algorithm where to start at, using the named parameter `root_vertex` that is also documented on that page.

You can do that like this:

```c++
std::vector<vertex_desc> pred_map(n);
auto p = boost::make_iterator_property_map(pred_map.begin(),
                                           boost::get(boost::vertex_index, G));
/*std::vector<int> distmap(n);
        auto d = boost::make_iterator_property_map(distmap.begin(),
                boost::get(boost::vertex_index, G));*/
boost::prim_minimum_spanning_tree(G, p,                      boost::root_vertex(tat)/*.distance_map(d)*/);
```

The commented-out part shows you how to use more than one named parameter.

### Custom Property Map (Vector)

```c++
#include <boost/property_map/property_map.hpp>
auto propmap = boost::make_iterator_property_map(myvec.begin(), boost::get(boost::vertex_index, G));
```

Writing to the newly created property map with the vertex as index will first look up the vertex in the second argument and because we specified a property map there that retrieves the `vertex_index`, it looks up the `vertex_index`. With `vecS` as vertex list in the graph, the `vertex_index` of a graph is often simply the `vertex_descriptor` but this here is how it's done correctly.

Then that index is used to perform an access to the vector `myvec` at that index. So modifying `propmap` will modify the vector.

Note: This would not work with edges that easily. Because edges don't have an `edge_index` by default. You would have to add that `boost::edge_index_t` property in the graph typedef and also set them yourself for each edge you add. Then you could use the same concept.

### Custom Property Map (Lambda)

```c++
// https://stackoverflow.com/questions/26979585/is-it-possible-to-have-several-edge-weight-property-maps-for-one-graph
// You can compose a property map in various ways. The simplest approach would seem something like:

#include <boost/property_map/function_property_map.hpp>
auto wmap1 = make_function_property_map<unsigned, float>([&weight_data](unsigned vertex_id) { return weight_data.at(vertex_id).weight1; });
```

We can even use it for writing to! E.g. with a vector backing it as storage:

```c++
std::vector<unsigned long> compmap(mbridges); // storage

// iterator property map (alternative)
auto propmap = boost::make_iterator_property_map(
    compmap.begin(), boost::get(boost::edge_index, OG));

// function property map
auto funmap = boost::make_function_property_map<edge_desc, unsigned long&>(
         [&OG, &compmap](edge_desc edge)->unsigned long&{
             // use edge_id as index in vector and return writeable reference
             size_t index = boost::get(boost::edge_index, OG)[edge];
                    return compmap[index];
                });
// use the map
unsigned long nbc = boost::biconnected_components(OG, funmap /* could also use propmap here */);

```

The second template argument is not necessary, because it can be inferred from the return type of the lambda. The lambda must have a return type that is a non-const reference, otherwise it cannot be written to - only read - and hence is not an lvalue. Final note: the lambda does not use an `&` in the return statement.

## How To CGAL

There are mostly three things we do with that library:

* Delaunay Triangulations
* Linear Programs
* Something with geometry and exact number types.

### Delaunay

#### Store Info

Copy paste the code for the typedefs from the last slide of the Delaunay slides. Then adapt the includes and typedefs so you can store a struct with each Vertex and / or each Face:

```c++
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
class VertexInfo{
    public:
        uint id;
};
class FaceInfo{};
// we want to store an index with each vertex
typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo,K>   Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo,K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>            Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                  Triangulation;
typedef Triangulation::Face_handle Fh;
typedef Triangulation::Vertex_handle Vh;
typedef K::Point_2 P;
```

#### Read Points

```c++
K::Point_2 p;
std::cin >> p;
```

This works. No need to first store both coordinates in ints.

#### Insert Points

There's a big speed difference between the different possibilities.

Ideally, you fill a vector `v` with all the points and let CGAL handle the ordering. (Instead of inserting one after the other. In that case you should maybe shuffle them to avoid adversarial inputs that decrease the runtime enormously!).

You [can](https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#ac5e9bc8adef80dc01a0b31c2d0234545) store the `VertexInfo` defined above with each Point in a Pair in the vector if you want to:

```c++
Triangulation tri;
std::vector<std::pair<Point_2, VertexInfo> > vec;
// [fill vec]
tri.insert(vec.begin(), vec.end());
```

#### Documentation

Relevant are

* [Delaunay_triangulation_2](https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Delaunay__triangulation__2.html)
* In the sidebar a bit below it [Triangulation_2](https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html). The Delaunay Triangulation has all those functions as well.
* [Triangulation_data_structure_2](https://doc.cgal.org/latest/TDS_2/classTriangulationDataStructure__2.html) with some weird concepts.

#### Faces

if you have a `Triangulation::Face_handle`, you can use the `->` operator to get the [neighbor](https://doc.cgal.org/latest/TDS_2/classTriangulationDataStructure__2_1_1Face.html#a096a50bd226daf09826eeacf54727f0e) face at an index between 0 and 2 (including).

That neighbour face at index $i$ is the one that is adjacent to edge $i$, which is opposite to the vertex $i$. You can get that [vertex](https://doc.cgal.org/latest/TDS_2/classTriangulationDataStructure__2_1_1Face.html#ac4d03671704cd164b279706a098fdf5e) handle with the index as well. `fh->vertex(1)`. Have a look at the `fh->cw()` and `fh->ccw()` method for juggling with those indices.

[circumcenter](https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#a4185c75ba2c5ec34181fdef8fa57401c) can be useful.

#### Edges

Edges are just a pair of a Face handle and an index. You can use `.first` and `.second` to get the face and the edge index $i$.

The triangulation has the method [mirror_edge](https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#ab97ce60b20674d0a7a4455e88c2eadb1) which gives you the same edge again, but from perspective of the neighbouring face.

tri.[segment](https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#a06f3967c92db0fe28368c31ff671c658) can be useful to compute the length if you don't want to use `CGAL::squared_distance` on the two endpoints for some reason.

#### Circulators

They are like iterators, but scary. E.g. `incident_faces` and `incident_vertices` are returning circulators.

* They can be null if you only have one Point in the triangulation
* They make it easy to do something wrong in the loop condition. Remember to increment them, always.
* They can give you something infinite.
* Use a `do...while` loop. *if* it is not null.
* Stop the loop once the circulator that you incremented is equal to the original circulator.

### Linear Programs

I think the slides suffice on this. I will have some common pitfalls regarding this though. Or if I never get to writing that section, you can find them in the `codingtricks.md` scattered all over the different tasks.

* doubles are generally faster than ints, counter-intuitively.
* Check out QP_Bland rule

### Geometry

* Geometry is slow. So it can pay off to avoid constructions.
* Use the worst possible kernel, it's the fastest as well.
* Lower numbers are generally faster than higher numbers.
* To avoid using square root, just work with the squares.
* Don't construct intersections directly - first check if they actually intersect. That's faster.

## C++

Just see the `codingtricks.md` for random but sometimes useful facts that I didn't know before (or maybe I did).

Things like how to sort a `priority_queue` with a lambda function.

## Common Pitfalls

TODO: collect them from all over `codingtricks.md`.

## Ideas For Exercises

If you're stuck, there are many good repos out there with full solutions. If you either don't get them or you want to not see code yet and only the idea, check out the latter part of `codingtricks.md`.