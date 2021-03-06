[![](https://jitpack.io/v/micycle1/JMedialAxis.svg)](https://jitpack.io/#micycle1/JMedialAxis)

# JMedialAxis
Medial Axes (2D skeletons) for JTS Geometry

## Overview

*JMedialAxis* produces medial axis **trees** from [JTS](https://en.wikipedia.org/wiki/JTS_Topology_Suite) geometries or collections of coordinates (as opposed to binary 2D images). A.k.a medial axis transform: the input shape can be approximate reconstruction of the original shape.

The library models the medial axis of a given geometry as a (rooted) tree of medial disks (sometimes called *Voronoi Balls*). The result is a medial axis that can be traversed recursively by starting at the largest disk and following the child nodes until they reach the boundary of the geometry.


The library is geared towards medial axis visualisation & animation and enables easy and powerful navigation of medial axes, and getting edges, disks and branches.

The implementation uses ideas presented in *Voronoi Ball Models for Computational Shape Applications*
by Roger C. Tam.

## Features

* JTS geometry as input
* For a given medial disk, easy access to its:
  * Descendants
  * Ancestors
  * Belonging branch
  * Parent disk
  * Child disk(s)
* Branch pruning, via:
  * Feature area pruning
* Supports geometries with holes (genus > 0)
* Output as:
  * List of edges
  * List of segments
  * Dissolved (simplified) JTS geometry
  * JTS LineMergeGraph

## Showcase
...
## Example
..
## TODO
* [ ] VDM simplification via clustering
* [ ] Option to model braches/segments as bezier curves
* [ ] Feature reconstruction
* [ ] Prevent cycles (in shapes with holes) from being pruned