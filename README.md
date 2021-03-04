# JMedialAxis
Medial Axes for JTS Geometry

## Overview

*JMedialAxis* produces medial axis **trees** from [JTS](https://en.wikipedia.org/wiki/JTS_Topology_Suite) geometries or collections of coordinates (as opposed to binary 2D images).

The library models the medial axis of a given geometry as a (rooted) tree of medial disks. The result is a medial axis that can be traversed recursively by starting at the largest
disk and following the child nodes until they reach the boundary of the geometry.

Medial disks reference a parent and have up to 3 children (disks with more than one child represent bifurcation of the axis). The root of the tree is the disk whose underlying triangle has the largest circumcircle; this is also the single disk that trifurcates.

The library is geared towards medial axis visualisation & animation and enables easy and powerful navigation of medial axes.

The implementation uses ideas introduced in *Voronoi Ball Models for Computational Shape Applications*
by Roger C. Tam.

## Features

* JTS Geometry as input
* For a given medial disk, easy access to its:
  * Descendants
  * Ancestors
  * Belonging segment
  * Parent disk
  * Children disks
* Branch/Segment pruning, via:
  * Feature area pruning
* Supports geometries with holes (genus > 0)

## Showcase
...
## Example
..
## TODO
* [ ] VDM simplification via clustering
* [ ] Option to model braches/segments as bezier curves 
* [ ] Prevent cycles (in shapes with holes) from being pruned