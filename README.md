# JMedialAxis
Medial Axes for JTS Geometry

## Overview

JMedialAxis produces medial axis **trees** from JTS geometries (or collections of coordinates).

The library models the medial axis of a given geometry as a (rooted) tree of medial disks. Medial disks reference a parent and have up to 3 children (disks with more than one child represent bifurcation of the axis). The root of the tree is the disk whose underlying triangle has the largest circumcircle; this is also the single disk that trifurcates.

The library is geared towards medial axis visualisation & animation and enables easy and powerful navigation of medial axes.

## Features

* JTS Geometry as input
* For a given medial disk, easy access to its:
  * Descendants
  * Ancestors
  * Belonging segment
  * Parent disk
  * Children disks
* Branch/Segment pruning
* Supports geometries with holes (genus > 0)