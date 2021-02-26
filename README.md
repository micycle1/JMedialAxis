# JMedialAxis
Medial Axes for JTS Geometry

## 

JMedialAxis produces medial axis trees from JTS geometries.

The library models medial axes of shapes as a (rooted) tree of medial disks. Medial disks reference a parent and have up to 3 children (the medial axis bifurcates/forks at disks with 2 children). The root is the disk with the largest circumcircle, which is also the single disk that trifurcates.

The library is geared towards medial axis visualisation & animation and enables easy and powerful navigation of medial axes.