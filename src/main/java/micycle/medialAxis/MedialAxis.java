package micycle.medialAxis;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.operation.linemerge.LineMergeGraph;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.HilbertSort;
import org.tinfour.utils.TriangleCollector;
import org.tinspin.index.kdtree.KDTree;

import micycle.medialAxis.Utility.LinearRingIterator;
import processing.core.PApplet;
import processing.core.PConstants;

/**
 * Based on ideas from
 * https://www.cs.ubc.ca/labs/imager/th/2003/Tam2004/Tam2004.pdf and Voronoi
 * Ball/Disk Models for Computational Shape Applications
 * 
 * <p>
 * Produces a directed acyclic graph of the MD / more specifically rooted tree
 * <p>
 * Medial Axis library models medial axes of shapes as a (rooted) tree of medial
 * disks. Medial disks reference a parent and have upto 3 children (the medial
 * axis bifurcates/forks at disks with 2 children). The root is the disk with
 * the largest circumcirle, also the single disk that trifurcates. All disks
 * have an indegree of 1.
 * <p>
 * Medial axes are...
 * <p>
 * For geometries of genus zero, the graph is naturally a tree. For objects with
 * holes, we break each cycle by imposing an appropriate breakpoint in the loop.
 * This is only for the purposes of traversing the graph without running into
 * infinite loops. In order to preserve the topology of the original object,
 * cycles are never pruned
 * <p>
 * Decompose the disk coordinates into curves.
 * <p>
 * If angle of successive line branches is the same, merge them
 * <p>
 * * One way to get a smoother initial axis (i.e. without pruning) is to apply
 * on a smoothed version of the shape.
 * 
 * @author Michael Carleton
 *
 */
public class MedialAxis {

	// TODO see https://www.youtube.com/watch?v=WiU0foXVSsM
	// also see https://lgg.epfl.ch/publications/2009/scale_axis_transform_2009.pdf

	// TODO convert disks into undirected graph (e.g. to find the distance from the
	// furthest node to every other node).
	// to enable traversal up and down the direction (And find loop connected loops)

	/*
	 * or have a map of bifurcations
	 */

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	private Geometry g;
	private KDTree<MedialDisk> kdTree;

	private Geometry dissolved;
	private LineMergeGraph lineMergeGraph;

	private final ArrayList<MedialDisk> voronoiDisks = new ArrayList<>();
	private MedialDisk rootNode; // the medial axis has a trifuraction at the root
	MedialDisk deepestNode;
	public MedialDisk furthestNode; // node that is furthest spatial distance away from root
	private List<MedialDisk> leaves;
	private List<MedialDisk> bifurcations; // birfurcating/forking nodes
	private List<Branch> branches;
	/**
	 * Gives a mapping for the edge 'parent->e' for a given edge e.
	 */
	private final HashMap<MedialDisk, Edge> edges; // a map of edge-tail -> edge for easier access in pruning methods

	public static boolean debug = false;

	private double minimumAxialGradient = Double.MAX_VALUE; // most negative axial gradient value.
	private double maximumAxialGradient = -Double.MAX_VALUE; // most positive axial gradient value.

//	public ArrayList<ArrayList<VD>> branches = new ArrayList<>(); // node sections between forks

	public MedialAxis(Coordinate[] coordinates) {
		// TODO support holes
		this(GEOM_FACTORY.createLinearRing(coordinates));
	}

//	public MedialAxis(List<Coordinate> coordinates) {
//		this(GEOM_FACTORY.createLinearRing(coordinates));
//	}

	private IndexedPointInAreaLocator pointLocator;

	/**
	 * The medial axis is computed for the geometry during construction.
	 * 
	 * @param g Geometry whose medial axis to compute. Can contain holes. Must not
	 *          be a multigeometry. Must be fairly dense.
	 */
	public MedialAxis(Geometry g) {
//		if (!g.isSimple()) {
//			throw new Exception("Invalid Geometry -- it is not simple (may self-intersect or have repeated points.");
//			return;
//		}

		/*
		 * Simplify if complex
		 */
		if (debug) {
			System.out.println("Input coordinates #: " + g.getCoordinates().length);
		}
//		pointLocator = new IndexedPointInAreaLocator(g); // ideally before densification
//		if (g.getCoordinates().length > 4000) {
//			Polygon poly = (Polygon) g;
//			if (poly.getNumInteriorRing() > 0) {
//				g = TopologyPreservingSimplifier.simplify(g, 1);
//			} else {
//				g = DouglasPeuckerSimplifier.simplify(g, 1);
//			}
//			if (debug) {
//				System.out.println("Simplifying shape. New coordinates #: " + g.getCoordinates().length);
//			}
//		}

		this.g = g;
		/*
		 * ideally would have adaptive densification that densifies more concave/narrow
		 * parts more
		 */
//		this.g = Densifier.densify(g, 10); // NOTE constant=10
		if (debug) {
			System.out.println("Densified coordinates #" + this.g.getCoordinates().length);
		}

		TC tc = prepareTin();
		// TODO clustering after tin generation

		// having done triangulation, build up directed tree of medial axis
		rootNode = new MedialDisk(null, 0, tc.largestDiskTriangle, tc.largestCircumcircle, 0);
		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.push(rootNode); // start tree with rootnode

		HashSet<SimpleTriangle> remaining = new HashSet<>(tc.triangles.size()); // use bit set?
		remaining.addAll(tc.triangles);
		remaining.remove(tc.largestDiskTriangle);

		edges = new HashMap<>();

		int depth = 1; // number of edges in the path from the root to the node // BF ID?
		int id = 1; // breadth-first ID
		double highestDistance = 0;

		while (!stack.isEmpty()) {

			MedialDisk live = stack.pop(); // FIFO (BFS)
			voronoiDisks.add(live);

			if (live.distance > highestDistance) {
				highestDistance = live.distance;
				furthestNode = live;
			}

			// get half-edge duals (edge shared with other triangles)
			final SimpleTriangle n1 = tc.map.get(live.t.getEdgeA().getDual());
			final SimpleTriangle n2 = tc.map.get(live.t.getEdgeB().getDual());
			final SimpleTriangle n3 = tc.map.get(live.t.getEdgeC().getDual());

			if (remaining.contains(n1)) {
				MedialDisk child = new MedialDisk(live, id++, n1, tc.cccMap.get(n1), depth);
				stack.add(child);
				live.addchild(child);
				remaining.remove(n1);
				edges.put(child, new Edge(live, child));
				child.distance = live.distance + distance(live.position, child.position);
			}
			if (remaining.contains(n2)) {
				MedialDisk child = new MedialDisk(live, id++, n2, tc.cccMap.get(n2), depth);
				stack.add(child);
				live.addchild(child);
				remaining.remove(n2);
				edges.put(child, new Edge(live, child));
				child.distance = live.distance + distance(live.position, child.position);
			}
			if (remaining.contains(n3)) {
				MedialDisk child = new MedialDisk(live, id++, n3, tc.cccMap.get(n3), depth);
				stack.add(child);
				live.addchild(child);
				remaining.remove(n3);
				edges.put(child, new Edge(live, child));
				child.distance = live.distance + distance(live.position, child.position);
				// TODO add in feature area here?
			}

			depth++;
		}

		deepestNode = voronoiDisks.get(voronoiDisks.size() - 1); // TODO check
//		calculateFeatureArea();
//		calcFeatureArea();
		if (debug) {
			System.out.println("Total Area: " + rootNode.featureArea);
		}
	}

	public List<MedialDisk> getDisks() {
		return voronoiDisks;
	}

	/**
	 * Returns the root disk -- the disk with largest radius and 3 children.
	 * 
	 * @return
	 */
	public MedialDisk getRoot() {
		return rootNode;
	}

	public MedialDisk nearestDisk(double x, double y) {
		// TODO use PhTree/Quadtree?
		// TODO use triangulation instead (and remove tinspin dependency)?
		// TODO shuffle before insertion?
		if (kdTree == null) { // Lazy initialisation
			kdTree = KDTree.create(2);
			voronoiDisks.forEach(disk -> kdTree.insert(new double[] { disk.position.x, disk.position.y }, disk));
		}
		return kdTree.nnQuery(new double[] { x, y }).value();
	}

	public MedialDisk nearestDisk(Coordinate coordinate) {
		return nearestDisk(coordinate.x, coordinate.y);
	}

	/**
	 * Returns the children of a disk into a linear array
	 * 
	 * @param parent contains chil out: contains children and the parent (at
	 *               position 0)
	 * @return
	 */
	public List<MedialDisk> getDescendants(MedialDisk parent) {

		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.add(parent);

		ArrayList<MedialDisk> out = new ArrayList<>();

		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop();
			out.add(live);
			live.children.forEach(child -> {
				stack.add(child);
			});
		}
		return out;
	}

	public List<MedialDisk> getDescendants(MedialDisk parent, int maxDepth) {

		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.add(parent);

		ArrayList<MedialDisk> out = new ArrayList<>();
		int depth = 0;

		while (!stack.isEmpty() && depth < maxDepth) {
			MedialDisk live = stack.pop();
			out.add(live);
			live.children.forEach(child -> {
				stack.add(child);
			});
			depth++;
		}
		return out;
	}

	/**
	 * The ancestors of a node d are the nodes on the path from d to the root. This
	 * includes node d.
	 * 
	 * <p>
	 * TODO another method that returns paths from fork VDs in the ancestors (Except
	 * rootnode)
	 * 
	 * @param child
	 * @return
	 */
	public List<MedialDisk> getAncestors(MedialDisk child) {

		ArrayList<MedialDisk> out = new ArrayList<>();

		MedialDisk live = child;
		while (live.parent != null) {
			out.add(live);
			live = live.parent;
		}
		return out;
	}

	/**
	 * Returns a JTS geometry where medial axis edges are dissolved into a set of
	 * maximal-length Linestrings. The output can be further simplfied with
	 * DouglasPeuckerSimplifier for example. * Returns a simplified medial axis
	 * graph leaving only maximal-length lines in which every unique segment appears
	 * once only. The output lines run between node vertices of the input, which are
	 * vertices which have either degree 1, or degree 3 or greater.
	 * 
	 * 
	 * @return
	 */
	public Geometry getDissolvedGeometry() {
		if (dissolved == null) { // Lazy initialisation
			LineDissolver ld = new LineDissolver();
			getEdges().forEach(e -> {
				ld.add(e.lineString);
			});
			dissolved = ld.getResult();
//			System.out.println(String.format("Dissolved %s axis edges into %s linestrings.", edges.size(),
//					dissolved.getNumGeometries()));
		}
		return dissolved;
	}

	/**
	 * Get the medial axis in the form of an undirected planar graph. The graph is
	 * based on the dissolved geometry.
	 * 
	 * @return
	 */
	public LineMergeGraph getLineMergeGraph() {
		if (lineMergeGraph == null) { // Lazy initialisation
			getDissolvedGeometry();
			lineMergeGraph = new LineMergeGraph();
			for (int i = 0; i < dissolved.getNumGeometries(); i++) {
				lineMergeGraph.addEdge((LineString) dissolved.getGeometryN(i));
			}
		}
		return lineMergeGraph;
	}

	private void asGraph() {
		// directed jGraphT graph
//		convert to ?cyclical? jGrapt graph // TOD// then org.jgrapht.alg.cycle.CycleDetector
	}

	/**
	 * End nodes / A VDs of degree 1/no children
	 * 
	 * @return
	 */
	public List<MedialDisk> getLeaves() {
		if (leaves == null) { // Lazy initialisation
			leaves = voronoiDisks.stream().filter(vd -> vd.degree == 0).collect(Collectors.toList());
		}
		return leaves; // TODO IMMUTABLE ARRAY, SORTED BY DISK DISTANCE/FEATUREAREA highest first
	}

	public MedialDisk getDeepestDisk() {
//		getLeaves().sort((a, b) -> Double.compare(b.id, a.id));
//		return leaves.get(0); // TODO copy array first
		return deepestNode;
	}

	public MedialDisk getFurthestDisk() {
//		getLeaves().sort((a, b) -> Double.compare(b.id, a.id));
//		return leaves.get(0); // TODO copy array first
		return furthestNode;
	}

	private MedialDisk getDeepestBFDisk() {
		getLeaves().sort((a, b) -> Double.compare(b.depthBF, a.depthBF));
		return leaves.get(0); // TODO copy array first
	}

	/**
	 * Nodes/medial disks with two descendent lineages (two children disks). Note
	 * the root node (a trifurcating node) is not included in the output.
	 * 
	 * @return
	 */
	public List<MedialDisk> getBifurcations() {
		if (bifurcations == null) { // Lazy initialisation
			bifurcations = voronoiDisks.stream().filter(vd -> vd.degree == 2).collect(Collectors.toList());
		}
		return bifurcations;
	}

	/**
	 * 
	 * @return list of edges (in breadth-first order from the root node) connecting
	 *         medial disks.
	 */
	public Collection<Edge> getEdges() {
		return edges.values();
	}

	/**
	 * 
	 * @return a map of edge-tail -> edge for easier access in pruning methods
	 */
	public Map<MedialDisk, Edge> getEdgeMap() {
		return edges;
	}

	/**
	 * Returns the edges comprising the path from the given disk to the root of the
	 * medial axis.
	 * 
	 * @param d
	 * @return
	 */
	public List<Edge> getEdgesToRoot(MedialDisk d) {
		return getAncestors(d).stream().map(a -> edges.get(a)).collect(Collectors.toList());
	}

	/**
	 * Aka features, aka branches Segment is a linear portion of medial disks.
	 * 
	 * @return List of {@link micycle.medialAxis.MedialAxis.Branch branches}
	 */
	public List<Branch> getBranches() {
		if (branches == null) { // Lazy initialisation
			branches = new ArrayList<>();

			Branch branch;
			ArrayList<MedialDisk> forks = new ArrayList<>(getBifurcations());
			forks.add(rootNode);

			for (MedialDisk disk : forks) {
				for (MedialDisk child : disk.children) {
					branch = new Branch(disk); // init branch with bifurcating node
					MedialDisk node = child;
					while (node.degree == 1) {
						branch.add(node);
						node = node.children.get(0); // get the only child node
					}
					branch.end(node); // end with axis leaf or next bifurcating node
					branches.add(branch);
				}
			}
		}
		return this.branches;
	}

	/**
	 * Returns a subset of the axis' edges; in this method edges are pruned by their
	 * axial gradient value. Edges with negative axial gradient values indicate a
	 * narrowing of the shape along that edge, and are more likely to be noise. The
	 * threshold takes into account this shape's minimum and maximum axial values
	 * when filtering.
	 * <p>
	 * Presently this method prunes based on global min and max axial values (rather
	 * than per branch).
	 * 
	 * @param threshold between 0...1, where 0 is no pruning and 1 is maximal
	 *                  pruning. the shape may be fully pruned before threshold = 1
	 * @return
	 */
	public List<Edge> getPrunedEdges(double threshold) {
		threshold = Math.min(1, Math.max(0, threshold)); // constrain 0...1
		final double mappedThreshold = minimumAxialGradient + (maximumAxialGradient - minimumAxialGradient) * (threshold);

		final ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.add(rootNode);

		final ArrayList<Edge> out = new ArrayList<>();

		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop();
			live.children.forEach(child -> {
				if (child.axialGradient >= mappedThreshold) {
					stack.add(child);
					out.add(edges.get(child));
				}
			});
		}
		return out;
	}

	/**
	 * Returns a subset of the axis' edges; in this method edges are pruned by their
	 * axial gradient value and axis distance from the root node.
	 * 
	 * @param axialThreshold    between 0...1, where 0 is no pruning and 1 is
	 *                          maximal pruning
	 * @param distanceThreshold between 0...1, where 0 is no pruning and 1 is
	 *                          maximal pruning
	 * @return
	 */
	public List<Edge> getPrunedEdges(double axialThreshold, double distanceThreshold) {
		axialThreshold = Math.min(1, Math.max(0, axialThreshold)); // constrain 0...1
		axialThreshold = axialThreshold * axialThreshold * axialThreshold; // make 0...1 more linear
		final double mappedAxial = minimumAxialGradient + (maximumAxialGradient - minimumAxialGradient) * axialThreshold;

		distanceThreshold = 1 - Math.min(1, Math.max(0, distanceThreshold)); // constrain 0...1
		final double mappedDistance = furthestNode.distance * distanceThreshold;

		final ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.add(rootNode);

		final ArrayList<Edge> out = new ArrayList<>();

		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop();
			live.children.forEach(child -> {
				if (child.axialGradient >= mappedAxial && child.distance <= mappedDistance) {
					stack.add(child);
					out.add(edges.get(child));
				}
			});
		}
		return out;
	}

	/**
	 * 
	 * @param axialThreshold
	 * @param distanceThreshold
	 * @param areaThreshold     between 0...1, where 0 is no pruning and 1 is
	 *                          maximal pruning for this feature
	 * @return
	 */
	public List<Edge> getPrunedEdges(double axialThreshold, double distanceThreshold, double areaThreshold) {
		axialThreshold = Math.min(1, Math.max(0, axialThreshold)); // constrain 0...1
		axialThreshold = axialThreshold * axialThreshold * axialThreshold; // make 0...1 more linear
		final double mappedAxial = minimumAxialGradient + (maximumAxialGradient - minimumAxialGradient) * axialThreshold;

		distanceThreshold = 1 - Math.min(1, Math.max(0, distanceThreshold)); // constrain 0...1
		final double mappedDistance = furthestNode.distance * distanceThreshold;

		calcFeatureArea();
		areaThreshold = Math.min(1, Math.max(0, areaThreshold)); // constrain 0...1
		areaThreshold = areaThreshold * areaThreshold * areaThreshold; // make 0...1 more linear
		final double mappedArea = areaThreshold * rootNode.featureArea;

		final ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.add(rootNode);

		final ArrayList<Edge> out = new ArrayList<>();

		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop();
			live.children.forEach(child -> {
				if (child.axialGradient >= mappedAxial && child.distance <= mappedDistance && child.featureArea >= mappedArea) {
					stack.add(child);
					out.add(edges.get(child));
				}
			});
		}
		return out;
	}

	/**
	 * Will call consumer with the next disk in a BF manner, starting at the root
	 * node.
	 * 
	 * @param consumer
	 */
	public void iteratorBF(Consumer<MedialDisk> consumer) { // accept starting node as parameter
		// TODO consumer integer arg, corresponding to recursion depth (Same for all
		// children) (or user can reference depthBF int?)
		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.add(rootNode);

		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop(); // pop first of queue
			consumer.accept(live);
			live.children.forEach(stack::add); // insert at end of queue
		}
	}

	public void iteratorDF(Consumer<MedialDisk> consumer) { // accept enum to determine DF via area, distance, bf id, etc
		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.add(rootNode);
		calcFeatureArea();
		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop(); // pop first of queue
			consumer.accept(live);
			live.children.sort((a, b) -> Double.compare(a.distance, b.distance)); // NOTE EXAMPLE
			live.children.forEach(stack::push); // insert at start of queue
		}
	}

	public void iteratorBfUp(Consumer<MedialDisk> consumer) {
		// start at leaves, add parent, cleanup leaves, continue
		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		Set<MedialDisk> seen = new HashSet<>(getDisks().size());
		getLeaves().forEach(stack::add); // sort leaves by farthest?
		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop(); // pop first of queue
			consumer.accept(live);
			if (live.parent != null && seen.add(live.parent)) {
				stack.add(live.parent); // insert at end of queue
			}
		}
	}

	public void iteratorDfUp(Consumer<MedialDisk> consumer) {
		// start at leaf, add parents until either root node, or next bifurcating node
		// when stack not empty, stop when parent has been seen, and pop next leaf.
		Set<MedialDisk> seen = new HashSet<>(getDisks().size());
		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		getLeaves().sort((a, b) -> Double.compare(a.distance, b.distance)); // NOTE
		getLeaves().forEach(stack::push);
		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop(); // pop first of queue
			consumer.accept(live);
			if (live.parent != null && seen.add(live.parent)) {
				stack.push(live.parent); // insert at start of queue
			}
		}

	}

	/**
	 * Is disk A a parent of disk B?
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	private boolean isParentOf(MedialDisk a, MedialDisk b) {
		return false; // TODO
	}

	private boolean isChildOf(MedialDisk a, MedialDisk b) {
		return false; // TODO
	}

	private boolean isRelated(MedialDisk a, MedialDisk b) {
		return false; // TODO
	}

	private void longestPath() {
		/*-
		 * TODO
		 * find the two disks that make up the longest path of the medial axis.
		 * the path must pass through the root node.
		 * To compute, find non-related (i.e. not child/parent) disk pairs and sum their distance values. 
		 * 		Pair with largest value denotes longest path.
		 * Group the disks into 3 groups, depending on which of the three root branches they belong to,
		 * 		then use pairs from distinct groups. 
		 */
	}

	/**
	 * Computes a random point that lies on the medial axis. Not to be confused with
	 * getting the center point from a random medial disk.
	 * 
	 * @return
	 */
	private Coordinate getRandomPoint() {
		// pick random edge, weighted by its length, then random point on that edge
		return null;
	}

	/**
	 * 90% of CPU time is here
	 * 
	 * @return
	 */
	private TC prepareTin() {
		IncrementalTin tin = new IncrementalTin(10);
//		tin.getNeighborhoodPointsCollector().

		// Triangulate geometry coordinates
		Coordinate[] coords = g.getCoordinates();
		final ArrayList<Vertex> vertices = new ArrayList<>();
		for (int i = 0; i < coords.length - 1; i++) {
			vertices.add(new Vertex(coords[i].x, coords[i].y, 0));
		}
		HilbertSort hs = new HilbertSort();
		hs.sort(vertices);
		tin.add(vertices, null); // insert point set; points are triangulated upon insertion
		computeConstraints(tin);

		TC tc = new TC(tin);
		TriangleCollector.visitSimpleTriangles(tin, tc);
		return tc;
	}

	private void computeConstraints(IIncrementalTin tin) {
		List<IConstraint> constraints = new ArrayList<>();
		for (int n = 0; n < g.getNumGeometries(); n++) {
			boolean hole = false;

			LinearRingIterator lri = new LinearRingIterator(g.getGeometryN(n));
			for (LinearRing ring : lri) {
				ArrayList<Vertex> points = new ArrayList<>();
				Coordinate[] c = ring.getCoordinates();
				if (c.length == 0) {
					continue;
				}
				if (Orientation.isCCW(c) && !hole) {
					for (int i = 0; i < c.length; i++) {
						points.add(new Vertex(c[i].x, c[i].y, 0));
					}
				} else {
					// holes are CW; region to keep lies to left the of constraints
					for (int i = c.length - 1; i >= 0; i--) {
						points.add(new Vertex(c[i].x, c[i].y, 0));
					}
				}
				constraints.add(new PolygonConstraint(points));
				hole = true; // all rings except the first are holes
			}
		}
		if (!constraints.isEmpty()) {
			tin.addConstraints(constraints, false); // true/false is negligible?
		}
	}

	private void calculateDepthFirstIndices() {
		/**
		 * calculate index of nodes visited When tree is walked in DFS manner.
		 */
	}

	/**
	 * @deprecated in favour of recursive method
	 */
	private void calculateFeatureArea() {

		/**
		 * Starting at leaves, sum each branch, stopping when reach bifurcating node;
		 * add bifurcation node to array and use that next time.
		 */

		List<MedialDisk> nodes = getLeaves();
		HashSet<MedialDisk> nextLeaves = new HashSet<>(); // what to begin with next iteration
		HashMap<MedialDisk, Integer> waitAtNodes = new HashMap<>();

		getBifurcations().forEach(n -> waitAtNodes.put(n, 0));

		while (!nodes.isEmpty()) {

			for (MedialDisk node : nodes) {

				MedialDisk current = node;

				while (current.degree < 2) { // FIXME what happens when start with bifurcating node
					current.parent.featureArea += current.featureArea;
					current = current.parent;
				}

				// current is bifurcating node

				final int childrenReached = waitAtNodes.merge(current, 1, Integer::sum);

				// exit when reach parent bifurcating node

				/**
				 * Only propagate from a bifurcating node one all its child paths have reached
				 * it.
				 */
				if (childrenReached == current.degree && current.parent != null) { // skip root node
					/**
					 * Don't add parent node until every child path to it has reached it (use
					 * counter per bifurcating node?)
					 */

					// skip (current.degree < 2) being false
					current.parent.featureArea += current.featureArea;
					current = current.parent;
					nextLeaves.add(current);
				}

			}
			nodes = new ArrayList<>(nextLeaves);
			nextLeaves.clear();
		}
	}

	public void calcFeatureArea() {
		if (rootNode.area == rootNode.featureArea) { // lazy compute
			recurseFeatureArea(rootNode);
		}
	}

	private static double recurseFeatureArea(MedialDisk node) {
		// TODO do as part of initial build?
		if (node.degree == 0) {
			return node.area;
		}
		for (MedialDisk child : node.children) {
			node.featureArea += recurseFeatureArea(child);
		}
		return node.featureArea;
	}

	public void drawVDM(PApplet p) {
		// looks different when drawing with DFS vs BFS
		p.pushStyle();
		p.strokeWeight(5);
		p.colorMode(PConstants.HSB, 1, 1, 1, 1);

		// limit using depth or vds index to grow skeleton

		final float depth = deepestNode.depthBF;

		iteratorBF(live -> {
			live.children.forEach(child -> {
				p.stroke((live.depthBF / depth) * .8f, 1, 1, 0.8f);
				p.line((float) live.position.x, (float) live.position.y, (float) child.position.x, (float) child.position.y);
			});
		});
		p.popStyle();
	}

	/**
	 * 
	 * @param p
	 * @param maxDepth inclusive
	 */
	public void drawVDM(PApplet p, int maxDepth) {

		// iterate by DFS for easier break

		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.push(rootNode);
		p.strokeWeight(5);
		p.colorMode(PConstants.HSB, 1, 1, 1, 1);

		// limit using depth or vds index to grow skeleton

		final float depth = deepestNode.depthBF;

		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop();

//			Vertex v = live.t.getVertexA();
//			p.beginShape();
////			p.noStroke();
//			p.strokeWeight(1);
//			p.stroke(live.depthBF / depth, 1, 1, 0.8f);
//			p.fill(live.depthBF / depth, 1, 1);
//			p.vertex((float) v.x, (float) v.y);
//			v = live.t.getVertexB();
//			p.vertex((float) v.x, (float) v.y);
//			v = live.t.getVertexC();
//			p.vertex((float) v.x, (float) v.y);
//			p.endShape();

			if (live.depthBF < maxDepth) {
				live.children.forEach(child -> {
					p.stroke(live.depthBF / depth, 1, 1, 0.8f);
					p.line((float) live.position.x, (float) live.position.y, (float) child.position.x, (float) child.position.y);
					stack.push(child); // DFS
				});
			}

		}
		p.colorMode(PConstants.RGB, 255, 255, 255, 255);
	}

	/**
	 * 
	 * @param p
	 * @param fraction 0...1
	 */
	public void drawVDM(PApplet p, double fraction) {
		drawVDM(p, (int) (fraction * deepestNode.depthBF));
	}

	/**
	 * Prune disks with feature area that is smaller in area than the given
	 * significance threshold. The significance value of a feature can be determined
	 * by summing the areas of all triangles in the subtree associated with that
	 * feature. Any subtree that has a value below the threshold is pruned. Each
	 * feature can be seen as being supported by the branches of the subtree, so
	 * when the subtree is pruned, the feature is eliminated.
	 * 
	 * @param p
	 * @param threshold percentage of the total area of the object 0...1 (where 0
	 *                  includes everything and 1 root node only)
	 */
	public void drawVDMPrune(PApplet p, double threshold) {

		// square threshold to make the effect more linear TODO different function?
		final double areaLimit = rootNode.featureArea * threshold * threshold * threshold;

		// iterate by DFS for easier break

		ArrayDeque<MedialDisk> stack = new ArrayDeque<>();
		stack.push(rootNode);
		p.strokeWeight(5);
		p.colorMode(PConstants.HSB, 1, 1, 1, 1);

		// limit using depth or vds index to grow skeleton

		final float depth = deepestNode.depthBF;

		while (!stack.isEmpty()) {
			MedialDisk live = stack.pop();

			live.children.forEach(child -> {
				if (child.featureArea > areaLimit) {
					p.stroke(live.depthBF / depth, 1, 1, 0.8f);
					p.line((float) live.position.x, (float) live.position.y, (float) child.position.x, (float) child.position.y);
					stack.push(child); // DFS
				}
			});
		}
		p.colorMode(PConstants.RGB, 255, 255, 255, 255);
	}

	private class TC implements Consumer<SimpleTriangle> {

		SimpleTriangle largestDiskTriangle = null; // triangle with largest circumcircle
		double[] largestCircumcircle = new double[3]; // largest inscribed circle
		HashMap<IQuadEdge, SimpleTriangle> map; // TODO use SparseArray?

		/**
		 * Map of triangles to their circumcircle [x,y,r]
		 */
		HashMap<SimpleTriangle, double[]> cccMap = new HashMap<>();
		ArrayList<Coordinate> coords = new ArrayList<>();
		HashMap<SimpleTriangle, Coordinate> coord = new HashMap<>();
		ArrayList<SimpleTriangle> triangles = new ArrayList<>();

		public TC(IncrementalTin tin) {
			map = new HashMap<>(tin.getEdges().size()); // *3?
		}

		@Override
		public void accept(SimpleTriangle t) {
			double[] circumcircle = circumcircle(t.getVertexA(), t.getVertexB(), t.getVertexC());

			Coordinate center = new Coordinate(circumcircle[0], circumcircle[1]);

			/**
			 * If geometry contains triangle circumcircle center then include it in medial
			 * axis
			 */

//			if (pointLocator.locate(center) == Location.INTERIOR) {
			IConstraint constraint = t.getContainingRegion();
			if (constraint != null && constraint.definesConstrainedRegion()) {
				// TODO use rasterised PShape?

				coords.add(center);
				coord.put(t, center);

				triangles.add(t);

				map.put(t.getEdgeA(), t);
				map.put(t.getEdgeB(), t);
				map.put(t.getEdgeC(), t);
				cccMap.put(t, circumcircle);

				if (circumcircle[2] > largestCircumcircle[2]) {
					largestDiskTriangle = t;
					largestCircumcircle = circumcircle;
				}
			}
		}
	}

	/**
	 * Voronoi Disk. medial disk? maximal disk? Also inscribed disks.
	 * 
	 * <p>
	 * A voronoi disk represents each node in the medial axis and has a
	 * corresponding triangle in the underlying Delaunay triangulation. A disk is
	 * the shape's maximal-inscribed circle located at the disk's center point -- a
	 * disk's radius is the distance to the nearest point on the shape edge from the
	 * disk center point.
	 * 
	 * @author Michael Carleton
	 *
	 */
	public class MedialDisk {

		/** The underlying delaunay triangle associated with this disk */
		public SimpleTriangle t;
		/** This disk's parent node. Null if root node */
		public MedialDisk parent;
		/**
		 * This disk's children nodes. Nodes have upto 3 children. Leaf nodes have 0
		 * children.
		 */
		public List<MedialDisk> children;
		public final int depthBF; // breadth-first depth from the root node, better for drawing. // TODO
		public int depthDF; // distance from the root, length of the path from n to the root TODO

//		public double[] circumcircle;
		/** The number of children / outdegree */
		public int degree = 0;
		/** The area of the underlying delaunay triangle associated with this disk */
		double area;
		/** The sum of triangle areas of this disk and all its descendants */
		public double featureArea;

		/** euclidean shortest path distance from root node circumcircle */
		public double distance = 0; // TODO

//		"branch"/"path" b // branch is the set of VDs between two nodes of degree > 1, or between a node of degree >1 and a leaf
		// branch should point to start and end

		boolean forkChild = false; // is direct child of a bifurcating disk OR "isbranchParent"
		boolean forkParent = false; // (is child a node of degree >1?)

		/**
		 * Measures the change in the width of the shape per unit length of the axis
		 * (edge segment). When positive, this part of the object widens as we progress
		 * along the axis branch; if it’s negative the part narrows. The gradient for a
		 * given disk is measured using its parent node (i.e. how the gradient changes
		 * as we move towards this disk), rather than child (how gradient changes as we
		 * move away from this disk).
		 */
		public final double axialGradient; // axial gradient between this and its parent. (rChild-rParent/d)

		/** Centerpoint of this disk */
		public final Coordinate position;
		/**
		 * The radius of the circumcircle of the disk's underlying triangle. Also
		 * equivalent to the distance between the disk center and the closest point on
		 * the shape
		 */
		public final double radius;

		public final int id; // BFS id from root node, unique for each disk (unlike depthDF)

		private int bifurcationDepth; // TODO, determines how many bifurcations ("hops") are between this disk and the
										// root

		/**
		 * the "kink" angle (0...PI) this disk makes with its parent and child node (if
		 * it has only one child). if the node has 0 or more than 1 child, this value is
		 * NaN
		 */
		private double angle = Double.NaN; // TODO initially unknown until axis is completed.

		MedialDisk(MedialDisk parent, int id, SimpleTriangle t, double[] circumcircle, int depthBF) {
			this.parent = parent;
			this.id = id;
			this.t = t;
			children = new ArrayList<>(3);
			this.depthBF = depthBF;
			area = t.getArea();
			featureArea = area;
			position = new Coordinate(circumcircle[0], circumcircle[1]);
			radius = circumcircle[2];
			if (parent != null) {
				axialGradient = (radius - parent.radius) / distance(position, parent.position);
				minimumAxialGradient = Math.min(minimumAxialGradient, axialGradient);
				maximumAxialGradient = Math.max(maximumAxialGradient, axialGradient);
			} else { // root node
				axialGradient = 0;
			}
		}

		/**
		 * Add a child disk to this disk and increment its degree.
		 */
		private void addchild(MedialDisk child) {
			children.add(child);
			degree++;
		}

		boolean isLeaf() {
			return degree == 0;
		}

		@Override
		public int hashCode() {
			// or cantour pairing, or (int) double.tolongbits
			return t.getEdgeA().hashCode(); // unique for a given tin
		}

		@Override
		public boolean equals(Object obj) {
			if (obj instanceof MedialDisk) {
				MedialDisk other = (MedialDisk) obj;
				return other.t.equals(t);
			}
			return false;
		}
	}

	/**
	 * Edges are straight-line segments that connect two adjacent disks.
	 */
	public class Edge {

		public final MedialDisk head;
		public final MedialDisk tail;
		int depth;
		public final double axialGradient;
		public final LineString lineString;

		Edge(MedialDisk head, MedialDisk tail) {
			this.head = head;
			this.tail = tail;
			lineString = GEOM_FACTORY.createLineString(new Coordinate[] { head.position, tail.position });
			depth = head.depthBF;
			axialGradient = tail.axialGradient;
		}

		/**
		 * Hashcode for this edge.
		 */
		@Override
		public int hashCode() {
			return head.t.getEdgeA().hashCode() * (int) (tail.position.y + 1); // TODO check
		}

		@Override
		public boolean equals(Object obj) {
			if (obj instanceof Edge) {
				return hashCode() == ((Edge) obj).hashCode();
			}
			return false;
		}
	}

	/**
	 * Branches are series of successive disks found in between. Branches terminate
	 * when a disk is a leaf, or when a disk bifurcates. branch/trunk branch
	 * 
	 * comparable based on fork degree?
	 *
	 */
	public class Branch {

		public final MedialDisk root; // root disk
		public MedialDisk leaf; // leaf disk
		/**
		 * Disks between root and leaf of this branch, in descending order from the
		 * branch's root.
		 */
		public List<MedialDisk> innerDisks;
		/**
		 * most branches have one sibling; 3 branches have 2 siblings (trifurcating
		 * origin)
		 */
		Branch sibling; // branch that shares parent disk
		// contains bezier interpolation
		int forkDegree; // how many forks are visited from the rootnode to the root of this branch
		double length; // sum of edge lengths (not count)
		public final List<Edge> edges;
		public LineString lineString; // TODO

		/**
		 * 
		 * @param root a tri/bifurcating disk
		 */
		Branch(MedialDisk root) {
			this.root = root;
			innerDisks = new ArrayList<>();
			edges = new ArrayList<>();
		}

		void add(MedialDisk disk) {
			if (innerDisks.isEmpty()) {
				edges.add(new Edge(root, disk));
			} else {
				edges.add(new Edge(innerDisks.get(innerDisks.size() - 1), disk));
			}
			innerDisks.add(disk);
		}

		/**
		 * Add the last (leaf) disk for the branch
		 * 
		 * @param disk
		 */
		void end(MedialDisk disk) {
			leaf = disk;
			if (innerDisks.isEmpty()) { // bifurcations that link directly to other bifucations
				edges.add(new Edge(root, disk));
			} else {
				edges.add(new Edge(innerDisks.get(innerDisks.size() - 1), disk));
			}
			Coordinate[] coords = new Coordinate[innerDisks.size() + 2];
			coords[0] = root.position;
			for (int i = 1; i < innerDisks.size(); i++) {
				coords[i] = innerDisks.get(i - 1).position;
			}
			coords[coords.length - 1] = leaf.position;
			lineString = GEOM_FACTORY.createLineString(coords);
		}

		void smooth(int smoothingType) {
			// TODO
		}

		/**
		 * Whether this branch terminates at a leaf node.
		 */
		public boolean terminates() {
			return root.isLeaf();
		}
	}

	/**
	 * @return double[3] of [x, y, radius]
	 */
	private static double[] circumcircle(Vertex a, Vertex b, Vertex c) {

		double D = (a.getX() - c.getX()) * (b.getY() - c.getY()) - (b.getX() - c.getX()) * (a.getY() - c.getY());
		double px = (((a.getX() - c.getX()) * (a.getX() + c.getX()) + (a.getY() - c.getY()) * (a.getY() + c.getY())) / 2
				* (b.getY() - c.getY())
				- ((b.getX() - c.getX()) * (b.getX() + c.getX()) + (b.getY() - c.getY()) * (b.getY() + c.getY())) / 2
						* (a.getY() - c.getY()))
				/ D;

		double py = (((b.getX() - c.getX()) * (b.getX() + c.getX()) + (b.getY() - c.getY()) * (b.getY() + c.getY())) / 2
				* (a.getX() - c.getX())
				- ((a.getX() - c.getX()) * (a.getX() + c.getX()) + (a.getY() - c.getY()) * (a.getY() + c.getY())) / 2
						* (b.getX() - c.getX()))
				/ D;

		double rs = (c.getX() - px) * (c.getX() - px) + (c.getY() - py) * (c.getY() - py);

		return new double[] { px, py, Math.sqrt(rs) };
	}

	private static double distance(Coordinate a, Coordinate b) {
		double deltaX = a.y - b.y;
		double deltaY = a.x - b.x;
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

}
