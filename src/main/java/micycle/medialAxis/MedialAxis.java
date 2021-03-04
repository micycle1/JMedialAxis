package micycle.medialAxis;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.operation.linemerge.LineMergeGraph;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;
import org.tinspin.index.kdtree.KDTree;

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
 * For geometries of genus zero, the graph is naturally a tree. For objects with
 * holes, we break each cycle by imposing an appropriate breakpoint in the loop.
 * This is only for the purposes of traversing the graph without running into
 * infinite loops. In order to preserve the topology of the original object,
 * cycles are never pruned
 * <p>
 * Decompose the disk coordinates into curves.
 * <p>
 * If angle of successive line segments is the same, merge them
 * <p>
 * * One way to get a smoother initial axis (i.e. without pruning) is to apply
 * on a smoothed version of the shape.
 * 
 * @author Michael Carleton
 *
 */
public class MedialAxis {

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(
			new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	private Geometry g;
//	private final PreparedGeometry cache;
	private KDTree<VD> kdTree;
//	private IncrementalTin tin;
//	private TC tc;
	private Geometry dissolved;

	ArrayList<VD> voronoiDisks = new ArrayList<>();
	public VD rootNode; // the medial axis has a trifuraction at the root
	VD deepestNode;
	private List<VD> leaves;
	private List<VD> bifurcations; // birfurcating/forking nodes
	private List<Segment> segments;

//	public ArrayList<ArrayList<VD>> segments = new ArrayList<>(); // node sections between forks

	public MedialAxis(Coordinate[] coordinates) {
		// TODO support holes
		this(GEOM_FACTORY.createLinearRing(coordinates));
	}

//	public MedialAxis(List<Coordinate> coordinates) {
//		this(GEOM_FACTORY.createLinearRing(coordinates));
//	}

	private final IndexedPointInAreaLocator pointLocator;

	/**
	 * The medial axis is computed for the geometry during construction.
	 * 
	 * @param g
	 */
	public MedialAxis(Geometry g) {
//		if (!g.isSimple()) {
//			throw new Exception("Invalid Geometry -- it is not simple (may self-intersect or have repeated points.");
//			return;
//		}

		/**
		 * Simplify if complex
		 */
//		System.out.println("before " + g.getCoordinates().length);
		if (g.getCoordinates().length > 2000) {
			Polygon poly = (Polygon) g;
			if (poly.getNumInteriorRing() > 0) {
				g = TopologyPreservingSimplifier.simplify(g, 1);
			} else {
				g = DouglasPeuckerSimplifier.simplify(g, 1);
			}
		}
//		System.out.println("w/simplify " + g.getCoordinates().length);

		pointLocator = new IndexedPointInAreaLocator(g);
		this.g = Densifier.densify(g, 10); // NOTE constant=10
//		System.out.println("after " + this.g.getCoordinates().length);

		TC tc = prepareTin();

		// having done triangulation, build up directed tree of medial axis
		rootNode = new VD(null, 0, tc.largestDiskTriangle, tc.largestCircumcircle, 0);
		deepestNode = rootNode;
		ArrayDeque<VD> stack = new ArrayDeque<>(15);
		stack.push(rootNode); // start tree with rootnode

		HashSet<SimpleTriangle> remaining = new HashSet<>();
		remaining.addAll(tc.triangles);
		remaining.remove(tc.largestDiskTriangle);

		int depth = 1; // number of edges in the path from the root to the node
		int id = 1; // breadth-first ID

		while (!stack.isEmpty()) {

			VD live = stack.pop(); // FIFO
			voronoiDisks.add(live);

			// get half-edge duals (shared with other triangles)
			final SimpleTriangle n1 = tc.map.get(live.t.getEdgeA().getDual());
			final SimpleTriangle n2 = tc.map.get(live.t.getEdgeB().getDual());
			final SimpleTriangle n3 = tc.map.get(live.t.getEdgeC().getDual());

			if (remaining.contains(n1)) {
				VD child = new VD(live, id++, n1, tc.cccMap.get(n1), depth);
				stack.add(child);
				live.addchild(child);
				remaining.remove(n1);
			}
			if (remaining.contains(n2)) {
				VD child = new VD(live, id++, n2, tc.cccMap.get(n2), depth);
				stack.add(child);
				live.addchild(child);
				remaining.remove(n2);
			}
			if (remaining.contains(n3)) {
				VD child = new VD(live, id++, n3, tc.cccMap.get(n3), depth);
				stack.add(child);
				live.addchild(child);
				remaining.remove(n3);
			}

			depth++;
		}
		deepestNode = voronoiDisks.get(voronoiDisks.size() - 1); // TODO check
//		calculateFeatureArea();
		calcFeatureArea();
//		System.out.println("total area: " + totalArea);
//		System.out.println(rootNode.featureArea);
	}

	public List<VD> getDisks() {
		return voronoiDisks;
	}

	public VD nearestDisk(double x, double y) {
		// TODO use cover tree?
		if (kdTree == null) { // Lazy initialisation
			kdTree = KDTree.create(2);
			voronoiDisks.forEach(disk -> kdTree.insert(new double[] { disk.position.x, disk.position.y }, disk));
		}
		return kdTree.nnQuery(new double[] { x, y }).value();
	}

	/**
	 * Returns the children of a disk into a linear array
	 * 
	 * @param parent contains chil out: contains children and the parent (at
	 *               position 0)
	 * @return
	 */
	public List<VD> getDescendants(VD parent) {

		ArrayDeque<VD> stack = new ArrayDeque<>();
		stack.add(parent);

		ArrayList<VD> out = new ArrayList<>();

		while (!stack.isEmpty()) {
			VD live = stack.pop();
			out.add(live);
			live.children.forEach(child -> {
				stack.add(child);
			});
		}
		return out;
	}

	public List<VD> getDescendants(VD parent, int maxDepth) {

		ArrayDeque<VD> stack = new ArrayDeque<>();
		stack.add(parent);

		ArrayList<VD> out = new ArrayList<>();
		int depth = 0;

		while (!stack.isEmpty() && depth < maxDepth) {
			VD live = stack.pop();
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
	public List<VD> getAncestors(VD child) {

		ArrayList<VD> out = new ArrayList<>();

		VD live = child;
		while (live.parent != null) {
			out.add(live);
			live = live.parent;
		}
		return out;
	}

	public void getDissolvedGeometry() {
		LineDissolver ld = new LineDissolver();
		// dissolved geom
	}

	public void getLineMergeGraph() {
		LineMergeGraph lmGraph = new LineMergeGraph();
	}

	/**
	 * End nodes / A VDs of degree 1/no children
	 * 
	 * @return
	 */
	public List<VD> getLeaves() {
		if (leaves == null) { // Lazy initialisation
			leaves = voronoiDisks.stream().filter(vd -> vd.degree == 0).collect(Collectors.toList());
		}
		return leaves;
	}

	/**
	 * Nodes/VDs with two descendent lineages (two children)
	 * 
	 * @return
	 */
	public List<VD> getBifurcations() {
		if (bifurcations == null) { // Lazy initialisation
			bifurcations = voronoiDisks.stream().filter(vd -> vd.degree == 2).collect(Collectors.toList());
		}
		return bifurcations;
	}

	public void getEdges() {
	}

	/**
	 * Segments are sets of points between two forking disks or one forking disk and
	 * a leaf. Segments include parent bifurcating node
	 * 
	 * @return
	 */
	public List<ArrayList<VD>> getSegments() {

		/**
		 * During construction we use stack.add() to naviagate nodes in a BFS manner.
		 * Here we use push() to do so so in a DFS manner.
		 */
		ArrayList<ArrayList<VD>> segments = new ArrayList<>();
		ArrayList<VD> segment;
		ArrayList<VD> forks = new ArrayList<>(getBifurcations());
		forks.add(rootNode);

		for (VD disk : forks) {

			for (VD child : disk.children) {
				segment = new ArrayList<>();
				segment.add(disk); // add bifurcating parent
				VD node = child;

				while (node.degree == 1) {
					segment.add(node);
					node = node.children.get(0);
				}
				segment.add(node);
				segments.add(segment);
			}
		}
		return segments;
	}

	public List<Segment> getSegmentsAsSegments() {

		// TODO test/use

		if (segments == null) { // Lazy initialisation
			segments = new ArrayList<>();

			/**
			 * During construction we use stack.add() to naviagate nodes in a BFS manner.
			 * Here we use push() to do so so in a DFS manner.
			 */
			Segment segment;
			ArrayList<VD> forks = new ArrayList<>(getBifurcations());
			forks.add(rootNode);

			for (VD disk : forks) {

				segment = new Segment(disk);
				for (VD child : disk.children) {

					segment.add(disk); // add bifurcating parent
					VD node = child;

					while (node.degree == 1) {
						segment.add(node);
						node = node.children.get(0);
					}
					segment.add(node);
					segments.add(segment);
				}
			}
		}

		return this.segments;
	}

	/**
	 * 90% of CPU time is here
	 * 
	 * @return
	 */
	private TC prepareTin() {
		IncrementalTin tin = new IncrementalTin();

		// Triangulate geometry coordinates
		Coordinate[] coords = g.getCoordinates();
		final ArrayList<Vertex> vertices = new ArrayList<>();
		for (int i = 0; i < coords.length; i++) {
			vertices.add(new Vertex(coords[i].x, coords[i].y, 0));
		}
		tin.add(vertices, null); // insert point set; points are triangulated upon insertion

		TC tc = new TC(tin);
		TriangleCollector.visitSimpleTriangles(tin, tc);

		return tc;
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
		 * Starting at leaves, sum each segment, stopping when reach bifurcating node;
		 * add bifurcation node to array and use that next time.
		 */

		List<VD> nodes = getLeaves();
		HashSet<VD> nextLeaves = new HashSet<>(); // what to begin with next iteration
		HashMap<VD, Integer> waitAtNodes = new HashMap<>();

		getBifurcations().forEach(n -> waitAtNodes.put(n, 0));

		while (!nodes.isEmpty()) {

			for (VD node : nodes) {

				VD current = node;

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
			nodes = new ArrayList<VD>(nextLeaves);
			nextLeaves.clear();
		}
	}

	private void calcFeatureArea() {
		if (rootNode.area == rootNode.featureArea) { // lazy compute
			recurseFeatureArea(rootNode);
		}
	}

	private static double recurseFeatureArea(VD node) {
		if (node.degree == 0) {
			return node.area;
		}
		for (VD child : node.children) {
			node.featureArea += recurseFeatureArea(child);
		}
		return node.featureArea;
	}

	public void drawVDM(PApplet p) {

		// looks different when drawing with DFS vs BFS

		ArrayDeque<VD> stack = new ArrayDeque<>();
		stack.push(rootNode);
		p.strokeWeight(5);
		p.colorMode(PConstants.HSB, 1, 1, 1, 1);

		// limit using depth or vds index to grow skeleton

		final float depth = deepestNode.depthBF;

		while (!stack.isEmpty()) {
			VD live = stack.pop();

			live.children.forEach(child -> {
				p.stroke((live.depthBF / depth) * .8f, 1, 1, 0.8f);
				p.line((float) live.position.x, (float) live.position.y, (float) child.position.x,
						(float) child.position.y);
				stack.push(child);
			});
		}
		p.colorMode(PConstants.RGB, 255, 255, 255, 255);
	}

	/**
	 * 
	 * @param p
	 * @param maxDepth inclusive
	 */
	public void drawVDM(PApplet p, int maxDepth) {

		// iterate by DFS for easier break

		ArrayDeque<VD> stack = new ArrayDeque<>();
		stack.push(rootNode);
		p.strokeWeight(5);
		p.colorMode(PConstants.HSB, 1, 1, 1, 1);

		// limit using depth or vds index to grow skeleton

		final float depth = deepestNode.depthBF;

		while (!stack.isEmpty()) {
			VD live = stack.pop();

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
					p.line((float) live.position.x, (float) live.position.y, (float) child.position.x,
							(float) child.position.y);
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
	 * @param threshold 0...1 (where 0 includes everything and 1 root node only)
	 */
	public void drawVDMPrune(PApplet p, double threshold) {

		// square threshold to make the effect more linear TODO different function?
		final double areaLimit = rootNode.featureArea * threshold * threshold * threshold;

		// iterate by DFS for easier break

		ArrayDeque<VD> stack = new ArrayDeque<>();
		stack.push(rootNode);
		p.strokeWeight(5);
		p.colorMode(PConstants.HSB, 1, 1, 1, 1);

		// limit using depth or vds index to grow skeleton

		final float depth = deepestNode.depthBF;

		while (!stack.isEmpty()) {
			VD live = stack.pop();

			live.children.forEach(child -> {
				if (child.featureArea > areaLimit) {
					p.stroke(live.depthBF / depth, 1, 1, 0.8f);
					p.line((float) live.position.x, (float) live.position.y, (float) child.position.x,
							(float) child.position.y);
					stack.push(child); // DFS
				}
			});
		}
		p.colorMode(PConstants.RGB, 255, 255, 255, 255);
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

			if (pointLocator.locate(center) == Location.INTERIOR) {
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
	 * corresponding triangle in the underlying Delaunay triangulation.
	 * 
	 * @author MCarleton
	 *
	 */
	public static class VD {

		/** The underlying delaunay triangle associated with this disk */
		public SimpleTriangle t;
		/** This disk's parent node. Null if root node */
		public VD parent;
		/**
		 * This disk's children nodes. Nodes have upto 3 children. Leaf nodes have 0
		 * children.
		 */
		ArrayList<VD> children;
		public final int depthBF; // breadth-first depth from the root node, better for drawing. // TODO
		public int depthDF; // distance from the root, length of the path from n to the root TODO

//		public double[] circumcircle;
		/** The number of children / outdegree */
		public int degree = 0;
		/** The area of the underlying delaunay triangle associated with this disk */
		double area;
		/** The sum of triangle areas of this disk and all its descendants */
		public double featureArea;

		VD segmentParent;

//		"segment"/"path" b // segment is the set of VDs between two nodes of degree > 1, or between a node of degree >1 and a leaf
		// segment should point to start and end

		boolean forkChild = false; // is direct child of a bifurcating disk OR "issegmentParent"
		boolean forkParent = false; // (is child a node of degree >1?)

		/**
		 * Measures the change in the width of the shape per unit length of the axis.If
		 * positive, this part of the object widens as we progress along the axis
		 * segment; if it’s negative the part narrows
		 */
		final double axialGradient; // axial gradient between this and its parent. (rChild-rParent/d)

		public final Coordinate position;
		/** The radius of the circumcircle of the disk's underlying triangle */
		public final double radius;

		final int id; // BFS id from root node, unique (unlike depthDF)

		VD(VD parent, int id, SimpleTriangle t, double[] circumcircle, int depthBF) {
			this.parent = parent;
			this.id = id;
			this.t = t;
			children = new ArrayList<>(3);
			this.depthBF = depthBF;
//			this.circumcircle = circumcircle;
			area = t.getArea();
			featureArea = area;
			position = new Coordinate(circumcircle[0], circumcircle[1]);
			radius = circumcircle[2];
			if (parent != null) {
				axialGradient = (radius - parent.radius) / distance(position, parent.position);
			} else {
				axialGradient = 0;
			}
		}

		/**
		 * Add a child disk to this disk and increment its degree.
		 */
		private void addchild(VD child) {
			children.add(child);
			degree++;
		}

		boolean isLeaf() {
			return children.size() == 0;
		}

		@Override
		public int hashCode() {
			// or cantour pairing, or (int) double.tolongbits
			return t.getEdgeA().hashCode();
		}
	}

	/**
	 * Edges model links between two disks. Return graphs based on merged edges
	 * 
	 * @author MCarleton
	 *
	 */
	public static class Edge {

		VD head;
		VD tail;
		int depth = head.depthBF;

		public Edge() {
			// TODO Auto-generated constructor stub
		}
	}

	/**
	 * segment/trunk segment
	 * 
	 * comparable based on fork degree?
	 *
	 */
	public static class Segment {
		VD root; // root disk
		VD leaf; // leaft disk
		/** In descending order from the segment root */
		private List<VD> innerDisks;
		LineString lineString;
		Segment sibling; // segment that shares parent disk
		// contains bezier interpolation
		int forkDegree; // how many forks are visited from the rootnode to the root of this segment
		double length; // sum of edge lengths (not count)

		public Segment(VD root) {
			this.root = root;
			innerDisks = new ArrayList<>();
//			length = root.length;
		}

		public void add(VD disk) {
			innerDisks.add(disk);
//			length++;
		}
	}

	private static double distance(Coordinate a, Coordinate b) {
		double deltaX = a.y - b.y;
		double deltaY = a.x - b.x;
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

}
