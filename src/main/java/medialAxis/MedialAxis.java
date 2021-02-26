package medialAxis;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
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
 * https://www.cs.ubc.ca/labs/imager/th/2003/Tam2004/Tam2004.pdf
 * 
 * Produces a directed acyclic graph of the MD / more specifically rooted tree
 * 
 * Medial Axis library models medial axes of shapes as a (rooted) tree of medial
 * disks. Medial disks reference a parent and have upto 3 children (the medial
 * axis bifurcates/forks at disks with 2 children). The root is the disk with
 * the largest circumcirle, also the single disk that trifurcates. All disks
 * have an indegree of 1.
 * 
 * For geometries of genus zero, the graph is naturally a tree. For objects with
 * holes, we break each cycle by imposing an appropriate breakpoint in the loop.
 * This is only for the purposes of traversing the graph without running into
 * infinite loops. In order to preserve the topology of the original object,
 * cycles are never pruned
 * 
 * Decompose the disk coordinates into curves.
 * 
 * If angle of successive line segments is the same, merge them
 * 
 * @author Michael Carleton
 *
 */
public class MedialAxis {

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(1));

	private Geometry g;
	private final PreparedGeometry cache;
	private KDTree<VD> kdTree;
//	private IncrementalTin tin;
//	private TC tc;
	private Geometry dissolved;

	ArrayList<VD> voronoiDisks = new ArrayList<>();
	VD rootNode; // the medial axis has a trifuraction at the root
	VD deepestNode;
	private List<VD> leaves;
	private List<VD> bifurcations; // birfurcating/forking nodes

	ArrayList<ArrayList<VD>> branches = new ArrayList<>(); // node sections between forks

	public MedialAxis(Coordinate[] coordinates) {
		this(GEOM_FACTORY.createLinearRing(coordinates));
	}

//	public MedialAxis(List<Coordinate> coordinates) {
//		this(GEOM_FACTORY.createLinearRing(coordinates));
//	}

	public MedialAxis(Geometry g) {
//		if (!g.isSimple()) {
//			throw new Exception("Invalid Geometry -- it is not simple (may self-intersect or have repeated points.");
//			return;
//		}

//		System.out.println("before " + g.getCoordinates().length);
		cache = PreparedGeometryFactory.prepare(g);
		this.g = Densifier.densify(g, 10).reverse(); // NOTE constant=10
//		System.out.println("after " + this.g.getCoordinates().length);

		TC tc = prepareTin();

		// build up directed tree of medial axis
		rootNode = new VD(null, tc.largestDiskTriangle, tc.largestCircumcircle, 0);
		deepestNode = rootNode;
		ArrayDeque<VD> stack = new ArrayDeque<>(15);
		stack.push(rootNode); // start tree with rootnode

		HashSet<SimpleTriangle> remaining = new HashSet<>();
		remaining.addAll(tc.triangles);
		remaining.remove(tc.largestDiskTriangle);

		int depth = 1; // number of edges in the path from the root to the node

		/**
		 * Branches are sets of disks between two nodes that are forks (or leaf)
		 */
		ArrayList<VD> branch = new ArrayList<>();
		VD lastBranchParent = rootNode;
		while (!stack.isEmpty()) {

			VD live = stack.pop();
			voronoiDisks.add(live);

			// get half-edge duals (shared with other triangles)
			final SimpleTriangle n1 = tc.map.get(live.t.getEdgeA().getDual());
			final SimpleTriangle n2 = tc.map.get(live.t.getEdgeB().getDual());
			final SimpleTriangle n3 = tc.map.get(live.t.getEdgeC().getDual());

			int children = 0;
			if (remaining.contains(n1)) {
				VD child = new VD(live, n1, tc.cccMap.get(n1), live.depthDF + 1);
//				child.branchParent = lastBranchParent;
				stack.push(child);
				live.addchild(child);
				remaining.remove(n1);
				children++;
			}
			if (remaining.contains(n2)) {
				VD child = new VD(live, n2, tc.cccMap.get(n2), live.depthDF + 1);
//				child.branchParent = lastBranchParent;
				stack.push(child);
				live.addchild(child);
				remaining.remove(n2);
				children++;
			}
			if (remaining.contains(n3)) {
				VD child = new VD(live, n3, tc.cccMap.get(n3), live.depthDF + 1);
//				child.branchParent = lastBranchParent;
				stack.push(child);
				live.addchild(child);
				remaining.remove(n3);
				children++;
			}

			if (live.forkChild == true) {
				if (branch.size() > 0) { // skip rootnode branch
					if (children > 1) {
						branch.add(live.children.get(0));
					}
					branches.add(new ArrayList<>(branch)); // branch.clone~?
				}
				branch = new ArrayList<>();
				branch.add(live.parent);
			}

			/**
			 * Means the VD we popped is a bifurcation node
			 */
			if (children > 1 || children == 0) {
//				System.out.println("Asdasd");
				lastBranchParent = live;
				for (VD child : live.children) {
					child.forkChild = true;
					// bool childOfFork
				}
			}

			if (children < 2) {
				branch.add(live);
				if (live.depthDF > deepestNode.depthDF) {
					deepestNode = live;
				}
			}
		}
	}

	public VD nearestDisk(double x, double y) {
		if (kdTree == null) { // Lazy initialisation
			kdTree = KDTree.create(2);
			voronoiDisks
					.forEach(disk -> kdTree.insert(new double[] { disk.circumcircle[0], disk.circumcircle[1] }, disk));
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
	 */
	public void getAncestors(VD child) {

		ArrayDeque<VD> stack = new ArrayDeque<>();
		stack.add(child);

		ArrayList<VD> out = new ArrayList<>();

		while (!stack.isEmpty()) {
			VD live = stack.pop();

			if (live.parent != null) {
			}
		}

	}

	public void getDissolvedGeometry() {
		// dissolved geom
	}

	public void getLineMergeGraph() {
		// LineMergeGraph
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

//	public HashMap<VD, List<VD>> getBranches() { // TODO rename
////		ArrayList<ArrayList<VD>> branches = new ArrayList<>();
////		ArrayList<VD> branch = new ArrayList<>();
//
//		HashMap<VD, List<VD>> branches = new HashMap<>();
//		getBifurcations().forEach(b -> branches.put(b, new ArrayList<VD>()));
//		branches.put(rootNode, new ArrayList<MedialAxis.VD>());
//
//		ArrayList<VD> x = (ArrayList<VD>) voronoiDisks.clone();
//		x.remove(rootNode);
//
//		x.forEach(vd -> {
//			branches.get(vd.branchParent).add(vd);
//		});
//
//		return branches;
//
//	}

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

	public void drawVDM(PApplet p) {

		// looks different when drawing with DFS vs BFS

		ArrayDeque<VD> stack = new ArrayDeque<>();
		stack.push(rootNode);
		p.strokeWeight(5);
		p.colorMode(PConstants.HSB, 1, 1, 1, 1);

		// limit using depth or vds index to grow skeleton

		final float depth = deepestNode.depthDF;

		while (!stack.isEmpty()) {
			VD live = stack.pop();

			live.children.forEach(child -> {
				p.stroke(live.depthDF / depth, 1, 1, 0.8f);
				p.line((float) live.circumcircle[0], (float) live.circumcircle[1], (float) child.circumcircle[0],
						(float) child.circumcircle[1]);
				stack.push(child);
			});
		}
		p.colorMode(PConstants.RGB, 255, 255, 255, 255);
	}

	void drawVDM(int maxDepth) {

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
			if (cache.covers(GEOM_FACTORY.createPoint(center))) {

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
	 * Voronoi Disk. medial disk? maximal disk?
	 * 
	 * <p>
	 * A voronoi disk represents each node in the medial axis and has a
	 * corresponding triangle in the Delaunay triangulation.
	 * 
	 * @author MCarleton
	 *
	 */
	public static class VD {

		public SimpleTriangle t;
		public VD parent; // null if root node
		ArrayList<VD> children;
		public final int depthDF; // distance from the root, length of the path from n to the root
		int depthBF; // breadth-first depth from the root node, better for drawing. // TODO
		public double[] circumcircle;
		double[] position;
		double radius;
		public int degree = 0; // number of children / outdegree
		double area; // area of underlying triangle only (compute once)
		double featureArea; // area of triangles associated with axis graph subtree of this disk |
							// "underlying area"?

		VD branchParent;

//		"branch"/"path" b // branch is the set of VDs between two nodes of degree > 1, or between a node of degree >1 and a leaf
		// branch should point to start and end

		boolean forkChild = false; // is direct child of a bifurcating disk OR "isBranchParent"
		boolean forkParent = false; // (is child a node of degree >1?)

		/**
		 * If positive, this part of the object widens as we progress along the axis
		 * segment; if it’s negative the part narrows
		 */
		double axialGradient; // axial gradient between this and its parent.

		VD(VD parent, SimpleTriangle t, double[] circumcircle, int depth) {
			this.parent = parent;
			this.t = t;
			children = new ArrayList<>(3);
			this.depthDF = depth;
			this.circumcircle = circumcircle;
			area = t.getArea();
			featureArea = area;
		}

		/**
		 * Add a child disk to this disk and increment its degree.
		 */
		private void addchild(VD child) {
			children.add(child);
			degree++;
		}

		boolean terminates() {
			// isendpoint
			return children.size() == 0;
		}

		@Override
		public int hashCode() {
			// or cantour pairing, or (int) double.tolongbits
			return t.getEdgeA().hashCode();
		}
	}

	public static class Edge {

		VD head;
		VD tail;
		int depth = head.depthDF;

		public Edge() {
			// TODO Auto-generated constructor stub
		}
	}

}
