import numpy as np
import math
from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional
import matplotlib.pyplot as plt

""" 2D """

@dataclass(slots=True)
class Vertex:
    x: float
    y: float
    edge: Optional["HalfEdge"] = field(default=None, repr=False)

    def __eq__(self, other):
        if isinstance(other, Vertex):
            return math.isclose(self.x, other.x) and math.isclose(self.y, other.y)
        
    def __add__(self, other):
        if isinstance(other, Vertex):
            return Vertex(self.x+other.x, self.y+other.y, edge=None)

    def __repr__(self):
        return f"Vertex(x={self.x:.4f},y={self.y:.4f})"

@dataclass
class HalfEdge:
    origin: Vertex = None
    twin: Optional["HalfEdge"] = field(default=None)
    next: Optional["HalfEdge"] = field(default=None)
    face: Optional["Face"] = field(default=None)

    def __post_init__(self):
        self.origin.edge = self

    @cached_property
    def length(self):
        assert self.origin != None
        assert self.next.origin != None
        lx = (self.origin.x-self.next.origin.x)**2
        ly = (self.origin.y-self.next.origin.y)**2
        return math.sqrt(lx+ly)

    def __eq__(self, other):
        if isinstance(other, HalfEdge):
            return self.origin == other.origin and self.next.origin == other.next.origin
        
    def __repr__(self):
        w =f"HalfEdge(origin={self.origin}, "
        w += (f"next={self.next.origin}" if self.next else "next=None") + ", "
        w += (f"twin={self.twin.origin, self.twin.next.origin}" if self.twin else "twin=None") + ", "
        w += f"face_centroid={self.face.centroid})"
        return w

@dataclass
class Face:
    half_edge: Optional[HalfEdge] = None # one of his half_edges
    parent: Optional["Face"] = field(default=None, repr=False)
    children: list["Face"] = field(default_factory=list, repr=False)
    level: int = 0

    def __repr__(self):
        w =f"""Face(points={self.half_edge.origin, self.half_edge.next.origin, self.half_edge.next.next.origin},
        Edges= {self.half_edge, self.half_edge.next, self.half_edge.next.next}
        has_parent={self.parent!=None},
        children={len(self.children)},
        level={self.level})"""
        return w

    @property
    def longest_half_edge(self):
        return max(self.edge_list, key=lambda he: he.length) # Hypotenuse edge

    @cached_property
    def centroid(self) -> Vertex:
        return Vertex(sum(he.origin.x for he in self.edge_list)/3, sum(he.origin.y for he in self.edge_list)/3)

    @property
    def edge_list(self):
        assert self.half_edge != None
        assert self.half_edge.next != None
        assert self.half_edge.next.next != None
        return [self.half_edge, self.half_edge.next, self.half_edge.next.next]

    @classmethod
    def make_Face(cls, v0:Vertex, v1:Vertex, v2:Vertex, level: int = 0):
        """
            Returns a face composed with parameter induced half edges.
            Remark: face.halfedge is v0 to v1
        """
        he = HalfEdge(v0, next=HalfEdge(v1, next=HalfEdge(v2)))
        he.next.next.next=he # close the triangle
        f = Face(he, level=level)
        he.face = he.next.face = he.next.next.face = f
        return f

    def refine(self, common_edge:Optional[HalfEdge]=None):
        """
            Refines face by hypotenuse bisect.
            If half-edge is given, it means we are looking for refining near a specific
            half-edge.
            Return (face, face), the two new faces created
        """
        hyp: HalfEdge = self.longest_half_edge

        # create mid point and new faces
        midpoint: Vertex = Vertex((hyp.origin.x+hyp.next.origin.x)/2,(hyp.origin.y+hyp.next.origin.y)/2)

        A: Vertex = hyp.next.next.origin # opposite of hypotenuse
        B: Vertex = hyp.next.origin # "Left" (not necessarily)
        C: Vertex = hyp.origin      # "Right"

        f0: Face = Face.make_Face(C, midpoint, A, level=self.level+1) # v0, v2 for clockwise
        f1: Face = Face.make_Face(midpoint, B, A, level=self.level+1) # Also, f0,f1.halfedge gives the newly created half edge (hypotenuse split)

        # Pair common edge
        f0.half_edge.next.twin=f1.half_edge.next.next
        f1.half_edge.next.next.twin=f0.half_edge.next

        # Pair new hypotenuse to old neighbor
        f0.half_edge.next.next.twin=hyp.next.next.twin
        if hyp.next.next.twin: hyp.next.next.twin.twin=f0.half_edge.next.next # Important: use .twin.twin ... (serious debugging to se this)

        # Pair new hypotenuse to old neighbor
        f1.half_edge.next.twin=hyp.next.twin
        if hyp.next.twin: hyp.next.twin.twin=f1.half_edge.next # Important: use .twin.twin ... (serious debugging to se this)

        self.children.extend([f0,f1]) # Congratulations
        f0.parent = f1.parent = self

        if not (common_edge and hyp == common_edge.twin) and hyp.twin is not None:
            neighbor = hyp.twin.face
            while neighbor.level <= self.level:
            # equal (or higher) levels garantees connectivity (non dangling vertices)
                neighbor0, neighbor1 = neighbor.refine(hyp)
                if hyp.twin in neighbor0.edge_list:
                    neighbor = neighbor0
                else:
                    neighbor = neighbor1
            neighbor0.half_edge.twin = f1.half_edge
            f1.half_edge.twin = neighbor0.half_edge
            neighbor1.half_edge.twin = f0.half_edge
            f0.half_edge.twin = neighbor1.half_edge
            # else: hypotenuse.twin.face is border so we dgaf
        return f0, f1


@dataclass(slots=True)
class Mesh:
    nx: int
    ny: int
    dx: float = field(init=False)
    dy: float = field(init=False)
    bounds: tuple[tuple[float,float], tuple[float,float]]
    vertices: list[Vertex] = field(default_factory=list)
    active_faces: list[Face] = field(default_factory=list)
    # half_edges: list[HalfEdge] = field(default_factory=list)
    faces: list[Face] = field(default_factory=list)

    def __post_init__(self):
        self._generate_mesh()

    def _generate_mesh(self):

        print("Mesh post init")

        self.dx = (self.bounds[0][1]-self.bounds[0][0])/(self.nx-1)
        self.dy = (self.bounds[1][1]-self.bounds[1][0])/(self.ny-1)

        # Creating vertices
        Y = np.linspace(*self.bounds[1], self.ny)
        X = np.linspace(*self.bounds[0], self.nx)
        for y in Y:
            for x in X:
                self.vertices.append(Vertex(x, y))
        print(f"{len(self.vertices)} vertices created.")

        edges={}

        def pair(he):
            key = ((he.origin.x,he.origin.y), (he.next.origin.x, he.next.origin.y))
            twin_key = ((he.next.origin.x, he.next.origin.y), (he.origin.x,he.origin.y))
            if twin_key in edges:
                twin_he = edges[twin_key]
                twin_he.twin = he
                he.twin = twin_he
            else:
                edges[key] = he


        # Creating half-edges and faces:
        for iy, y in enumerate(Y[:-1]):
            for ix, x in enumerate(X[:-1]):

                _id = ix + self.nx*iy

                bottom_left = self.vertices[_id] # (0,0)
                bottom_right = self.vertices[_id+1] # (1,0)
                top_left = self.vertices[_id+self.nx] # (0,1)
                top_right = self.vertices[_id+1+self.nx] # (1,1)

                f0 = Face.make_Face(top_left, bottom_right, bottom_left)
                f1 = Face.make_Face(bottom_right, top_left, top_right)

                f0.half_edge.twin= f1.half_edge
                f1.half_edge.twin= f0.half_edge

                for face in [f0,f1]:
                    for he in face.edge_list:
                        pair(he)
                self.faces.extend([f0,f1])

        self.active_faces = self.faces
        # print(f"{len(self.half_edges)} half edges created.")
        print(f"{len(self.faces)} faces created.")
        print(f"{self.nx} by {self.ny} mesh built")
    
    def refine(self, face):
        assert face in self.active_faces, f"Face: {face} not present in mesh's faces"
        face.refine()
        self.get_active_faces()

    def get_active_faces(self):
        self.active_faces = []
        def find_leaves(face):
            if len(face.children):
                yield from find_leaves(face.children[0])
                yield from find_leaves(face.children[1])
            else: yield face

        for face in self.faces:
            self.active_faces.extend(list(find_leaves(face)))