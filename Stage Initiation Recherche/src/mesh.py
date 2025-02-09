import numpy as np
import math
from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional

""" 2D """

@dataclass(slots=True)
class Vertex:
    x: float
    y: float
    edge: Optional["HalfEdge"] = None

@dataclass(slots=True)
class HalfEdge:
    origin: Optional[Vertex] = None
    hypothenuse: bool
    twin: Optional["HalfEdge"] = None
    next: Optional["HalfEdge"] = None
    face: Optional["Face"] = None

    @cached_property
    def length(self):
        lx = (self.origin.x-self.next.origin.x)**2
        ly = (self.origin.y-self.next.origin.y)**2
        return math.sqrt(lx+ly)

@dataclass(slots=True)
class Face:
    half_edge: Optional[HalfEdge] = None # one of his half_edges
    parent: Optional["Face"] = None
    children: list["Face"] = field(default_factory=list)

    @cached_property
    def centroid(self) -> Vertex:
        return Vertex(sum([he.origin.x for he in self.edge_list])/3, sum([he.origin.y for he in [he1,he2,he3]])/3)

    @property
    def edge_list(self):
        assert self.half_edge != None
        assert self.half_edge.next != None
        assert self.half_edge.next.next != None
        return [self.half_edge, self.half_edge.next, self.half_edge.next.next]
    

@dataclass(slots=True)
class Mesh:
    nx: int
    ny: int
    dx: float = field(init=False)
    dy: float = field(init=False)
    bounds: tuple[tuple[float,float], tuple[float,float]]
    vertices: list[Vertex] = field(default_factory=list)
    half_edges: list[HalfEdge] = field(default_factory=list)
    faces: list[Face] = field(default_factory=list)

    def __post_init__(self):

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

        # Creating half-edges and faces:
        previous_row_top = []
        for iy, y in enumerate(Y[:-1]):
            previous_face = None
            for ix, x in enumerate(X[:-1]):

                id = ix + self.nx*iy

                v0 = self.vertices[id] # (0,0)
                v1 = self.vertices[id+1] # (0,1)
                v2 = self.vertices[id+self.nx] # (1,0)
                v3 = self.vertices[id+1+self.nx] # (1,1)

                he0 = HalfEdge(v2, next=HalfEdge(v1, next=HalfEdge(v0)))
                he0.next.next.next=he0
                he1 = HalfEdge(v1, next=HalfEdge(v2, next=HalfEdge(v3)))
                he1.next.next.next=he1

                he0.twin=he1
                he1.twin=he0
                self.half_edges.extend([he0, he0.next, he0.next.next,he1,he1.next,he1.next.next])

                f0 = Face(he0)
                f1 = Face(he1)

                he0.face = f0
                he1.face = f1

                self.faces.extend([f0,f1])
                previous_row_top.append(f1)

                if previous_face == None:
                    previous_face = f1
                else:
                    f0.half_edge.next.next.twin = previous_face.half_edge.next.next
                    previous_face.half_edge.next.next.twin = f0.half_edge.next.next

                if ix!=0:
                    down_face = previous_row_top.pop()
                    f0.half_edge.twin = down_face.half_edge.next
                    down_face.half_edge.next.twin = he0

        print(f"{len(self.half_edges)} half edges created.")
        print(f"{len(self.faces)} faces created.")
        print(f"{self.nx} by {self.ny} mesh built")
    
    # [href=https://www.math.uci.edu/~chenlong/ifemdoc/afem/bisectdoc.html] not really used in fact since triangles are rectangle
    def refine_face(self, face: Face):
        """Refines face by hypotenuse bisect"""
        hes = face.edge_list
        hypotenuse = min(hes, key=lambda he: he.length)
        face.children.append()
        