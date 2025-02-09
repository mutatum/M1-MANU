import numpy as np
import math
from typing import Final
from dataclasses import field, dataclass

""" 2D """

@dataclass(slots=True)
class Point:
    id: int
    x: float
    y: float

@dataclass(slots=True)
class Face: # A line essentially
    id: int
    vertices: tuple[int, int] # Ids
    length: float # same here, length could be calculated based on triangle's order

    # def __post_init__(self):
    #     self.center = (self.vertices[1]+self.vertices[0])/2

@dataclass(slots=True)
class Triangle: # Cell, really. Triangles will be isosceles
    id: int
    faces: tuple[int, int, int] # Hypothenuse is last ? anyway, computable by edge length
    centroid: Point
    area: float # area is really = base_area / 2**order
    order: int = 0 # Number of times Shrunk, this is for connecting properly with neighbors
               # For instance, an adjacent side of order n must either be for its neighbor
               # (1) an adjacent side and the neighbor be a triangle of same order
               # (2) the neighbor's hypothenuse and neighbor be a triangle of order-1
               # Rephrased:
               # an edge that is a hypothenuse on either side implies that the sides share the same order
               # an edge that is adjacent on either side implies that the neighbors share the same order
               # if the edge is adjacent on one side and hypothenuse on the other, the hypothenuse side must be 1 order higher

class Mesh:

    def __init__(self, shape: tuple[int, int], bounds: tuple[tuple[int,int], ...]):
        self.shape = shape
        self.bounds = bounds # No incorrect input checking for now | Rectangular domain
        self.dx = (bounds[0][1]-bounds[0][0])/(shape[0]-1)
        self.dy = (bounds[1][1]-bounds[1][0])/(shape[1]-1)
        self.points : dict[int: Point] = {}
        self.centroids : dict[int: Point] = {}
        self.faces : dict[int: Face] = {}
        self.cells : dict[int: Triangle]={} # id : Cell
        self.connectivity: dict[int: tuple[int,int,int]]={} # Cell id to Neighbor Cells id
        self.boundary_cells = list[Triangle]

        point_id = 0
        face_id = 0
        cell_id = 0

        prev_column_points = []
        prev_column_faces = []
        for ix in range(self.shape[0]-1):
            for iy in range(self.shape[1]-1):

                if ix==0:
                    pUL = Point(point_id,self.bounds[0][0]+self.dx*ix, self.bounds[1][0]+self.dy*(iy+1))
                    self.points[point_id] = pUL
                    point_id+=1
                    if iy == 0:
                        pDL = Point(point_id, self.bounds[0][0]+self.dx*ix, self.bounds[1][0]+self.dy*iy)
                        self.points[point_id] = pDL
                        point_id+=1
                        pDR = Point(point_id, self.bounds[0][0]+self.dx*(ix+1), self.bounds[1][0]+self.dy*iy)
                        self.points[point_id] = pDR
                        point_id+=1
                        prev_column_points.append(pDR)
                        bottom = Face(face_id, [pDL.id, pDR.id], length=self.dx)
                        self.faces[face_id] = bottom
                        face_id+=1
                        left = Face(face_id, [pDL.id, pUL.id], length=self.dy)
                        self.faces[face_id] = left
                        face_id+=1
                    else:
                        left = Face(face_id, [pDL.id, pUL.id], length=self.dy)
                        self.faces[face_id] = left
                        face_id+=1
                        pDL = pUL
                        pDR = pUR
                        bottom = top
                else:
                    if iy==0:
                        pDR = Point(point_id, self.bounds[0][0]+self.dx*(ix+1), self.bounds[1][0]+self.dy*iy)
                        self.points[point_id] = pDR
                        point_id+=1
                        prev_column_points.append(pDR)
                        pDL = prev_column_points.pop(0)
                        bottom = Face(face_id, [pDL.id, pDR.id], length=self.dx)
                        self.faces[face_id] = bottom
                        face_id+=1
                    else:
                        pDR = pUR
                        bottom = top
                    pUL = prev_column_points.pop(0)
                    left = prev_column_faces.pop(0)

                pUR = Point(point_id, self.bounds[0][0]+self.dx*(ix+1), self.bounds[1][0]+self.dy*(iy+1))
                self.points[point_id] = pUR
                point_id+=1
                prev_column_points.append(pUR) # appending points to Y column

                right = Face(face_id, [pDR.id, pUR.id], length=self.dy)
                self.faces[face_id] = right
                prev_column_faces.append(right)
                face_id +=1

                middle = Face(face_id, [pUL.id, pDR.id], length=np.sqrt(self.dx**2 + self.dy**2))
                self.faces[face_id] = middle
                face_id += 1

                top = Face(face_id, [pUL.id, pUR.id], length=self.dx)
                self.faces[face_id] = top
                face_id += 1

                sqrt2 = math.sqrt(2)
                left_centroid = Point(point_id, pDL.x + sqrt2 * self.dx/2, pDL.y + sqrt2 * self.dy/2)
                self.centroids[point_id] = left_centroid
                point_id+=1

                self.cells[cell_id] = Triangle(cell_id, [bottom.id, left.id, middle.id], left_centroid, self.dx*self.dy/2)
                cell_id += 1
                right_centroid= Point(point_id, pUR.x - sqrt2 * self.dx/2, pUR.y - sqrt2 * self.dy/2)
                self.centroids[point_id] = right_centroid
                point_id+=1
                self.cells[cell_id] = Triangle(cell_id, [top.id, right.id, middle.id], right_centroid, self.dx*self.dy/2)
                cell_id += 1

                
                # if ix == 0 or iy == 0 or ix==self.shape[0]-1 or iy==self.shape[1]-1:
                #     self.boundary_cells.append(id)
                # self.cells[id: Triangle(id, )]
