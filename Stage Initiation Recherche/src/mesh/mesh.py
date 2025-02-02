import numpy as np
from typing import Final
from dataclasses import field, dataclass


@dataclass(slots=True)
class Point:
    x:float
    y:float
    z:float

class Face:
    def __init__(self, face_id, vertex_ids, cell_ids) -> None:
        self.face_id = face_id
        self.vertex_ids = vertex_ids
        self.cell_ids = cell_ids
        self.area = None
        self.normal = None

@dataclass(slots=True)
class Cell: # regular square/cube cell for now
    centroid: tuple[int, ...] # coordinates (x,y,z) or (x,y) or (x)
    vertex_id: list[int]
    face_ids: list[int]
    volume: float


@dataclass(slots=True)
class Mesh:
    nx: int
    ny: int
    nz: int
    vertices: dict[int, Point]
    cells: dict[int,Cell]={}

    def __post_init__(self) -> None:
        for i in range(nx):
            for j in range(self.ny)
                for 
