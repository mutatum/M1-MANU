# %%
import matplotlib.pyplot as plt
import math
import mesh  # your mesh module with Mesh, Face, Vertex, HalfEdge

def draw_face(face,color):
    he = face.half_edge
    vertices = []
    start = he
    while True:
        vertices.append((he.origin.x, he.origin.y))
        he = he.next
        if he == start:
            break
    vertices.append(vertices[0])  # close the loop
    xs, ys = zip(*vertices)
    plt.plot(xs, ys, '-', color=color)  # draw triangle boundary in black

def draw_edge(he,color):
    vertices = []
    vertices.append((he.origin.x, he.origin.y))
    vertices.append((he.next.origin.x, he.next.origin.y))
    xs, ys = zip(*vertices)
    plt.plot(xs, ys, '-', color=color)  # draw triangle boundary in black

def plot_mesh(mesh: mesh.Mesh, arrows:bool = False):
    """
    Visualize the mesh:
      - Draw each triangle (face) with a dashed boundary.
      - Label the face at its centroid.
      - Draw half-edge arrows with an inward offset.
        Twin edges are colored oppositely: if a half-edge has a twin,
        it is 'red' if id(he) < id(he.twin) and 'blue' otherwise.
    """
    n = max(5, len(mesh.active_faces)/2-5)
    plt.figure(figsize=(n,n))
    
    for i, face in enumerate(mesh.active_faces):
        # Get face vertices by traversing the half-edge cycle.
        # if face.half_edge.origin.x==face.half_edge.next.origin.x or face.half_edge:
        #     print("FARTS", face)
        he = face.half_edge
        vertices = []
        start = he
        while True:
            vertices.append((he.origin.x, he.origin.y))
            he = he.next
            if he == start:
                break
        vertices.append(vertices[0])  # close the loop for plotting.
        xs, ys = zip(*vertices)
        plt.plot(xs, ys, 'k--', alpha=0.6)
        
        # Label the face at its centroid.
        
        # For each half-edge, draw an arrow.
        if arrows:
            color = ['red', 'blue', 'green'][face.level % 3]
            plt.text(face.centroid.x, face.centroid.y, f"{face.level}", 
                    color=color, fontsize=12, ha='center', va='center')
            for he in face.edge_list:
                # Compute the midpoint of the half-edge.
                x_start, y_start = he.origin.x, he.origin.y
                x_end, y_end = he.next.origin.x, he.next.origin.y
                a,b = .73, .27
                xm, ym = (x_start+ x_end)/2, (y_start + y_end)/2
                xm1, ym1 = (x_start*a+ x_end*b), (y_start*a + y_end*b)
                
                # Compute an inward offset: from midpoint towards face centroid.
                cx, cy = face.centroid.x, face.centroid.y
                nx, ny = cx - xm, cy - ym
                norm = math.hypot(nx, ny)
                if norm != 0:
                    nx, ny = nx / norm, ny / norm
                offset = 0.011  # adjust offset magnitude as needed.
                xm_offset, ym_offset = xm1 + nx * offset, ym1 + ny * offset
                
                # Arrow vector along half-edge (scaled for clarity).
                dx, dy = (x_end - x_start) * 0.45, (y_end - y_start) * 0.45
                
                plt.arrow(xm_offset, ym_offset, dx, dy, head_width=0.013, head_length=norm*.3,
                        fc=color, ec=color, length_includes_head=True)
    
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Mesh Visualization: Faces and Inward-Offset Half-Edges')
    plt.show()


m = mesh.Mesh(7,7, ((0,1), (0,1)))
# vx = [v.x for v in m.vertices]
# vy = [v.y for v in m.vertices]
# plt.scatter(vx,vy)
# plt.grid()
# for f in m.faces:
    # draw_face(f)

# plot_mesh(m,arrows=True)
m.refine(m.active_faces[9*9//2]);
# plot_mesh(m,arrows=True)
m.refine(m.active_faces[9*9//2]);
# plot_mesh(m,arrows=True)
for i in range(5):
    m.refine(m.active_faces[2]);
    plot_mesh(m,arrows=True)