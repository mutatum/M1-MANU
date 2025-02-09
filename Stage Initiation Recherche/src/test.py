# %%
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import mesh

m = mesh.Mesh((6,6), ((0,1),(0,1)))
print(m.points)
print(m.points.keys())
print(m.points.values())
plt.figure(figsize=(6,6))
plt.scatter([p.x for p in m.points.values()],[p.y for p in m.points.values()], color="r")
plt.scatter([p.x for p in m.centroids.values()],[p.y for p in m.centroids.values()], color="y")
lc = LineCollection([[[m.points[id].x,m.points[id].y] for id in f.vertices] for f in m.faces.values()], colors='k', linewidth=2)
ax = plt.gca()
ax.add_collection(lc)
ax.autoscale()
print(m.faces.values())


# %%
