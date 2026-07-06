import numpy as np

# Cylinder parameters
radius = 0.25      # radius in meters
height = 1.2      # full tank height in meters
n_points = 36     # number of points around circumference
x_center = 5.0    # middle of tank length
y_center = 2.0    # center of tank width

# Generate cylinder vertices
triangles = []
angles = np.linspace(0, 2*np.pi, n_points, endpoint=False)

# Bottom and top circle centers
bottom_center = np.array([x_center, y_center, 0])
top_center = np.array([x_center, y_center, height])

for i in range(n_points):
    a1 = angles[i]
    a2 = angles[(i+1) % n_points]
    
    p1 = np.array([x_center + radius*np.cos(a1), y_center + radius*np.sin(a1), 0])
    p2 = np.array([x_center + radius*np.cos(a2), y_center + radius*np.sin(a2), 0])
    p3 = np.array([x_center + radius*np.cos(a1), y_center + radius*np.sin(a1), height])
    p4 = np.array([x_center + radius*np.cos(a2), y_center + radius*np.sin(a2), height])
    
    triangles.append((p1, p2, p3))
    triangles.append((p2, p4, p3))
    triangles.append((bottom_center, p2, p1))
    triangles.append((top_center, p3, p4))

# Write STL file
with open('cylinder.stl', 'w') as f:
    f.write('solid cylinder\n')
    for tri in triangles:
        v1, v2, v3 = tri
        normal = np.cross(v2-v1, v3-v1)
        norm = np.linalg.norm(normal)
        if norm > 0:
            normal = normal/norm
        f.write(f'  facet normal {normal[0]:.6f} {normal[1]:.6f} {normal[2]:.6f}\n')
        f.write(f'    outer loop\n')
        f.write(f'      vertex {v1[0]:.6f} {v1[1]:.6f} {v1[2]:.6f}\n')
        f.write(f'      vertex {v2[0]:.6f} {v2[1]:.6f} {v2[2]:.6f}\n')
        f.write(f'      vertex {v3[0]:.6f} {v3[1]:.6f} {v3[2]:.6f}\n')
        f.write(f'    endloop\n')
        f.write(f'  endfacet\n')
    f.write('endsolid cylinder\n')

print("cylinder.stl created successfully!")
