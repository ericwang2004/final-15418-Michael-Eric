n_verts = 100
n_edges = n_verts * (n_verts - 1) // 2

with open("k{}.txt".format(n_verts), "w") as f:
    f.write(f"{n_verts} {n_edges}\n")
    for i in range(n_verts):
        for j in range(n_verts):
            if i != j:
                f.write(f"{i} {j}\n")