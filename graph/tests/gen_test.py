import random

n_verts = 100
n_edges = n_verts * (n_verts - 1) // 2

with open(f"k{n_verts}.txt", "w") as f:
    f.write(f"{n_verts} {n_edges}\n")
    pairs = []
    for i in range(n_verts):
        for j in range(n_verts):
            if i != j:
                pairs.append((i, j))
    random.shuffle(pairs)
    for i, j in pairs:
        f.write(f"{i} {j}\n")