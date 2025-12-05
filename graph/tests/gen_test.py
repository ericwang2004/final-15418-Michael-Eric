import random

n_verts = 500



# n_edges = n_verts * (n_verts - 1) // 2
# with open(f"k{n_verts}.txt", "w") as f:
#     f.write(f"{n_verts} {n_edges}\n")
#     pairs = []
#     for i in range(n_verts):
#         for j in range(n_verts):
#             if i != j:
#                 pairs.append((i, j))
#     random.shuffle(pairs)
#     for i, j in pairs:
#         f.write(f"{i} {j}\n")

n_edges = 0
pairs = []
for i in range(n_verts):
    for j in range(i, n_verts):
        if i!= j and random.random() < 1:
            w = random.randint(1, 100)
            pairs.append((i,j,w))
            pairs.append((j,i,w))
            n_edges += 1
random.shuffle(pairs)
with open(f"random_weights_{n_verts}.txt", "w") as f:
    f.write(f"{n_verts} {n_edges}\n")
    for i, j, w in pairs:
        f.write(f"{i} {j} {w}\n")