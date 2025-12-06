import random
import argparse

# n_verts = 500

# creates a complete graph with n_verts without weights
def complete_graph_no_weights(n_verts):
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

# creates a complete graph with n_verts with weights
def complete_graph_weights(n_verts):
    n_edges = n_verts * (n_verts - 1) // 2
    with open(f"weights_k{n_verts}.txt", "w") as f:
        f.write(f"{n_verts} {n_edges}\n")
        pairs = []
        for i in range(n_verts):
            for j in range(i+1,n_verts):
                if i != j:
                    w = random.randint(1, 100)
                    pairs.append((i, j, w))
                    pairs.append((j, i, w))
        random.shuffle(pairs)
        for i, j, w in pairs:
            f.write(f"{i} {j} {w}\n")

# creates a n_verts graph with probability of adding an edge being prob
def random_graph_with_weights(n_verts, prob):
    n_edges = 0
    pairs = []
    for i in range(n_verts):
        for j in range(i+1, n_verts):
            if random.random() < prob:
                w = random.randint(1, 100)
                pairs.append((i,j,w))
                pairs.append((j,i,w))
                n_edges += 1
    random.shuffle(pairs)
    with open(f"random_weights_{n_verts}.txt", "w") as f:
        f.write(f"{n_verts} {n_edges}\n")
        for i, j, w in pairs:
            f.write(f"{i} {j} {w}\n")

# creates a n_verts graph with probability of adding an edge being prob
# eforces that vertex i has edge with i+1 so that the graph is connected
def random_connected_graph_with_weights(n_verts, prob):
    n_edges = 0
    pairs = []
    for i in range(n_verts):
        for j in range(i+1, n_verts):
            if j == i+1:
                w = random.randint(1,100)
                pairs.append((i,j,w))
                pairs.append((j,i,w))
                n_edges += 1
            else:
                if random.random() < prob:
                    w = random.randint(1, 100)
                    pairs.append((i,j,w))
                    pairs.append((j,i,w))
                    n_edges += 1
    random.shuffle(pairs)
    with open(f"random{prob}_connected_weights_{n_verts}.txt", "w") as f:
        f.write(f"{n_verts} {n_edges}\n")
        for i, j, w in pairs:
            f.write(f"{i} {j} {w}\n")

# To run:
# python3 gen_test.py -v <number vertices> -p <probability of edge if needed, default = 0.5> -t <type of graph>
# type of graph is of following:
# complete_no_weights, complete_weights, random_with_weights, random_connected_with_weights
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', type=int, 
                        help="number of vertices")
    parser.add_argument('-p', type=float, 
                        default = 0.5,
                        help="probability of making an edge")
    parser.add_argument('-t', type=str, 
                        help="type of graph is one of the following: complete_no_weights, complete_weights, random_with_weights, random_connected_with_weights")
    args = parser.parse_args()
    n_verts = args.v
    prob = args.p

    if args.t == "complete_no_weights":
        complete_graph_no_weights(n_verts)
        print(f"generated graph stored with name: k{n_verts}.txt")
    elif args.t == "complete_weights":
        complete_graph_weights(n_verts)
        print(f"generated graph stored with name: weights_k{n_verts}.txt")
    elif args.t == "random_with_weights":
        random_graph_with_weights(n_verts, prob)
        print(f"generated graph stored with name: random_weights_{n_verts}.txt")
    elif args.t == "random_connected_with_weights":
        random_connected_graph_with_weights(n_verts, prob)
        print(f"generated graph stored with name: random{prob}_connected_weights_{n_verts}.txt")
    else:
        print("graph type must be one of: complete_no_weights, complete_weights, random_with_weights, random_connected_with_weights")