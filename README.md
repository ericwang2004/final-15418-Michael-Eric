# 15-418 Proposal: Star Contraction

Eric Wang · Michael Zhou  
November 2025

## URL
https://github.com/ericwang2004/final-15418-Michael-Eric/

## Summary
The goal of our project is to implement **star contraction**, which can be seen as a template for many parallel graph algorithms, and apply it to **Borůvka’s minimum spanning tree (MST)** algorithm. The platform we use is **C++/OpenMP**; however, if time permits, we would also like to implement these in **CUDA**. We will then compare and contrast the implementations.

## Background
A fruitful source of problems amenable to parallelism is graph algorithms. Consider the problem:

CONN : {G | G is a graph} → ℕ  
G ↦ (# of connected components of G)

DFS/BFS are sequential. To parallelize CONN, we use graph contraction:

Let G = G(0) be the input graph, and let G₁, …, Gₙ be connected subgraphs of G(0) that partition V(G(0)). Contract each Gᵢ into a single vertex to obtain G(1). Recurse until reaching G(k) with no edges. Return |V(G(k))|.

Graph contraction is a divide‑and‑conquer method. Each thread contracts a subgraph, and results are combined. The key data structure enabling this is **array‑based sequences**, so we must implement the **sequence interface** (or at least the required subset).

## Challenges
### 1. Sufficiently Random Star Partition
The partition step must yield parallel speedup. In **star contraction**, each component is a star (a K₁,r). Each vertex v flips a coin:
- If heads → v becomes a **center**.
- If tails → v looks at neighbors:
  - If any flipped heads → v chooses one **randomly** and becomes a **satellite**.
  - If none → v becomes a center.

Randomness here is crucial. For example, in Kₙ, a naïve partition could produce one huge star and many singletons, destroying parallelism.

### 2. Sequence Filter
The **filter** operation (returning the subsequence of elements satisfying a predicate) should have span O(log n) for constant‑time predicates. Achieving close‑to‑ideal speedup is challenging.

### 3. Graph Representation
Choices include adjacency matrix, adjacency sequence, and edge sequence. The chosen representation affects:
- supported operations,
- span of algorithms,
- memory locality.

We may need to experiment to determine the best representation.

## Resources
Primary reference: the **15‑210 textbook**, especially chapters on sequence data structures, graph contraction, and MSTs.  
Link: https://www.cs.cmu.edu/~15210/docs/book.pdf

We will use **OpenMP**, and possibly **CUDA**.

## Goals and Deliverables
### Plan to Achieve
- **50% goal** — Implement the subset of the sequence interface needed for the project.
- **75% goal** — Implement star partition and graph connectivity with good OpenMP speedup.
- **100% goal** — Implement Borůvka’s MST.

### Hope to Achieve
- **125% goal** — CUDA implementation of sequences, star partition, and connectivity.
- **150% goal** — CUDA implementation of Borůvka’s MST.

### Demo
At the poster session, we plan to show speedup graphs for:
- sequence operations (likely **filter**),
- graph connectivity,
- Borůvka’s algorithm.

If CUDA implementations are completed, we will compare those as well.

## Platform Choice
OpenMP’s **shared‑address‑space model** is appropriate because threads will need to contend for parts of the graph. Parallelism naturally aligns with `parallel for` loops.

CUDA is also a reasonable option (especially for sequence operations), but using shared memory for graph operations may be tricky.

## Schedule
| Week of (starting Monday) | Implementation goals |
|---------------------------|-----------------------|
| **November 17** | Sequence operations: inject and filter |
| **November 24** | Star partition |
| **December 1** | Borůvka’s MST (C++/OpenMP) |
