#include <iostream>
#include <fstream>
#include <cassert>
#include <random>
#include <ranges>
#include <vector>
#include <chrono>
#include <unistd.h>
#include <omp.h>
#include <algorithm>
#include <set>

using namespace std;

/**
 * Triple is the data structure for triples
 * defined == and != for eqaulity/inequality between triples
 * used for data like (vertex1 : size_t, vertex2 : size_t, weight : int)
 */
struct triple {
    size_t fst;
    size_t snd;
    int data;
    bool operator==(const triple& other) const {
        return fst == other.fst &&
               snd == other.snd &&
               data == other.data;
    }
    bool operator!=(const triple& other) const {
        return !(*this == other);
    }
};

/**
 * Lexicographical compare function for triples - determining inequality using
 * elements from left to right
 */
bool cmp(const triple& a, const triple& b) {
    if (a.fst != b.fst) {
        return a.fst < b.fst;
    }
    if (a.snd != b.snd) {
        return a.snd < b.snd;
    }
    return a.data < b.data;
}

/**
 * maps a vector with entries of type T in place
 */
template <typename T, typename F>
void map_inplace(vector<T> &data, const F &f) {
    #pragma omp parallel for
    for (size_t i = 0; i < data.size(); i++) {
        data[i] = f(data[i]);
    }
}

/**
 * input: integer vectors data and sums
 * computes prefix sum of data and stores in sums
 */
void prefix_sum(const vector<int>& data, vector<int>& sums) {
    size_t n = data.size();

    sums[0] = 0;
    int prefix = 0;
    #pragma omp parallel for reduction(inscan, +:prefix)
    for (size_t i = 0; i < n; i++) {
        prefix += data[i];
        #pragma omp scan inclusive(prefix)
        sums[i + 1] = prefix;
    }
}

/**
 * input: data:vector of type T, flags: int vector, psum : int vector, out: vector of type T
 * prefix sum of flags (stored in psum), is index of next available location in out to store a value
 * e.g. flags = [0,1,1,0,1] has psum = [0,0,1,2,2] so store entry 1 at 0, entry 2 at 3 and entry 4 at 2
 * returns length of output vector (number of elements left after filtering)
 */
template <typename T>
size_t generalizedFilter(const vector<T> &data, const vector<int> &flags, vector<int> &psum, vector<T> &out) {
    size_t N = data.size();
    prefix_sum(flags, psum);
    #pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
        if (flags[i])
            out[psum[i]] = data[i];
    }
    return psum[N];
}

/**
 * input: data:vector of type T, flags: int vector, psum : int vector
 * returns out = vector of type T, which is the result of filtering
 */
template <typename T>
vector<T> generalizedFilter(const vector<T> &data, const vector<int> &flags) {
    size_t N = data.size();
    vector<T> out(N);
    vector<int> psum(N + 1);
    size_t newN = generalizedFilter(data, flags, psum, out);
    out.resize(newN);
    return out;
}

/**
 * input: data:vector of type T, flags: int vector, psum : int vector, out : vector of type T, pred : a predicate function
 * Filters the data vector using pred and stores result in out
 * returns the length of out (the filtered vector)
 */
template <typename T, typename tau>
size_t filter(const vector<T> &data, vector<int> &flags, vector<int> &psum, vector<T> &out, tau &pred) {
    size_t total = 0;
    size_t n = data.size();
    #pragma omp parallel for reduction(+:total)
    for (size_t i = 0; i < n; i++) {
        if (pred(data[i])) {
            flags[i] = 1;
            total++;
        } else {
            flags[i] = 0;
        }
    }
    prefix_sum(flags, psum);
    #pragma omp parallel for
    for (size_t i = 0; i < n; i++) {
        if (flags[i])
            out[psum[i]] = data[i];
    }
    return total;
}

/**
 * input: data:vector of type T, flags: int vector, psum : int vector, out : vector of type T, pred : a predicate function
 * Filters the data vector using pred
 * returns filtered data in a vector of type T
 */
template <typename T, typename tau>
vector<T> filter(const vector<T> &data, tau &pred) {
    size_t N = data.size();
    vector<int> flags(N);
    vector<int> psum(N + 1);
    vector<T> out(N);
    size_t newN = filter(data, flags, psum, out, pred);
    out.resize(newN);
    return out;
}

/**
 * inputs: flips:bool vector, prob:double between [0,1]
 * stores result of coin flips with probability prob in flips
 */
void flipCoins(vector<bool> &flips, double prob) {
    #pragma omp parallel
    {
        // thread local rng
        random_device rd;
        mt19937 gen(rd() + omp_get_thread_num());
        bernoulli_distribution coin(prob);

        // flip a coin for each entry of flips
        #pragma omp for
        for (size_t i = 0; i < flips.size(); i++) {
            flips[i] = coin(gen);
        }
    }
}

/**
 * input: mapping:size_t vector
 * the mapping vector shows where each vertex should be mapped
 * For example mapping[i] = j means that vertex i should be mapped to vertex j
 * InitializeMapping maps every vertex to itself - representing that no contractions have occurred yet
 */
void initializeMapping(vector<size_t> &mapping) {
    #pragma omp parallel for
    for (size_t i = 0; i < mapping.size(); i++) {
        mapping[i] = i;
    }
}

/**
 * inputs: graph:vector of type triple (edge list), flips:vector of bools (vector of coin flips for each vertex), 
 * mapping:vector of size_t storing where each vector should contract to
 * If a vertex is heads, it is a "center" of a star, and if tails, it is a potential satellite vertex
 * updates mapping vector based on the coin flips for each vertex
 */
void starPartition(const vector<triple> &graph, const vector<bool> &flips, vector<size_t> &mapping) {
    // assume mapping is initialized to {0, 1, ..., n - 1}
    // i.e. each vertex initially maps to itself
    #pragma omp parallel for
    for (size_t e = 0; e < graph.size(); e++) {
        // loop through the edges of the graph
        // int e = myvector[e];
        size_t u = graph[e].fst;
        size_t v = graph[e].snd;
        if (!flips[u] && flips[v]) {
            // u |-> v if u flipped tails, v flipped heads, and they're adjacent
            // ties are broken arbitrarily
            mapping[u] = v;
        }
    }
}

/**
 * input: mapping:vector of size_t
 * output: size_t value of the number of centers in the mapping
 */
size_t countCenters(const vector<size_t> &mapping) {
    size_t total = 0;
    #pragma omp parallel for reduction(+:total)
    for (size_t i = 0; i < mapping.size(); i++) {
        if (mapping[i] == i) {
            total++;
        }
    }
    return total;
}

/**
 * inputs: mapping:vector of size_t, max_label:size_t
 * max_label stores the largest label number (or equivalently the number of unique 
 * label numbers since the label numbers are consecutive) after applying the mapping
 * outputs: size_t vector containing the update vertex labels after using the mapping
 */
vector<size_t> changeLabels(const vector<size_t> &mapping, size_t &max_label) {
    vector<size_t> res(mapping.size());
    size_t count = 0;
    for (size_t i = 0; i < mapping.size(); i++) {
        if (mapping[i] == i) {
            res[i] = count;
            count++;
        }
    }
    max_label = count;
    for (size_t i = 0; i < mapping.size(); i++) {
        res[i] = res[mapping[i]];
    }
    return res;
}

/**
 * inputs: graph:vector of triples (edge list representation)
 * If there are duplicate edges e.g. (1,2,w), (1,2,w) then only keep one of those edges
 * returns: graph without duplicates
 */
vector<triple> removeDuplicates(vector<triple> &graph) {
    size_t N = graph.size();
    // auto cmp = [&](const triple &e1, const triple &e2) {
    //     if (e1.fst < e2.fst) {
    //         return true;
    //     }
    //     if (e1.fst > e2.fst) {
    //         return false;
    //     }
    //     return e1.snd < e2.snd;
    // };
    sort(graph.begin(), graph.end(), cmp);
    // now the edges are sorted and we wish to remove duplicates and self-edges
    vector<int> flags(N);
    if (graph[0].fst != graph[0].snd) {
        flags[0] = 1;
    }
    #pragma omp parallel for
    for (size_t i = 1; i < N; i++) {
        triple e = graph[i];
        if (e.fst != e.snd && !(e.fst == graph[i - 1].fst && e.snd == graph[i - 1].snd)) {
            flags[i] = 1;
        }
    }
    return generalizedFilter(graph, flags);
}

/**
 * inputs: graph:vector of triple, mapping:vector of size_t, new_num_verts:size_t
 * new_num_verts is updated to how many vertices there are after applying the mapping to the graph
 * output: vector of triple of the new edge list representing the new contracted graph
 */
vector<triple> quotient(vector<triple> &graph, const vector<size_t> &mapping, size_t &new_num_verts) {
    vector<size_t> newLabels = changeLabels(mapping, new_num_verts);
    auto f = [&](const triple &e) -> triple {
        return triple { newLabels[e.fst], newLabels[e.snd], e.data };
    };
    map_inplace(graph, f);
    return removeDuplicates(graph);
}

/**
 * inputs: graph:vector of triple, verts:size_t
 * graph is an edge list, and verts is the number of vertices
 * output: size_t value of how many connected components there are
 */
size_t countComponents(vector<triple> &graph, size_t verts) {
    // cout << "\nreached loop again with " << verts << "vertices \n";
    // for (size_t i = 1; i < graph.size(); i++) {
    //     cout << "edge" << graph[i].first << "," << graph[i].second << "\n";
    // }
    if (graph.size() == 0) {
        return verts;
    }
    vector<size_t> mapping(verts);
    vector<bool> flips(verts);
    initializeMapping(mapping);
    flipCoins(flips, 0.5);
    starPartition(graph, flips, mapping);
    size_t new_num_verts;
    vector<triple> new_graph = quotient(graph, mapping, new_num_verts);
    return countComponents(new_graph, new_num_verts);
}

// (vertex1, vertex2, weight)
/**
 * inputs: graph:vector of triple, num_verts:size_t
 * Find all the vertex bridges (lightest edge out of each vertex)
 * result is a vector of triple where result[v] should be lightest edge out of vertex v
 * outputs: result:vector of triple representing edge list of vertex bridges
 */
vector<triple> vertexBridges (const vector<triple> &graph, size_t &num_verts){
    vector<triple> result(num_verts);

    # pragma omp parallel for
    for (size_t i = 0; i < num_verts; i++){
        result[i] = {i, i, numeric_limits<int>::max()};
    }

    #pragma omp parallel for
    for (size_t v = 0; v < num_verts; v++){
        size_t min_vx = v;
        int min_weight = numeric_limits<int>::max();
        for (auto &e : graph){
            if (e.fst == v){
                if (min_weight > e.data){
                    min_weight = e.data;
                    min_vx = e.snd;
                }
            }
            result[v].snd = min_vx;
            result[v].data = min_weight;
        }
    }
    // ! important !
    // result[u] = (u, v, w) means
    // vertex u's min weight edge is (u, v) with weight w
    return result;
}

/**
 * inputs: bridges:vector of triple, flips:vector of bool, mapping:size_t vector
 * Updates mapping by going though the lightest edge out of each vertex and 
 * using coin flips to determine if the vertex is a satellite and if it should 
 * be contracted to a neighbour (which must be a center)
 */
void bridgeStarPartition(const vector<triple> &bridges, const vector<bool> &flips, vector<size_t> &mapping) {
    // assume mapping is initialized to {0, 1, ..., n - 1}
    // i.e. each vertex initially maps to itself
    // note that this function only cares about the bridges, not the entire graph
    #pragma omp parallel for
    for (size_t e = 0; e < bridges.size(); e++) {
        // bridges[e] = (u, v, weight)
        const size_t u = bridges[e].fst;
        const size_t v = bridges[e].snd;
        if (!flips[u] && flips[v]) {
            mapping[u] = v;
        }
        if (!flips[v] && flips[u]) {
            mapping[v] = u;
        }
    }

    // size_t total = 0;
    // #pragma omp parallel for reduction(+:total)
    // for (size_t e = 0; e < bridges.size(); e++) {
    //     const auto& [u, v, weight] = bridges[e];
    //     if ((mapping[u] == v || mapping[v] == u) && (u < v || (v < u && bridges[v].snd != u))) {
    //         total += weight;
    //     }
    // }
    // return total;
}

/**
 * inputs: bridges:vector of triple, mapping:size_t vector
 * output: size_t value of the sum of weights of the contracted edges
 */
size_t sumContracted(const vector<triple> &bridges, const vector<size_t> &mapping) {
    size_t total = 0;
    #pragma omp parallel for reduction(+:total)
    for (size_t u = 0; u < mapping.size(); u++) {
        size_t v = mapping[u];
        if (u != v) {
            if (bridges[u].snd == v) {
                total += bridges[u].data;
            } else if (bridges[v].snd == u) {
                total += bridges[v].data;
            }
        }
    }
    return total;
}

/**
 * Boruvka's algorithm
 * inputs: graph:vector of triple, n_verts:triple
 * output: size_t value representing the weight of the MST
 */
size_t boruvka(vector<triple> &graph, size_t &n_verts) {
    if (graph.size() == 0) {
        return 0;
    }
    vector<size_t> mapping(n_verts);
    initializeMapping(mapping);
    vector<bool> flips(n_verts);
    flipCoins(flips, 0.5);
    vector<triple> bridges = vertexBridges(graph, n_verts);
    bridgeStarPartition(bridges, flips, mapping);
    size_t contractedWeights = sumContracted(bridges, mapping);
    size_t new_n_verts;
    vector<triple> newGraph = quotient(graph, mapping, new_n_verts);
    return contractedWeights + boruvka(newGraph, new_n_verts);
}

/* (vertex1, (vertex2, weight)) */
/**
 * inputs: input_filename:string, n:size_t, m_size_t
 * reads a graph from a text file where it is formatted as
 * n_vertices n_edges
 * vertex1 vertex1' weight1
 * vertex2 vertex2' weight2
 * ...
 * Note that graphs should be undirected so if edge (v1,v2) present so should (v2,v1)
 * outputs: vector of triple representing an edge list corresponding to input text file
 */
vector<triple> readGraph(const string &input_filename, size_t &n, size_t &m) {
    ifstream fin(input_filename);
    size_t n_verts, n_edges;
    fin >> n_verts >> n_edges;
    vector<triple> graph(2 * n_edges);
    for (auto &edge : graph) {
        fin >> edge.fst >> edge.snd >> edge.data;
    }
    n = n_verts;
    m = n_edges;
    return graph;
}

/**
 * after running make, test by using
 * ./graph -f <filename> -n <number threads>
 */
int main(int argc, char *argv[]) {
    string input_filename;
    int num_threads = 0;
    int opt;
    while ((opt = getopt(argc, argv, "f:n:")) != -1) {
        switch (opt) {
            case 'f':
                input_filename = optarg;
                break;
            case 'n':
                num_threads = atoi(optarg);
                break;
            default:
                cerr << "blah blah blah usage\n";
                exit(EXIT_FAILURE);
        }
    }
    if (empty(input_filename) || num_threads <= 0) {
        cerr << "blah blah blah usage\n";
        exit(EXIT_FAILURE);
    }

    cout << "Number of threads: " << num_threads << '\n';
    cout << "Input file: " << input_filename << '\n';

    // ifstream fin(input_filename);
    size_t n_verts, n_edges;
    // fin >> n_verts >> n_edges;
    // vector<pair<int, int>> graph(2 * n_edges);
    // for (auto &edge : graph) {
    //     fin >> edge.first >> edge.second;
    // }

    vector<triple> graph = readGraph(input_filename, n_verts, n_edges);

    omp_set_num_threads(num_threads);

    // vector<size_t> mapping(n_verts);
    // initializeMapping(mapping);
    // vector<bool> flips(n_verts);
    // flipCoins(flips, 0.5);
    // cout << "Coin flips: ";
    // for (bool b : flips) {
    //     cout << b << " ";
    // }

    const auto compute_start = chrono::steady_clock::now();
    size_t mst_weight = boruvka(graph, n_verts);
    cout << "weight of mst: " << mst_weight << "\n";

    // vector<triple> bridges = vertexBridges(graph, n_verts);
    // cout << "\nbridges: ";
    // for (triple &e : bridges) {
    //     cout << "(" << e.fst << ", " << e.snd << ", " << e.data << ") ";
    // }
    // bridgeStarPartition(bridges, flips, mapping);
    // size_t contractedWeights = sumContracted(bridges, mapping);
    // cout << "\nmapping: ";
    // for (size_t i = 0; i < mapping.size(); i++) {
    //     cout << "(" << i << ", " << mapping[i] << ") ";
    // }
    // size_t new_n_verts;
    // vector<triple> newGraph = quotient(graph, mapping, new_n_verts);
    // cout << "\nnew graph: ";
    // for (triple &e : newGraph) {
    //     cout << "(" << e.fst << ", " << e.snd << ", " << e.data << ") ";
    // }
    // cout << "\nsum of contracted weights: " << contractedWeights << "\n";

    // int count = countComponents(graph, n_verts);
    // cout << "\nnumber of components: " << count << "\n";

    // initializeMapping(mapping);
    // flipCoins(flips, 0.5);
    // starPartition(graph, flips, mapping);
    // cout << "printing coin flips: ";
    // for (bool b : flips) {
    //     cout << b << " ";
    // }
    // cout << "\nprinting mapping: ";
    // for (int i : mapping) {
    //     cout << i << " ";
    // }
    // int new_num_verts;
    // vector<int> res = changeLabels(mapping, new_num_verts);
    // cout << "\nprinting changeLabels(mapping): ";
    // for (int i : res) {
    //     cout << i << " ";
    // }
    // auto new_graph = quotient(graph, mapping, new_num_verts);
    // cout << "\nnew number of vertices: " << new_num_verts << "\n";
    // cout << "\nprinting quotient graph: ";
    // for (auto &e : new_graph) {
    //     cout << "(" << e.first << ", " << e.second << ") ";
    // }
    const double compute_time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - compute_start).count();
    cout << "Computation time (sec): " << compute_time << "\n";
}