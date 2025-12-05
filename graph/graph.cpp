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

bool cmp(const triple& a, const triple& b) {
    if (a.fst != b.fst) {
        return a.fst < b.fst;
    }
    if (a.snd != b.snd) {
        return a.snd < b.snd;
    }
    return a.data < b.data;
}

template <typename T, typename F>
void map_inplace(vector<T> &data, const F &f) {
    #pragma omp parallel for
    for (size_t i = 0; i < data.size(); i++) {
        data[i] = f(data[i]);
    }
}

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

template <typename T>
vector<T> generalizedFilter(const vector<T> &data, const vector<int> &flags) {
    size_t N = data.size();
    vector<T> out(N);
    vector<int> psum(N + 1);
    size_t newN = generalizedFilter(data, flags, psum, out);
    out.resize(newN);
    return out;
}

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

void initializeMapping(vector<size_t> &mapping) {
    #pragma omp parallel for
    for (size_t i = 0; i < mapping.size(); i++) {
        mapping[i] = i;
    }
}

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

vector<triple> quotient(vector<triple> &graph, const vector<size_t> &mapping, size_t &new_num_verts) {
    vector<size_t> newLabels = changeLabels(mapping, new_num_verts);
    auto f = [&](const triple &e) -> triple {
        return triple { newLabels[e.fst], newLabels[e.snd], e.data };
    };
    map_inplace(graph, f);
    return removeDuplicates(graph);
}

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

/* return the sum of the weights of the contracted edges */
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