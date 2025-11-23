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

void initializeMapping(vector<int> &mapping) {
    #pragma omp parallel for
    for (size_t i = 0; i < mapping.size(); i++) {
        mapping[i] = i;
    }
}

void starPartition(const vector<pair<int, int>> &graph, const vector<bool> &flips, vector<int> &mapping) {
    // assume mapping is initialized to {0, 1, ..., n - 1}
    // i.e. each vertex initially maps to itself
    #pragma omp parallel for
    for (size_t e = 0; e < graph.size(); e++) {
        // loop through the edges of the graph
        // int e = myvector[e];
        int u = graph[e].first;
        int v = graph[e].second;
        if (!flips[u] && flips[v]) {
            // u |-> v if u flipped tails, v flipped heads, and they're adjacent
            // ties are broken arbitrarily
            mapping[u] = v;
        }
    }
}

int countCenters(const vector<int> &mapping) {
    int total = 0;
    #pragma omp parallel for reduction(+:total)
    for (size_t i = 0; i < mapping.size(); i++) {
        if (mapping[i] == (int)i) {
            total++;
        }
    }
    return total;
}

vector<int> changeLabels(const vector<int> &mapping, int &max_label) {
    vector<int> res(mapping.size());
    int count = 0;
    for (size_t i = 0; i < mapping.size(); i++) {
        if (mapping[i] == (int)i) {
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

vector<pair<int, int>> removeDuplicates(vector<pair<int, int>> &graph) {
    size_t N = graph.size();
    auto cmp = [&](const pair<int, int> &e1, const pair<int, int> &e2) {
        if (e1.first < e2.first) {
            return true;
        }
        if (e1.first > e2.first) {
            return false;
        }
        return e1.second < e2.second;
    };
    sort(graph.begin(), graph.end(), cmp);
    // now the edges are sorted and we wish to remove duplicates and self-edges
    vector<int> flags(N);
    if (graph[0].first != graph[0].second) {
        flags[0] = 1;
    }
    #pragma omp parallel for
    for (size_t i = 1; i < N; i++) {
        pair<int, int> e = graph[i];
        if (e.first != e.second && e != graph[i - 1]) {
            flags[i] = 1;
        }
    }
    return generalizedFilter(graph, flags);
}

vector<pair<int, int>> quotient(vector<pair<int, int>> &graph, const vector<int> &mapping, int &new_num_verts) {
    vector<int> newLabels = changeLabels(mapping, new_num_verts);
    auto f = [&](const pair<int, int> &e) {
        return pair<int, int>(newLabels[e.first], newLabels[e.second]);
    };
    map_inplace(graph, f);
    return removeDuplicates(graph);
}

int countComponents(vector<pair<int, int>> &graph, int verts) {
    // cout << "\nreached loop again with " << verts << "vertices \n";
    // for (size_t i = 1; i < graph.size(); i++) {
    //     cout << "edge" << graph[i].first << "," << graph[i].second << "\n";
    // }
    if (graph.size() == 0) {
        return verts;
    }
    vector<int> mapping(verts);
    vector<bool> flips(verts);
    initializeMapping(mapping);
    flipCoins(flips, 0.5);
    starPartition(graph, flips, mapping);
    int new_num_verts;
    vector<pair<int, int>> new_graph = quotient(graph, mapping, new_num_verts);
    return countComponents(new_graph, new_num_verts);
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

    ifstream fin(input_filename);
    size_t n_verts, n_edges;
    fin >> n_verts >> n_edges;
    vector<pair<int, int>> graph(2 * n_edges);
    for (auto &edge : graph) {
        fin >> edge.first >> edge.second;
    }

    omp_set_num_threads(num_threads);

    vector<int> mapping(n_verts);
    vector<bool> flips(n_verts);

    const auto compute_start = chrono::steady_clock::now();

    int count = countComponents(graph, n_verts);
    cout << "\nnumber of components: " << count << "\n";

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