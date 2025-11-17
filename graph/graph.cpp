#include <iostream>
#include <fstream>
#include <cassert>
#include <random>
#include <ranges>
#include <vector>
#include <chrono>
#include <unistd.h>
#include <omp.h>

using namespace std;

bool is_even(int x) { 
    return x % 2 == 0; 
}

void prefix_sum(const std::vector<int>& data, std::vector<int>& sums) {
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

int filter_parallel(const vector<int> &data, vector<int> &flags, vector<int> &psum, vector<int> &out) {
    int total = 0;
    size_t n = data.size();
    #pragma omp parallel for reduction(+:total)
    for (size_t i = 0; i < n; i++) {
        if (data[i] % 2 == 0) {
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
        int u = graph[e].first;
        int v = graph[e].second;
        if (!flips[u] && flips[v]) {
            // u |-> v if u flipped tails, v flipped heads, and they're adjacent
            // ties are broken arbitrarily
            mapping[u] = v;
        }
    }
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
    initializeMapping(mapping);
    flipCoins(flips, 0.5);
    cout << "printing coin flips: ";
    for (bool b : flips) {
        cout << b << " ";
    }
    cout << "\n";
    starPartition(graph, flips, mapping);
    const double compute_time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - compute_start).count();
    cout << "Computation time (sec): " << compute_time << "\n";
    for (int x : mapping) {
        cout << x << " ";
    }
    cout << "\n";
}