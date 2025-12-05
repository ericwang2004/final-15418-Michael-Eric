#include <iostream>
#include <fstream>
#include <cassert>
#include <random>
#include <ranges>
#include <vector>
#include <chrono>
#include <unistd.h>
#include <algorithm>
#include <set>
#include <queue>

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

class compare_weight{
    public:
    bool operator()(triple x, triple y){
        // return x.data < y.data; 
        // Since the queue is maxheap and we want smallest weight first,
        // we have to flip the results of comparison
        return x.data > y.data;
    }
};

vector<triple> readGraph(const string &input_filename, size_t &n, size_t &m, bool hasweights) {
    ifstream fin(input_filename);
    size_t n_verts, n_edges;
    fin >> n_verts >> n_edges;
    vector<triple> graph(2 * n_edges);
    if (hasweights){
        for (auto &edge : graph) {
            fin >> edge.fst >> edge.snd >> edge.data;
        }
    } else {
        for (auto &edge : graph) {
            fin >> edge.fst >> edge.snd;
            // default weight is 1
            edge.data = 1;
        }
    }
    n = n_verts;
    m = n_edges;
    return graph;
}

// using https://www.cs.cmu.edu/afs/cs/academic/class/15210-f12/www/recis/rec10.pdf
size_t prim_MST(vector<triple> &graph, size_t &n_vertices){
    vector<bool> seen(n_vertices, false);
    int graph_matrix[n_vertices][n_vertices];
    for (size_t i = 0; i < n_vertices; i++){
        for (size_t j = 0; j < n_vertices; j++){
            graph_matrix[i][j] = -1;
        }
    }

    for (auto &edge : graph) {
        graph_matrix[edge.fst][edge.snd] = edge.data;
    }
    size_t weight = 0;
    seen[0] = true;
    // the weight stored in tuple will be -weight since priority queue is max heap
    // a priority queue allows us to find the lightest edge from the seen vertices
    // to the unseen vertices the quickest
    priority_queue<triple, vector<triple>, compare_weight> frontier_edges;
    for (size_t i = 0; i < n_vertices; i++){
        if (graph_matrix[0][i] != -1){
            triple curr = {0,i,graph_matrix[0][i]};
            frontier_edges.push(curr);
        }
    }
    // print out priority queue to test
    // while (!frontier_edges.empty()){
    //     triple curr;
    //     curr = frontier_edges.top();
    //     cout << "(" << curr.fst << "," << curr.snd << "," << curr.data << "),";
    //     frontier_edges.pop();
    // }
    int vertices_left = n_vertices-1; //since one vertex has been added already
    while (vertices_left > 0){
        if (frontier_edges.empty()){
            cout << "the graph is not connected \n";
            exit(EXIT_FAILURE);
        }
        triple curr_edge = frontier_edges.top();
        frontier_edges.pop();
        if (seen[curr_edge.fst] && !seen[curr_edge.snd]){
            seen[curr_edge.snd] = true;
            weight += curr_edge.data;
            vertices_left -= 1;
            // add new frontier edges into the queue
            for (size_t i = 0; i < n_vertices; i++){
                if (graph_matrix[curr_edge.snd][i] != -1){
                    triple curr = {curr_edge.snd,i,graph_matrix[curr_edge.snd][i]};
                    frontier_edges.push(curr);
                }
            }
        }
    }
    return weight;
}

int main(int argc, char *argv[]) {
    string input_filename;
    int opt;
    while ((opt = getopt(argc, argv, "f:")) != -1) {
        switch (opt) {
            case 'f':
                input_filename = optarg;
                break;
            default:
                cerr << "blah blah blah usage1\n";
                exit(EXIT_FAILURE);
        }
    }
    if (empty(input_filename)) {
        cerr << "blah blah blah usage2\n";
        exit(EXIT_FAILURE);
    }

    cout << "Input file: " << input_filename << '\n';

    size_t n_verts, n_edges;
    bool hasweights = true;
    vector<triple> graph = readGraph(input_filename, n_verts, n_edges, hasweights);

    // int graph_matrix[n_verts][n_verts] = {numeric_limits<int>::max()};
    // for (auto &edge : graph) {
    //     graph_matrix[edge.fst][edge.snd] = edge.data;
    // }


    const auto compute_start = chrono::steady_clock::now();
    size_t mst_weight = prim_MST(graph, n_verts);
    cout << "weight of mst: " << mst_weight << "\n";

    const double compute_time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - compute_start).count();
    cout << "Computation time (sec): " << compute_time << "\n";
}