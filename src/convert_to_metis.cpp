// File: convert_to_metis.cpp
// Compile: g++ -std=c++17 -O2 -o convert_to_metis convert_to_metis.cpp
// Run: ./convert_to_metis HiggsCombined.edgelist higgs_metis.graph



// this was made to process HiggsCombined.edgelist for metis processing

// File: edgelist_to_metis.cpp
// Compile: g++ -std=c++17 -O2 -o convert convert_to_metis.cpp
// Usage:   ./convert HiggsCombined.edgelist higgs.graph

#include <bits/stdc++.h>
using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0]
             << " <input.edgelist> <output.graph>\n";
        return 1;
    }
    ifstream fin(argv[1]);
    if (!fin) {
        cerr << "Error: cannot open " << argv[1] << "\n";
        return 1;
    }

    // Step 1: Read edges, build undirected adjacency
    unordered_map<int, unordered_set<int>> adj;
    int u, v; double w;
    int max_node = 0;
    while (fin >> u >> v >> w) {
        if (u == v) continue;        // skip self-loops
        adj[u].insert(v);
        adj[v].insert(u);
        max_node = max({max_node, u, v});
    }
    fin.close();

    int n = max_node;
    // Ensure every node 1..n appears, even if isolated
    for (int i = 1; i <= n; ++i)
        adj.try_emplace(i);

    // Step 2: Count undirected edges once
    long long m = 0;
    for (auto &p : adj) {
        m += p.second.size();
    }
    m /= 2;

    // Step 3: Write METIS .graph
    ofstream fout(argv[2]);
    if (!fout) {
        cerr << "Error: cannot write to " << argv[2] << "\n";
        return 1;
    }
    fout << n << " " << m << "\n";
    for (int i = 1; i <= n; ++i) {
        auto &nbrs = adj[i];
        vector<int> list(nbrs.begin(), nbrs.end());
        sort(list.begin(), list.end());
        for (int x : list) fout << x << " ";
        fout << "\n";
    }
    fout.close();

    cout << "Converted to METIS format: " << n << " nodes, " << m << " edges\n";
    return 0;
}

