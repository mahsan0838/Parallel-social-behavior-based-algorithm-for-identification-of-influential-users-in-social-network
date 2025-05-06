/*
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>

using Node = int;

// Simple adjacency list
using AdjList = std::unordered_map<Node, std::vector<Node>>;

// Load the nodes assigned to this rank
std::unordered_set<Node> loadPartition(const std::string& fname) {
    std::unordered_set<Node> S;
    std::ifstream f(fname);
    Node u;
    while (f >> u) S.insert(u);
    return S;
}

// Extract subgraph: edges where src ∈ assigned
AdjList extractSubgraph(const std::unordered_set<Node>& assigned, const std::string& gle) {
    AdjList G;
    std::ifstream f(gle);
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream is(line);
        Node u, v; int l, s, c;
        if (!(is >> u >> v >> l >> s >> c)) continue;
        if (assigned.count(u)) {
            G[u].push_back(v);
        }
    }
    return G;
}

// Phase1: Compute influence power per node = out‐degree (parallel)
std::vector<std::pair<Node,int>> phase1(const AdjList& G,
                                        const std::unordered_set<Node>& assigned,
                                        int K) 
{
    // Copy into vector for OpenMP
    std::vector<Node> nodes(assigned.begin(), assigned.end());
    std::vector<std::pair<Node,int>> scores;
    scores.reserve(nodes.size());

    #pragma omp parallel
    {
        std::vector<std::pair<Node,int>> local;
        #pragma omp for nowait
        for (int i = 0; i < (int)nodes.size(); ++i) {
            Node u = nodes[i];
            int sc = G.count(u) ? G.at(u).size() : 0;
            local.emplace_back(u, sc);
        }
        #pragma omp critical
        {
            scores.insert(scores.end(), local.begin(), local.end());
        }
    }

    std::sort(scores.begin(), scores.end(),
              [](auto &a, auto &b) { return a.second > b.second; });
    if ((int)scores.size() > K)
        scores.resize(K);
    return scores;
}

// Phase2: Build influence‐BFS tree from each candidate
struct BFSInfo {
    int size;
    double avgDist;
};
BFSInfo buildTree(const AdjList& G, Node seed) {
    std::queue<std::pair<Node,int>> q;
    std::unordered_set<Node> vis;
    q.push({seed, 0});
    vis.insert(seed);
    long long sumd = 0;

    while (!q.empty()) {
        auto [u, d] = q.front(); q.pop();
        sumd += d;

        // safe adjacency lookup
        auto it = G.find(u);
        if (it == G.end()) continue;  // no outgoing edges

        for (Node v : it->second) {
            if (!vis.count(v)) {
                vis.insert(v);
                q.push({v, d + 1});
            }
        }
    }

    int S = vis.size();
    double A = S ? double(sumd) / S : 1e9;
    return {S, A};
}

// Phase2: Select best seed among candidates
Node phase2(const AdjList& G, const std::vector<Node>& cands) {
    Node best = cands[0];
    BFSInfo bi = buildTree(G, best);
    for (size_t i = 1; i < cands.size(); ++i) {
        Node cand = cands[i];
        auto cur = buildTree(G, cand);
        if (cur.size > bi.size ||
            (cur.size == bi.size && cur.avgDist < bi.avgDist)) {
            best = cand;
            bi = cur;
        }
    }
    return best;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, sz;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &sz);

    if (argc < 5) {
        if (!rank)
            std::cerr << "Usage: mpirun -np P ./psaiim <global.edgelist> <prefix> <K1> <K2>\n";
        MPI_Finalize();
        return 1;
    }

    std::string globalFile = argv[1];
    std::string prefix     = argv[2];
    int K1 = std::stoi(argv[3]);
    int K2 = std::stoi(argv[4]);  // for future use

    // Load this rank's partition and its subgraph
    auto assigned = loadPartition(prefix + std::to_string(rank) + ".txt");
    auto G        = extractSubgraph(assigned, globalFile);

    // Phase1: top‐K1 by out‐degree
    auto topK = phase1(G, assigned, K1);

    // Prepare candidates for Phase2
    std::vector<Node> cands;
    cands.reserve(topK.size());
    for (auto &p : topK) cands.push_back(p.first);

    // Phase2: select best seed among candidates
    Node seed = phase2(G, cands);
    std::cout << "Rank " << rank << " selected seed " << seed << "\n";

    MPI_Finalize();
    return 0;
}
*/

// mpic++ -fopenmp -O3 -o psaiim mpicode.cpp
// mpirun -np 3 ./psaiim out.edgelist partition_ 50 5




#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>

using Node = int;
// Adjacency list
using AdjList = std::unordered_map<Node, std::vector<Node>>;

//— Load this rank’s nodes from partition_<rank>.txt —
std::unordered_set<Node> loadPartition(const std::string& fname) {
    std::unordered_set<Node> S;
    std::ifstream f(fname);
    Node u;
    while (f >> u) S.insert(u);
    return S;
}

//— Extract local subgraph: keep edges u→v if u ∈ assigned —
AdjList extractSubgraph(const std::unordered_set<Node>& assigned,
                        const std::string& gle) 
{
    AdjList G;
    std::ifstream f(gle);
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream is(line);
        Node u, v; int l, s, c;
        if (!(is >> u >> v >> l >> s >> c)) continue;
        if (assigned.count(u)) {
            G[u].push_back(v);
        }
    }
    return G;
}

//— Phase 1: Top-K1 by out-degree (parallel via OpenMP) —
std::vector<std::pair<Node,int>> phase1(const AdjList& G,
                                        const std::unordered_set<Node>& assigned,
                                        int K1) 
{
    std::vector<Node> nodes(assigned.begin(), assigned.end());
    std::vector<std::pair<Node,int>> scores;
    scores.reserve(nodes.size());

    #pragma omp parallel
    {
        std::vector<std::pair<Node,int>> local;
        #pragma omp for nowait
        for (int i = 0; i < (int)nodes.size(); ++i) {
            Node u = nodes[i];
            int sc = G.count(u) ? G.at(u).size() : 0;
            local.emplace_back(u, sc);
        }
        #pragma omp critical
        scores.insert(scores.end(), local.begin(), local.end());
    }

    std::sort(scores.begin(), scores.end(),
              [](auto& a, auto& b){ return a.second > b.second; });
    if ((int)scores.size() > K1)
        scores.resize(K1);
    return scores;
}

//— Phase 2: Build influence‐BFS tree from seed and return (size, avgDist) —
struct BFSInfo { int size; double avgDist; };
BFSInfo buildTree(const AdjList& G, Node seed) {
    std::queue<std::pair<Node,int>> q;
    std::unordered_set<Node> vis;
    q.push({seed,0});
    vis.insert(seed);
    long long sumd = 0;
    while (!q.empty()) {
        auto [u,d] = q.front(); q.pop();
        sumd += d;
        auto it = G.find(u);
        if (it == G.end()) continue;
        for (Node v : it->second) {
            if (!vis.count(v)) {
                vis.insert(v);
                q.push({v, d+1});
            }
        }
    }
    int S = vis.size();
    double A = S ? double(sumd)/S : 1e9;
    return {S, A};
}

//— Phase 2: Among candidates, pick best by (max size, then min avgDist) —
Node phase2(const AdjList& G, const std::vector<Node>& cands) {
    Node best = cands[0];
    BFSInfo bi = buildTree(G, best);
    for (size_t i = 1; i < cands.size(); ++i) {
        BFSInfo cur = buildTree(G, cands[i]);
        if (cur.size > bi.size || (cur.size==bi.size && cur.avgDist < bi.avgDist)) {
            best = cands[i];
            bi = cur;
        }
    }
    return best;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, sz;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &sz);

    if (argc < 5) {
        if (!rank)
            std::cerr << "Usage: mpirun -np <P> ./psaiim "
                         "<global.edgelist> <partition_prefix> <K1> <K2>\n";
        MPI_Finalize();
        return 1;
    }

    std::string globalFile  = argv[1];
    std::string prefix      = argv[2];
    int K1 = std::stoi(argv[3]);
    // int K2 = std::stoi(argv[4]); // reserved for later

    // 1) Load this rank’s assigned node set
    auto assigned = loadPartition(prefix + std::to_string(rank) + ".txt");

    // 2) Extract local subgraph from global edgelist
    auto G = extractSubgraph(assigned, globalFile);

    // 3) Phase 1: pick top-K1 by out-degree
    auto topK = phase1(G, assigned, K1);

    // 4) Phase 2: from these, pick best seed via BFS influence
    std::vector<Node> cands;
    cands.reserve(topK.size());
    for (auto &p : topK) cands.push_back(p.first);
    Node localSeed = phase2(G, cands);

    // 5) Global coordination: gather all localSeeds to rank 0
    std::vector<Node> allSeeds;
    if (rank == 0) allSeeds.resize(sz);
    MPI_Gather(&localSeed, 1, MPI_INT,
               allSeeds.data(), 1, MPI_INT,
               0, MPI_COMM_WORLD);

    // 6) Rank 0 prints the combined seed list
    if (rank == 0) {
        std::cout << "=== Global Seed List ===\n";
        for (int r = 0; r < sz; ++r) {
            std::cout << " Rank " << r
                      << " → Seed " << allSeeds[r] << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}

