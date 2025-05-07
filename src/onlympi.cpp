

// mpicxx -o ompi onlympi.cpp
// time mpirun -np 3 ./ompi out.edgelist partition_ 10 0.1 0.5 0.8 0.85 100 50


#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <random>
#include <limits>
#include <cmath>

using Node = int;
using AdjList = std::unordered_map<Node, std::vector<std::pair<Node, double>>>;

// Load partition_<rank>.txt
std::unordered_set<Node> loadPartition(const std::string& fn, int rank) {
    std::unordered_set<Node> S;
    std::ifstream f(fn);
    if (!f.is_open()) {
        std::cerr << "[Rank " << rank << "] Error: Cannot open partition file " << fn << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (Node u; f >> u; ) S.insert(u);
    if (S.empty()) {
        std::cerr << "[Rank " << rank << "] Error: Partition file " << fn << " is empty or incorrectly formatted\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::cout << "[Rank " << rank << "] Loaded partition " << fn << " (" << S.size() << " nodes)\n";
    return S;
}

// Extract edges if u or v ∈ assigned. Collect ghosts.
AdjList extractSubgraph(
    const std::unordered_set<Node>& assigned,
    std::unordered_set<Node>& ghosts,
    const std::string& globalFile,
    double α, double β, double γ,
    int rank
) {
    AdjList G;
    std::ifstream f(globalFile);
    if (!f.is_open()) {
        std::cerr << "[Rank " << rank << "] Error: Cannot open global edge file " << globalFile << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::string line;
    size_t kept = 0;
    while (std::getline(f, line)) {
        std::istringstream is(line);
        Node u, v; int likes, shares, comments;
        if (!(is >> u >> v >> likes >> shares >> comments)) continue;
        bool inU = assigned.count(u), inV = assigned.count(v);
        if (!inU && !inV) continue;
        if (!inU) ghosts.insert(u);
        if (!inV) ghosts.insert(v);
        double w = α * likes + β * shares + γ * comments;
        G[u].emplace_back(v, w);
        ++kept;
    }
    if (kept == 0) {
        std::cerr << "[Rank " << rank << "] Warning: No edges extracted for subgraph (assigned nodes: " << assigned.size() << ")\n";
    }
    std::cout << "[Rank " << rank << "] Extracted subgraph: " << kept
              << " edges, " << ghosts.size() << " ghost nodes\n";
    return G;
}

// Phase 1: Weighted PageRank on full G, but only keep top K1 *assigned* scores
std::vector<std::pair<Node, double>> phase1_pagerank(
    const AdjList& G,
    const std::unordered_set<Node>& assigned,
    int K1,
    double d, int maxIter, double tol,
    int rank
) {
    std::cout << "[Rank " << rank << "] Starting PageRank (d=" << d << ", maxIter=" << maxIter << ")…\n";
    // Build node list & index
    std::vector<Node> nodes; nodes.reserve(G.size());
    for (auto& kv : G) nodes.push_back(kv.first);
    int N = nodes.size();
    if (N == 0) {
        std::cerr << "[Rank " << rank << "] Warning: Graph is empty, skipping PageRank\n";
        return std::vector<std::pair<Node, double>>();
    }
    std::unordered_map<Node, int> idx; idx.reserve(N);
    for (int i = 0; i < N; ++i) idx[nodes[i]] = i;

    // Build reverse adjacency + out-weight sums
    std::vector<double> outW(N, 0);
    std::vector<std::vector<std::pair<int, double>>> inAdj(N);
    for (int i = 0; i < N; ++i) {
        for (auto& e : G.at(nodes[i])) {
            auto it = idx.find(e.first);
            if (it == idx.end()) continue;
            int j = it->second;
            outW[i] += e.second;
            inAdj[j].emplace_back(i, e.second);
        }
    }

    // PageRank iterations
    std::vector<double> PR(N, 1.0 / N), next(N);
    for (int it = 0; it < maxIter; ++it) {
        double diff = 0;
        for (int i = 0; i < N; ++i) {
            double sum = 0;
            for (auto& p : inAdj[i]) {
                int j = p.first; double w = p.second;
                if (outW[j] > 0) sum += PR[j] * (w / outW[j]);
            }
            next[i] = (1 - d) / N + d * sum;
            diff += std::fabs(next[i] - PR[i]);
        }
        PR.swap(next);
        if (diff < tol) {
            std::cout << "[Rank " << rank << "] PageRank converged at iter " << it << "\n";
            break;
        }
        if (it == maxIter - 1)
            std::cout << "[Rank " << rank << "] PageRank reached maxIter\n";
    }

    // Collect only *assigned* nodes and pick top K1
    std::vector<std::pair<Node, double>> scored;
    scored.reserve(assigned.size());
    for (int i = 0; i < N; ++i) {
        if (assigned.count(nodes[i]))
            scored.emplace_back(nodes[i], PR[i]);
    }
    std::partial_sort(
        scored.begin(),
        scored.begin() + std::min(K1, (int)scored.size()),
        scored.end(),
        [](auto& a, auto& b) { return a.second > b.second; }
    );
    if ((int)scored.size() > K1) scored.resize(K1);

    std::cout << "[Rank " << rank << "] Top " << K1 << " seeds by PageRank:\n";
    for (auto& pr : scored)
        std::cout << "   Node " << pr.first << " → PR=" << pr.second << "\n";

    return scored;
}

// BFS-based influence zone to measure spread size & avg dist
struct BFSInfo { int size; double avgDist; };
BFSInfo buildTree(const AdjList& G, Node seed) {
    std::queue<std::pair<Node, int>> q;
    std::unordered_set<Node> vis;
    q.push({seed, 0}); vis.insert(seed);
    long long sumd = 0;
    while (!q.empty()) {
        auto [u, d] = q.front(); q.pop();
        sumd += d;
        auto it = G.find(u);
        if (it != G.end()) {
            for (auto& p : it->second) {
                Node v = p.first;
                if (!vis.count(v)) {
                    vis.insert(v);
                    q.push({v, d + 1});
                }
            }
        }
    }
    int S = vis.size();
    double A = S ? double(sumd) / S : 1e9;
    return {S, A};
}

// Simple greedy on small candidate set: pick K1 maximizing BFS size
std::vector<Node> greedyChoose(
    const AdjList& G,
    const std::vector<Node>& cands,
    int K1,
    int rank
) {
    std::cout << "[Rank " << rank << "] Running greedy over " << cands.size() << " candidates…\n";
    std::vector<Node> chosen;
    std::unordered_set<Node> used;
    for (int k = 0; k < K1 && (int)used.size() < (int)cands.size(); ++k) {
        Node bestN = -1; int bestSz = -1; double bestAvg = 0;
        for (auto c : cands) {
            if (used.count(c)) continue;
            auto info = buildTree(G, c);
            if (info.size > bestSz) {
                bestSz = info.size;
                bestAvg = info.avgDist;
                bestN = c;
            }
        }
        if (bestN < 0) break;
        chosen.push_back(bestN);
        used.insert(bestN);
        std::cout << "   Chosen " << bestN << " (spread=" << bestSz
                  << ", avgDist=" << bestAvg << ")\n";
    }
    return chosen;
}

// IC simulation (Independent Cascade) returns total infected count in subgraph
int simulateIC(
    const AdjList& G,
    const std::vector<Node>& seeds,
    int trials,
    std::mt19937_64& rnd,
    int rank
) {
    std::cout << "[Rank " << rank << "] Running IC simulation (" << trials << " trials)…\n";
    // Build adjacency lookup
    std::unordered_map<Node, const std::vector<std::pair<Node, double>>*> adjPtr;
    for (auto& kv : G) adjPtr[kv.first] = &kv.second;

    int total = 0;
    for (int t = 0; t < trials; ++t) {
        std::unordered_set<Node> activated(seeds.begin(), seeds.end());
        std::queue<Node> q;
        for (auto s : seeds) q.push(s);
        std::mt19937_64 gen(rnd());
        std::uniform_real_distribution<double> dist(0, 1);
        while (!q.empty()) {
            Node u = q.front(); q.pop();
            auto it = adjPtr.find(u);
            if (it == adjPtr.end()) continue;
            for (auto& p : *it->second) {
                Node v = p.first; double w = p.second;
                if (!activated.count(v) && dist(gen) < w) {
                    activated.insert(v);
                    q.push(v);
                }
            }
        }
        total += activated.size();
    }
    double avg = double(total) / trials;
    std::cout << "[Rank " << rank << "] IC local average infected = " << avg << "\n";
    return total;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, sz; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &sz);

    if (argc != 10 && !rank) {
        std::cerr << "Usage: mpirun -np P ./psaiim "
                     "out.edgelist part_pref K1 α β γ d maxIter ICtrials\n";
    }
    if (argc != 10) {
        MPI_Finalize();
        return 0;
    }

    std::string globalFile = argv[1];
    std::string prefix = argv[2];
    int K1 = std::stoi(argv[3]);
    double α = std::atof(argv[4]),
           β = std::atof(argv[5]),
           γ = std::atof(argv[6]),
           d = std::atof(argv[7]);
    int maxIter = std::stoi(argv[8]),
        ICtrials = std::stoi(argv[9]);
    double tol = 1e-4;

    // Validate damping factor
    if (d <= 0 || d >= 1) {
        if (!rank) {
            std::cerr << "Error: Damping factor d (" << d << ") must be between 0 and 1\n";
        }
        MPI_Finalize();
        return 1;
    }

    // 1) Load partition & build local+ghost subgraph
    auto assigned = loadPartition(prefix + std::to_string(rank) + ".txt", rank);
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize to ensure all ranks report

    std::unordered_set<Node> ghosts;
    auto G = extractSubgraph(assigned, ghosts, globalFile, α, β, γ, rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // 2) Phase 1: local top-K₁ candidates by weighted-PageRank
    auto localTopK = phase1_pagerank(G, assigned, K1, d, maxIter, tol, rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // 3) Gather *all* candidates at rank 0
    int cnt = localTopK.size();
    std::vector<int> counts(sz);
    MPI_Gather(&cnt, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<Node> flat(cnt);
    for (int i = 0; i < cnt; ++i) flat[i] = localTopK[i].first;

    std::vector<int> displs, recvbuf;
    if (rank == 0) {
        displs.resize(sz);
        int off = 0;
        for (int r = 0; r < sz; ++r) {
            displs[r] = off; off += counts[r];
        }
        recvbuf.resize(off);
    }
    MPI_Gatherv(
        flat.data(), cnt, MPI_INT,
        recvbuf.data(), counts.data(), displs.data(), MPI_INT,
        0, MPI_COMM_WORLD
    );
    MPI_Barrier(MPI_COMM_WORLD);

    // 4) Rank 0 does greedy to pick final K₁ seeds
    std::vector<Node> finalSeeds;
    if (rank == 0) {
        finalSeeds = greedyChoose(G, recvbuf, K1, rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // 5) Broadcast finalSeeds
    int M = finalSeeds.size();
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    finalSeeds.resize(M);
    MPI_Bcast(finalSeeds.data(), M, MPI_INT, 0, MPI_COMM_WORLD);

    // 6) Simulate IC on the *global* seed set in each local subgraph
    std::mt19937_64 rnd(rank + 1234567);
    int localSpread = simulateIC(G, finalSeeds, ICtrials, rnd, rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // 7) Reduce spreads to rank 0
    int totalSpread = 0;
    MPI_Reduce(&localSpread, &totalSpread, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!rank) {
        std::cout << "\n=== Final Global Seeds ===\n";
        for (int i = 0; i < M; ++i)
            std::cout << " Seed " << i << " → " << finalSeeds[i] << "\n";
        std::cout << "Estimated total infected (avg over "
                  << ICtrials << " trials): "
                  << double(totalSpread) / ICtrials << "\n";
    }

    MPI_Finalize();
    return 0;
}



