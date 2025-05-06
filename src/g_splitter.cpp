// File: psa_iim_partition.cpp
// Compile with: g++ -std=c++17 -O2 -o pp gptp.cpp
// Usage: ./pp HiggsCombined.edgelist 3

#include <bits/stdc++.h>
using namespace std;

// Aliases
using VI = vector<int>;
using VVI = vector<VI>;
using ll = long long;

// Globals for graph
vector<vector<int>> graph, rev_graph;
vector<bool> visited;
stack<int> finish_stack;
VVI scc_list;
vector<int> node_to_scc;

// Iterative Kosaraju Phase 1: fill finish_stack
void kosaraju_phase1(int N) {
    visited.assign(N, false);
    for (int start = 0; start < N; ++start) {
        if (visited[start]) continue;
        // stack of (node, next-child-index)
        stack<pair<int,int>> st;
        st.push({start, 0});
        while (!st.empty()) {
            auto &top = st.top();
            int u = top.first, &idx = top.second;
            if (!visited[u]) visited[u] = true;
            if (idx < (int)graph[u].size()) {
                int v = graph[u][idx++];
                if (!visited[v]) st.push({v, 0});
            } else {
                finish_stack.push(u);
                st.pop();
            }
        }
    }
}

// Iterative Kosaraju Phase 2: build SCCs
void kosaraju_phase2(int N) {
    visited.assign(N, false);
    node_to_scc.assign(N, -1);
    while (!finish_stack.empty()) {
        int start = finish_stack.top();
        finish_stack.pop();
        if (visited[start]) continue;
        // new component
        VI comp;
        stack<int> st;
        st.push(start);
        while (!st.empty()) {
            int u = st.top(); st.pop();
            if (visited[u]) continue;
            visited[u] = true;
            node_to_scc[u] = scc_list.size();
            comp.push_back(u);
            for (int v : rev_graph[u]) {
                if (!visited[v]) st.push(v);
            }
        }
        scc_list.push_back(move(comp));
    }
}

// Build the SCC meta-graph (DAG)
void build_scc_graph(VVI &scc_graph) {
    int S = scc_list.size();
    scc_graph.assign(S, {});
    unordered_set<ll> seen;
    for (int u = 0; u < (int)graph.size(); ++u) {
        int su = node_to_scc[u];
        for (int v : graph[u]) {
            int sv = node_to_scc[v];
            if (su != sv) {
                ll code = ((ll)su << 32) | (unsigned long long)sv;
                if (!seen.count(code)) {
                    scc_graph[su].push_back(sv);
                    seen.insert(code);
                }
            }
        }
    }
}

// Topological sort on DAG
void topo_sort(const VVI &dag, VI &order) {
    int S = dag.size();
    VI indeg(S, 0);
    for (int u = 0; u < S; ++u)
        for (int v : dag[u])
            indeg[v]++;
    queue<int> q;
    for (int i = 0; i < S; ++i)
        if (indeg[i] == 0) q.push(i);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        order.push_back(u);
        for (int v : dag[u])
            if (--indeg[v] == 0)
                q.push(v);
    }
}

// Extract CACs from SCC-DAG
void extract_cacs(const VVI &dag, vector<VI> &cac_list) {
    int S = dag.size();
    // Compute in-degrees
    VI indeg(S, 0);
    for (int u = 0; u < S; ++u)
        for (int v : dag[u])
            indeg[v]++;
    vector<bool> used(S, false);
    VI order;
    topo_sort(dag, order);
    for (int u : order) {
        if (used[u]) continue;
        // build one CAC
        VI cac;
        stack<int> st; 
        st.push(u);
        while (!st.empty()) {
            int x = st.top(); st.pop();
            if (used[x]) continue;
            used[x] = true;
            cac.push_back(x);
            // only traverse children with single in-edge
            for (int v : dag[x]) {
                if (!used[v] && indeg[v] == 1)
                    st.push(v);
            }
        }
        cac_list.push_back(move(cac));
    }
}

// Balanced partitioning of CACs by total node count
void balanced_partition(const vector<VI> &cac_list,
                        vector<vector<int>> &partitions,
                        int k) {
    partitions.assign(k, {});
    vector<ll> part_size(k, 0);
    // Sort CACs descending by size
    vector<int> idx(cac_list.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) {
        return cac_list[a].size() > cac_list[b].size();
    });
    // Greedy balance
    for (int id : idx) {
        // find partition with min size
        int p = min_element(part_size.begin(), part_size.end()) - part_size.begin();
        // assign all nodes of each SCC in this CAC to p
        for (int scc_id : cac_list[id]) {
            for (int node : scc_list[scc_id]) {
                partitions[p].push_back(node);
            }
        }
        part_size[p] += cac_list[id].size(); 
    }
}

// Main
int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0]
             << " <combined.edgelist> <num_partitions>\n";
        return 1;
    }
    string infile = argv[1];
    int k = stoi(argv[2]);

    // 1) Read edges, find max node
    vector<tuple<int,int,double>> edges;
    edges.reserve(50000000);
    ifstream fin(infile);
    if (!fin) {
        cerr << "Cannot open " << infile << "\n";
        return 1;
    }
    int u, v; double w;
    int max_node = 0;
    while (fin >> u >> v >> w) {
        edges.emplace_back(u,v,w);
        max_node = max({max_node, u, v});
    }
    fin.close();

    int N = max_node + 1;
    graph.assign(N, {});
    rev_graph.assign(N, {});

    // 2) Build graph (ignore weights for partitioning)
    for (auto &e : edges) {
        tie(u,v,w) = e;
        graph[u].push_back(v);
        rev_graph[v].push_back(u);
    }

    // 3) Compute SCCs
    kosaraju_phase1(N);
    kosaraju_phase2(N);

    // 4) Build SCC-DAG
    VVI scc_graph;
    build_scc_graph(scc_graph);

    // 5) Extract CACs
    vector<VI> cac_list;
    extract_cacs(scc_graph, cac_list);

    // 6) Balanced partition CACs
    vector<vector<int>> partitions;
    balanced_partition(cac_list, partitions, k);

    // 7) Write partitions
    for (int i = 0; i < k; ++i) {
        ofstream fout("partition_" + to_string(i) + ".txt");
        for (int node : partitions[i])
            fout << node << "\n";
    }

    cout << "Done: " << k << " balanced partitions written.\n";
    return 0;
}

