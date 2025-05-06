
/*
    performance Analysis of sequential code: 
    NO_OF_EDGES        execution time (s)
    10                0.01s user 0.00s system 94% cpu 0.006 total
    100               0.01s user 0.00s system 105% cpu 0.010 total
    1000              0.05s user 0.01s system 101% cpu 0.057 total
    10000             0.53s user 0.00s system 107% cpu 0.494 total
    100000            6.92s user 0.01s system 102% cpu 6.772 total
    150000            33.97s user 0.08s system 99% cpu 34.342 total
    200000            86.87s user 0.27s system 98% cpu 1:28.07 total
    250000            170.33s user 0.84s system 99% cpu 2:51.55 total


*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <numeric>
#include <limits>
#include <utility>
#include <cmath>
#define VERBOSE 0
#define NO_OF_EDGES 250000
using namespace std;

// Custom hash function for std::pair<int, int>
struct PairHash {
    size_t operator()(const pair<int, int>& p) const {
        // Combine hashes of first and second using a simple but effective method
        size_t h1 = hash<int>{}(p.first);
        size_t h2 = hash<int>{}(p.second);
        return h1 ^ (h2 << 1); // Shift and XOR to combine
    }
};

// Structure to represent a social network interaction edge
struct Edge {
    int src;
    int dest;
    int weight1; // Value of liking post
    int weight2; // Value of sharing post
    int weight3; // Value of commenting post
};

// Structure to represent a vertex in the graph
struct Vertex {
    int id;
    int index = -1;      // DFS index
    int lowlink = -1;    // Lowlink value for SCC
    int level = -1;      // Level in the hierarchy
    int depth = -1;      // Depth in DFS tree
    char type = '\0';    // 's' for SCC or 'c' for CAC
    int comp = -1;       // Component ID
    bool initialized = false;
    bool onStack = false; // Flag to track if vertex is on stack
    vector<int> neighbors; // Stores indices of neighboring vertices
    vector<int> incoming_edges; // Stores indices of vertices with edges to this vertex
};

class Graph {
private:
    vector<Vertex> vertices;
    unordered_map<int, int> idToIndex; // Maps vertex ID to index in vertices
    stack<int> dfsStack;
    int index = 0;
    int nextCompId = 0;
    vector<double> influenceScores;
    vector<int> seedCandidates;
    vector<int> seeds;
    vector<vector<int>> bfsTrees; // Cache for BFS trees

    // Initialize vertex state
    void resetVertexStates() {
        for (auto& vertex : vertices) {
            vertex.index = -1;
            vertex.lowlink = -1;
            vertex.level = -1;
            vertex.depth = -1;
            vertex.type = '\0';
            vertex.comp = -1;
            vertex.initialized = false;
            vertex.onStack = false;
        }
        index = 0;
        nextCompId = 0;
        while (!dfsStack.empty()) dfsStack.pop();
    }

    // Discover function for DFS
    void discover(int v) {
        Vertex& vertex = vertices[v];
        vertex.index = index;
        vertex.lowlink = index;
        vertex.level = 1;
        vertex.depth = 1;
        index++;
        dfsStack.push(v);
        vertex.onStack = true;
        vertex.initialized = true;
    }

    // Explore function for DFS
    void explore(int v) {
        Vertex& vVertex = vertices[v];
        for (int w_idx : vVertex.neighbors) {
            if (w_idx < 0 || w_idx >= vertices.size()) continue; // Safety check
            Vertex& wVertex = vertices[w_idx];
            if (!wVertex.initialized) {
                discover(w_idx);
                explore(w_idx);
                vVertex.lowlink = min(vVertex.lowlink, vertices[w_idx].lowlink);
            } else if (wVertex.onStack) {
                vVertex.lowlink = min(vVertex.lowlink, wVertex.index);
            }
        }
        // Only call finish if v is the root of an SCC/CAC
        if (vVertex.lowlink == vVertex.index) {
            finish(v);
        }
    }

    // Finish function to assign SCC/CAC
    void finish(int v) {
        Vertex& vVertex = vertices[v];
        int compId = nextCompId++;
        vector<int> componentMembers;
        int size = 0;
        int w_idx;

        do {
            w_idx = dfsStack.top();
            dfsStack.pop();
            Vertex& wVertex = vertices[w_idx];
            wVertex.onStack = false;
            wVertex.comp = compId;
            componentMembers.push_back(w_idx);
            size++;
        } while (w_idx != v);

        char componentType = (size > 1) ? 's' : 'c';
        for (int member_idx : componentMembers) {
            vertices[member_idx].type = componentType;
        }
    }

    // Compute followers (in-degree) and followees (out-degree)
    void computeFollowersFollowees(vector<int>& followers, vector<int>& followees) const {
        followers.assign(vertices.size(), 0);
        followees.assign(vertices.size(), 0);
        for (size_t i = 0; i < vertices.size(); ++i) {
            followees[i] = vertices[i].neighbors.size();
            followers[i] = vertices[i].incoming_edges.size();
        }
    }

    // Calculate influence probability psi(u,v) using edge weights
    double calculatePsi(int u_idx, int v_idx) {
        // Example: Check if there's a direct edge from u to v
        const auto& neighbors = vertices[u_idx].neighbors;
        bool connected = find(neighbors.begin(), neighbors.end(), v_idx) != neighbors.end();

        if (connected) {
            // Basic example: return 1.0 / followees[u_idx] if connected, else 0
            // Or use edge weights if available and relevant
             size_t out_degree = vertices[u_idx].neighbors.size();
             return (out_degree > 0) ? (1.0 / out_degree) : 0.0; // Simplified influence distribution
        }
        return 0.0; // Return 0 if not connected or based on your formula
    }

    // Build BFS tree for influence calculation
    void buildInfluenceBFSTree(int start_node_idx, vector<int>& tree) const {
        tree.clear();
        if (start_node_idx < 0 || start_node_idx >= vertices.size()) return;

        queue<int> q;
        vector<bool> visited(vertices.size(), false);
        q.push(start_node_idx);
        visited[start_node_idx] = true;
        tree.push_back(start_node_idx);

        while (!q.empty()) {
            int current_idx = q.front();
            q.pop();
            for (int neighbor_idx : vertices[current_idx].neighbors) {
                if (neighbor_idx >= 0 && neighbor_idx < vertices.size() && !visited[neighbor_idx]) {
                    visited[neighbor_idx] = true;
                    q.push(neighbor_idx);
                    tree.push_back(neighbor_idx);
                }
            }
        }
    }

public:
    // Constructor with edge validation
    Graph(const vector<Edge>& edges) {
        unordered_set<pair<int, int>, PairHash> edgeSet;
        unordered_map<int, int> tempIdMap;
        int currentIdx = 0;

        // First pass: Identify unique vertices
        for (const auto& edge : edges) {
            if (edge.src < 0 || edge.dest < 0) {
                cerr << "Warning: Skipping edge with negative ID (" << edge.src << " -> " << edge.dest << ")." << endl;
                continue;
            }
            if (edge.weight1 < 0 || edge.weight2 < 0 || edge.weight3 < 0) {
                cerr << "Warning: Edge (" << edge.src << " -> " << edge.dest << ") has negative weights." << endl;
            }
            tempIdMap.emplace(edge.src, 0);
            tempIdMap.emplace(edge.dest, 0);
        }

        // Assign indices to vertices
        for (auto& pair : tempIdMap) {
            pair.second = currentIdx++;
        }

        vertices.resize(tempIdMap.size());
        for (const auto& pair : tempIdMap) {
            idToIndex[pair.first] = pair.second;
            vertices[pair.second].id = pair.first;
        }

        // Second pass: Add edges
        for (const auto& edge : edges) {
            if (idToIndex.count(edge.src) == 0 || idToIndex.count(edge.dest) == 0) continue;
            int srcIdx = idToIndex[edge.src];
            int destIdx = idToIndex[edge.dest];
            auto edgePair = make_pair(srcIdx, destIdx);
            if (edgeSet.insert(edgePair).second) {
                vertices[srcIdx].neighbors.push_back(destIdx);
                vertices[destIdx].incoming_edges.push_back(srcIdx);
            }
        }

        // Reserve space for BFS trees cache
        bfsTrees.resize(vertices.size());
    }

    // SCC/CAC partitioning
    void sccCacPartitioning() {
        if (vertices.empty()) {
            cerr << "Warning: Graph is empty. No SCC/CAC partitioning performed." << endl;
            return;
        }

        resetVertexStates();

        // Process all vertices, including isolated ones
        for (size_t v = 0; v < vertices.size(); ++v) {
            if (!vertices[v].initialized) {
                if (vertices[v].neighbors.empty() && vertices[v].incoming_edges.empty()) {
                    // Isolated vertex: treat as CAC
                    vertices[v].comp = nextCompId++;
                    vertices[v].type = 'c';
                    vertices[v].level = 1;
                    vertices[v].initialized = true;
                } else {
                    discover(v);
                    explore(v);
                }
            }
        }

        // Build component graph for level assignment
        vector<vector<int>> componentAdj(nextCompId);
        vector<int> componentInDegree(nextCompId, 0);
        unordered_set<pair<int, int>, PairHash> edgeExists;

        for (size_t u_idx = 0; u_idx < vertices.size(); ++u_idx) {
            int compU = vertices[u_idx].comp;
            if (compU < 0) continue;
            for (int v_idx : vertices[u_idx].neighbors) {
                int compV = vertices[v_idx].comp;
                if (compV < 0 || compU == compV) continue;
                auto edge = make_pair(compU, compV);
                if (edgeExists.insert(edge).second) {
                    componentAdj[compU].push_back(compV);
                    componentInDegree[compV]++;
                }
            }
        }

        // Topological sort for levels
        vector<int> componentLevels(nextCompId, 1);
        queue<int> q;
        for (int i = 0; i < nextCompId; ++i) {
            if (componentInDegree[i] == 0) q.push(i);
        }

        while (!q.empty()) {
            int currentComp = q.front();
            q.pop();
            for (int neighborComp : componentAdj[currentComp]) {
                componentLevels[neighborComp] = max(componentLevels[neighborComp], componentLevels[currentComp] + 1);
                if (--componentInDegree[neighborComp] == 0) {
                    q.push(neighborComp);
                }
            }
        }

        for (size_t v = 0; v < vertices.size(); ++v) {
            if (vertices[v].comp >= 0) {
                vertices[v].level = componentLevels[vertices[v].comp];
            }
        }
    }

    // Print SCC/CAC results
    void printResults() const {
        cout << "SCC/CAC Partitioning Results:" << endl;
        cout << "------------------------------" << endl;

        map<int, vector<int>> components;
        map<int, char> componentTypes;
        map<int, int> componentLevels;
        int sccCount = 0, cacCount = 0;

        for (const auto& vertex : vertices) {
            if (vertex.comp >= 0) {
                components[vertex.comp].push_back(vertex.id);
                componentTypes[vertex.comp] = vertex.type;
                componentLevels[vertex.comp] = vertex.level;
                if (vertex.type == 's') sccCount++;
                else cacCount++;
            }
        }

        if(VERBOSE) {

        for (const auto& pair : components) {
            int compId = pair.first;
            const auto& compVertices = pair.second;
            string typeStr = componentTypes[compId] == 's' ? "SCC" : "CAC";
                cout << "Component " << compId << " (" << typeStr << ") - Level: "
                     << componentLevels[compId] << endl;
                       cout << "Vertices: ";
            
            vector<int> sortedVertices = compVertices;
            sort(sortedVertices.begin(), sortedVertices.end());
            for (size_t i = 0; i < sortedVertices.size(); ++i) {
                    cout << sortedVertices[i] << (i < sortedVertices.size() - 1 ? ", " : "");
            
            cout << endl << endl;
        }
        }
        }

        cout << "Summary:" << endl;
        cout << "Total components: " << nextCompId << endl;
        cout << "Non-empty components: " << components.size() << endl;
        cout << "SCC count: " << sccCount << endl;
        cout << "CAC count: " << cacCount << endl;
    }

    void parallelInfluencePowerMeasure(const vector<string>& actions, const vector<Edge>& edges) {
        if (vertices.empty()) {
            cerr << "Warning: Graph is empty. No influence scores calculated." << endl;
            return;
        }
        influenceScores.assign(vertices.size(), 0.0);
        vector<int> followers, followees;
        computeFollowersFollowees(followers, followees);
        double d = 0.85;
        for (size_t u = 0; u < vertices.size(); ++u) {
            double sum_psi = 0.0;
            for (int v_idx : vertices[u].neighbors) {
                sum_psi += calculatePsi(u, v_idx);
            }
            double F_u = vertices.size() > 0 ? static_cast<double>(followers[u]) / vertices.size() : 0.0;
            influenceScores[u] = (1.0 - d) * F_u + d * sum_psi;
            if (influenceScores[u] < 0 || !isfinite(influenceScores[u])) {
                influenceScores[u] = 0.0;
            }
        }
    }

    // Seed candidates selection
    void seedCandidatesSelection() {
        seedCandidates.clear();
        if (vertices.empty() || influenceScores.empty()) {
            cerr << "Warning: No seed candidates selected. Graph or influence scores empty." << endl;
            return;
        }

        vector<pair<double, int>> nodeIPs;
        bfsTrees.assign(vertices.size(), vector<int>()); // Clear cache

        for (size_t i = 0; i < vertices.size(); ++i) {
            if (influenceScores[i] <= 0) continue;
            double I_P_v = 0.0;
            buildInfluenceBFSTree(i, bfsTrees[i]);
            for (int node_idx : bfsTrees[i]) {
                if (node_idx >= 0 && node_idx < influenceScores.size()) {
                    I_P_v += influenceScores[node_idx];
                }
            }
            nodeIPs.emplace_back(I_P_v, i);
        }

        sort(nodeIPs.rbegin(), nodeIPs.rend());
        size_t numCandidates = max<size_t>(1, vertices.size() / 5);
        numCandidates = min(numCandidates, nodeIPs.size());

        for (size_t i = 0; i < numCandidates; ++i) {
            seedCandidates.push_back(nodeIPs[i].second);
        }
        sort(seedCandidates.begin(), seedCandidates.end());

        if (seedCandidates.empty() && !vertices.empty()) {
            double max_L = -1.0;
            int best_cand = -1;
            for (size_t i = 0; i < influenceScores.size(); ++i) {
                if (influenceScores[i] > max_L) {
                    max_L = influenceScores[i];
                    best_cand = i;
                }
            }
            if (best_cand >= 0) {
                seedCandidates.push_back(best_cand);
            }
        }
    }

    // Seed selection
    vector<int> seedSelection(int k) {
        seeds.clear();
        if (k <= 0 || vertices.empty() || influenceScores.empty()) {
            cerr << "Warning: Invalid k or empty graph/scores. No seeds selected." << endl;
            return seeds;
        }

        k = min(k, static_cast<int>(vertices.size()));
        vector<int> nodesPool = seedCandidates.empty() ? vector<int>(vertices.size()) : seedCandidates;
        if (nodesPool.empty()) {
            iota(nodesPool.begin(), nodesPool.end(), 0);
        }
        k = min(k, static_cast<int>(nodesPool.size()));

        set<int> availableNodes(nodesPool.begin(), nodesPool.end());
        set<int> selectedSeedSet;

        for (int i = 0; i < k && !availableNodes.empty(); ++i) {
            int best_v_idx = -1;
            double minRank = numeric_limits<double>::max();

            for (int v_idx : availableNodes) {
                if (v_idx < 0 || v_idx >= vertices.size()) continue;
                const auto& tree = bfsTrees[v_idx]; // Use cached BFS tree
                double currentRank = 0.0;
                for (int u_idx : tree) {
                    if (u_idx >= 0 && u_idx < influenceScores.size() &&
                        selectedSeedSet.count(u_idx) == 0) {
                        currentRank += influenceScores[u_idx];
                    }
                }
                double effectiveRank = currentRank - static_cast<double>(v_idx) * 1e-9;
                if (effectiveRank < minRank) {
                    minRank = effectiveRank;
                    best_v_idx = v_idx;
                }
            }

            if (best_v_idx >= 0) {
                seeds.push_back(best_v_idx);
                selectedSeedSet.insert(best_v_idx);
                availableNodes.erase(best_v_idx);
            } else {
                cerr << "Warning: Could not select seed in iteration " << i + 1 << "." << endl;
                break;
            }
        }
        return seeds;
    }

    // Print influence results
    void printInfluenceResults() const {
        cout << "\nInfluence Calculation Results:" << endl;
        cout << "------------------------------" << endl;

        if (influenceScores.empty()) {
            cout << "Influence scores not calculated." << endl;
            return;
        }

        vector<pair<double, int>> scoreIdPairs;
        for (size_t i = 0; i < vertices.size(); ++i) {
            if (i < influenceScores.size()) {
                scoreIdPairs.emplace_back(influenceScores[i], vertices[i].id);
            }
        }
        sort(scoreIdPairs.rbegin(), scoreIdPairs.rend());

        cout << "Influence Scores (Sorted Descending):" << endl;
        for (const auto& pair : scoreIdPairs) {
            if (isfinite(pair.first) && pair.first >= 0) {
                if(VERBOSE) {
                    cout << "Vertex " << pair.second << ": " << pair.first << endl;
                }
            } else {
                cout << "Vertex " << pair.second << ": Invalid Score" << endl;
            }
        }

        if (VERBOSE) { 
        cout << "\nSeed Candidates: ";
        if (seedCandidates.empty()) {
            cout << "None.";
        } else {
            for (size_t i = 0; i < seedCandidates.size(); ++i) {
                if (seedCandidates[i] >= 0 && seedCandidates[i] < vertices.size()) {
                    cout << vertices[seedCandidates[i]].id << (i == seedCandidates.size() - 1 ? "" : ", ");
                }
            }
        }
        cout << endl;
    }
        cout << "Selected Seeds (k=" << seeds.size() << "): ";
        if (seeds.empty()) {
            cout << "None.";
        } else {
            for (size_t i = 0; i < seeds.size(); ++i) {
                if (seeds[i] >= 0 && seeds[i] < vertices.size()) {
                    cout << vertices[seeds[i]].id << (i == seeds.size() - 1 ? "" : ", ");
                }
            }
        }
        cout << endl;
    }
};

vector<Edge> readEdgeList(const string& filename) {
    vector<Edge> edges; // Vector to store the edges
    edges.reserve(NO_OF_EDGES); // Preallocate space for up to 1000 edges
    ifstream file(filename); // Open the file for reading
    if (!file.is_open()) { // Check if the file was opened successfully
        cerr << "Error opening file: " << filename << endl; // Print an error message to the console
        return edges; // Return an empty vector if the file could not be opened
    }

    string line; // Variable to store each line read from the file
    size_t line_count = 0; // Counter for valid edges
    while (getline(file, line) && line_count < NO_OF_EDGES) { // Read up to 1000 valid lines
        istringstream iss(line); // Create a string stream from the line
        Edge edge; // Create an Edge object to store the edge data
        // Parse the line to extract the source, destination, and weights of the edge
        if (iss >> edge.src >> edge.dest >> edge.weight1 >> edge.weight2 >> edge.weight3) {
            edges.push_back(edge); // Add the Edge object to the vector of edges
            line_count++; // Increment counter only for valid edges
        }
    }
    file.close(); // Close the file after reading
    cerr << "Read " << edges.size() << " edges from " << filename << endl; // Log number of edges read
    return edges; // Return the vector containing all edges
}

int main() {

    string filename = "out.edgelist";
    vector<Edge> edges = readEdgeList(filename);

    // vector<Edge> edges = {
    //     // SCC 1 (Users 1-10): Tightly connected group of close friends
    //     {1, 2, 5, 0, 0},    {2, 1, 3, 0, 0},    // Bidirectional connection between 1 & 2
    //     {2, 3, 0, 4, 0},    {3, 2, 0, 0, 6},    // Bidirectional connection between 2 & 3
    //     {3, 4, 7, 0, 0},    {4, 3, 0, 2, 0},    // Bidirectional connection between 3 & 4
    //     {4, 5, 0, 0, 3},    {5, 4, 4, 0, 0},    // Bidirectional connection between 4 & 5
    //     {5, 6, 0, 5, 0},    {6, 5, 0, 0, 6},    // Bidirectional connection between 5 & 6
    //     {6, 7, 3, 0, 0},    {7, 6, 0, 4, 0},    // Bidirectional connection between 6 & 7
    //     {7, 8, 0, 0, 5},    {8, 7, 6, 0, 0},    // Bidirectional connection between 7 & 8
    //     {8, 9, 0, 7, 0},    {9, 8, 0, 0, 4},    // Bidirectional connection between 8 & 9
    //     {9, 10, 5, 0, 0},   {10, 9, 0, 6, 0},   // Bidirectional connection between 9 & 10
    //     {10, 1, 0, 0, 7},   {1, 10, 8, 0, 0},   // Bidirectional connection between 10 & 1
    //     {1, 5, 0, 5, 0},    {5, 1, 0, 0, 6},    // Extra connections within SCC 1
    //     {2, 7, 4, 0, 0},    {7, 2, 0, 3, 0},
    //     {3, 8, 0, 0, 5},    {8, 3, 6, 0, 0},
    //     {4, 9, 0, 7, 0},    {9, 4, 0, 0, 4},
        
    //     // SCC 2 (Users 11-20): Another tightly connected group of friends
    //     {11, 12, 6, 0, 0},  {12, 11, 0, 5, 0},  // Bidirectional connection between 11 & 12
    //     {12, 13, 0, 0, 7},  {13, 12, 8, 0, 0},  // Bidirectional connection between 12 & 13
    //     {13, 14, 0, 6, 0},  {14, 13, 0, 0, 5},  // Bidirectional connection between 13 & 14
    //     {14, 15, 7, 0, 0},  {15, 14, 0, 8, 0},  // Bidirectional connection between 14 & 15
    //     {15, 16, 0, 0, 6},  {16, 15, 5, 0, 0},  // Bidirectional connection between 15 & 16
    //     {16, 17, 0, 7, 0},  {17, 16, 0, 0, 8},  // Bidirectional connection between 16 & 17
    //     {17, 18, 6, 0, 0},  {18, 17, 0, 5, 0},  // Bidirectional connection between 17 & 18
    //     {18, 19, 0, 0, 7},  {19, 18, 8, 0, 0},  // Bidirectional connection between 18 & 19
    //     {19, 20, 0, 6, 0},  {20, 19, 0, 0, 5},  // Bidirectional connection between 19 & 20
    //     {20, 11, 7, 0, 0},  {11, 20, 0, 8, 0},  // Bidirectional connection between 20 & 11
    //     {11, 16, 0, 0, 6},  {16, 11, 5, 0, 0},  // Extra connections within SCC 2
    //     {12, 17, 0, 7, 0},  {17, 12, 0, 0, 8},
    //     {13, 18, 6, 0, 0},  {18, 13, 0, 5, 0},
    //     {14, 19, 0, 0, 7},  {19, 14, 8, 0, 0},
        
    //     // SCC 3 (Users 21-30): Work colleagues
    //     {21, 22, 3, 0, 0},  {22, 21, 0, 4, 0},  // Bidirectional connection between 21 & 22
    //     {22, 23, 0, 0, 5},  {23, 22, 6, 0, 0},  // Bidirectional connection between 22 & 23
    //     {23, 24, 0, 7, 0},  {24, 23, 0, 0, 3},  // Bidirectional connection between 23 & 24
    //     {24, 25, 4, 0, 0},  {25, 24, 0, 5, 0},  // Bidirectional connection between 24 & 25
    //     {25, 26, 0, 0, 6},  {26, 25, 7, 0, 0},  // Bidirectional connection between 25 & 26
    //     {26, 27, 0, 3, 0},  {27, 26, 0, 0, 4},  // Bidirectional connection between 26 & 27
    //     {27, 28, 5, 0, 0},  {28, 27, 0, 6, 0},  // Bidirectional connection between 27 & 28
    //     {28, 29, 0, 0, 7},  {29, 28, 3, 0, 0},  // Bidirectional connection between 28 & 29
    //     {29, 30, 0, 4, 0},  {30, 29, 0, 0, 5},  // Bidirectional connection between 29 & 30
    //     {30, 21, 6, 0, 0},  {21, 30, 0, 7, 0},  // Bidirectional connection between 30 & 21
        
    //     // SCC 4 (Users 31-40): Sports team members
    //     {31, 32, 4, 0, 0},  {32, 31, 0, 5, 0},  // Bidirectional connection between 31 & 32
    //     {32, 33, 0, 0, 6},  {33, 32, 7, 0, 0},  // Bidirectional connection between 32 & 33
    //     {33, 34, 0, 3, 0},  {34, 33, 0, 0, 4},  // Bidirectional connection between 33 & 34
    //     {34, 35, 5, 0, 0},  {35, 34, 0, 6, 0},  // Bidirectional connection between 34 & 35
    //     {35, 36, 0, 0, 7},  {36, 35, 3, 0, 0},  // Bidirectional connection between 35 & 36
    //     {36, 37, 0, 4, 0},  {37, 36, 0, 0, 5},  // Bidirectional connection between 36 & 37
    //     {37, 38, 6, 0, 0},  {38, 37, 0, 7, 0},  // Bidirectional connection between 37 & 38
    //     {38, 39, 0, 0, 3},  {39, 38, 4, 0, 0},  // Bidirectional connection between 38 & 39
    //     {39, 40, 0, 5, 0},  {40, 39, 0, 0, 6},  // Bidirectional connection between 39 & 40
    //     {40, 31, 7, 0, 0},  {31, 40, 0, 3, 0},  // Bidirectional connection between 40 & 31
        
    //     // SCC 5 (Users 41-50): Gaming community
    //     {41, 42, 5, 0, 0},  {42, 41, 0, 6, 0},  // Bidirectional connection between 41 & 42
    //     {42, 43, 0, 0, 7},  {43, 42, 3, 0, 0},  // Bidirectional connection between 42 & 43
    //     {43, 44, 0, 4, 0},  {44, 43, 0, 0, 5},  // Bidirectional connection between 43 & 44
    //     {44, 45, 6, 0, 0},  {45, 44, 0, 7, 0},  // Bidirectional connection between 44 & 45
    //     {45, 46, 0, 0, 3},  {46, 45, 4, 0, 0},  // Bidirectional connection between 45 & 46
    //     {46, 47, 0, 5, 0},  {47, 46, 0, 0, 6},  // Bidirectional connection between 46 & 47
    //     {47, 48, 7, 0, 0},  {48, 47, 0, 3, 0},  // Bidirectional connection between 47 & 48
    //     {48, 49, 0, 0, 4},  {49, 48, 5, 0, 0},  // Bidirectional connection between 48 & 49
    //     {49, 50, 0, 6, 0},  {50, 49, 0, 0, 7},  // Bidirectional connection between 49 & 50
    //     {50, 41, 3, 0, 0},  {41, 50, 0, 4, 0},  // Bidirectional connection between 50 & 41
        
    //     // SCC 6 (Users 51-60): Book club members
    //     {51, 52, 4, 0, 0},  {52, 51, 0, 5, 0},  // Bidirectional connection between 51 & 52
    //     {52, 53, 0, 0, 6},  {53, 52, 7, 0, 0},  // Bidirectional connection between 52 & 53
    //     {53, 54, 0, 3, 0},  {54, 53, 0, 0, 4},  // Bidirectional connection between 53 & 54
    //     {54, 55, 5, 0, 0},  {55, 54, 0, 6, 0},  // Bidirectional connection between 54 & 55
    //     {55, 56, 0, 0, 7},  {56, 55, 3, 0, 0},  // Bidirectional connection between 55 & 56
    //     {56, 57, 0, 4, 0},  {57, 56, 0, 0, 5},  // Bidirectional connection between 56 & 57
    //     {57, 58, 6, 0, 0},  {58, 57, 0, 7, 0},  // Bidirectional connection between 57 & 58
    //     {58, 59, 0, 0, 3},  {59, 58, 4, 0, 0},  // Bidirectional connection between 58 & 59
    //     {59, 60, 0, 5, 0},  {60, 59, 0, 0, 6},  // Bidirectional connection between 59 & 60
    //     {60, 51, 7, 0, 0},  {51, 60, 0, 3, 0},  // Bidirectional connection between 60 & 51
        
    //     // CAC 1 (Users 61-70): Celebrity (User 61) and followers (one-way connections)
    //     {61, 62, 8, 0, 0},  // Celebrity likes follower's post
    //     {61, 63, 0, 7, 0},  // Celebrity shares follower's post
    //     {61, 64, 0, 0, 9},  // Celebrity comments on follower's post
    //     {62, 61, 9, 8, 7},  // Follower interacts with celebrity
    //     {63, 61, 8, 7, 9},  // Follower interacts with celebrity
    //     {64, 61, 7, 9, 8},  // Follower interacts with celebrity
    //     {65, 61, 9, 7, 8},  // Follower interacts with celebrity
    //     {66, 61, 8, 9, 7},  // Follower interacts with celebrity
    //     {67, 61, 7, 8, 9},  // Follower interacts with celebrity
    //     {68, 61, 9, 9, 9},  // Follower interacts with celebrity
    //     {69, 61, 8, 8, 8},  // Follower interacts with celebrity
    //     {70, 61, 7, 7, 7},  // Follower interacts with celebrity
        
    //     // CAC 2 (Users 71-80): News outlet (User 71) and readers (one-way connections)
    //     {71, 72, 6, 8, 7},  // News outlet interacts with a reader
    //     {71, 73, 7, 6, 8},  // News outlet interacts with a reader
    //     {71, 74, 8, 7, 6},  // News outlet interacts with a reader
    //     {72, 71, 9, 0, 0},  // Reader likes news outlet's post
    //     {73, 71, 0, 9, 0},  // Reader shares news outlet's post
    //     {74, 71, 0, 0, 9},  // Reader comments on news outlet's post
    //     {75, 71, 8, 0, 0},  // Reader likes news outlet's post
    //     {76, 71, 0, 8, 0},  // Reader shares news outlet's post
    //     {77, 71, 0, 0, 8},  // Reader comments on news outlet's post
    //     {78, 71, 7, 0, 0},  // Reader likes news outlet's post
    //     {79, 71, 0, 7, 0},  // Reader shares news outlet's post
    //     {80, 71, 0, 0, 7},  // Reader comments on news outlet's post
        
    //     // CAC 3 (Users 81-90): Isolated users with minimal interaction
    //     {81, 82, 3, 0, 0},  // One-way interaction
    //     {83, 84, 0, 4, 0},  // One-way interaction
    //     {85, 86, 0, 0, 5},  // One-way interaction
    //     {87, 88, 6, 0, 0},  // One-way interaction
    //     {89, 90, 0, 7, 0},  // One-way interaction
        
    //     // Connections between SCCs (creating a directed acyclic graph of SCCs)
    //     {10, 11, 5, 5, 5},  // Connection from SCC 1 to SCC 2
    //     {20, 21, 6, 6, 6},  // Connection from SCC 2 to SCC 3
    //     {30, 31, 7, 7, 7},  // Connection from SCC 3 to SCC 4
    //     {40, 41, 8, 8, 8},  // Connection from SCC 4 to SCC 5
    //     {50, 51, 9, 9, 9},  // Connection from SCC 5 to SCC 6
        
    //     // Connections to CACs
    //     {10, 61, 4, 4, 4},  // Connection from SCC 1 to CAC 1 (Celebrity)
    //     {20, 71, 5, 5, 5},  // Connection from SCC 2 to CAC 2 (News outlet)
    //     {30, 81, 6, 6, 6},  // Connection from SCC 3 to CAC 3 (Isolated users)
        
    //     // Users 91-100: Completely isolated users (each forms a separate CAC)
    //     {91, 91, 1, 0, 0},  // Self-interaction
    //     {92, 92, 0, 1, 0},  // Self-interaction
    //     {93, 93, 0, 0, 1},  // Self-interaction
    //     {94, 94, 1, 1, 0},  // Self-interaction
    //     {95, 95, 0, 1, 1},  // Self-interaction
    //     {96, 96, 1, 0, 1},  // Self-interaction
    //     {97, 97, 1, 1, 1},  // Self-interaction
    //     {98, 98, 2, 2, 2},  // Self-interaction
    //     {99, 99, 3, 3, 3},  // Self-interaction
    //     {100, 100, 4, 4, 4} // Self-interaction
    // };
    
    // Create graph and run the algorithm
    Graph graph(edges);
    graph.sccCacPartitioning();
    graph.printResults();
    
    // --- Influence Power Measurement and Seed Selection ---
    vector<string> actions = {"like", "share", "comment"}; // Example actions (currently unused)
    vector<int> interest = {1, 2, 3}; // Example interest factors (currently unused)

    // 1. Calculate Influence Scores
    graph.parallelInfluencePowerMeasure(actions, edges);

    // 2. Select Seed Candidates (Optional intermediate step)
    graph.seedCandidatesSelection();

    // 3. Select Final Seeds (e.g., k=3)
    int k_seeds = 3;
    vector<int> final_seeds_indices = graph.seedSelection(k_seeds);

    // 4. Print Influence and Seed Results
    graph.printInfluenceResults();

    return 0;
}
