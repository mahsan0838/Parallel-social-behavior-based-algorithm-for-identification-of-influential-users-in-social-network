/*
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <tuple>

void add_edges_with_weight(const std::string& filename, double weight, std::ofstream& outFile, std::unordered_set<std::tuple<int, int, double>>& edgeSet, bool directed = true) {
    std::ifstream inFile(filename);
    int src, dst;
    while (inFile >> src >> dst) {
        if (filename.find("social") == std::string::npos) inFile >> std::ws; // ignore third column if not social
        std::tuple<int, int, double> edge = std::make_tuple(src, dst, weight);
        if (edgeSet.find(edge) == edgeSet.end()) {
            edgeSet.insert(edge);
            outFile << src << " " << dst << " " << weight << "\n";
            if (!directed) {
                // Write the reverse edge for undirected graphs
                edgeSet.insert(std::make_tuple(dst, src, weight));
                outFile << dst << " " << src << " " << weight << "\n";
            }
        }
    }
    inFile.close();
}

int main() {
    std::ofstream outFile("HiggsCombined.edgelist");
    std::unordered_set<std::tuple<int, int, double>> edgeSet;

    add_edges_with_weight("higgs-retweet_network.edgelist", 3.0, outFile, edgeSet);
    add_edges_with_weight("higgs-reply_network.edgelist", 2.0, outFile, edgeSet);
    add_edges_with_weight("higgs-mention_network.edgelist", 1.0, outFile, edgeSet);
    add_edges_with_weight("higgs-social_network.edgelist", 0.5, outFile, edgeSet, false); // undirected

    outFile.close();
    std::cout << "Combined edge list with weights written to HiggsCombined.edgelist\n";
    return 0;
}

*/
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <tuple>
#include <string>
#include <functional>

// Custom hash function for tuple<int, int, double>
struct TupleHash {
    std::size_t operator()(const std::tuple<int, int, double>& t) const {
        auto h1 = std::hash<int>{}(std::get<0>(t));
        auto h2 = std::hash<int>{}(std::get<1>(t));
        auto h3 = std::hash<double>{}(std::get<2>(t));
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

// Function to add edges from a file with a specific weight
void add_edges_with_weight(
    const std::string& filename,
    double weight,
    std::ofstream& outFile,
    std::unordered_set<std::tuple<int, int, double>, TupleHash>& edgeSet,
    bool directed = true
) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open " << filename << "\n";
        return;
    }

    int src, dst;
    while (inFile >> src >> dst) {
        if (filename.find("social") == std::string::npos) {
            // ignore third column if it exists (retweet, mention, reply files)
            inFile >> std::ws;
        }

        if (!directed && src > dst) std::swap(src, dst); // Normalize undirected edges

        std::tuple<int, int, double> edge = std::make_tuple(src, dst, weight);
        if (edgeSet.find(edge) == edgeSet.end()) {
            edgeSet.insert(edge);
            outFile << src << " " << dst << " " << weight << "\n";
            if (!directed) {
                // Add reverse edge too
                edgeSet.insert(std::make_tuple(dst, src, weight));
                outFile << dst << " " << src << " " << weight << "\n";
            }
        }
    }

    inFile.close();
}

int main() {
    std::ofstream outFile("HiggsCombined.edgelist");
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open HiggsCombined.edgelist for writing\n";
        return 1;
    }

    std::unordered_set<std::tuple<int, int, double>, TupleHash> edgeSet;

    // Add edges from multiple files with corresponding weights
    add_edges_with_weight("higgs-retweet_network.edgelist", 3.0, outFile, edgeSet);
    add_edges_with_weight("higgs-reply_network.edgelist",   2.0, outFile, edgeSet);
    add_edges_with_weight("higgs-mention_network.edgelist", 1.0, outFile, edgeSet);
    add_edges_with_weight("higgs-social_network.edgelist",  0.5, outFile, edgeSet, false); // undirected

    outFile.close();
    std::cout << "Combined edge list with weights written to HiggsCombined.edgelist\n";
    return 0;
}


