// g++ -std=c++11   -I./metis-5.1.0/include   -I./metis-5.1.0/build/Linux-x86_64/include   partition_graph2.cpp   ./metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a   -o pg

// ./pg higgs.graph 3

/*
#include <iostream>
#include <fstream>
#include <vector>
#include <metis.h>

void writePartitions(const std::vector<idx_t>& part, idx_t nparts) {
    std::vector<std::ofstream> files(nparts);
    
    // Open one output file for each partition
    for (idx_t i = 0; i < nparts; ++i) {
        files[i].open("partition_" + std::to_string(i) + ".txt");
        if (!files[i].is_open()) {
            std::cerr << "Error opening output file for partition " << i << "\n";
            exit(1);
        }
    }

    // Write each node to the file of its assigned partition
    for (size_t node = 0; node < part.size(); ++node) {
        idx_t p = part[node];
        files[p] << node << "\n";
    }

    // Close files
    for (idx_t i = 0; i < nparts; ++i) {
        files[i].close();
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./partition_graph <graph_file> <nparts>\n";
        return 1;
    }

    std::string filename = argv[1];
    idx_t nparts = std::stoi(argv[2]);

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Unable to open file: " << filename << "\n";
        return 1;
    }

    idx_t nvtxs, nedges;
    infile >> nvtxs >> nedges;

    std::vector<idx_t> xadj(nvtxs + 1);
    std::vector<idx_t> adjncy(2 * nedges); // Undirected graph: 2 edges per connection

    for (idx_t i = 0; i <= nvtxs; ++i) {
        infile >> xadj[i];
    }

    for (idx_t i = 0; i < 2 * nedges; ++i) {
        infile >> adjncy[i];
    }

    infile.close();

    std::vector<idx_t> part(nvtxs); // Output partition vector
    idx_t objval;

    int result = METIS_PartGraphKway(
        &nvtxs,
        nullptr,          // ncon: number of balancing constraints (default 1)
        xadj.data(),
        adjncy.data(),
        nullptr, nullptr, // vwgt and vsize (optional vertex weights)
        nullptr,          // adjwgt (optional edge weights)
        &nparts,
        nullptr, nullptr, // tpwgts and ubvec (balancing)
        nullptr,          // options (use defaults)
        &objval,
        part.data()
    );

    if (result != METIS_OK) {
        std::cerr << "METIS_PartGraphKway failed with code " << result << "\n";
        return 1;
    }

    std::cout << "Partitioning successful. Writing to files...\n";
    writePartitions(part, nparts);

    std::cout << "Done. Objective value: " << objval << "\n";
    return 0;
}

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <metis.h>

void writePartitions(const std::vector<idx_t>& part, idx_t nparts) {
    std::vector<std::ofstream> files(nparts);

    for (idx_t i = 0; i < nparts; ++i) {
        files[i].open("partition_" + std::to_string(i) + ".txt");
        if (!files[i]) {
            std::cerr << "Error opening file for partition " << i << "\n";
            exit(1);
        }
    }

    for (size_t i = 0; i < part.size(); ++i) {
        files[part[i]] << i << "\n";
    }

    for (auto& f : files) {
        f.close();
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./pg <graph_file> <nparts>\n";
        return 1;
    }

    std::string filename = argv[1];
    idx_t nparts = std::stoi(argv[2]);

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Unable to open file: " << filename << "\n";
        return 1;
    }

    idx_t nvtxs, nedges;
    infile >> nvtxs >> nedges;
    infile.ignore(); // skip to the next line

    std::vector<idx_t> xadj(nvtxs + 1, 0);
    std::vector<idx_t> adjncy;
    std::string line;
    idx_t edge_counter = 0;

    for (idx_t i = 0; i < nvtxs; ++i) {
        if (!std::getline(infile, line)) {
            std::cerr << "Error reading adjacency list for vertex " << i << "\n";
            return 1;
        }

        std::istringstream iss(line);
        idx_t neighbor;
        while (iss >> neighbor) {
            adjncy.push_back(neighbor - 1); // convert to 0-based indexing
            edge_counter++;
        }
        xadj[i + 1] = adjncy.size(); // running total
    }

    infile.close();

    if (adjncy.size() != static_cast<size_t>(2 * nedges)) {
        std::cerr << "Warning: Expected " << 2 * nedges << " edges, but got " << adjncy.size() << ".\n";
    }

    std::vector<idx_t> part(nvtxs);
    idx_t objval;
    idx_t ncon = 1; // REQUIRED for METIS

    int status = METIS_PartGraphKway(
        &nvtxs,
        &ncon,
        xadj.data(),
        adjncy.data(),
        nullptr, nullptr, nullptr,
        &nparts,
        nullptr, nullptr, nullptr,
        &objval,
        part.data()
    );

    if (status != METIS_OK) {
        std::cerr << "METIS_PartGraphKway failed with code " << status << "\n";
        return 1;
    }

    std::cout << "Partitioning successful. Objective value: " << objval << "\n";
    writePartitions(part, nparts);
    std::cout << "Partition files written successfully.\n";

    return 0;
}

