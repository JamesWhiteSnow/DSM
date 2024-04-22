#include <vector>

#include "types.h"
#include "graph.h"
#include "matching.h"

class Graphflow : public matching
{
public:
    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> order_csrs_;
    std::vector<std::vector<uint>> order_offs_;

    Graphflow(Dynamic_Graph& query_graph, Dynamic_Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo);
    ~Graphflow() override {};

    void Preprocessing() override;
    void InitialMatching() override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

    void GenerateMatchingOrder();
    void FindMatches(uint order_index, uint depth,
        std::vector<uint> m, size_t& num_results);
};
