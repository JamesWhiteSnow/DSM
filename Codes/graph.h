#pragma once
#include <queue>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "types.h"
#include "utils.h"
#include "config.h"

class Dynamic_Graph
{
protected:
    uint edge_count_;
    uint vlabel_count_;
    uint elabel_count_;
    std::vector<std::vector<uint>> neighbors_;
    std::vector<std::vector<uint>> elabels_;

public:
    std::queue<InsertUnit> updates_;
    std::vector<uint> vlabels_;
    std::vector < std::vector<uint>> edges_;

public:
    Dynamic_Graph();

    virtual uint NumVertices() const { return vlabels_.size(); }
    virtual uint NumEdges() const { return edge_count_; }
    uint NumVLabels() const { return vlabel_count_; }
    uint NumELabels() const { return elabel_count_; }
    uint GetDiameter() const;

    void AddVertex(uint id, uint label);
    void RemoveVertex(uint id);
    void AddEdge(uint v1, uint v2, uint label);
    void RemoveEdge(uint v1, uint v2);

    uint GetVertexLabel(uint u) const;
    const std::vector<uint>& GetNeighbors(uint v) const;
    const std::vector<uint>& GetNeighborLabels(uint v) const;
    uint GetDegree(uint v) const;
    std::tuple<uint, uint, uint> GetEdgeLabel(uint v1, uint v2) const;

    void LoadFromFile(const std::string &path);
    void LoadUpdateStream(const std::string &path);
    void PrintMetaData() const;
};


class Static_Graph {
private:
    bool enable_label_offset_;

    ui vertices_count_;
    ui edges_count_;
    ui labels_count_;
    ui max_degree_;
    ui max_label_frequency_;

    ui* offsets_;
    VertexID* neighbors_;
    LabelID* labels_;
    ui* reverse_index_offsets_;
    ui* reverse_index_;

    int* core_table_;
    ui core_length_;

    std::unordered_map<LabelID, ui> labels_frequency_;

#if OPTIMIZED_LABELED_GRAPH == 1
    ui* labels_offsets_;
    std::unordered_map<LabelID, ui>* nlf_;
#endif

public:
    void BuildReverseIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
    void BuildNLF();
    void BuildLabelOffset();
#endif

public:
    Static_Graph();

    Static_Graph(const bool enable_label_offset) {
        enable_label_offset_ = enable_label_offset;

        vertices_count_ = 0;
        edges_count_ = 0;
        labels_count_ = 0;
        max_degree_ = 0;
        max_label_frequency_ = 0;
        core_length_ = 0;

        offsets_ = NULL;
        neighbors_ = NULL;
        labels_ = NULL;
        reverse_index_offsets_ = NULL;
        reverse_index_ = NULL;
        core_table_ = NULL;
        labels_frequency_.clear();

#if OPTIMIZED_LABELED_GRAPH == 1
        labels_offsets_ = NULL;
        nlf_ = NULL;
#endif
    }

    ~Static_Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] labels_;
        delete[] reverse_index_offsets_;
        delete[] reverse_index_;
        delete[] core_table_;
#if OPTIMIZED_LABELED_GRAPH == 1
        delete[] labels_offsets_;
        delete[] nlf_;
#endif
    }

public:
    void loadGraphFromFile(const std::string& file_path);
    void loadGraphFromFileCompressed(const std::string& degree_path, const std::string& edge_path,
        const std::string& label_path);
    void storeComparessedGraph(const std::string& degree_path, const std::string& edge_path,
        const std::string& label_path);
    void printGraphMetaData();
    void loadGraphFromDynamicGraph(Dynamic_Graph dynamic_graph);
public:
    const ui getLabelsCount() const {
        return labels_count_;
    }

    const ui getVerticesCount() const {
        return vertices_count_;
    }

    const ui getEdgesCount() const {
        return edges_count_;
    }

    const ui getGraphMaxDegree() const {
        return max_degree_;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const ui getVertexDegree(const VertexID id) const {
        return offsets_[id + 1] - offsets_[id];
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }

    const ui getCoreValue(const VertexID id) const {
        return core_table_[id];
    }

    const ui get2CoreSize() const {
        return core_length_;
    }
    const LabelID getVertexLabel(const VertexID id) const {
        return labels_[id];
    }

    const ui* getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }


    const ui* getVerticesByLabel(const LabelID id, ui& count) const {
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

#if OPTIMIZED_LABELED_GRAPH == 1
    const ui* getNeighborsByLabel(const VertexID id, const LabelID label, ui& count) const {
        ui offset = id * labels_count_ + label;
        count = labels_offsets_[offset + 1] - labels_offsets_[offset];
        return neighbors_ + labels_offsets_[offset];
    }

    const std::unordered_map<LabelID, ui>* getVertexNLF(const VertexID id) const {
        return nlf_ + id;
    }

    bool checkEdgeExistence(const VertexID u, const VertexID v, const LabelID u_label) const {
        ui count = 0;
        const VertexID* neighbors = getNeighborsByLabel(v, u_label, count);
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }
#endif

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors = getVertexNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    void buildCoreTable();
};