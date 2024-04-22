#pragma once
#define UTILS_TYPES

#include <chrono>
#include <climits>
#include <cstdint>
#include <functional>
#include <stdlib.h>

#define NOT_EXIST UINT_MAX
#define UNMATCHED UINT_MAX

#define Get_Time() std::chrono::high_resolution_clock::now()
#define Duration(start) std::chrono::duration_cast<\
    std::chrono::microseconds>(Get_Time() - start).count()/(float)1000
#define Print_Time(str, start) std::cout << str << Duration(start) << \
    "ms" << std::endl

typedef unsigned int uint;
typedef unsigned int ui;

typedef unsigned int VertexID;
typedef unsigned int LabelID;

struct InsertUnit {
    char type; 
    bool is_add;
    uint id1;  
    uint id2;  
    uint label; 
    InsertUnit(char type_arg, bool is_add_arg, uint id1_arg, uint id2_arg, uint label_arg)
        : type(type_arg), is_add(is_add_arg), id1(id1_arg), id2(id2_arg), label(label_arg) {}
};

template <typename T>
inline void hash_combine(std::size_t& seed, const T& val) {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <typename T> inline void hash_val(std::size_t& seed, const T& val) {
    hash_combine(seed, val);
}
template <typename T, typename... Types>
inline void hash_val(std::size_t& seed, const T& val, const Types &... args) {
    hash_combine(seed, val);
    hash_val(seed, args...);
}

template <typename... Types>
inline std::size_t hash_val(const Types &... args) {
    std::size_t seed = 0;
    hash_val(seed, args...);
    return seed;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        return hash_val(p.first, p.second);
    }
};

enum MatchingIndexType {
    VertexCentric = 0,
    EdgeCentric = 1
};

class TreeNode {
public:
    VertexID id_;
    VertexID parent_;
    int level_;
    int under_level_count_;
    int children_count_;
    int bn_count_;
    int fn_count_;
    VertexID* under_level_;
    VertexID* children_;
    VertexID* bn_;
    VertexID* fn_;
    size_t estimated_embeddings_num_;
public:
    TreeNode() {
        id_ = 0;
        under_level_ = NULL;
        bn_ = NULL;
        fn_ = NULL;
        children_ = NULL;
        parent_ = 0;
        level_ = 0;
        under_level_count_ = 0;
        children_count_ = 0;
        bn_count_ = 0;
        fn_count_ = 0;
        estimated_embeddings_num_ = 0;
    }

    ~TreeNode() {
        delete[] under_level_;
        delete[] bn_;
        delete[] fn_;
        delete[] children_;
    }

    void initialize(const int size) {
        under_level_ = new VertexID[size];
        bn_ = new VertexID[size];
        fn_ = new VertexID[size];
        children_ = new VertexID[size];
    }
};

class Edges {
public:
    int* offset_;
    int* edge_;
    int vertex_count_;
    int edge_count_;
    int max_degree_;
public:
    Edges() {
        offset_ = NULL;
        edge_ = NULL;
        vertex_count_ = 0;
        edge_count_ = 0;
        max_degree_ = 0;
    }

    ~Edges() {
        delete[] offset_;
        delete[] edge_;
    }
};