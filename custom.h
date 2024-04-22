#include "graph.h"
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <vector>
#include <chrono>
#include <random>
#include <numeric>
#include <chrono>
#include <limits>
#include <omp.h>
#include <unordered_map>
#include <iomanip>
#include "graphflow.h"

#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)
#define NANOSECTOMSEC(elapsed_time) ((elapsed_time)/(double)1000000)

using namespace std;
using namespace chrono;

ui MAX_LIMIT = UINT_MAX;

ui dimension_x = 2;
ui dimension_y = dimension_x;
ui dimension_z = 2 * dimension_x;
ui dimension_o = 2 * dimension_x;
ui dimension_o_prime = 2 * dimension_x;
ui dimension_mbr = 2 * dimension_o_prime;

double alpha_value = 0.01;
double beta_value = 10;

ui cell_num = 5;

ui degree_group_num = 3;

class Data_Vertex {
public:
	ui degree;
	ui label;

	vector<ui> storage_pos;

	vector<double> keys;

	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> o;
	vector<double> o_prime;

	vector<vector<double>> neighbor_x;
	vector<vector<double>> neighbor_x_sort;

	vector<vector<double>> MBRs;

	unordered_map<ui, ui> neighbors;
};

class Query_Vertex {
public:
	ui label;
	ui degree;
	double key;

	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> o;
	vector<double> o_prime;

	vector<ui> candidates;

	unordered_map<ui, ui> neighbors;
	unordered_map<ui, ui> candidates_hash;
};

class Point {
public:
	ui vertex_id;
	double key;

	vector<double> upper_bound;
	vector<vector<double>> MBRs;

	bool operator<(const Point& other) const {
		return key > other.key;
	}
};

class Cell {
public:
	double key;

	vector<double> bounds;
	vector<Point> points;

	bool operator<(const Cell& other) const {
		return key > other.key;
	}
};

class Grid_Index {
public:
	ui deg_lb;
	ui deg_up;

	vector<Cell> cell_list;
};

class DAS3 {
public:
	vector<ui> degree_index_map;
	vector<Grid_Index> index_list;
};

vector<double> gen_VDE_x(const ui& vertex_label, const ui& dimension) {
	ui seed = vertex_label;
	mt19937 psr_gen(seed);
	uniform_real_distribution<double> dis(0.0, 1.0);

	vector<double> embedding_x(dimension);
	for (ui i = 0; i < dimension; i++) {
		embedding_x[i] = dis(psr_gen);
	}
	return embedding_x;
}

vector<double> gen_VDE_z(const ui& vertex_label, const ui& dimension) {
	ui seed = vertex_label*2;
	mt19937 psr_gen(seed);
	uniform_real_distribution<double> dis(0.0, 1.0);

	vector<double> embedding_z(dimension);
	for (ui i = 0; i < dimension; i++) {
		embedding_z[i] = dis(psr_gen);
	}

	double sum = std::accumulate(embedding_z.begin(), embedding_z.end(), 0.0);
	for (ui i = 0; i < dimension; i++) {
		embedding_z[i] = embedding_z[i] / sum;
	}

	return embedding_z;
}

vector<double> gen_VDE_z_zipf(const ui& vertex_label, const ui& dimension, const vector<double>& zipf) {
	ui seed = vertex_label;
	mt19937 psr_gen(seed);
	uniform_int_distribution<ui> dis(0,zipf.size());

	vector<double> embedding_z(dimension);
	for (ui i = 0; i < dimension; i++) {
		embedding_z[i] = zipf[dis(psr_gen)];
	}

	return embedding_z;
}


vector<Data_Vertex> gen_VDE_data(const Dynamic_Graph& initial_graph,const vector<double>& zipf) {
	vector<Data_Vertex> data_vertices(initial_graph.NumVertices());
	//embedding_x and embedding_z
	for (ui i = 0; i < initial_graph.NumVertices(); i++) {
		data_vertices[i].label=initial_graph.GetVertexLabel(i);
		data_vertices[i].degree = initial_graph.GetDegree(i);
		data_vertices[i].x = gen_VDE_z_zipf(initial_graph.GetVertexLabel(i), dimension_x,zipf);
		data_vertices[i].z = gen_VDE_z(initial_graph.GetVertexLabel(i), dimension_z);
	}
	for (ui i = 0; i < initial_graph.NumVertices(); i++) {
		//embedding_y
		vector<ui> neighbors = initial_graph.GetNeighbors(i);
		data_vertices[i].neighbor_x.resize(neighbors.size());
		data_vertices[i].y = vector<double>(dimension_y, 0);
		for (ui j = 0; j < neighbors.size(); j++) {
			data_vertices[i].neighbors.insert(make_pair(neighbors[j], neighbors[j]));
			data_vertices[i].neighbor_x[j] = data_vertices[neighbors[j]].x;
			for (ui k = 0; k < dimension_y; k++) {
				data_vertices[i].y[k] += data_vertices[i].neighbor_x[j][k];
			}
		}
		//embedding_o
		data_vertices[i].o = vector<double>();
		for (ui j = 0; j < dimension_x; j++) {
			data_vertices[i].o.push_back(data_vertices[i].x[j]);
		}
		for (ui j = 0; j < dimension_y; j++) {
			data_vertices[i].o.push_back(data_vertices[i].y[j]);
		}
		//embedding_o_prime
		data_vertices[i].o_prime = vector<double>();
		for (ui j = 0; j < dimension_o; j++) {
			data_vertices[i].o_prime.push_back(alpha_value * data_vertices[i].o[j] + beta_value * data_vertices[i].z[j]);
		}
	}

	//sort neighbors' embeddings
	for (ui i = 0; i < initial_graph.NumVertices(); i++) {
		if (data_vertices[i].neighbor_x.size() != 0) {
			data_vertices[i].neighbor_x_sort.resize(data_vertices[i].neighbor_x.size());
			for (ui col = 0; col < dimension_x; col++) {
				vector<double> column;
				for (ui row = 0; row < data_vertices[i].neighbor_x.size(); row++) {
					column.push_back(data_vertices[i].neighbor_x[row][col]);
				}
				sort(column.begin(), column.end());
				for (ui row = 0; row < data_vertices[i].neighbor_x.size(); row++) {
					data_vertices[i].neighbor_x_sort[row].push_back(column[row]);
				}
			}
		}
	}

	//calculate MBRs
	for (ui vid = 0; vid < initial_graph.NumVertices(); vid++) {
		data_vertices[vid].MBRs.resize(initial_graph.GetDegree(vid) + 1);
		for (ui deg = 0; deg <= initial_graph.GetDegree(vid); deg++) {
			vector<double> mbr_min = data_vertices[vid].x;
			vector<double> mbr_max = data_vertices[vid].x;
			for (ui dim = 0; dim < dimension_x; dim++) {
				double min_value = 0;
				double max_value = 0;
				for (ui i = 0; i < deg; i++) {
					min_value += data_vertices[vid].neighbor_x_sort[i][dim];
				}
				for (ui i = data_vertices[vid].neighbor_x_sort.size() - deg; i < data_vertices[vid].neighbor_x_sort.size(); i++) {
					max_value += data_vertices[vid].neighbor_x_sort[i][dim];
				}
				mbr_min.push_back(min_value);
				mbr_max.push_back(max_value);
			}
			for (ui i = 0; i < dimension_o_prime; i++) {
				data_vertices[vid].MBRs[deg].push_back(alpha_value * mbr_min[i] + beta_value * data_vertices[vid].z[i]);
				data_vertices[vid].MBRs[deg].push_back(alpha_value * mbr_max[i] + beta_value * data_vertices[vid].z[i]);
			}
		}
	}
	return data_vertices;
}

void edge_insertion_update(vector<Data_Vertex>& data_vertices, DAS3& das3, ui v1, ui v2, double& update_cost) {
	auto start = chrono::high_resolution_clock::now();
	//update vertex information
	data_vertices[v1].degree += 1;
	data_vertices[v2].degree += 1;

	//update embedding
	data_vertices[v1].neighbor_x.push_back(data_vertices[v2].x);
	data_vertices[v2].neighbor_x.push_back(data_vertices[v1].x);
	for (ui i = 0; i < dimension_x; i++) {
		data_vertices[v1].y[i] += data_vertices[v2].x[i];
		data_vertices[v2].y[i] += data_vertices[v1].x[i];
	}
	for (ui i = dimension_x; i < dimension_o; i++) {
		data_vertices[v1].o[i] = data_vertices[v1].y[i];
		data_vertices[v2].o[i] = data_vertices[v2].y[i];
	}
	for (ui i = dimension_x; i < dimension_o; i++) {
		data_vertices[v1].o_prime[i] = (alpha_value * data_vertices[v1].o[i] + beta_value * data_vertices[v1].z[i]);
		data_vertices[v2].o_prime[i] = (alpha_value * data_vertices[v2].o[i] + beta_value * data_vertices[v2].z[i]);
	}

	//sort neighbor's embeddings
	data_vertices[v1].neighbor_x_sort.resize(data_vertices[v1].neighbor_x.size());
	for (ui col = 0; col < dimension_x; col++) {
		vector<double> column;
		for (ui row = 0; row < data_vertices[v1].neighbor_x.size(); row++) {
			column.push_back(data_vertices[v1].neighbor_x[row][col]);
		}
		sort(column.begin(), column.end());
		for (ui row = 0; row < data_vertices[v1].neighbor_x.size(); row++) {
			data_vertices[v1].neighbor_x_sort[row].push_back(column[row]);
		}
	}
	data_vertices[v2].neighbor_x_sort.resize(data_vertices[v2].neighbor_x.size());
	for (ui col = 0; col < dimension_x; col++) {
		vector<double> column;
		for (ui row = 0; row < data_vertices[v2].neighbor_x.size(); row++) {
			column.push_back(data_vertices[v2].neighbor_x[row][col]);
		}
		sort(column.begin(), column.end());
		for (ui row = 0; row < data_vertices[v2].neighbor_x.size(); row++) {
			data_vertices[v2].neighbor_x_sort[row].push_back(column[row]);
		}
	}

	//Update MBRs
	data_vertices[v1].MBRs.resize(data_vertices[v1].degree + 1);
	for (ui deg = 0; deg <= data_vertices[v1].degree; deg++) {
		vector<double> mbr_min = data_vertices[v1].x;
		vector<double> mbr_max = data_vertices[v1].x;
		for (ui dim = 0; dim < dimension_x; dim++) {
			double min_value = 0;
			double max_value = 0;
			for (ui i = 0; i < deg; i++) {
				min_value += data_vertices[v1].neighbor_x_sort[i][dim];
			}
			for (ui i = data_vertices[v1].neighbor_x_sort.size() - deg; i < data_vertices[v1].neighbor_x_sort.size(); i++) {
				max_value += data_vertices[v1].neighbor_x_sort[i][dim];
			}
			mbr_min.push_back(min_value);
			mbr_max.push_back(max_value);
		}
		for (ui i = 0; i < dimension_o_prime; i++) {
			data_vertices[v1].MBRs[deg].push_back(alpha_value * mbr_min[i] + beta_value * data_vertices[v1].z[i]);
			data_vertices[v1].MBRs[deg].push_back(alpha_value * mbr_max[i] + beta_value * data_vertices[v1].z[i]);
		}
	}
	data_vertices[v2].MBRs.resize(data_vertices[v2].degree + 1);
	for (ui deg = 0; deg <= data_vertices[v2].degree; deg++) {
		vector<double> mbr_min = data_vertices[v2].x;
		vector<double> mbr_max = data_vertices[v2].x;
		for (ui dim = 0; dim < dimension_x; dim++) {
			double min_value = 0;
			double max_value = 0;
			for (ui i = 0; i < deg; i++) {
				min_value += data_vertices[v2].neighbor_x_sort[i][dim];
			}
			for (ui i = data_vertices[v2].neighbor_x_sort.size() - deg; i < data_vertices[v2].neighbor_x_sort.size(); i++) {
				max_value += data_vertices[v2].neighbor_x_sort[i][dim];
			}
			mbr_min.push_back(min_value);
			mbr_max.push_back(max_value);
		}
		for (ui i = 0; i < dimension_o_prime; i++) {
			data_vertices[v2].MBRs[deg].push_back(alpha_value * mbr_min[i] + beta_value * data_vertices[v2].z[i]);
			data_vertices[v2].MBRs[deg].push_back(alpha_value * mbr_max[i] + beta_value * data_vertices[v2].z[i]);
		}
	}

	auto end = chrono::high_resolution_clock::now();
	update_cost += double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

void edge_insertion_refine(Graphflow& mm, const vector<Query_Vertex>& query_vertices, vector<Data_Vertex>& data_vertices, ui v1, ui v2, ui label, double& continuous_cost) {
	vector<uint> m(mm.query_.NumVertices(), UNMATCHED);
	size_t num_results = 0ul;

	auto start = chrono::high_resolution_clock::now();
	for (uint i = 0; i < mm.query_.NumEdges(); i++)
	{
		uint u1 = mm.order_vs_[i][0], u2 = mm.order_vs_[i][1];
		auto temp_q_labels = mm.query_.GetEdgeLabel(u1, u2);

		if (query_vertices[u1].degree <= data_vertices[v1].degree
			&& query_vertices[u2].degree <= data_vertices[v2].degree
			&& query_vertices[u1].label==data_vertices[v1].label
			&& query_vertices[u2].label==data_vertices[v2].label) {
			ui j;
			for (j = 0; j < dimension_o; j++) {
				if (!(query_vertices[u1].o_prime[j] >= data_vertices[v1].MBRs[query_vertices[u1].degree][2 * j]
					&& query_vertices[u1].o_prime[j] <= data_vertices[v1].MBRs[query_vertices[u1].degree][2 * j + 1]
					&& query_vertices[u2].o_prime[j] >= data_vertices[v2].MBRs[query_vertices[u2].degree][2 * j]
					&& query_vertices[u2].o_prime[j] <= data_vertices[v2].MBRs[query_vertices[u2].degree][2 * j + 1])) {
					break;
				}
			}
			if (j == dimension_o)
			{
				m[u1] = v1;
				m[u2] = v2;
				mm.visited_[v1] = true;
				mm.visited_[v2] = true;

				mm.FindMatches(i, 2, m, num_results);

				mm.visited_[v1] = false;
				mm.visited_[v2] = false;
				m[u1] = UNMATCHED;
				m[u2] = UNMATCHED;

				if (num_results >= mm.max_num_results_) goto END_ENUMERATION;
			}
		}

		if (query_vertices[u1].degree <= data_vertices[v2].degree
			&& query_vertices[u2].degree <= data_vertices[v1].degree
			&& query_vertices[u1].label==data_vertices[v2].label
			&& query_vertices[u2].label==data_vertices[v1].label) {
			ui j;
			for (j = 0; j < dimension_o; j++) {
				if (!(query_vertices[u1].o_prime[j] >= data_vertices[v2].MBRs[query_vertices[u1].degree][2 * j]
					&& query_vertices[u1].o_prime[j] <= data_vertices[v2].MBRs[query_vertices[u1].degree][2 * j + 1]
					&& query_vertices[u2].o_prime[j] >= data_vertices[v1].MBRs[query_vertices[u2].degree][2 * j]
					&& query_vertices[u2].o_prime[j] <= data_vertices[v1].MBRs[query_vertices[u2].degree][2 * j + 1])) {
					break;
				}
			}
			if (j == dimension_o)
			{
				m[u1] = v2;
				m[u2] = v1;
				mm.visited_[v2] = true;
				mm.visited_[v1] = true;

				mm.FindMatches(i, 2, m, num_results);

				mm.visited_[v2] = false;
				mm.visited_[v1] = false;
				m[u1] = UNMATCHED;
				m[u2] = UNMATCHED;

				if (num_results >= mm.max_num_results_) goto END_ENUMERATION;
			}
		}
	}
END_ENUMERATION:
	auto end = chrono::high_resolution_clock::now();
	mm.num_positive_results_ += num_results;
	continuous_cost += double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

void edge_deletion_update(vector<Data_Vertex>& data_vertices, DAS3& das3, ui v1, ui v2, double& update_cost) {
	auto start = chrono::high_resolution_clock::now();
	data_vertices[v1].degree -= 1;
	data_vertices[v2].degree -= 1;
	data_vertices[v1].neighbors.erase(v2);
	data_vertices[v2].neighbors.erase(v1);

	vector<vector<double>>::iterator it = find(data_vertices[v1].neighbor_x.begin(), data_vertices[v1].neighbor_x.end(), data_vertices[v2].x);
	data_vertices[v1].neighbor_x.erase(it);
	it = find(data_vertices[v2].neighbor_x.begin(), data_vertices[v2].neighbor_x.end(), data_vertices[v1].x);
	data_vertices[v2].neighbor_x.erase(it);
	for (ui i = 0; i < dimension_x; i++) {
		data_vertices[v1].y[i] -= data_vertices[v2].x[i];
		data_vertices[v2].y[i] -= data_vertices[v1].x[i];
	}
	for (ui i = dimension_x; i < dimension_o; i++) {
		data_vertices[v1].o[i] = data_vertices[v1].y[i];
		data_vertices[v2].o[i] = data_vertices[v2].y[i];
	}
	for (ui i = dimension_x; i < dimension_o; i++) {
		data_vertices[v1].o_prime[i] = (alpha_value * data_vertices[v1].o[i] + beta_value * data_vertices[v1].z[i]);
		data_vertices[v2].o_prime[i] = (alpha_value * data_vertices[v2].o[i] + beta_value * data_vertices[v2].z[i]);
	}

	data_vertices[v1].neighbor_x_sort.resize(data_vertices[v1].neighbor_x.size());
	for (ui col = 0; col < dimension_x; col++) {
		vector<double> column;
		for (ui row = 0; row < data_vertices[v1].neighbor_x.size(); row++) {
			column.push_back(data_vertices[v1].neighbor_x[row][col]);
		}
		sort(column.begin(), column.end());
		for (ui row = 0; row < data_vertices[v1].neighbor_x.size(); row++) {
			data_vertices[v1].neighbor_x_sort[row].push_back(column[row]);
		}
	}
	data_vertices[v2].neighbor_x_sort.resize(data_vertices[v2].neighbor_x.size());
	for (ui col = 0; col < dimension_x; col++) {
		vector<double> column;
		for (ui row = 0; row < data_vertices[v2].neighbor_x.size(); row++) {
			column.push_back(data_vertices[v2].neighbor_x[row][col]);
		}
		sort(column.begin(), column.end());
		for (ui row = 0; row < data_vertices[v2].neighbor_x.size(); row++) {
			data_vertices[v2].neighbor_x_sort[row].push_back(column[row]);
		}
	}

	//Update MBRs
	data_vertices[v1].MBRs.resize(data_vertices[v1].degree + 1);
	for (ui deg = 0; deg <= data_vertices[v1].degree; deg++) {
		vector<double> mbr_min = data_vertices[v1].x;
		vector<double> mbr_max = data_vertices[v1].x;
		for (ui dim = 0; dim < dimension_x; dim++) {
			double min_value = 0;
			double max_value = 0;
			for (ui i = 0; i < deg; i++) {
				min_value += data_vertices[v1].neighbor_x_sort[i][dim];
			}
			for (ui i = data_vertices[v1].neighbor_x_sort.size() - deg; i < data_vertices[v1].neighbor_x_sort.size(); i++) {
				max_value += data_vertices[v1].neighbor_x_sort[i][dim];
			}
			mbr_min.push_back(min_value);
			mbr_max.push_back(max_value);
		}
		for (ui i = 0; i < dimension_o_prime; i++) {
			data_vertices[v1].MBRs[deg].push_back(alpha_value * mbr_min[i] + beta_value * data_vertices[v1].z[i]);
			data_vertices[v1].MBRs[deg].push_back(alpha_value * mbr_max[i] + beta_value * data_vertices[v1].z[i]);
		}
	}
	data_vertices[v2].MBRs.resize(data_vertices[v2].degree + 1);
	for (ui deg = 0; deg <= data_vertices[v2].degree; deg++) {
		vector<double> mbr_min = data_vertices[v2].x;
		vector<double> mbr_max = data_vertices[v2].x;
		for (ui dim = 0; dim < dimension_x; dim++) {
			double min_value = 0;
			double max_value = 0;
			for (ui i = 0; i < deg; i++) {
				min_value += data_vertices[v2].neighbor_x_sort[i][dim];
			}
			for (ui i = data_vertices[v2].neighbor_x_sort.size() - deg; i < data_vertices[v2].neighbor_x_sort.size(); i++) {
				max_value += data_vertices[v2].neighbor_x_sort[i][dim];
			}
			mbr_min.push_back(min_value);
			mbr_max.push_back(max_value);
		}
		for (ui i = 0; i < dimension_o_prime; i++) {
			data_vertices[v2].MBRs[deg].push_back(alpha_value * mbr_min[i] + beta_value * data_vertices[v2].z[i]);
			data_vertices[v2].MBRs[deg].push_back(alpha_value * mbr_max[i] + beta_value * data_vertices[v2].z[i]);
		}
	}

	auto end = chrono::high_resolution_clock::now();
	update_cost += double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());

}

void edge_deletion_refine(Graphflow& mm, const vector<Query_Vertex>& query_vertices, vector<Data_Vertex>& data_vertices, ui v1, ui v2, ui label, double& continuous_cost) {
	std::vector<uint> m(mm.query_.NumVertices(), UNMATCHED);
	size_t num_results = 0ul;

	auto start = chrono::high_resolution_clock::now();
	for (uint i = 0; i < mm.query_.NumEdges(); i++)
	{
		uint u1 = mm.order_vs_[i][0], u2 = mm.order_vs_[i][1];
		auto temp_q_labels = mm.query_.GetEdgeLabel(u1, u2);

		if (query_vertices[u1].degree <= data_vertices[v1].degree
			&& query_vertices[u2].degree <= data_vertices[v2].degree
			&& query_vertices[u1].label==data_vertices[v1].label
			&& query_vertices[u2].label==data_vertices[v2].label) {
			ui j;
			for (j = 0; j < dimension_o; j++) {
				if (!(query_vertices[u1].o_prime[j] >= data_vertices[v1].MBRs[query_vertices[u1].degree][2 * j]
					&& query_vertices[u1].o_prime[j] <= data_vertices[v1].MBRs[query_vertices[u1].degree][2 * j + 1]
					&& query_vertices[u2].o_prime[j] >= data_vertices[v2].MBRs[query_vertices[u2].degree][2 * j]
					&& query_vertices[u2].o_prime[j] <= data_vertices[v2].MBRs[query_vertices[u2].degree][2 * j + 1])) {
					break;
				}
			}
			if (j == dimension_o)
			{
				m[u1] = v1;
				m[u2] = v2;
				mm.visited_[v1] = true;
				mm.visited_[v2] = true;

				mm.FindMatches(i, 2, m, num_results);

				mm.visited_[v1] = false;
				mm.visited_[v2] = false;
				m[u1] = UNMATCHED;
				m[u2] = UNMATCHED;

				if (num_results >= mm.max_num_results_) goto END_ENUMERATION;
			}
		}
		if (query_vertices[u1].degree <= data_vertices[v2].degree
			&& query_vertices[u2].degree <= data_vertices[v1].degree
			&& query_vertices[u1].label==data_vertices[v2].label
			&& query_vertices[u2].label==data_vertices[v1].label) {
			ui j;
			for (j = 0; j < dimension_o; j++) {
				if (!(query_vertices[u1].o_prime[j] >= data_vertices[v2].MBRs[query_vertices[u1].degree][2 * j]
					&& query_vertices[u1].o_prime[j] <= data_vertices[v2].MBRs[query_vertices[u1].degree][2 * j + 1]
					&& query_vertices[u2].o_prime[j] >= data_vertices[v1].MBRs[query_vertices[u2].degree][2 * j]
					&& query_vertices[u2].o_prime[j] <= data_vertices[v1].MBRs[query_vertices[u2].degree][2 * j + 1])) {
					break;
				}
			}
			if (j == dimension_o)
			{
				m[u1] = v2;
				m[u2] = v1;
				mm.visited_[v2] = true;
				mm.visited_[v1] = true;

				mm.FindMatches(i, 2, m, num_results);

				mm.visited_[v2] = false;
				mm.visited_[v1] = false;
				m[u1] = UNMATCHED;
				m[u2] = UNMATCHED;

				if (num_results >= mm.max_num_results_) goto END_ENUMERATION;
			}
		}
	}
END_ENUMERATION:
	auto end = chrono::high_resolution_clock::now();
	mm.num_negative_results_ += num_results;

	continuous_cost += double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

DAS3 build_index(const Dynamic_Graph& initial_graph, vector<Data_Vertex>& data_vertices) {
	ui max_degree = 0;
	for (ui i = 0; i < initial_graph.NumVertices(); i++) {
		if (initial_graph.GetDegree(i) > max_degree) {
			max_degree = initial_graph.GetDegree(i);
		}
	}
	vector<ui> degree_vertices_num(max_degree + 1, 0);
	vector<vector<ui>> degree_vertices_id(max_degree + 1);
	for (ui i = 0; i < initial_graph.NumVertices(); i++) {
		degree_vertices_num[initial_graph.GetDegree(i)]++;
		degree_vertices_id[initial_graph.GetDegree(i)].push_back(i);
	}

	DAS3 das3;
	das3.degree_index_map.resize(max_degree + 1);
	das3.index_list.resize(1);

	das3.index_list[0].deg_lb = 0;
	das3.index_list[0].deg_up = 0;
	das3.degree_index_map[0] = 0;

	ui avg_vertice_num = initial_graph.NumVertices() / degree_group_num;
	ui current_group = 1;
	ui current_vertex_num = 0;
	ui start_deg = 1;
	for (ui i = 1; i <= max_degree; i++) {
		current_vertex_num += degree_vertices_num[i];
		das3.degree_index_map[i] = current_group;
		if (current_vertex_num >= avg_vertice_num || i == max_degree) {
			Grid_Index index;
			index.deg_lb = start_deg;
			index.deg_up = i;

			das3.index_list.push_back(index);
			current_group++;
			start_deg = i + 1;
			current_vertex_num = 0;
			if (current_group > degree_group_num) {
				das3.index_list[degree_group_num].deg_up = max_degree;
				for (ui j = i; j < max_degree; j++) {
					das3.degree_index_map[i] = degree_group_num;
				}
				break;
			}
		}
	}

	vector<vector<double>> ranges(das3.index_list.size());
	vector<vector<double>> min_values(das3.index_list.size());
	vector<vector<double>> max_values(das3.index_list.size());
	for (ui index_id = 1; index_id < das3.index_list.size(); index_id++) {
		ranges[index_id] = vector<double>(dimension_o_prime, 0);
		min_values[index_id] = vector<double>(dimension_o_prime, numeric_limits<double>::max());
		max_values[index_id] = vector<double>(dimension_o_prime, 0);
	}
	for (ui vid = 0; vid < initial_graph.NumVertices(); vid++) {
		if (initial_graph.GetDegree(vid) != 0) {
			for (ui deg = 1; deg <= initial_graph.GetDegree(vid); deg++) {
				for (ui i = 0; i < dimension_o_prime; i++) {
					if (min_values[das3.degree_index_map[deg]][i] > data_vertices[vid].MBRs[deg][2 * i]) {
						min_values[das3.degree_index_map[deg]][i] = data_vertices[vid].MBRs[deg][2 * i];
					}
					if (max_values[das3.degree_index_map[deg]][i] < data_vertices[vid].MBRs[deg][2 * i + 1]) {
						max_values[das3.degree_index_map[deg]][i] = data_vertices[vid].MBRs[deg][2 * i + 1];
					}
				}
			}
		}
	}

	for (ui index_id = 1; index_id < das3.index_list.size(); index_id++) {
		for (ui i = 0; i < dimension_o; i++) {
			ranges[index_id][i] = (max_values[index_id][i] - min_values[index_id][i]) / (1.0 * cell_num);
		}
		ui sum_cell = pow(cell_num, dimension_o_prime);
		for (ui i = 0; i < sum_cell; i++) {
			vector<ui> dim;
			double factor = 1;
			for (ui j = 0; j < dimension_o; j++) {
				dim.push_back(ui(i / factor) % cell_num);
				factor *= cell_num;
			}

			Cell cell;

			for (ui j = 0; j < dimension_o; j++) {
				cell.bounds.push_back(min_values[index_id][j] + dim[j] * ranges[index_id][j]);
				cell.bounds.push_back(min_values[index_id][j] + (dim[j] + 1) * ranges[index_id][j]);
			}

			cell.key = 0;
			for (ui j = 0; j < dimension_o; j++) {
				cell.key += pow(cell.bounds[2 * j + 1], 2.0);
			}

			das3.index_list[index_id].cell_list.push_back(cell);
		}

		sort(das3.index_list[index_id].cell_list.begin(), das3.index_list[index_id].cell_list.end());
	}

	for (ui vid = 0; vid < initial_graph.NumVertices(); vid++) {
		data_vertices[vid].storage_pos.push_back(0);
		data_vertices[vid].keys.push_back(0.0);
		if (data_vertices[vid].degree != 0) {
			for (ui index_id = 1; index_id < das3.index_list.size(); index_id++) {
				ui min_deg, max_deg;
				if (data_vertices[vid].degree < das3.index_list[index_id].deg_lb) {
					break;
				}
				else if (data_vertices[vid].degree >= das3.index_list[index_id].deg_up) {
					min_deg = das3.index_list[index_id].deg_lb;
					max_deg = das3.index_list[index_id].deg_up;
				}
				else {
					min_deg = das3.index_list[index_id].deg_lb;
					max_deg = data_vertices[vid].degree;
				}
				Point point;
				point.vertex_id = vid;
				for (ui i = 0; i < dimension_o_prime; i++) {
					point.upper_bound.push_back(data_vertices[vid].MBRs[max_deg][2 * i + 1]);
				}
				for (ui deg = min_deg; deg <= max_deg; deg++) {
					point.MBRs.push_back(data_vertices[vid].MBRs[deg]);
				}

				point.key = 0;
				for (ui i = 0; i < dimension_o; i++) {
					point.key += pow(point.upper_bound[i], 2.0);
				}

				for (ui i = 0; i < das3.index_list[index_id].cell_list.size(); i++) {
					ui j;
					for (j = 0; j < dimension_o_prime; j++) {
						if (!(point.upper_bound[j] >= das3.index_list[index_id].cell_list[i].bounds[2 * j] && point.upper_bound[j] <= das3.index_list[index_id].cell_list[i].bounds[2 * j + 1])) {
							break;
						}
					}
					if (j == dimension_o_prime) {
						das3.index_list[index_id].cell_list[i].points.push_back(point);
						data_vertices[vid].storage_pos.push_back(i);
						data_vertices[vid].keys.push_back(point.key);
						break;
					}
				}
			}
		}
	}

	for (ui i = 0; i < das3.index_list.size(); i++) {
		for (ui j = 0; j < das3.index_list[i].cell_list.size(); j++) {
			sort(das3.index_list[i].cell_list[j].points.begin(), das3.index_list[i].cell_list[j].points.end());
		}
	}

	return das3;
}

vector<Query_Vertex> gen_VDE_query(const Dynamic_Graph& query_graph,const vector<double>& zipf) {
	vector<Query_Vertex> query_vertices(query_graph.NumVertices());
	for (ui i = 0; i < query_graph.NumVertices(); i++) {
		query_vertices[i].label=query_graph.GetVertexLabel(i);
		query_vertices[i].degree = query_graph.GetDegree(i);
		query_vertices[i].x = gen_VDE_z_zipf(query_graph.GetVertexLabel(i), dimension_x,zipf);
		query_vertices[i].z = gen_VDE_z(query_graph.GetVertexLabel(i), dimension_z);
	}
	for (ui i = 0; i < query_graph.NumVertices(); i++) {
		vector<ui> neighbors = query_graph.GetNeighbors(i);
		query_vertices[i].y = vector<double>(dimension_y, 0);
		for (ui j = 0; j < neighbors.size(); j++) {
			query_vertices[i].neighbors.insert(make_pair(neighbors[j], neighbors[j]));
			for (ui k = 0; k < dimension_y; k++) {
				query_vertices[i].y[k] += query_vertices[neighbors[j]].x[k];
			}
		}
		query_vertices[i].o = vector<double>();
		for (ui j = 0; j < dimension_x; j++) {
			query_vertices[i].o.push_back(query_vertices[i].x[j]);
		}
		for (ui j = 0; j < dimension_y; j++) {
			query_vertices[i].o.push_back(query_vertices[i].y[j]);
		}
		query_vertices[i].o_prime = vector<double>();
		for (ui j = 0; j < dimension_o; j++) {
			query_vertices[i].o_prime.push_back(alpha_value * query_vertices[i].o[j] + beta_value * query_vertices[i].z[j]);
		}

		query_vertices[i].key = 0;
		for (ui j = 0; j < dimension_o; j++) {
			query_vertices[i].key += pow(query_vertices[i].o_prime[j], 2.0);
		}
	}
	return query_vertices;
}

double candidateSearchStatic(const DAS3& das3, const vector<Data_Vertex>& data_vertices, vector<Query_Vertex>& query_vertices) {
	vector<double> times(query_vertices.size());
#pragma omp parallel for
	for (int qvid = 0; qvid < query_vertices.size(); qvid++) {
		auto start = chrono::high_resolution_clock::now();

		ui dgid = das3.degree_index_map[query_vertices[qvid].degree];

		for (ui cell_id = 0; cell_id < das3.index_list[dgid].cell_list.size(); cell_id++) {
			if (das3.index_list[dgid].cell_list[cell_id].points.size() == 0) {
				continue;
			}

			if (query_vertices[qvid].key > das3.index_list[dgid].cell_list[cell_id].key) {
				break;
			}

			ui degree_id = query_vertices[qvid].degree - das3.index_list[dgid].deg_lb;
			for (ui pid = 0; pid < das3.index_list[dgid].cell_list[cell_id].points.size(); pid++) {
				if (query_vertices[qvid].key > das3.index_list[dgid].cell_list[cell_id].points[pid].key) {
					break;
				}

				ui j;
				for (j = 0; j < dimension_o; j++) {
					if (das3.index_list[dgid].cell_list[cell_id].points[pid].upper_bound[j] < query_vertices[qvid].o_prime[j]) {
						break;
					}
				}
				if (j == dimension_o) {
					if (das3.index_list[dgid].cell_list[cell_id].points[pid].MBRs.size() > degree_id) {
						ui k;
						for (k = 0; k < dimension_o; k++) {
							if (!(das3.index_list[dgid].cell_list[cell_id].points[pid].MBRs[degree_id][2 * k] <= query_vertices[qvid].o_prime[k] &&
								das3.index_list[dgid].cell_list[cell_id].points[pid].MBRs[degree_id][2 * k + 1] >= query_vertices[qvid].o_prime[k])) {
								break;
							}
						}

						if (k == dimension_o) {
							query_vertices[qvid].candidates.push_back(das3.index_list[dgid].cell_list[cell_id].points[pid].vertex_id);
						}
					}
				}
			}
		}
		auto end = chrono::high_resolution_clock::now();
		times[qvid]=double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
	}

	return *max_element(times.begin(),times.end());
}

VertexID selectStartVertex(const Static_Graph* query_graph, ui* candidates_count) {
	ui start_vertex = 0;

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
		VertexID cur_vertex = i;

		if (candidates_count[cur_vertex] < candidates_count[start_vertex]) {
			start_vertex = cur_vertex;
		}
		else if (candidates_count[cur_vertex] == candidates_count[start_vertex]
			&& query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(start_vertex)) {
			start_vertex = cur_vertex;
		}
	}

	return start_vertex;
}

void updateValidVertices(const Static_Graph* query_graph, VertexID query_vertex, std::vector<bool>& visited,
	std::vector<bool>& adjacent) {
	visited[query_vertex] = true;
	ui nbr_cnt;
	const ui* nbrs = query_graph->getVertexNeighbors(query_vertex, nbr_cnt);

	for (ui i = 0; i < nbr_cnt; ++i) {
		ui nbr = nbrs[i];
		adjacent[nbr] = true;
	}
}

double generateQueryPlan(const Static_Graph* data_graph, const Static_Graph* query_graph, ui* candidates_count, ui*& order, ui*& pivot) {
	std::vector<bool> visited_vertices(query_graph->getVerticesCount(), false);
	std::vector<bool> adjacent_vertices(query_graph->getVerticesCount(), false);
	order = new ui[query_graph->getVerticesCount()];
	pivot = new ui[query_graph->getVerticesCount()];

	auto start = chrono::high_resolution_clock::now();
	VertexID start_vertex = selectStartVertex(query_graph, candidates_count);
	order[0] = start_vertex;
	updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
		VertexID next_vertex;
		ui min_value = data_graph->getVerticesCount() + 1;
		for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
			VertexID cur_vertex = j;

			if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
				if (candidates_count[cur_vertex] < min_value) {
					min_value = candidates_count[cur_vertex];
					next_vertex = cur_vertex;
				}
				else if (candidates_count[cur_vertex] == min_value && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(next_vertex)) {
					next_vertex = cur_vertex;
				}
			}
		}
		updateValidVertices(query_graph, next_vertex, visited_vertices, adjacent_vertices);
		order[i] = next_vertex;
	}

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
		VertexID u = order[i];
		for (ui j = 0; j < i; ++j) {
			VertexID cur_vertex = order[j];
			if (query_graph->checkEdgeExistence(u, cur_vertex)) {
				pivot[i] = cur_vertex;
				break;
			}
		}
	}
	auto end = chrono::high_resolution_clock::now();
	return double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

bool checkEdgeExistence(const vector<Data_Vertex>& data_vertices, VertexID u, VertexID v) {
	if (data_vertices[u].degree < data_vertices[v].degree) {
		std::swap(u, v);
	}

	if (data_vertices[v].neighbors.find(u) != data_vertices[v].neighbors.end()) {
		return true;
	}

	return false;
}

void generateBN(const Static_Graph* query_graph, ui* order, ui* pivot, ui**& bn, ui*& bn_count) {
	ui query_vertices_num = query_graph->getVerticesCount();
	bn_count = new ui[query_vertices_num];
	std::fill(bn_count, bn_count + query_vertices_num, 0);
	bn = new ui * [query_vertices_num];
	for (ui i = 0; i < query_vertices_num; ++i) {
		bn[i] = new ui[query_vertices_num];
	}

	std::vector<bool> visited_vertices(query_vertices_num, false);
	visited_vertices[order[0]] = true;
	for (ui i = 1; i < query_vertices_num; ++i) {
		VertexID vertex = order[i];

		ui nbrs_cnt;
		const ui* nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
		for (ui j = 0; j < nbrs_cnt; ++j) {
			VertexID nbr = nbrs[j];

			if (visited_vertices[nbr] && nbr != pivot[i]) {
				bn[i][bn_count[i]++] = nbr;
			}
		}

		visited_vertices[vertex] = true;
	}
}

void generateValidCandidates(const Static_Graph* query_graph, const Static_Graph* data_graph, const vector<Data_Vertex>& data_vertices, ui depth, vector<ui>& embedding,
	ui* idx_count, ui** valid_candidate, bool* visited_vertices, ui** bn, ui* bn_cnt, ui* order, ui* pivot) {
	VertexID u = order[depth];
	LabelID u_label = query_graph->getVertexLabel(u);
	ui u_degree = query_graph->getVertexDegree(u);

	idx_count[depth] = 0;

	VertexID p = embedding[pivot[depth]];
	ui nbr_cnt;
	const VertexID* nbrs = data_graph->getVertexNeighbors(p, nbr_cnt);

	for (ui i = 0; i < nbr_cnt; ++i) {
		VertexID v = nbrs[i];

		if (!visited_vertices[v] && u_label == data_graph->getVertexLabel(v) &&
			u_degree <= data_graph->getVertexDegree(v)) {
			bool valid = true;

			for (ui j = 0; j < bn_cnt[depth]; ++j) {
				VertexID u_nbr = bn[depth][j];
				VertexID u_nbr_v = embedding[u_nbr];

				if(!checkEdgeExistence(data_vertices,v, u_nbr_v)){
					valid = false;
					break;
				}
			}

			if (valid) {
				valid_candidate[depth][idx_count[depth]++] = v;
			}
		}
	}
}

double enumerateCandidates(const Static_Graph* data_graph, const Static_Graph* query_graph, const vector<Data_Vertex>& data_vertices, ui** candidates,
	ui* candidates_count, ui* order, ui* pivot, size_t output_limit_num) {
	size_t embedding_cnt = 0;
	int cur_depth = 0;
	int max_depth = query_graph->getVerticesCount();
	VertexID start_vertex = order[0];

	ui** bn;
	ui* bn_count;

	ui* idx;
	ui* idx_count;
	VertexID** valid_candidate;
	bool* visited_vertices;

	idx = new ui[max_depth];
	idx_count = new ui[max_depth];
	vector<ui> embedding(max_depth);
	visited_vertices = new bool[data_graph->getVerticesCount()];
	std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
	valid_candidate = new ui * [max_depth];

	ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
	for (ui i = 0; i < max_depth; ++i) {
		valid_candidate[i] = new VertexID[max_candidate_count];
	}

	idx[cur_depth] = 0;
	idx_count[cur_depth] = candidates_count[start_vertex];
	std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
		valid_candidate[cur_depth]);

	auto start = chrono::high_resolution_clock::now();

	generateBN(query_graph, order, pivot, bn, bn_count);

	while (true) {
		while (idx[cur_depth] < idx_count[cur_depth]) {
			VertexID u = order[cur_depth];
			VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
			embedding[u] = v;
			visited_vertices[v] = true;
			idx[cur_depth] += 1;

			if (cur_depth == max_depth - 1) {
				embedding_cnt += 1;

				visited_vertices[v] = false;
				if (embedding_cnt >= output_limit_num) {
					goto EXIT;
				}
			}
			else {
				cur_depth += 1;
				idx[cur_depth] = 0;
				generateValidCandidates(query_graph, data_graph, data_vertices, cur_depth, embedding, idx_count, valid_candidate,
					visited_vertices, bn, bn_count, order, pivot);
			}
		}

		cur_depth -= 1;
		if (cur_depth < 0)
			break;
		else
			visited_vertices[embedding[order[cur_depth]]] = false;
	}

EXIT:
	auto end = chrono::high_resolution_clock::now();
	delete[] bn_count;
	delete[] idx;
	delete[] idx_count;
	delete[] visited_vertices;
	for (ui i = 0; i < max_depth; ++i) {
		delete[] bn[i];
		delete[] valid_candidate[i];
	}

	delete[] bn;
	delete[] valid_candidate;

	return double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

double refineCandidates(const Static_Graph* static_data_graph, const Dynamic_Graph& query_graph, const vector<Data_Vertex>& data_vertices,
	const vector<Query_Vertex>& query_vertices) {
	Static_Graph* static_query_graph = new Static_Graph(true);
	static_query_graph->loadGraphFromDynamicGraph(query_graph);

	ui** candidates = new ui * [query_vertices.size()];
	ui* candidates_count = new ui[query_vertices.size()];
	for (ui i = 0; i < query_vertices.size(); i++) {
		candidates[i] = new ui[query_vertices[i].candidates.size()];
		candidates_count[i] = query_vertices[i].candidates.size();
		for (ui j = 0; j < query_vertices[i].candidates.size(); j++) {
			candidates[i][j] = query_vertices[i].candidates[j];
		}
	}
	ui* matching_order = new ui[query_vertices.size()];
	ui* pivots = new ui[query_vertices.size()];

	double time_ns = 0;

	time_ns += generateQueryPlan(static_data_graph, static_query_graph, candidates_count, matching_order, pivots);

	size_t output_limit = MAX_LIMIT;

	time_ns += enumerateCandidates(static_data_graph, static_query_graph, data_vertices, candidates, candidates_count, matching_order, pivots, output_limit);

	delete static_query_graph, candidates_count, matching_order, pivots;
	for (ui i = 0; i < query_vertices.size(); i++) {
		delete candidates[i];
	}

	return time_ns;
}

string to_string_with_precision(float value, int precision) {
    ostringstream out;
    out << fixed << setprecision(precision) << value;
    return out.str();
}