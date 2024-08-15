#include "custom.h"
#include "CLI11.hpp"
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
#include <omp.h>
#include <chrono>
#include <cstdio>

using namespace std;
using namespace chrono;

int main(int argc, char *argv[])
{
	uint max_num_results = MAX_LIMIT, time_limit = UINT_MAX, initial_time_limit = UINT_MAX;
	bool print_prep = false, print_enum = false, homo = false, report_initial = false;

	ui query_num = 100;

	ui edge_deletion = 0;

	string data_graph_name = "";

	CLI::App app{"App description"};
	app.add_option("-d,--data", data_graph_name, "initial data graph path")->required();
	app.add_option("-e,--deletion", edge_deletion, "edge_deletion");
	app.add_option("-g,--group", degree_group_num, "degree group num");
	app.add_option("-c,--cell", cell_num, "cell num");
	app.add_option("-x,--dimx", dimension_x, "embedding dimension");
	app.add_option("-a,--alpha", alpha_value, "alpha value");
	app.add_option("-b,--beta", beta_value, "beta value");
	app.add_option("-m,--max", MAX_LIMIT, "max limit");

	CLI11_PARSE(app, argc, argv);

	dimension_y = dimension_x;
	dimension_z = 2 * dimension_x;
	dimension_o = 2 * dimension_x;
	dimension_o_prime = 2 * dimension_x;
	dimension_mbr = 2 * dimension_o_prime;

	string dataset_path = "Datasets/" + data_graph_name + "/";

	string initial_path = dataset_path + "initial_graph.graph";
	string data_graph_path = dataset_path + "data_graph.graph";
	string insertion_path = dataset_path + "insertion_graph.graph";
	string deletion_path = dataset_path + "deletion_graph.graph";

	vector<double> zipf;
	string filename = "Datasets/zipf.txt";
	ifstream fin(filename);
	ui zipf_num;
	fin >> zipf_num;
	for (ui i = 0; i < zipf_num; i++)
	{
		double value;
		fin >> value;
		zipf.push_back(value);
	}
	fin.close();

	Dynamic_Graph initial_graph{};
	if (edge_deletion == 0)
	{
		initial_graph.LoadFromFile(initial_path);
	}
	else
	{
		initial_graph.LoadFromFile(data_graph_path);
	}

	cout << data_graph_name << endl;
	initial_graph.PrintMetaData();
	auto start = chrono::high_resolution_clock::now();
	vector<Data_Vertex> data_vertices = gen_VDE_data(initial_graph, zipf);
	DAS3 das3 = build_index(initial_graph, data_vertices);
	auto end = chrono::high_resolution_clock::now();
	double Continuous_Offline_Cost = double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());

	Static_Graph *static_data_graph = new Static_Graph(true);
	static_data_graph->loadGraphFromDynamicGraph(initial_graph);

	vector<Dynamic_Graph> query_graphs(query_num);
	vector<Graphflow> mms;
	vector<vector<Query_Vertex>> query_vertices(query_num);

	vector<double> time_ns_preprocessing(query_num, 0);
	vector<double> time_ns_search_static(query_num, 0);
	vector<double> time_ns_refinement_static(query_num, 0);
	vector<double> time_ns_search(query_num, 0);
	vector<double> time_ns_refinement(query_num, 0);
	vector<double> initial_cost(query_num, 0);
	vector<double> continuous_cost(query_num, 0);

	double Initial_Cost = 0, Continuous_Cost = 0, Update_Cost = 0;

	double candidates_num = 0;
	double vertex_sum = 0;

	for (ui qid = 0; qid < query_num; qid++)
	{
		string query_path = dataset_path + "query_graphs/" + "query_graph-" + to_string(qid) + ".graph";
		query_graphs[qid].LoadFromFile(query_path);
		mms.push_back(Graphflow(query_graphs[qid], initial_graph, max_num_results, print_prep, print_enum, homo));

		start = chrono::high_resolution_clock::now();
		mms[qid].Preprocessing();
		query_vertices[qid] = gen_VDE_query(query_graphs[qid], zipf);
		end = chrono::high_resolution_clock::now();
		time_ns_preprocessing[qid] = double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());

		time_ns_search[qid] = candidateSearchStatic(das3, data_vertices, query_vertices[qid]);

		time_ns_refinement[qid] = refineCandidates(static_data_graph, query_graphs[qid], data_vertices, query_vertices[qid]);

		initial_cost[qid] = time_ns_preprocessing[qid] + time_ns_search[qid] + time_ns_refinement[qid];
	}
	delete static_data_graph;

	if (edge_deletion == 0)
	{
		initial_graph.LoadUpdateStream(insertion_path);
	}
	else
	{
		initial_graph.LoadUpdateStream(deletion_path);
	}

	ui num_e_updates = 0;
	while (!initial_graph.updates_.empty())
	{
		if (num_e_updates % 10000 == 0)
		{
			cout << "Update Edge Number: " << num_e_updates << endl;	
		}

		InsertUnit insert = initial_graph.updates_.front();
		initial_graph.updates_.pop();

		if (insert.type == 'e' && insert.is_add)
		{
			edge_insertion_update(data_vertices, das3, insert.id1, insert.id2, Update_Cost);
			initial_graph.AddEdge(insert.id1, insert.id2, insert.label);
			for (ui qid = 0; qid < query_num; qid++)
			{
				edge_insertion_refine(mms[qid], query_vertices[qid], data_vertices, insert.id1, insert.id2, insert.label, continuous_cost[qid]);
			}
			num_e_updates++;
		}
		else if (insert.type == 'e' && !insert.is_add)
		{
			edge_deletion_update(data_vertices, das3, insert.id1, insert.id2, Update_Cost);
			for (ui qid = 0; qid < query_num; qid++)
			{
				edge_deletion_refine(mms[qid], query_vertices[qid], data_vertices, insert.id1, insert.id2, insert.label, continuous_cost[qid]);
			}
			initial_graph.RemoveEdge(insert.id1, insert.id2);
			num_e_updates++;
		}
	}
	cout << "Update Edge Number: " << num_e_updates << endl;

	for (ui qid = 0; qid < query_num; qid++)
	{
		Initial_Cost += initial_cost[qid];
		Continuous_Cost += continuous_cost[qid];
	}

	cout << "Query ID | Query Cost (ms)" << endl;
	for (ui qid = 0; qid < query_num; qid++) {
		cout << qid << " | " << NANOSECTOMSEC(initial_cost[qid]+continuous_cost[qid]) << endl;
	}
	cout << "Average Cost (ms): " << NANOSECTOMSEC((Initial_Cost + Continuous_Cost + Update_Cost) / query_num) << endl;

	return 0;
}
