# DSM

In this project, we open-source the source code of our S-DSM/C-DSM approaches and upload the technical report of our paper "Dynamic Subgraph Matching via Cost-Model-based Vertex Dominance Embeddings (Technical Report)".

## Check List

Readme File ✔

Technical Report ✔

Source Code ✔

Real/Synthetic Graph Data Sets ✔

## Source Code
We will introduce how to reproduce the results of our experiments, where all parameters are set by their default values.

1. Download the whole project.

2. Under the root directory of the online directory, execute the following commands to compile the source code.

```
mkdir build
cd build
cmake ..
make
```

3. Execute the following command to run the experiment over the Yeast dataset.

```
./build/DSM -d Yeast
```

4. You can also try other datasets or the case of edge deletion following the table below.

| Parameter of Command Line | Description | Optional Values |
| :-----------: | ----------- | ----------- |
| -d | data graph name | Syn-Uni, Syn-Gau, Syn-Zipf, Yeast, HPRD, DBLP, YouTube, USPatents|
| -e | edge deletion | **0** for edge insertion, 1 for edge deletion |

## Real/Synthetic Graph Data Sets

The real-world and default synthetic datasets used in our paper are stored in the datasets directory. As some synthetic datasets are large, we do not upload them. You can easily generate them by following the instructions in our paper.

The statistics of these data sets are summarized as follows:

| Data Sets | # of Vertices | # of Edges | # of Distinct Labels | Average Degree |
| :-----------: | :-----------: | :-----------: | :-----------: | :-----------: |
| Syn-Uni | 50,000 | 125,000 | 15 | 5 |
| Syn-Gau | 50,000 | 125,000 | 15 | 5 |
| Syn-Zipf | 50,000 | 125,000 | 15 | 5 |
| Yeast | 3,112 | 12,519 | 71 | 8 |
| HPRD | 9,460 | 34,998 | 307 | 7.4 |
| DBLP | 317,080 | 1,049,866 | 15 | 6.6 |
| YouTube | 1,134,890 | 2,987,624 | 25 | 5.3 |
| USPatents | 3,774,768 | 16,518,947 | 20 | 8.8 |

Note that, due to the upload size limitations, we compress the "data_graph.graph" and "initial_graph.graph" of USPatents data set. Thus, we need to decompress these two files as follows before using them:

```
cd Datasets/USPatents
tar -xzvf data_graph.tar.gz
tar -xzvf initial_graph.tar.gz
```
