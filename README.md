# DSM

In this project, we open-source the source code of our S-DSM/C-DSM approaches.

The real-world and default synthetic datasets used in our paper are stored in the datasets directory. As some synthetic datasets are large, we do not upload them. You can easily generate them by following the instructions in our paper.

On Git Hub, we will introduce how to reproduce the results of our experiments.

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