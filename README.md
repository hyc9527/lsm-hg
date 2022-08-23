# A latent space model of general hypergraphs

R implementation of examples in chapter $4$

### Folders
- `code`: source code
- `data`: Starwar hypergraph dataset
- `demo`:  R implementation of three experiemnts in chapter $4$



### Requirements
The project is written in R 4.1.3 with several R packages. The file `requirments.R` contains a list of R packages requried for demo. To install these packages, run:
``` shell
    Rscript requirments.R
```



### Usages

We documents experiments in jupyter notebook files in folde `/demo/`.
To run demo, one need jupyter notebook, R kenerl and several R packages installed.

```shell
    cd demo
```

- Experiment 1: simulation study on dependence between shape parameter and cluster coefficient
    - run `demo/demo_cluster_coef.ipynb`

- Experiment 2: demo on single run of mcem on uniform/general hypergraph
    - run `demo/demo_mcem_general_hg.ipynb`

- Experiment 3: analytics on Star War hypergraphs
    - run `demo/demo_analysis_starwar.ipynb`



