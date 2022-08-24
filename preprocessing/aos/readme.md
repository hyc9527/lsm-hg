# Data preprocessing : AOS coauthorship 2021

Basic steps:

1. [Web scraping](#web-scraping)
2. [Extracting hypergraph](#extracting-hypergraph)

### Web scrabing

A python spider is in `./crawler_aos`, which  is to scrawl the coauthorship data from the AOS website.


### Extracting hypergraph
Extract and clean hypergraph from the spider output; Save the hypergraph in a format of author-by-paper adjacency matrix `outout/aos_adjMat.csv` alone with author list `output/aos_charactor_lst.csv`.


### Usage

To run the spider, extract hypergraph and save results in `./output/`.
```shell
./main.sh
```

