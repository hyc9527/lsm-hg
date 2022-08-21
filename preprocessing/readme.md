# Data preprocessing : Star War IV script

Basic steps:

1. [Web scraping](#web-scraping)
2. [Extracting hypergraph](#extracting-hypergraph)

### Web scrabing

Scrawl the movie script "Star War IV: New Hope" from [IMDB website](https://imsdb.com/scripts/Star-Wars-A-New-Hope.html). A node-by-edge adjacency matrix is  saved in `spider_output.csv`.



### Extracting hypergraph
Extract and clean hypergraph from the spider output; Save the hypergraph in a format of actor-by-scene adjacency matrix `outout/starwar4_adjMat.csv` alone with charactor list `output/test_starwars4_charactor_lst.csv`.





### Usage

```shell
./extract_hyperedge_from_spider_output.py
```

