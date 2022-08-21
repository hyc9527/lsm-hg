# Data preprocessing : Star War IV script

Basic steps:

1. [Web scraping](#web-scraping)
2. [Extracting hypergraph](#extracting-hypergraph)

### Web scrabing

A python spider is in `./my_starwar_crawler`, which  is to scrawl the movie script "Star War IV: New Hope" from the [IMDB website](https://imsdb.com/scripts/Star-Wars-A-New-Hope.html).


### Extracting hypergraph
Extract and clean hypergraph from the spider output; Save the hypergraph in a format of actor-by-scene adjacency matrix `outout/starwar4_adjMat.csv` alone with charactor list `output/starwars4_charactor_lst.csv`.


### Usage

To run the spider, extract hypergraph and save in `output`.
```shell
./main.sh
```

