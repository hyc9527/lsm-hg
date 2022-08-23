#!/bin/bash

# The main folder where `main.sh` sits
data_preprocessing_folder="YOUR_LOCAL_PATH_TO_THE_PREPROCESSING_FOLDER"

# The spider folder for web scraping
spider_folder="$data_preprocessing_folder/my_starwar_crawler/my_starwar_crawler"

# The spider output file
spider_output_file="$data_preprocessing_folder/output/spider_output.csv"

# helper to extract adjacney matrix, author list from spider output
helper_extractor="extract_hyperedge_from_spider_output.py"


# web scraping for movie script
cd $spider_folder
scrapy crawl 'new_hope' -o $spider_output_file                   # save spider output csv file in preprocessing folder
echo "#### Spider output file is saved in $spider_output_file"

# extract hypergraph from the movie script
cd $data_preprocessing_folder
echo "#### Extract author-paper adj matrix and authorlist from spider output. "
python3 $helper_extractor
echo '#### Done.'
