#!/bin/bash

# The current folder where output csv files(adjacney matrix, author list, original data from spider) sit at
data_preprocessing_folder="/Users/yichen/Desktop/p1/lsm-hg/preprocessing/aos"
#data_preprocessing_folder="YOUR_LOCAL_PREPROCESSING_FOLDER"

# python script to extract adjacney matrix, author list from original data from spider output
helper_extractor="extract_hyperedge_from_spider_output.py"

# run spider for aos dataset
## data size (issue date from when to when) is specified in the spider script
spider_folder="/Users/yichen/Desktop/p1/my-github-local/preprocessing/aos/crawler_aos/crawler_aos"

## spider_output_file="${data_preprocessing_folder}/output_spider.csv"
spider_output_file="${data_preprocessing_folder}/output/output_spider.csv"


# Rull shell
cd $spider_folder
scrapy crawl 'aos-2021' -o $spider_output_file
echo "#### Spider output file is saved in $spider_output_file"
cd $data_preprocessing_folder
echo "#### Extract author-paper adj matrix and authorlist from spider output. "
python3 $helper_extractor
echo '#### Done.'
