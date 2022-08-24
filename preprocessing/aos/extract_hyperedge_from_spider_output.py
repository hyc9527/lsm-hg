#!/usr/bin/env python
# pyright: reportUnboundVariable=false, reportGeneralTypeIssues=false, reportMissingImports=false

import numpy as np
import re

# Convert spider output(a list of paper coauthorship) to author-paper adjacency mat and author list

# Input:
spider_output_path='./output/output_spider.csv'


# Output:
output_path_adjMat = './output/aos_adjMat.csv'
output_path_auth =  './output/aos_author_lst.csv'


# Main
author_count = {}  # record authours occurences for all lines
coauthorship_lst = []
for data in [spider_output_path]:
    fh = open(data, 'r')
    for line in fh.readlines()[1:]:  # skip header line
        line = line.strip('\n')
        authors = re.findall(r'"([^"]*)"', line)
        if authors:  # to deal with a speicial case in csv file, single author without "", returns []
            one_paper = authors[0].split(',')
            if len(one_paper) > 1:  # skip singe authored paper
                coauthorship_lst.append(one_paper)
                for auth in one_paper:
                    if auth in author_count:
                        author_count[auth] += 1
                    else:
                        author_count[auth] = 1
    fh.close()
    author_lst = list(author_count.keys())
    author_num = len(author_lst)
    paper_num = len(coauthorship_lst)
    author_paper_adj = np.zeros((author_num, paper_num))
    for j in range(paper_num):
        for auth in coauthorship_lst[j]:
            which_auth = author_lst.index(auth)
            author_paper_adj[which_auth, j] = 1

if True:
    np.savetxt(output_path_adjMat, author_paper_adj, delimiter=',', fmt='%i')
    print(f"*** Save adjacency matrinx in ${output_path_adjMat}")

if True:
    np.savetxt(output_path_auth, author_lst, delimiter=' ', fmt='%s')
    print(f'*** Charactor list is saved in ${output_path_auth}')




#   with open(output_file_adj, "w") as f: # "w" => creaate a new file eavh time and add data to it
#        np.savetxt(f, author_paper_adj, delimiter=',')
#    with open(output_file_auth,"w") as f:
#        np.savetxt(f, author_lst, delimiter=",", fmt='%s')







