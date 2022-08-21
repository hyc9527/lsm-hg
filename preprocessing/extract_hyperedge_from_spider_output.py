#!/usr/bin/env python
# pyright: reportUnboundVariable=false, reportGeneralTypeIssues=false, reportMissingImports=false
import numpy as np

'''
    Extract adjacency matrix and charactors list from a movie script scrawled by python spider.
'''

# Input:  a movie script from python spider.
spider_output_path='./output/spider_output.csv'

# Output: adj mat and charactor list
output_path_character_lst='./output/starwars4_charactor_lst.csv'
output_path_adjmat='./output/starwars4_adjMat.csv'




# 0. Extract hypergraph from spider output
for csv_file in [spider_output_path]:
    fh = open(csv_file, 'r')
    all_scene=[]
    one_scene=[]
    count_scene = 0
    clapper_board_for_one_scene = False
    for line in fh.readlines()[1:2984]:
        line = line.strip('"')
        line = line.strip()
        if not line:                                                               # when it is empty line, skip
            continue
        if ('INT.' in line) or ('EXT.' in line) or ('END' in line):                # when 'INT.' or 'EXT.' happens
            if not clapper_board_for_one_scene:
                clapper_board_for_one_scene = not clapper_board_for_one_scene
                continue
            else:
                all_scene.append(one_scene) #  one valid hyperedge is finished
                count_scene += 1
                one_scene = []
                clapper_board_for_one_scene = not clapper_board_for_one_scene
                continue
        if not one_scene:
            one_scene = [line] # add first actor to one scene
        else:
            if line not in one_scene:
                one_scene.append(line) # add second different actor if it exists
    print("*** Total number of scenes: {}".format(count_scene))
    fh.close()

#print(all_scene)

# 1. clean up hypergraph

## 1.1 find a list of charactors of top degree
hyperedges = [edge for edge in  all_scene if len(edge) > 2]  # extract edge and hyperedges
actor_dict = {} # histogram on all actor
actor_hg  = [] # keep edges between main actors only, whose names show up in charcter_lst. for instance "IMPERIAL OFFICER" is removed
for edge in hyperedges:
    new_edge = []
    for actor in edge:
        #if actor in character_lst:
        if actor:
            new_edge.append(actor)
            if actor in actor_dict:
                actor_dict[actor] += 1
            else:
                actor_dict[actor] = 1
charactor_lst = [key for (key,val) in actor_dict.items() if val >2]
#print(f'*** Main charactors before cleaning are : {charactor_lst}')
#print(actor_dict)
charactor_lst.remove('TROOPER')
charactor_lst.remove('OFFICER')
print(f'*** Main charactors after cleaning are : {charactor_lst}')



if True:
    #output_path_adjmat = 'starwars4_charactor_lst.csv'
    #output_character_lst_path = './output/test_starwars4_charactor_lst.csv'
    np.savetxt(output_path_character_lst, charactor_lst, delimiter=' ', fmt='%s')
    print(f'*** Charactor list is saved in {output_path_character_lst}')


##  1.2 Clean hypergraph by keep ONLY actors in charactor list
actor_dict = {} # histogram on all actor
actor_hg  = [] # keep edges between actors only, whose names show up in charcter_lst. for instance "IMPERIAL OFFICER" is removed
for edge in hyperedges:
    new_edge = []
    for actor in edge:
        if actor in charactor_lst:
            new_edge.append(actor)
            if actor in actor_dict:
                actor_dict[actor] += 1
            else:
                actor_dict[actor] = 1
    if len(new_edge) > 1:
        actor_hg.append(new_edge)
print("*** After cleaning minor charactors, total number of scenes in the episode is {}".format(len(actor_hg)))
# convert hypergraph to adjacency matrix.
actor_lst = list(actor_dict.keys())
actor_num = len(actor_lst)
scene_num = len(actor_hg)

actor_scene_adj = np.zeros((actor_num, scene_num))
for j in range(scene_num):
    for auth in actor_hg[j]:
        which_auth = actor_lst.index(auth)
        actor_scene_adj[which_auth, j] = 1


if True:
    np.savetxt(output_path_adjmat, actor_scene_adj, delimiter=',', fmt='%i')
print(f"*** Save adjacency matrinx in {output_path_adjmat}")

