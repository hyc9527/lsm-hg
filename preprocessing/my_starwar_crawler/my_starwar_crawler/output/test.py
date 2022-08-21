import numpy as np
import pandas as pd

# Extract a list of lists of chractors from Episode script(s)
for data in ['new_hope.csv']:
#for data in ['test.csv']:
    fh = open(data, 'r')
    all_scene=[]
    one_scene=[]
    count_scene = 0
    clapper_board_for_one_scene = False
    for line in fh.readlines():
        line = line.strip()
        if not line:                                                    # when it is empty line, skip
            continue
        if ('INT. ' in line) or ('EXT. ' in line) or ('END' in line):                # when 'INT.' or 'EXT.' happens
            if not clapper_board_for_one_scene:
                clapper_board_for_one_scene = not clapper_board_for_one_scene
                continue
            else:
                all_scene.append(one_scene)
                count_scene += 1
                #print("all_scene is {}".format(all_scene))
                one_scene = []
                clapper_board_for_one_scene = not clapper_board_for_one_scene
                continue
        if not one_scene:
            one_scene = [line]
        else:
            if line not in one_scene:
                one_scene.append(line)
    print("total number of scenes in the episode is {}".format(count_scene))
    fh.close()


# Extract full chractor hypergraph
full_charater_hypergraph = [edge for edge in  all_scene if len(edge) > 2]  # extract edge and full_charater_hypergraph
print("total number hyperedges in full_charater_hypergraph :  {}".format(len(full_charater_hypergraph)))




### Wrong way to clean hyperedge !!!
# clean full_charater_hypergraph by keeping ONLY actors in the charactor list and removing those not
# the full charactor list by gabasova :
# https://github.com/evelinag/StarWars-social-network/tree/master/data/characters.csv
# alias list by gabasova:
# https://github.com/evelinag/StarWars-social-network/blob/master/data/aliases.csv
# df_character = pd.read_csv('characters.csv')
# character_lst = [x[0] for x in df_character.values]
# actor_dict = {} # histogram on all actor
# actor_hg  = [] # keep edges between actors only, whose names show up in charcter_lst. for instance "IMPERIAL OFFICER" is removed
# for edge in full_charater_hypergraph:
#     new_edge = []
#     for actor in edge:
#         if actor in character_lst:
#             new_edge.append(actor)
#             if actor in actor_dict:
#                 actor_dict[actor] += 1
#             else:
#                 actor_dict[actor] = 1
#         if not new_edge:
#             actor_hg.append(new_edge)


character_lst = []

actor_dict = {} # histogram on all actor
for edge in full_charater_hypergraph:
    for actor in edge:
        if actor in actor_dict:
            actor_dict[actor] += 1
        else:
            actor_dict[actor] = 1
major_actor_lst = [key for (key,val) in actor_dict.items() if val >2] # keep top nine actors
major_actor_lst.remove('TROOPER')
major_actor_lst.remove('OFFICER')
print(major_actor_lst) 



# create a hypergraph of actors on major_actor_lst
major_actor_dict = {}
major_actor_hg  = [] # keep edges between actors only, whose names show up in charcter_lst. for instance "IMPERIAL OFFICER" is removed
for edge in full_charater_hypergraph:
    new_edge = []
    for actor in edge:
        if actor in major_actor_lst:
            new_edge.append(actor)
    if len(new_edge)>1:
        major_actor_hg.append(new_edge)
        for x in new_edge:
                if x in major_actor_dict:
                    major_actor_dict[x] += 1
                else:
                    major_actor_dict[x] = 1

major_actor_dict
major_actor_hg



# create a adj matrix: actors by scenes

actor_num = len(major_actor_lst)
scene_num = len(major_actor_hg)
actor_scene_adj = np.zeros((actor_num, scene_num))
for j in range(scene_num):
    for auth in major_actor_hg[j]:
        which_auth = major_actor_lst.index(auth)
        actor_scene_adj[which_auth, j] = 1
actor_scene_adj.shape

np.sum(actor_scene_adj, axis = 0) # num of actors per scene
np.sum(actor_scene_adj, axis = 1) # actors' total shoots in the episode






file_name = data.replace('.csv', '_adjMat.csv')
np.savetxt(file_name, actor_scene_adj, delimiter =',')
file_name = data.replace('.csv', '_actor_lst.csv')
np.savetxt(file_name, major_actor_lst, delimiter=",", fmt='%s')


