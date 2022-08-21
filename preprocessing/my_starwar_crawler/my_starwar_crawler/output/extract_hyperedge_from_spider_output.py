#import re,sys
import pandas as pd
import numpy as np



for csv_file in ['new_hope.csv']:
#for csv_file in ['test.csv']:
    fh = open(csv_file, 'r')
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
    
actor_hg
df_character = pd.read_csv('characters.csv')
character_lst = [x[0] for x in df_character.values]
actor_dict = {} # histogram on all actor
actor_hg  = [] # keep edges between actors only, whose names show up in charcter_lst. for instance "IMPERIAL OFFICER" is removed
for edge in hyperedges:
    new_edge = []
    for actor in edge:
        if actor in character_lst:
            new_edge.append(actor)
            if actor in actor_dict:
                actor_dict[actor] += 1
            else:
                actor_dict[actor] = 1
        if not new_edge:
            actor_hg.append(new_edge)
major_charactor_lst = [key for (key,val) in actor_dict.items() if val >2] # keep top nine actors
print(major_charactor_lst)
hyperedges = [edge for edge in  all_scene if len(edge) > 2]  # extract edge and hyperedges
    print("total number of hyperedges in the episode is {}".format(len(hyperedges))) 
    df_character = pd.read_csv('characters.csv')
    r_charactor_lst

    character_lst = [x[0] for x in df_character.values]
    # create main charater hypergraph
    actor_dict = {}
    for edge in hyperedges:
        print(edge)




#############
#fh = open('new_hope.csv', 'r')
#line = fh.readlines()[2]
#print(line)
#fh.close()

############
# adjM = np.load('tt_adjMat.npy')
# adjM.shape
# col_sum = np.sum(adjM, axis=0)
# row_sum = np.sum(adjM, axis=1)
