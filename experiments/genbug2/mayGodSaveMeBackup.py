from itertools import combinations
import math
nodes=[]
dist=lambda x,y:math.sqrt((nodes[x][0]-nodes[y][0])**2+(nodes[x][1]-nodes[y][1])**2)
combinations = []

lst=list(range(170))
def condition_func(l):
    if len(l)<=1: return True
    for i in l[:-1]:
        if dist(i,l[-1])<1: 
            return False
    return True

def readMap():
    with open("../Instances/roadmaps/sparse/map.graph") as f:
        node_num=int(f.readline())
        for i in range(node_num):
            l=f.readline()
            nodes.append(list(map(float,l.split(" "))))

def generate_combinations(r):
    generate_combinations_recursive(lst,r, [], combinations)

def generate_combinations_recursive(lst, r, current_combination, combinations):
    if r == 0:
        combinations.append(tuple(current_combination))
        #print(tuple(current_combination))
        return

    for i in range(len(lst)):
        current_combination.append(lst[i])
        if condition_func(current_combination):
            generate_combinations_recursive(lst[i+1:], r-1, current_combination, combinations)
        current_combination.pop()

if __name__=="__main__":
    r=4
    readMap()
    generate_combinations(r)
    print(len(combinations))
