# Hierarchical Pathfinding Algorithm
Path planning for autonomous driving in real scenarios, large-scale crowds or video games is a challenge to solve, due to the inherent dynamism in the scenes and vast search space, it is often necessary to reprocess the path or make an additional process to return to the original one after avoiding a collision, making the process computationally expensive. In this sense, the use of a hierarchical approach offers different advantages: 

- it divides the search space reducing the problem complexity
- we can obtain a partial path faster which allows to a character start moving and compute the rest of the path asynchronously
- we can reprocess only a part if necessary with different levels of abstraction.

Typically the approach employed is to construct a graph from an occupancy matrix and then to traverse the graph using graph-algorithms such as A*, Dijkstra or Breadth-First Search. All of the incumbent methods are very mature. However, they consume relatively high amounts of memory, are slow due to the number of pointer manipulations and are not readily amenable to hardware acceleration for real-time operation. Even enhancements such as visibility graphs do nothing to change the fundamental nature of the graph algorithms and require significant levels of pre-processing which are not feasible to implement in real-time.

This algorithm operates directly on the grid or bitmap rather than first converting it to a graph making the implementation very hardware friendly. It is domain-independent, and it could be adapted to dynamic scenes. The method makes use of a hierarchical representation which incrementally builds the equivalent of a visibility graph without pre-processing. The hierarchical approach reduces the search space, identifying the areas of interest in high levels of abstraction, allows the tree to be pruned as it is constructed using a heuristic based on density and distance. In this way, the algorithm gives priority to the shortest paths in areas with fewer obstacles.

## Hierarchical Data Structures
To reduce the search space, we create a hierarchical structure. We scale the original map creating auxiliary maps with different levels of abstraction. We reduce the resolution of the map while maintaining connectivity and density information in the abstraction. Every level is four times smaller than the previous one, e.g., if the size of the original map is 256x256, next level is 64x64, next one 16x16, the map with the smallest resolution is 4x4. 

Next figure shows the abstraction for the different levels. Dots represent available cells, and a different character represents occupied cells.

![image](https://github.com/vhpvmx/HPV/blob/main/Imgs/mapL1.png)
![image](https://github.com/vhpvmx/HPV/blob/main/Imgs/mapL2.png)
![image](https://github.com/vhpvmx/HPV/blob/main/Imgs/mapL3.png)
![image](https://github.com/vhpvmx/HPV/blob/main/Imgs/originalMap.png)


Each cell of a certain level (except the last one) represents a block of 16 cells (4x4) in the next level. The criterion to determine if a block should be considered occupied or free is: if at least one of the 16 cells that comprise it is free, then the whole block is deemed to be free. Otherwise, it would not be possible to define a path that would reach that free cell.

We storage the maps in blocks of 4x4, for level one there is only one block, for level 2 there are 16 blocks, for level 3 there are 256 blocks, for level 4 there are 4096 blocks. This distribution is an essential part of the algorithm because the main idea is to go from one block to the next one in each level, using as guide the path of the previous level.

Next figure shows an example how it the algorithm finds the path in different levels.

![image](https://github.com/vhpvmx/HPV/blob/main/Imgs/pathLevels.png)

