// HPV.cpp - HierarchicalPathfindingVola

//
// Name       : HPV.cpp
// Author     : Hugo Perez and David Moloney
// Company    : Movidius Ltd. (Intel), d'Olier St., Dublin 2, Ireland
// Date       : Oct 2017
// Purpose    : Hierarchical Pathfinding Algorithm
// Vola Level 1 4x4
// Vola Level 2 16x16
// Vola Level 3 64x64
// Vola Level 4 256x256
// Vola Level 5 1024x1024
// Arrays to store density
// For L4 256x256
// L1[0][16] <0,4096>
// L2[1][256] <0,256>
// L3[2][4096] <0,16>
// For L5 1024x1024
// L1[0][16] <0,65536>
// L2[1][256] <0,4096>
// L3[2][4096] <0,256>
// L4[3][65536] <0,16>
// JSON file structure
// {"Experiment": 1 "start": [100,118] "goal": [96,120]}
// consecutive number, level, x, y, action(select, eval)
// (0, 1, 1,1, select)


#include <iostream>
#include <vector>
#include <algorithm>    // std::reverse(vector)
#include <stdio.h>
#include "ScenarioLoader.h"
#include <boost/timer/timer.hpp>
#include <fstream>

using namespace std;

//#define DEBUG
//#define DEBUG_L2
//#define DEBUG_L3
//#define PREPROCESS
#define ROWS 4
#define COLS 4
#define ROWSxCOLS 16
#define MAX_NUM_REP 10

#ifdef DEBUG_L3
ofstream ofile("results.txt");
ofstream ofs("reprocess_stats.txt");
#endif

typedef unsigned long long u64;

int total_path_length = 0;
int vola_level = 0;
vector<int> failed_exp;
int *reprocess;
int exp_number;
float **density;
vector<bool> **g_grid;
u64** g_conn;
u64** g_conn_orig;
int *width;
int *max_density;
int *g_start;
int *g_goal;
vector<int> *g_paths;

//#ifdef PREPROCESS
vector<int> ****G_LUT;
vector<int> ****G_LUT_orig;
vector<bool> g_reprocess;
vector<int> g_reprocess_counter;
//#endif

// LUT to convert a path to a connectivity vector
// This LUT helps to check status of every bit which represent connectivity between two cells, left-up-right-down
// e.g. 1011 this value indicates connectivity to left, right and down
static u64 conn_dir[4][2] = {
		{8,2},	// 1000, 0010
		{2,8},	// 0010, 1000
		{4,1},	// 0100, 0001
		{1,4}	// 0001, 0100
};

// define directions to use in grid-search
#define DIRECTIONS 4 // 0-3 Manhattan, 4-7 Diagonal

// row/col transformations to apply to test new points on grid
static int tests[8] = {
  1, // Manhattan: E (right)
  -1, // Manhattan: W (left)
  ROWS, // Manhattan: S (down)
  -ROWS, // Manhattan: N (up)
  //
  -ROWS+1, // Diagonal : NE
  -ROWS-1, // Diagonal : NW
  ROWS+1, // Diagonal : SE
  ROWS-1  // Diagonal : SW
};

// Main functions
bool comparePaths(vector<int> l_paths[], int test_level);
bool findPathxLevel(vector<int> l_paths[], int startb, int goalb, int widtha, int widthb, vector<bool> l_grid[], int &level,
						bool b_reprocess, vector<int> *LUT[][ROWSxCOLS][ROWSxCOLS]);

//bool findPathxLevel(vector<int> l_paths[], int startb, int goalb, int widtha, int widthb, vector<bool> l_grid[], int level, vector<int> LUT[][ROWSxCOLS][ROWSxCOLS]);
bool findGlobalPath(int initial_level, int final_level, vector<int> l_paths[], bool b_reprocess);
bool findPartialPath(int i, int idx1, int idx2, int idx3, int idx4, vector<int> path_la, vector<int> &path_lb, int widthb, int widtha,
		vector<bool> l_grid[], int level, float &local_density, bool &conn_next, bool b_idx4 );
bool findPartialPathLast(int i, int idx1, int idx2, vector<int> path_la, vector<int> &path_lb, int widthb, int widtha,
		vector<bool> l_grid[], int level);
bool findPathL1(vector<int> l_paths[]);
bool processPartialGoals(int &idx2, int &idx3, int &idx4, int diff1, int diff2, int local_row, int local_col, int k, int widtha, int l_level,
		int block1, vector<bool>l_grid1, vector<bool>l_grid2, bool b_idx4);
void LoadMap(const char *fname, std::vector<bool> &map, int &width, int &height);
void create_auxiliary_maps(vector<bool> map);
bool pathValid(vector<bool> l_grid, int current, vector<bool> visited);
float getDensity(vector<int> tmp_path, int local_level, int offset);
u64 convert_conn_vector(vector<int> tmp_path);
bool findShortestPaths(vector<bool> l_grid, int src, int dst, vector<int> &tmp_path, vector<bool> visited, vector<int> &path, int offset, int local_level, int block);
void callFindShortestPaths(vector<bool> l_grid, int src, int dst, vector<int> &path, int offset, int local_level, int block);
void setup_connectivity();
void preprocess_paths();
void restore_preprocess_paths();

//Auxiliar functions
void printBlock(vector<bool> l_grid);
void printVector(vector<int> path);

/************************************************************************************************************/
/* check if possible moves on grid are legal and stay within bounds, don't cross themselves and are unique	*/
/************************************************************************************************************/
bool pathValid(vector<bool> l_grid, int current, vector<bool> visited) {
  if (current < 0 || current > ROWSxCOLS)    return false; // path invalid if row outside bounds
  if (l_grid[current] == false)   return false; // path invalid if grid is blocked
  if (visited[current] == true) return false; // path invalid if grid location already visited
  return true; // path is valid and path-finding can continue
}

/************************************************************************************************************/
/* Calculates "density": free cells between the distance to the square */
/************************************************************************************************************/
float getDensity(vector<int> tmp_path, int local_level, int offset){
	float local_density=0;
	int gv;

	for(int j = 0; j < tmp_path.size(); j++){
		if (local_level){  											//level 1 dont need conversion from local to global coords
			gv = offset + (tmp_path[j]/ROWS) * width[local_level] + tmp_path[j]%ROWS;
			local_density += density[local_level][gv];
		}else
			local_density += density[local_level][tmp_path[j]];
	}
	local_density = local_density/(tmp_path.size() * tmp_path.size());

	return local_density;
}

/************************************************************************************************************/
/* Calculates "density": free cells between the distance to the square */
/************************************************************************************************************/
float getDensity2(vector<int>path, int local_level){
	float local_density=0;

	//printVector(path);
	for(int j = 0; j < path.size(); j++)
		local_density += density[local_level][path[j]];

	local_density = local_density / (path.size() * path.size()) ;
	//cout << " density: " << local_density << endl;
	return local_density;
}

/************************************************************************************************************/
/* This function review the connectivity between cells inside a block		*/
/************************************************************************************************************/
u64 convert_conn_vector(vector<int> tmp_path){
	int diff=0;
	int lut_diff[] = {1, -1, ROWS, -ROWS};
	bool diff_pos_found = false;
	int diff_pos=0;
	u64 tmp_conn;
	u64 conn_array[tmp_path.size()];
	u64 conn_vector = 0;
	u64 connect = 0;

	if (tmp_path.size() > 1){
		for (int i=0; i < tmp_path.size(); i++){
			diff_pos = 0;
			diff_pos_found = false;
			if (i < tmp_path.size()-1)
			diff = tmp_path[i] - tmp_path[i+1];

			while(!diff_pos_found ){
				if (lut_diff[diff_pos] == diff)
					diff_pos_found = true;
				else
					diff_pos++;
			}
			if (i == 0)
				conn_array[i] = conn_dir[diff_pos][0] << tmp_path[i]*4;
			else if (i == tmp_path.size()-1)
				conn_array[i] = conn_dir[diff_pos][1] << tmp_path[i]*4;
			else{
				connect = conn_dir[diff_pos][0] | tmp_conn;
				conn_array[i] = connect << tmp_path[i]*4;
			}
			tmp_conn = conn_dir[diff_pos][1];
			conn_vector |= conn_array[i];
		}
	}else
		conn_vector = 0xFFFFFFFFFFFFFFFF;
	//cout << "conn_vector: " << conn_vector << endl;
	return conn_vector;
}

/************************************************************************************************************/
/* Recursive function to find the shortest paths in NxN grid												*/
/************************************************************************************************************/
bool findShortestPaths(vector<bool> l_grid, int src, int dst, vector<int> &tmp_path, vector<bool> visited, vector<int> &path, int offset, int local_level, int block) {
	bool assign_path = false;
	if (src == dst) {
		if (local_level == vola_level - 1 ){														// if is the last level we just compare distance
			if ( path.empty() ){
				assign_path = true;
			}else if ( tmp_path.size() < path.size() ){
				assign_path = true;
			}
		} else {
			u64 tmp_conn_vector = convert_conn_vector(tmp_path);
			if ( (tmp_conn_vector & g_conn[local_level][block]) == tmp_conn_vector){				//check connectivity
				if ( path.empty() ){
					assign_path = true;
				}else if ( getDensity(tmp_path, local_level, offset) >= getDensity2(path, local_level) ){			// choose the path less dense
					assign_path = true;
				}
			}
		}

		if (assign_path){
			path.clear();
			if (local_level){  											//distinction from Level 1, level 1 dont need conversion from local to global coords
				for (int j = 0; j < tmp_path.size(); j++) {
					int gv = offset + (tmp_path[j]/ROWS) * width[local_level] + tmp_path[j]%ROWS;
					path.push_back(gv);
				}
			}else
				path = tmp_path;
		}
	}
	// if you havent reached the destination keep exploring grid
	for (int i = 0; i < DIRECTIONS; i++) {
		int ptrn;
		ptrn = src + tests[i];
		if ((src == 3 or src == 7 or src == 11 or src == 15) and i == 0) ptrn = 10000; 							//there is no right neighborh
		if ((src == 0 or src == 4 or src == 8 or src == 12) and i == 1) ptrn = 10000; 							//there is no left neighborh
		if ((src == 12 or src == 13 or src == 14 or src == 15) and i == 2) ptrn = 10000; 						//there is no down neighborh
		if ((src == 0 or src == 1 or src == 2 or src == 3) and i == 3) ptrn = 10000; 							//there is no up neighborh

		if (pathValid(l_grid, ptrn, visited)) { 																// check if path is valid - within bounds, no obstacle and non-cyclic
			visited[ptrn] = true;  																				// mark grid location as already visited
			tmp_path.push_back(ptrn);               															// add cell location to path
			if (!findShortestPaths(l_grid, ptrn, dst, tmp_path, visited, path, offset, local_level, block)) { 	// restart if you didn't find a valid path
				visited[ptrn] = false; 																			// reinitialise visited to false - haven't visited anything yet
				tmp_path.pop_back();                   															// remove the last grid-element you tried to add to the path
			}
		}
	}
	return false; // terminate grid search for this branch
} // findShortestPaths()


/****************************************************************************************************************************************************/
/* Auxiliary function to avoid repetitive operations when looking for the shortest path					*/
/****************************************************************************************************************************************************/
void callFindShortestPaths(vector<bool> l_grid, int src, int dst, vector<int> &path, int offset, int local_level, int block){
    vector<int> tmp_path;
    vector<bool> visited;
    int processed = 0;

    if (l_grid.at(src) and l_grid.at(dst)){
    	if (src == dst){
			if (local_level){  											//distinction from Level 1, level 1 dont need conversion from local to global coords
				int gv = offset + (src/ROWS) * width[local_level] + src%ROWS;
				path.push_back(gv);
			}else
				path.push_back(src);
    	}
    	else{
			visited.resize(ROWSxCOLS);																	// initialise visited points on grid
			tmp_path.begin();                     														// initialise path
			tmp_path.push_back(src);              														// add start coordinate to end of path vector
			visited.at(src) = true; 																	// mark source as visited
			findShortestPaths(l_grid, src, dst, tmp_path, visited, path, offset, local_level, block); 					// find the shortest path from src -> dst
    	}
    }
}

/****************************************************************************************************************************************************/
/* Find a partial path in the last block																											*/
/****************************************************************************************************************************************************/
bool findPartialPathLast(int i, int idx1, int idx2, vector<int> path_la, vector<int> &path_lb, int widthb, int widtha,
		vector<bool> l_grid[], int level){
	int index1, index2;
	bool reverse_path = false;
	int upper_corner_index;
	vector<int>tmp_vector;

	if (idx1 >= idx2){
		index1 = idx1;
		index2 = idx2;
	}else{
		index1 = idx2;
		index2 = idx1;
		reverse_path = true;
	}
	//cout << "findPartialPathLast -- level: " << level << " block: " << path_la[i] << " g_conn[level][block]: " << g_conn[level][path_la[i]] << endl;

	upper_corner_index = (path_la[i]/widtha * widthb * ROWS) + (path_la[i]%widtha * ROWS);			// offset to calculate global coords

#ifdef PREPROCESS
	if ( g_reprocess.at(level)){
		//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@ findPartialPathLast Experiment: " << exp_number << " Reprocessing path @@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		G_LUT[level][path_la[i]][index1][index2].clear();
		callFindShortestPaths(l_grid[path_la[i]], index1, index2, G_LUT[level][path_la[i]][index1][index2], upper_corner_index, level, path_la[i]);
	}
	tmp_vector = G_LUT[level][path_la[i]][index1][index2];
	g_reprocess.at(level) = false;
#else
	callFindShortestPaths(l_grid[path_la[i]], index1, index2, tmp_vector, upper_corner_index, level, path_la[i]);
#endif



#ifdef DEBUG_L2
	cout << "Current Block: " << path_la[i] << endl;
	printBlock(l_grid[path_la[i]]);

	cout << "partial path: " << endl;
	printVector(tmp_vector);

#endif

	// Check if there is a path inside the current block
	if ( !tmp_vector.empty() ){
		if (reverse_path){
			std::reverse(tmp_vector.begin(), tmp_vector.end());
		}

#ifdef DEBUG_L3
		for(int j = 0; j < tmp_vector.size(); j++){
			ofile << "(" << path_lb.size()+j << ", " << level+1 << ", " << tmp_vector[j]%widthb << "," << tmp_vector[j]/widthb << ", select)" << endl;
		}
#endif
		path_lb.insert(path_lb.end(), tmp_vector.begin(), tmp_vector.end());
		return true;
	}else{
#ifdef DEBUG_L3
		if (reverse_path){
			std::reverse(tmp_vector.begin(), tmp_vector.end());
		}

		for(int j = 0; j < tmp_vector.size(); j++){
			ofile << "(" << path_lb.size()+j << ", " << level+1 << ", " << tmp_vector[j]%widthb << "," << tmp_vector[j]/widthb << ", eval)" << endl;
		}
#endif
		return false;
	}
}

/****************************************************************************************************************************************************/
/* Find a partial path																					*/
/****************************************************************************************************************************************************/
bool findPartialPath(int i, int idx1, int idx2, int idx3, int idx4, vector<int> path_la, vector<int> &path_lb,
						int widthb, int widtha, vector<bool> l_grid[], int level, float &local_density, bool &conn_next,
						bool b_idx4){

	int index1, index2, index3, index4;
	bool reverse_path = false;
	bool reverse_path_next = false;
	int upper_corner_index;
	vector<int>tmp_vector, tmp_vector_next;
	float tmp_density=0;

	if (idx1 >= idx2){
		index1 = idx1;
		index2 = idx2;
	}else{
		index1 = idx2;
		index2 = idx1;
		reverse_path = true;
	}

	//cout << "findPartialPath -- level: " << level << " block: " << path_la[i] << " g_conn[level][block]: " << g_conn[level][path_la[i]] << endl;

	upper_corner_index = (path_la[i]/widtha * widthb * ROWS) + (path_la[i]%widtha * ROWS);			// offset to calculate global coords in the LUT
#ifdef PREPROCESS
	if ( g_reprocess.at(level) and i >= g_reprocess_counter.at(level) ){
		//cout << "counter: "<< i << " g_reprocess_counter.at(level): " << g_reprocess_counter.at(level) << " G_LUT[level][path_la[i]][index1][index2].empty(): " << G_LUT[level][path_la[i]][index1][index2].empty() << " g_reprocess.at(level): " << g_reprocess.at(level) << endl;
		//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@ Experiment: " << exp_number << " Reprocessing path @@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		G_LUT[level][path_la[i]][index1][index2].clear();
		callFindShortestPaths(l_grid[path_la[i]], index1, index2, G_LUT[level][path_la[i]][index1][index2], upper_corner_index, level, path_la[i]);
	}
	tmp_vector = G_LUT[level][path_la[i]][index1][index2];
#else
	callFindShortestPaths(l_grid[path_la[i]], index1, index2, tmp_vector, upper_corner_index, level, path_la[i]);
#endif

	if ( !tmp_vector.empty() ) conn_next = true;

	if (idx3 >= idx4){
		index3 = idx3;
		index4 = idx4;
	}else{
		index3 = idx4;
		index4 = idx3;
		reverse_path_next = true;
	}

	upper_corner_index = (path_la[i+1]/widtha * widthb * ROWS) + (path_la[i+1]%widtha * ROWS);			// offset to calculate global coords in the LUT
#ifdef PREPROCESS
	if ( g_reprocess.at(level) and i + 1 >= g_reprocess_counter.at(level)){
		//cout << "counter: "<< i << " g_reprocess_counter.at(level): " << g_reprocess_counter.at(level) << " G_LUT[level][path_la[i+1]][index3][index4].empty(): " << G_LUT[level][path_la[i+1]][index3][index4].empty() << " g_reprocess.at(level): " << g_reprocess.at(level) << endl;
		//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@ Next Experiment: " << exp_number << " Reprocessing path @@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		G_LUT[level][path_la[i+1]][index3][index4].clear();
		callFindShortestPaths(l_grid[path_la[i+1]], index3, index4, G_LUT[level][path_la[i+1]][index3][index4], upper_corner_index, level, path_la[i+1]);
	}
	tmp_vector_next = G_LUT[level][path_la[i+1]][index3][index4];
#else
	callFindShortestPaths(l_grid[path_la[i+1]], index3, index4, tmp_vector_next, upper_corner_index, level, path_la[i+1]);
#endif

#ifdef DEBUG_L2
	cout << "Current Block: " << path_la[i] << endl;
	printBlock(l_grid[path_la[i]]);

	if ( i < path_la.size() - 1 ){
		cout << "Next Block: " << path_la[i+1] << endl;
		printBlock(l_grid[path_la[i+1]]);
	}
#endif

	// Check if there is a path inside the current block and the next one
	if ( !tmp_vector.empty() and !tmp_vector_next.empty() ){
		if (level < vola_level-1 ){
			//Calculate density
			for(int j = 0; j < tmp_vector.size(); j++){
				tmp_density += density[level][tmp_vector[j]];
				//cout << "level: " << level << " cell: " << tmp_vector[j] << " density: " << density[level][tmp_vector[j]] << endl;
			}

			for(int j = 0; j < tmp_vector_next.size(); j++){
				tmp_density += density[level][tmp_vector_next[j]];
				//cout << "level: " << level << " cell: " << tmp_vector_next[j] << " density: " << density[level][tmp_vector_next[j]] << endl;
			}

			tmp_density = tmp_density / ( ( tmp_vector.size()  +  tmp_vector_next.size() ) * ( tmp_vector.size()  +  tmp_vector_next.size() ) );
		}else{				//for last level we calculate distance instead of density
			tmp_density = 10000 - (tmp_vector.size() + tmp_vector_next.size());
		}

#ifdef DEBUG_L2
		cout << "level: " << level << " density: " << tmp_density << " partial path: " << endl;
		printVector(tmp_vector);
		printVector(tmp_vector_next);
#endif

		 if (tmp_density > local_density ){
			if (reverse_path){
				std::reverse(tmp_vector.begin(), tmp_vector.end());
			}

			local_density = tmp_density;

			path_lb.clear();
			path_lb.insert(path_lb.end(), tmp_vector.begin(),  tmp_vector.end());
			//if it is the penultime block it has the info about the last one too
			if(b_idx4){
				if (reverse_path_next){
					std::reverse(tmp_vector_next.begin(), tmp_vector_next.end());
				}
				path_lb.insert(path_lb.end(), tmp_vector_next.begin(),  tmp_vector_next.end());
			}

#ifdef DEBUG_L3
			for(int j = 0; j < tmp_vector.size(); j++){
				ofile << "(" << path_lb.size()+j << ", " << level+1 << ", " << tmp_vector[j]%widthb << "," << tmp_vector[j]/widthb << ", select)" << endl;
			}
#endif
			return true;
		 }else{
#ifdef DEBUG_L3
			if (reverse_path){
				std::reverse(tmp_vector.begin(), tmp_vector.end());
			}

			for(int j = 0; j < tmp_vector.size(); j++){
				ofile << "(" << path_lb.size()+j << ", " << level+1 << ", " << tmp_vector[j]%widthb << "," << tmp_vector[j]/widthb << ", eval)" << endl;
			}
#endif
		}
	}
		return false;
}

/****************************************************************************************************************************************************/
/* Update connectivity
/****************************************************************************************************************************************************/
void update_connectivity(int level, int cell1, int cell2, int &block1, int &block2, int &pos1, int &pos2, bool onOff, int &local_cell1, int &local_cell2){
	int diff_pos1=0, diff_pos2=0;
	int g_cell1 = cell1;
	int g_cell2 = cell2;

	if ( cell1 - cell2 == 1 ) diff_pos1 = 3;
	else if ( cell1 - cell2 == -1 ) diff_pos1 = 1;
	else if ( cell1 - cell2 > 1 ) diff_pos1 = 2;
	else if ( cell1 - cell2 < -1 ) diff_pos1 = 0;

	if ( cell2 - cell1 == 1 ) diff_pos2 = 3;
	else if ( cell2 - cell1 == -1 ) diff_pos2 = 1;
	else if ( cell2 - cell1 > 1 ) diff_pos2 = 2;
	else if ( cell2 - cell1 < -1 ) diff_pos2 = 0;

	if (level == 1){
		block1 = 0;
		block2 = 0;
	}else{
		block1 = ( ( ( cell1/width[level-1])/ROWS ) * width[level-2] ) + ( (cell1%width[level-1]) /ROWS );
		cell1 = ( cell1/width[level-1] - ( block1/width[level-2] * ROWS) ) * ROWS + ( cell1%width[level-1] - ( block1%width[level-2] * ROWS) );
		block2 = ( ( ( cell2/width[level-1])/ROWS ) * width[level-2] ) + ( (cell2%width[level-1]) /ROWS );
		cell2 = ( cell2/width[level-1] - ( block2/width[level-2] * ROWS) ) * ROWS + ( cell2%width[level-1] - ( block2%width[level-2] * ROWS) );
	}

	local_cell1 = cell1;
	local_cell2 = cell2;

#ifdef DEBUG_L2
	//cout << "update connectivity -- level: " << level << " block1: " << block1 << " global cell: " << g_cell1 << " local cell: " << cell1 << " diff_pos1: " << diff_pos1 << endl;
	//cout << "update connectivity -- level: " << level << " block2: " << block2 << " global cell: " << g_cell2 << " local cell: " << cell2 << " diff_pos2: " << diff_pos2 << endl;

	//cout << "before update connectivity g_conn[level-1][block1]: " << g_conn[level-1][block1] << endl;
	//cout << "before update connectivity g_conn[level-1][block2]: " << g_conn[level-1][block2] << endl;
#endif

	pos1 = (cell1*4)+diff_pos1;
	pos2 = (cell2*4)+diff_pos2;
	if(onOff){
		//enable connectivity
		g_conn[level-1][block1] |= (1UL << pos1);
		g_conn[level-1][block2] |= (1UL << pos2);
	}else{
		//disable connectivity
		g_conn[level-1][block1] &= ~(1UL << pos1);
		g_conn[level-1][block2] &= ~(1UL << pos2);

	}

#ifdef DEBUG_L2
	//cout << "after update connectivity g_conn[level-1][block1]: " << g_conn[level-1][block1] << endl;
	//cout << "after update connectivity g_conn[level-1][block2]: " << g_conn[level-1][block2] << endl;
#endif
}

/****************************************************************************************************************************/
/* This function finds an especific part of the path 	*/
/****************************************************************************************************************************/
bool findSecondaryPartialPath(vector<int> l_paths[], int &counter, int &level, vector<int> &old_path, bool conn_next){
	int block1, block2, pos1, pos2;
	int cells[level];
	int base_level = 0;
	int test_level = 2;
	int local_cell1, local_cell2;


	for (int p = 0; p <= level; p++){
#ifdef DEBUG_L2
		printVector(l_paths[p]);
#endif
		g_reprocess.at(p) = false;
	}


	if ( !l_paths[level].empty()){
		cells[level] = l_paths[level][l_paths[level].size() - 1];
		//cout << "cell: " << cells[level] << " level: " << level << endl;
		for ( int c = level; c > 0; c--){
			cells[c-1] = ( ( (cells[c]/width[c])/ROWS ) * width[c-1] ) + ( (cells[c]%width[c]) /ROWS );
			//cout << "cell: " << cells[c-1] << " level: " << c-1 << endl;
		}
	}//else cout << "the path in level: " << level << " is empty" << endl;


	bool b_onOff = false;
	if ( counter == l_paths[level-1].size() - 1 )				// if it did not find a path in the last block
		update_connectivity(level, old_path[0], old_path[1], block1, block2, pos1, pos2, b_onOff, local_cell1, local_cell2 );
	else if ( counter == l_paths[level-1].size() - 2 )
		update_connectivity(level, old_path[1], old_path[2], block1, block2, pos1, pos2, b_onOff, local_cell1, local_cell2 );
	else{
		if (conn_next)
			update_connectivity(level, old_path[2], old_path[3], block1, block2, pos1, pos2, b_onOff, local_cell1, local_cell2 );
		else
			update_connectivity(level, old_path[1], old_path[2], block1, block2, pos1, pos2, b_onOff, local_cell1, local_cell2 );
	}

#ifdef PREPROCESS
	if (level >= 2){
		g_reprocess.at(level-1) = true;
		//cout << "level-2: " << level-2 << " cells[level-2]: " << cells[level-2] << endl;
		//printVector(l_paths[level-2]);
		if ( !l_paths[level].empty()){
			vector <int>::iterator it = find (l_paths[level-2].begin(), l_paths[level-2].end(), cells[level-2]);

			if (it != l_paths[level-2].end()){
			   int nPosition = distance(l_paths[level-2].begin(), it);
			   //cout << "Value "<< *it;
			   //cout << " found in the vector at position: " << nPosition << endl;
			   g_reprocess_counter.at(level-1) = nPosition;
			}
		}else
			g_reprocess_counter.at(level-1) = 0;
	}
#endif


	// delete current paths
	for (int p = base_level; p <= level; p++)
		l_paths[p].erase( l_paths[p].begin(), l_paths[p].end() );

	if (level == 1){
		findPathL1(l_paths);

		if ( !l_paths[0].empty() ){
			counter = -1;
			return true;
		}else
			return false;
	}else{
		if ( comparePaths(l_paths, test_level ) ){
			if (level >= test_level){
				if ( findGlobalPath(test_level, level, l_paths, true ) ){
					counter = -1;
					return true;
				}else
					return false;
			}else{
#ifdef DEBUG_L2
				cout << "level: " << level << endl;
#endif
				level = test_level;
				counter = 100000;
				return true;
			}
		}else
			return false;
	}
}

/****************************************************************************************************************************/
/* Recursive call to find partial path																						*/
/****************************************************************************************************************************/
bool findPathxLevel(vector<int> l_paths[], int startb, int goalb, int widtha, int widthb, vector<bool> l_grid[], int &level,
						bool b_reprocess){
	int diff1, diff2;
	int idx1=0, idx2=0, idx3=0, idx4=0;
	int startb_local, goalb_local, local_row, local_col;
	int right_limit, left_limit, up_limit, down_limit;
	vector<int> old_path, tmp_path_lb;
	bool b_partial_path = true;
	float local_density = 0;
	vector<int> partial_start, partial_goal;
	bool conn_next = false;

	partial_start.resize(widthb);
	partial_goal.resize(widthb);

	//cout << "l_paths[level-1].size(): " << l_paths[level-1].size()  << endl;
	for(int i = 0; i < l_paths[level-1].size() ; i++){
		if ( i > 0)
			idx1 = partial_start[i];
		else{
			//if it is the first iteration we include the real start for this level
			startb_local = ( startb/widthb - ( l_paths[level-1][0]/widtha * ROWS) ) * ROWS + ( startb%widthb - ( l_paths[level-1][0]%widtha * ROWS) );
			idx1 = startb_local;
			//cout << "startb_local: " << startb_local << " startb: " << startb << endl;
		}

		if ( i == l_paths[level-1].size() - 1 ){		// we define the real goal for the last block
			goalb_local = ( ( goalb/widthb - ( (l_paths[level-1][l_paths[level-1].size()-1]/widtha) * ROWS) ) * ROWS ) + ( goalb%widthb - ( (l_paths[level-1][l_paths[level-1].size()-1]%widtha) * ROWS) );
			idx2 = goalb_local;
		}else{
			diff1 = l_paths[level-1][i] - l_paths[level-1][i+1];			//diff1 defines direction based on to next block
			if (i+2 >= l_paths[level-1].size()){
				diff2 = diff1;
			}else{
				diff2 = l_paths[level-1][i+1] - l_paths[level-1][i+2];	//diff2 defines direction for the next iteration based on the next next block
			}

			local_row = idx1/ROWS;
			local_col = idx1%ROWS;
			right_limit = ( (local_row + 1) * ROWS ) -1;
			left_limit = local_row * ROWS;
			up_limit = local_col;
			down_limit = (ROWS-1) * ROWS + local_col;
		}

		bool b_idx4 = false;

		if ( ( i < l_paths[level-1].size() - 1 ) and b_partial_path ){
			b_partial_path = false;
			local_density = 0;

			if (i+2 == l_paths[level-1].size()){
				// if it is the penultimate block, we set the final real goal for the next level
				goalb_local = ( goalb/widthb - ( l_paths[level-1][l_paths[level-1].size()-1]/widtha * ROWS) ) * ROWS + ( goalb%widthb - ( l_paths[level-1][l_paths[level-1].size()-1]%widtha * ROWS) );
				idx4 = goalb_local;
				b_idx4 = true;
			}

			conn_next = false;
			for (int k = 0; k < ROWS; k++){
				if ( processPartialGoals(idx2, idx3, idx4, diff1, diff2, local_row, local_col, k, widtha, level,
						l_paths[level-1][i], l_grid[l_paths[level-1][i]], l_grid[l_paths[level-1][i+1]], b_idx4 ) ) {

#ifdef DEBUG_L2
					cout << "level: " << level << " partial path -- Block: " << l_paths[level-1][i] << " next block: " << l_paths[level-1][i+1] << " partial_start: " << idx1
						<< " partial_goal: " << idx2 << " next start: " << idx3 << " next goal: " << idx4 << " diff1: " << diff1 << " diff2: " << diff2 << endl;
#endif

					if ( findPartialPath( i, idx1, idx2, idx3, idx4, l_paths[level-1], tmp_path_lb, widthb, widtha, l_grid, level, local_density, conn_next, b_idx4) ){
						b_partial_path = true;
						partial_goal[i] = idx2;
						partial_start[i+1] = idx3;
					}
				}
			}

			if (b_partial_path){ // for connectivity
				l_paths[level].insert(l_paths[level].end(), tmp_path_lb.begin(),  tmp_path_lb.end());
				if(b_idx4){
#ifdef DEBUG_L2
					cout << "Level: " << level << " g_reprocess.at(level): " << g_reprocess.at(level) << " g_reprocess_counter.at(level): " << g_reprocess_counter.at(level) << endl;
					cout << "############################## Final path: ######################################################" << endl;
					printVector(l_paths[level]);
#endif
					g_reprocess.at(level) = false;
					break;	//finish the for is not necessary to calculate the path for the last block because it is already done
				}
			}
		}
		if ( ( i == l_paths[level-1].size() - 1 ) and b_partial_path  ){
			if ( findPartialPathLast( i, idx1, idx2, l_paths[level-1], l_paths[level], widthb, widtha, l_grid, level) ) {		// Finding the path for the last block
				b_partial_path = true;
	#ifdef DEBUG_L2
				cout << "Final  Block: " << l_paths[level-1][l_paths[level-1].size() - 1] << " partial start: " << idx1 << " partial goal: " << idx2 << endl;
				cout << "size path: " << l_paths[level].size() << endl;
				// Print detailed path
				printVector(l_paths[level]);
	#endif
			}else
				b_partial_path = false;
		}

		if (!b_partial_path){
#ifdef DEBUG_L2
				cout << " counter: " << i << " Reprocessing l_paths[level-1][i-1] " << l_paths[level-1][i-1] << " l_paths[level-1][i] " << l_paths[level-1][i] << " l_paths[level-1][i+1]: " << l_paths[level-1][i+1] << " l_paths[level-1][i+2]: " << l_paths[level-1][i+2] << endl;
#endif
			if (b_reprocess and ( reprocess[level] < MAX_NUM_REP ) ){
				reprocess[level] += 1;	//general counter for reprocess in N Level, for statistics

				old_path = { l_paths[level-1][i-1], l_paths[level-1][i], l_paths[level-1][i+1], l_paths[level-1][i+2]};

				if ( findSecondaryPartialPath(l_paths, i, level, old_path, conn_next) )
					b_partial_path = true;
				else //	if reprocess
					return false;
			}
			else{
				//cout << "------------------ No reprocess ------------------" << endl;
				return false;
			}
		}
	} // for partial path
	return true;
}


/************************************************************************************************/
/* Recursive function to find the path in the last level (global path)							*/
/************************************************************************************************/
bool findGlobalPath(int initial_level, int final_level, vector<int> l_paths[], bool b_reprocess){
	for ( int c = initial_level; c < final_level; c++ ){

	#ifdef DEBUG_L2
		cout << "+++++++++++++++++ FINDING A PATH IN LEVEL " << c+1 << " ++++++++++++++++++++" << endl;
	#endif

		if ( c == 0){
			if (findPathL1(l_paths))
				continue;
			else
				return false;
		}

		if ( findPathxLevel( l_paths, g_start[c], g_goal[c], width[c-1], width[c], g_grid[c], c, b_reprocess) ){
			if (c == final_level-1)
					return true;
			else
				continue;
		}else
			return false;
	} // for levels
	return true;
}

/************************************************************************************************/
/* main function
/************************************************************************************************/
int main(int argc, char *argv[]) {
	vector<bool> map;
    int width_ll, height_ll;			//width last level
	char mapName[255];
	char scenarioName[255];
	int lut_vl[10] = {0, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144};		// lut to identify vola level
	bool vola_level_found = false;
	int counter_nac=0;		// counter for cases where start or goal are not available (nac:not available cell)

    if (argc != 2) {
   	    cout << "run as: executable <mapFile>" << endl;
   	    cout << "example: tests brc997d.map" << endl;
   	    cout << "All scenarios should be in the scenarios folder" << endl;
   	    cout << "All maps should be in the maps folder, the path is defined in its correspondant scenario file" << endl;
   	    cout << "example for scenarios/testL1.map.scen the route to the map is: maps/dao/testL1.map" << endl;
   	    return 1;
   	}

    boost::timer::cpu_timer timer;

	///////////////////////////////////////////////////////////////////////////////////
	// Load scenarios to solve -- benchmark
	// Bucket	map	map width	map height	start x-coordinate	start y-coordinate	goal x-coordinate	goal y-coordinate	optimal length
	///////////////////////////////////////////////////////////////////////////////////
    sprintf(scenarioName, "scenarios/%s.scen", argv[1]);
	ScenarioLoader *sl = new ScenarioLoader(scenarioName);
	sl->GetNthExperiment(0).GetMapName(mapName);
    LoadMap(mapName, map, width_ll, height_ll);

    // Identifying vola level
    while(!vola_level_found){
    	if (lut_vl[vola_level] == width_ll)
    		vola_level_found = true;
    	else
    		vola_level++;
    }

	g_conn = new u64*[vola_level];
	g_conn_orig = new u64*[vola_level];
	g_grid = new vector<bool>*[vola_level];
	density = new float*[vola_level-1];
	reprocess = new int[vola_level]();
	width = new int[vola_level]();
	max_density = new int[vola_level-1]();
	g_start = new int[vola_level]();
	g_goal = new int[vola_level]();
	g_paths = new vector<int>[vola_level]();

	//no se estan resolviendo todos los casos: 256 345 365 406 411
	//probare escribiendo los valores originales de los path despues de cada experimento
	//de la misma manera que se hace con la connectividad

#ifdef PREPROCESS
	G_LUT = new vector<int>***[vola_level];
	G_LUT_orig = new vector<int>***[vola_level];
#endif

	int tam = 1;
	int width_level = ROWS;
	for(int i = 0; i < vola_level; i++){
		g_conn[i] = new u64[tam]();
		g_conn_orig[i] = new u64[tam]();
		g_grid[i] = new vector<bool>[tam];

#ifdef PREPROCESS
		G_LUT[i] = new vector<int>**[tam];
		G_LUT_orig[i] = new vector<int>**[tam];
		for(int j = 0; j < tam; j++){
			G_LUT[i][j] = new vector<int>*[ROWSxCOLS];
			G_LUT_orig[i][j] = new vector<int>*[ROWSxCOLS];
			for(int k = 0; k < ROWSxCOLS; k++){
				G_LUT[i][j][k] = new vector<int>[ROWSxCOLS];
				G_LUT_orig[i][j][k] = new vector<int>[ROWSxCOLS];
			}
		}
#endif
		width[i] = width_level;
		tam *= ROWSxCOLS;
		width_level *= ROWS;
	}

	tam = ROWSxCOLS;
	for(int i = 0; i < vola_level-1; i++){
		density[i] = new float[tam]();
		max_density[vola_level - 2 - i] = tam;
		tam *= ROWSxCOLS;
	}

    cout << "==================== Loaded scenario and map, vola level: " << vola_level << "====================" << endl;
    ///////////////////////////////////////////////////////////////////////////////////
    // Create auxiliary maps for each level
    ///////////////////////////////////////////////////////////////////////////////////
    create_auxiliary_maps(map);

    ///////////////////////////////////////////////////////////////////////////////////
    // Set up connectivity
    ///////////////////////////////////////////////////////////////////////////////////
    setup_connectivity();

#ifdef PREPROCESS
    ///////////////////////////////////////////////////////////////////////////////////
    // Preprocess the shortest paths and save in LUTs for all levels except Level 1
    ///////////////////////////////////////////////////////////////////////////////////
    preprocess_paths();
#endif


    boost::timer::cpu_times elapsed_boost_1 = timer.elapsed();
	cout << " PREPROCESS CPU TIME: " << ((elapsed_boost_1.user + elapsed_boost_1.system) / 1e9) << " seconds"
		 << " PREPROCESS WALLCLOCK TIME: " << (elapsed_boost_1.wall / 1e9) << " seconds"
		 << endl;

	///////////////////////////////////////////////////////////////////////////////////
	// Solving experiments
	///////////////////////////////////////////////////////////////////////////////////
	int solved = 1;
	int local_src, local_dst;
	boost::timer::cpu_timer timer2;

	for (int x = 0; x < sl->GetNumExperiments(); x++)
    {
		bool b_connectivity = false;
		exp_number = x + 1;
		g_reprocess.clear();
		g_reprocess_counter.clear();

		//Identify starts and goals in all levels
		for ( int c=0; c < vola_level; c++ ){
			g_paths[c].clear();
			reprocess[c] = 0;
			g_reprocess.push_back(false);
			g_reprocess_counter.push_back(100000);

			if ( c == vola_level -1){
				g_start[c] = sl->GetNthExperiment(x).GetStartX() + sl->GetNthExperiment(x).GetStartY() * width[c];
				g_goal[c] = sl->GetNthExperiment(x).GetGoalX() + sl->GetNthExperiment(x).GetGoalY() * width[c];
			}else{
				g_start[c] = (sl->GetNthExperiment(x).GetStartY() / width[vola_level - c - 2] ) * width[c] + ( sl->GetNthExperiment(x).GetStartX() / width[vola_level - c - 2]);
				g_goal[c] = (sl->GetNthExperiment(x).GetGoalY() / width[vola_level - c - 2] ) * width[c] + ( sl->GetNthExperiment(x).GetGoalX() / width[vola_level - c - 2]);
			}
		}

#ifdef DEBUG_L3
		ofile << "{\"Experiment\": " << x + 1 << " \"start\": [" << sl->GetNthExperiment(x).GetStartX() << "," << sl->GetNthExperiment(x).GetStartY() << "] \"goal\": [" << sl->GetNthExperiment(x).GetGoalX() << "," << sl->GetNthExperiment(x).GetGoalY() << "] }" << endl;
#endif

#ifdef DEBUG_L2
		cout << "========================================== Experiment: "<< x + 1 << "  --  Solved: " << solved << " ==========================================" << endl;
		cout << "start_x_ll: " << sl->GetNthExperiment(x).GetStartX() << " start_y_ll: " << sl->GetNthExperiment(x).GetStartY() << endl;
		cout << "goal_x_ll: " << sl->GetNthExperiment(x).GetGoalX() << " goal_y_ll: " << sl->GetNthExperiment(x).GetGoalY() << endl;
		for (int i=0; i<vola_level; i++)
			cout << "Level: " << i+1 << " start: " << g_start[i] << " goal: " << g_goal[i] << endl;
#endif

		local_src = ( g_start[vola_level-1]/width[vola_level-1] - ( g_start[vola_level-2]/width[vola_level-2] * ROWS) ) * ROWS + ( g_start[vola_level-1]%width[vola_level-1] - ( g_start[vola_level-2]%width[vola_level-2] * ROWS) );
		local_dst = ( g_goal[vola_level-1]/width[vola_level-1] - ( g_goal[vola_level-2]/width[vola_level-2] * ROWS) ) * ROWS + ( g_goal[vola_level-1]%width[vola_level-1] - ( g_goal[vola_level-2]%width[vola_level-2] * ROWS) );

		bool b_l1_path = false;
		int test_level = 2;

		if ( comparePaths(g_paths, test_level ) )
			b_l1_path = findGlobalPath(test_level, vola_level, g_paths, true);

		tam = 1;
		for(int i = 0; i < vola_level-1; i++){
			// Restore connectivity
			memcpy( g_conn[i], g_conn_orig[i], sizeof(u64) * tam);
			tam *= ROWSxCOLS;
		}

#ifdef PREPROCESS
		//restore_preprocess_paths();
#endif


#ifdef DEBUG_L2
		cout << endl << "Map: " << mapName << " List of failed experiments:" << failed_exp.size() << endl;
		printVector(failed_exp);
#endif

#ifdef DEBUG_L3
		for (int m = 0; m < vola_level; m++)
			ofs << "Map: " << mapName << " :Experiment: " << exp_number << " :Level: " << m << " :reprocess: " << reprocess[m] << endl;
#endif

		if (b_l1_path){
#ifdef DEBUG_L2
			cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SOLVED!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
#endif
			total_path_length += g_paths[vola_level-1].size();
			solved++;
			continue;
		}else{
#ifdef DEBUG_L2
			cout << "########################################################## FAIL!!! ############################################################" << endl;
#endif
			failed_exp.push_back(x+1);
		}

    }// for experiments

	boost::timer::cpu_times elapsed_boost_2 = timer2.elapsed();

	cout << "Map: " << mapName << " Solved: " << solved << " of: "<< sl->GetNumExperiments() << ", %solved: " << (((float)solved/sl->GetNumExperiments()) * 100) <<  endl;

	cout << " PATH CPU TIME: " << ((elapsed_boost_2.user + elapsed_boost_2.system) / 1e9) << " seconds"
		          << " PATH WALLCLOCK TIME: " << (elapsed_boost_2.wall / 1e9) << " seconds"
		          << endl;

	cout << " AVG PATH CPU TIME: " << ((elapsed_boost_2.user + elapsed_boost_2.system) / 1e9) / sl->GetNumExperiments()  << " seconds"
		          << " AVG PATH WALLCLOCK TIME: " << (elapsed_boost_2.wall / 1e9) / sl->GetNumExperiments()<< " seconds"
		          << " AVG PATH LENGTH: " << total_path_length/solved
		          << endl;

	cout << " List of failed experiments:" << endl;
	printVector(failed_exp);


	//==clean==
	for(int i = 0; i < vola_level; i++){
		delete[] g_conn[i];
		delete[] g_conn_orig[i];
		delete[] g_grid[i];
//#ifdef PREPROCESS
	//	delete[] G_LUT[i];
	}
	delete[] g_conn;
	delete[] g_conn_orig;
	delete[] g_grid;

	for(int i = 0; i < vola_level-1; i++){
		delete[] density[i];
	}
	delete[] density;

	delete[] reprocess;
	delete[] width;
	delete[] max_density;
	delete[] g_start;
	delete[] g_goal;
	delete[] g_paths;

	while (1);
} // main()

/****************************************************************************************************************************************************/
/* Compare two path until a level that you can indicate as a parameter																				*/
/****************************************************************************************************************************************************/
bool comparePaths(vector<int> l_paths[], int test_level){
	vector<int> tmp_path1[test_level];
	vector<int> tmp_path2[test_level];
	float tmp_density1 = 0, tmp_density2 = 0;
	int block1, block2, pos1, pos2;
	int cell_avoid;
	float dense_cell=10000000;
	bool b_onOff = false;
	int local_cell1, local_cell2;

	if ( findGlobalPath(0, test_level, tmp_path1, true) ){												//AL1-N
		tmp_density1 = getDensity2(tmp_path1[test_level-1], test_level-1);
#ifdef DEBUG_L2
		cout << "tmp_density1: " << tmp_density1 << endl;
#endif
		if (tmp_path1[0].size() == 1){
			for(int a = 0; a < test_level; a++)
				l_paths[a] = tmp_path1[a];
			return true;
		}else{
			if (tmp_path1[0].size() == 2)
				cell_avoid = 1;
			else{
				//Find the most dense cell
				for(int j = 0; j < tmp_path1[0].size(); j++){
					if( j !=0 and j != tmp_path1[0].size() - 1  ){
						if ( density[0][tmp_path1[0][j]] < dense_cell){
							dense_cell = density[0][tmp_path1[0][j]];
							cell_avoid = j;
						}
					}
				}
			}
			// change connectivity
			update_connectivity(1, tmp_path1[0][cell_avoid-1], tmp_path1[0][cell_avoid], block1, block2, pos1, pos2, b_onOff, local_cell1, local_cell2 );

			if ( findGlobalPath(0, test_level, tmp_path2, false) )										//BL1-N
				tmp_density2 = getDensity2(tmp_path2[test_level-1], test_level-1);
			else tmp_density2 = 0;
#ifdef DEBUG_L2
			cout << "tmp_density2: " << tmp_density2 << endl;
#endif

			// restore connectivity
			g_conn[0][block1] |= (1UL << pos1);
			g_conn[0][block2] |= (1UL << pos2);

#ifdef DEBUG_L2
			cout << "restore connectivity g_conn[0][0]: " << g_conn[0][0] << endl;
			cout << "tmp_density1: " << tmp_density1 <<  " tmp_density2: " << tmp_density2 << endl;
#endif

			if ( tmp_density1 > tmp_density2 ){										// ALN vs BLN
				for(int a = 0; a < test_level; a++){
					//printVector(tmp_path1[a]);
					l_paths[a] = tmp_path1[a];
				}
			}else {
				for(int a = 0; a < test_level; a++)
					l_paths[a] = tmp_path2[a];
			}
			return true;
		}
	}else
		return false;
}/*

/****************************************************************************************************************************************************/
/* Find a path in level 1																					*/
/****************************************************************************************************************************************************/
bool findPathL1(vector<int> l_paths[]){
	int index1_l1, index2_l1;
	bool reverse_path_l1 = false;
	if (g_start[0] >= g_goal[0]){
		index1_l1 = g_start[0];
		index2_l1 = g_goal[0];
		reverse_path_l1 = false;
	}else{
		index1_l1 = g_goal[0];
		index2_l1 = g_start[0];
		reverse_path_l1 = true;
	}

	callFindShortestPaths(g_grid[0][0], index1_l1, index2_l1, l_paths[0], 0, 0, 0);					//AL1

	if ( !l_paths[0].empty()){
		if (reverse_path_l1){
			std::reverse(l_paths[0].begin(), l_paths[0].end() );
		}
#ifdef DEBUG_L2
		printVector(l_paths[0]);
		cout << "Density: " << getDensity2(l_paths[0], 0) << endl;
#endif
		return true;
	}else return false;
}

/************************************************************************************************************/
/* Load the original map																					*/
/************************************************************************************************************/
void LoadMap(const char *fname, std::vector<bool> &map, int &width, int &height)
{
	FILE *f;
	f = fopen(fname, "r");
	if (f)
    {
		fscanf(f, "type octile\nheight %d\nwidth %d\nmap\n", &height, &width);
		map.resize(height*width);
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				char c;
				do {
					fscanf(f, "%c", &c);
				} while (isspace(c));
				map[y*width+x] = (c == '.' || c == 'G' || c == 'S');
			}
		}
		fclose(f);
    }
}

/************************************************************************************************/
/* Create auxiliary maps for each level															*/
/************************************************************************************************/
void create_auxiliary_maps(vector<bool> map){
	int offset_x = 0;
    int offset_y = 0;
    int index_la = 0;	// upper level
    int index_lb = 0;	// lower level
    float bit_value = 0;
    int cell = 0;
    int offset=0;
    int free_cell_threshold = 1;

#ifdef DEBUG_L3
	FILE *pFile[10];
	char maps[10][255] = {"map_l1.txt", "map_l2.txt", "map_l3.txt", "map_l4.txt", "map_l5.txt", "map_l6.txt", "map_l7.txt", "map_l8.txt", "map_l9.txt", "map_l10.txt"};

	for (int k=0; k < vola_level-1; k++){
		pFile[k] = fopen (maps[k],"w");
		if (pFile[k]==NULL){
		    cout << "error opening file: " << maps[k] << endl;
		    return ;
		}
	}
#endif

    cout << "==================== Initializing Last Level: " << vola_level << "								====================" << endl;

    for(int itr = 0; itr < width[vola_level - 2]; itr++){
       	for(int it = 0; it < width[vola_level - 2]; it++){
        	index_la = itr*width[vola_level - 2] + it;
        	bit_value = 0;
        	//cout << "index_la: " << index_la << endl;
        	g_grid[vola_level - 1][index_la].resize(ROWSxCOLS);
    		for (int row = 0 + offset_y, row_g=0; row < ROWS + offset_y; row++, row_g++) {
    			for (int col = 0 + offset_x, col_g=0 ; col < COLS + offset_x; col++, col_g++) {
    				g_grid[vola_level - 1][index_la].at(col_g + row_g * ROWS) = map.at(row*width[vola_level - 1]+col);
    				bit_value += g_grid[vola_level - 1][index_la].at(col_g + row_g * ROWS) ? 1 : 0;
    				//cout << "index_l5: " << index_la << " idx: " << col_g + row_g * ROWS << " bit_value: " << bit_value << endl;
    			}
    		}
    		//cout << "before normalizing -- density: " << bit_value << " max_density[vola_level-2]: " << max_density[vola_level-2] << endl;

    		density[vola_level-2][index_la] = bit_value / max_density[vola_level-2];
    		//cout << "density: " << density[vola_level -2][index_la] << endl;
    		offset_x += 4;
    	}
    	offset_x = 0;
    	offset_y += 4;
    }

    for (int i = (vola_level - 3); i >= 0; i--){
        cout << "==================== Initializing Level: " << i+2 << "								====================" << endl;
        bit_value=0;
        offset_x = 0;
        offset_y = 0;

        for(int itr = 0; itr < width[i]; itr++){
           	for(int it = 0; it < width[i]; it++){
            	index_lb = itr*width[i] + it;
            	g_grid[i+1][index_lb].resize(ROWSxCOLS);
    			//cout << "block_l4: " << index_lb << endl;
    			for (int row=0; row<ROWS; row++) {
    				for (int col=0; col<COLS; col++) {
    					cell = col + row * ROWS;
    					index_la =  col + offset_x + ( (row + offset_y) * width[i+1]);
    					//cout << "index_l4: " << index_la << " cell_l4: " << cell;
    					for (int row2=0; row2<ROWS; row2++) {
    						for (int col2=0; col2<COLS; col2++) {
    							bit_value += g_grid[i+2][index_la][col2 + row2 * ROWS] ? 1 : 0;
    						}
    					}
    					//cout << " bit_value: " << bit_value;

    					if (i == (vola_level - 3) )
    						density[i][index_lb] += bit_value;
    					else
    						density[i][index_lb] += density[i+1][index_la];

    					g_grid[i+1][index_lb][cell] = bit_value >= free_cell_threshold ? true : false;
    					bit_value = 0;
    					//cout << endl;
    				}
    			}
    			//cout << "before normalizing -- density: " << density[i][index_lb] << " max_density[i]: " << max_density[i] << endl;
    			if (i == (vola_level - 3) )
    				density[i][index_lb] /= max_density[i];
    			else
    				density[i][index_lb] /= ROWSxCOLS;
    			//cout << "density: " << density[i][index_lb] << endl;
    			offset_x += 4;
    		}
        	offset_x = 0;
        	offset_y += 4;
        }

#ifdef DEBUG_L3
		//Save to a file the map to check that is ok
        offset=0;
		for (int j=0; j< width[i]; j++) {
			for (int row=0; row<ROWS; row++) {
				for(int it = 0; it < width[i]; it++){
					for (int col=0; col<COLS; col++) {
						if(g_grid[i+1][it + offset][col + row * ROWS])
							fputc ( '.' , pFile[i+1] );
						else
							fputc ( '@' , pFile[i+1] );
					}
				}
				fputc ( '\n' , pFile[i+1] );
			}
			offset += width[i];
		}
#endif
    }

    cout << "==================== Initializing Level: 1								====================" << endl;
    bit_value = 0;

    g_grid[0][0].resize(ROWSxCOLS);

	for (int row=0; row<ROWS; row++) {
		for (int col=0; col<COLS; col++) {
			cell = col + row * ROWS;
			//cout << " cell_l1: " << cell;
			for (int row2=0; row2<ROWS; row2++) {
				for (int col2=0; col2<COLS; col2++) {
					bit_value += g_grid[1][cell][col2 + row2 * ROWS] ? 1 : 0;
				}
			}
			//cout << " bit_value: " << bit_value;
			g_grid[0][0][cell] = bit_value >= free_cell_threshold ? true : false;
			bit_value = 0;
			//cout << endl;
		}
	}

#ifdef DEBUG_L3
	//Save to a file the map to check that is ok
   	for (int row=0; row<ROWS; row++) {
   		for (int col=0; col<COLS; col++) {
   			if(g_grid[0][0][col + row * ROWS])
   			   	fputc ( '.' , pFile[0] );
			else
				fputc ( '@' , pFile[0] );
    	}
   		fputc ( '\n' , pFile[0] );
    }

	for (int k=0; k < vola_level-1; k++){
		fclose ( pFile[k] );
	}
#endif
}

/************************************************************************************************/
/* Setup connectivity vector for each level														*/
/************************************************************************************************/
void setup_connectivity(){
	int block_b, offset_row=0, offset_col=0, pos1, pos2, block_a1, block_a2;
	bool b_onOff = true;
	u64	r_bit1 = 0, r_bit2 = 0, r_bit3 = 0, r_bit4 = 0;
	u64 c_bit1 = 0, c_bit2 = 0, c_bit3 = 0, c_bit4 = 0;
	int local_cell1, local_cell2;

	r_bit1 |= (1UL << 13), r_bit2 |= (1UL << 29), r_bit3 |= (1UL << 45), r_bit4 |= (1UL << 61);
	c_bit1 |= (1UL << 48), c_bit2 |= (1UL << 52), c_bit3 |= (1UL << 56), c_bit4 |= (1UL << 60);

	//int local_level = vola_level - 1;
	for(int local_level = vola_level - 1; local_level > 0 ; local_level--){	// 3, 2, 1
		for (int row = 0; row < width[local_level-1]; row++ ){ //0 -> 63 for L3
			for ( int col = 0, col2 = 0; col < width[local_level-1]-1; col++, col2 += width[local_level-1]){ //0 -> 4 for L3
				block_b = offset_row + col;	// 0 - 4096
				//cout << "local_level: " << local_level << " block_b: " << block_b << endl;
				if (local_level == vola_level -1){
					if ( (g_grid[local_level][block_b].at(3) and g_grid[local_level][block_b+1].at(0)) or
						(g_grid[local_level][block_b].at(7) and g_grid[local_level][block_b+1].at(4)) or
						(g_grid[local_level][block_b].at(11) and g_grid[local_level][block_b+1].at(8)) or
						(g_grid[local_level][block_b].at(15) and g_grid[local_level][block_b+1].at(12)) ){

						//printBlock(g_grid[local_level][block_b]);
						//printBlock(g_grid[local_level][block_b+1]);
						update_connectivity(local_level, block_b, block_b+1, block_a1, block_a2, pos1, pos2, b_onOff, local_cell1, local_cell2 );
					}
				}else{
					if ( (g_conn[local_level][block_b] & r_bit1 ) == r_bit1 or
							(g_conn[local_level][block_b] & r_bit2 ) == r_bit2 or
							(g_conn[local_level][block_b] & r_bit3 ) == r_bit3 or
							(g_conn[local_level][block_b] & r_bit4 ) == r_bit4 )

							update_connectivity(local_level, block_b, block_b+1, block_a1, block_a2, pos1, pos2, b_onOff, local_cell1, local_cell2 );
				}

				block_b = row + col2;	// 0 - 4096
				//cout << "Updating cols local_level: " << local_level << " block_b: " << block_b << " block_b+width[local_level-1]:" << block_b+width[local_level-1] <<endl;
				if (local_level == vola_level -1){
					if ( (g_grid[local_level][block_b].at(12) and g_grid[local_level][block_b+width[local_level-1]].at(0)) or
						(g_grid[local_level][block_b].at(13) and g_grid[local_level][block_b+width[local_level-1]].at(1)) or
						(g_grid[local_level][block_b].at(14) and g_grid[local_level][block_b+width[local_level-1]].at(2)) or
						(g_grid[local_level][block_b].at(15) and g_grid[local_level][block_b+width[local_level-1]].at(3)) ){

						//printBlock(g_grid[local_level][block_b]);
						//printBlock(g_grid[local_level][block_b+width[local_level-1]]);
						update_connectivity(local_level, block_b, block_b+width[local_level-1], block_a1, block_a2, pos1, pos2, b_onOff, local_cell1, local_cell2 );
					}
				}else{
					//no esta actualizando la conectividad en la ultima columna en el ultimo nivel
					if ( (g_conn[local_level][block_b] & c_bit1 ) == c_bit1 or
							(g_conn[local_level][block_b] & c_bit2 ) == c_bit2 or
							(g_conn[local_level][block_b] & c_bit3 ) == c_bit3 or
							(g_conn[local_level][block_b] & c_bit4 ) == c_bit4 )

							update_connectivity(local_level, block_b, block_b+width[local_level-1], block_a1, block_a2, pos1, pos2, b_onOff, local_cell1, local_cell2 );
				}
			}
			offset_row += width[local_level-1];
			//cout << "offset_row: " << offset_row << endl;
		}
		offset_row = 0;
	}
	//cout << "c_bit1: " << c_bit1 << " c_bit2: " << c_bit2 << " c_bit3: " << c_bit3 << " c_bit4: " << c_bit4 << endl;
	//cout << "r_bit1: " << r_bit1 << " r_bit2: " << r_bit2 << " r_bit3: " << r_bit3 << " r_bit4: " << r_bit4 << endl;
	int tam = 1;
	for(int i = 0; i < vola_level-1; i++){
		// Backup connectivity
		memcpy( g_conn_orig[i], g_conn[i], sizeof(u64) * tam);
		tam *= ROWSxCOLS;
	}
}

// LUT to look for connectivity
static int rr_lut[4][4] = { 	// right-right lut
		{3, 7, 11, 15},
		{7, 3, 11, 15},
		{11, 7, 15, 3},
		{15, 11, 7, 3}
};

static int ll_lut[4][4] = { 	// left-left lut
		{0, 4, 8, 12},
		{4, 0, 8, 12},
		{8, 4, 12, 0},
		{12, 8, 4, 0}
};

static int uu_lut[4][4] = { 	// up-up lut
		{0, 1, 2, 3},
		{1, 0, 2, 3},
		{2, 1, 3, 0},
		{3, 2, 1, 0}
};

static int dd_lut[4][4] = { 	// down-down lut
		{12, 13, 14, 15},
		{13, 12, 14, 15},
		{14, 13, 15, 12},
		{15, 14, 13, 12}
};


/************************************************************************************************/
/* Process partial goals 																		*/
/************************************************************************************************/
bool processPartialGoals(int &idx2, int &idx3, int &idx4, int diff1, int diff2, int local_row, int local_col, int k, int widtha, int l_level,
		int block1, vector<bool>l_grid1, vector<bool>l_grid2, bool b_idx4){
	int left=1, right=-1, up=widtha, down=-widtha, ne=widtha-1, nw=widtha+1, se=-widtha-1, sw=-widtha+1;
	int next_row, next_col, pos1;
	u64 bit=0;

	// LEFT
	if (diff1 == left and diff2 == left){
		idx2 = ll_lut[local_row][k];
		idx3 = idx2 + 3;
		next_row = idx3/ROWS;
		pos1 = 3;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = ll_lut[next_row][l];
				if(l_grid2[idx4])
					break;
			}
		}
	}
	else if (diff1 == left and diff2 == up){
		idx2 = (k * ROWS);
		idx3 = (ROWS-1) + (k * ROWS);
		pos1 = 3;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = (ROWS-1) - l;
				if(l_grid2[idx4])
					break;
			}
		}
	}
	else if (diff1 == left and diff2 == down){
		idx2 = (ROWS-1)*ROWS - (k * ROWS);
		idx3 = ROWS*ROWS -1 - (k * ROWS);
		pos1 = 3;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = ROWS*ROWS - l - 1;
				if(l_grid2[idx4])
					break;
			}
		}
	}
	// RIGHT
	else if (diff1 == right and diff2 == right){
		idx2 = rr_lut[local_row][k];
		idx3 = idx2 - 3;
		next_row = idx3/ROWS;
		pos1 = 1;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = rr_lut[next_row][l];
				if(l_grid2[idx4])
				break;
			}
		}
	}
	else if (diff1 == right and diff2 == up){
		idx2 = (ROWS-1) + (k * ROWS);
		idx3 = (k * ROWS);
		pos1 = 1;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = l;
				if(l_grid2[idx4])
					break;
			}
		}
	}
	else if (diff1 == right and diff2 == down){
		idx2 = ROWS*ROWS -1 - (k * ROWS);
		idx3 = (ROWS-1)*ROWS - (k * ROWS);
		pos1 = 1;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = (ROWS-1)*ROWS + l;
				if(l_grid2[idx4])
					break;
			}
		}
	}
	//UP
	else if (diff1 == up and diff2 == up){
		idx2 = uu_lut[local_col][k];
		idx3 = idx2 + 12;
		next_col = idx3%ROWS;
		pos1 = 2;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = uu_lut[next_col][l];
				if(l_grid2[idx4])
					break;
			}
		}
	}
	else if (diff1 == up and diff2 == left){
		idx2 = k;
		idx3 = (ROWS-1)*ROWS + k;
		pos1 = 2;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = (ROWS-1)*ROWS - (l * ROWS);
				if(l_grid2[idx4])
					break;
			}
		}
	}
	else if (diff1 == up and diff2 == right){
		idx2 = (ROWS-1) - k;
		idx3 = (ROWS*ROWS)-1 - k;
		pos1 = 2;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = (ROWS*ROWS)-1 - (l * ROWS);
				if(l_grid2[idx4])
					break;
			}
		}
	}
	// DOWN
	else if (diff1 == down and diff2 == down){
		idx2 = dd_lut[local_col][k];
		idx3 = idx2 - 12;
		next_col = idx3%ROWS;
		pos1 = 0;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = dd_lut[next_col][l];
				if(l_grid2[idx4])
					break;
			}
		}
	}
	else if (diff1 == down and diff2 == left){
		idx2 = (ROWS-1)*ROWS + k;
		idx3 = k;
		pos1 = 0;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = l * ROWS;
				if(l_grid2[idx4])
					break;
			}
		}
	}
	else if (diff1 == down and diff2 == right){
		idx2 = (ROWS*ROWS)-1 - k;
		idx3 = (ROWS-1) - k;
		pos1 = 0;
		if(!b_idx4){
			for (int l=0; l<ROWS; l++){
				idx4 = (ROWS-1) + (l * ROWS);
				if(l_grid2[idx4])
					break;
			}
		}
	}

	if ( l_grid1[idx2] and l_grid2[idx3] ){
		if (l_level == vola_level -1)
			return true;
		else{
			bit |= (1UL << idx2*4+pos1);
#ifdef DEBUG_L2
			cout << "bit: " << bit << " g_conn[l_level][block1]: " << g_conn[l_level][block1] << " g_conn[l_level][block1] and bit: " << (g_conn[l_level][block1] & bit) <<" level: " << l_level << " block: " << block1 << " idx2: " << idx2 << endl;
#endif
			if ( (g_conn[l_level][block1] & bit ) == bit)
				return true;
			else
				return false;
		}
	}else
		return false;
}


/************************************************************************************************************/
/* Auxiliar function to display a 4x4 block																	*/
/************************************************************************************************************/
void printBlock(vector<bool>l_grid){
	for (int row=0; row<ROWS; row++) {
		for (int col=0; col<COLS; col++) {
			if(l_grid[col + row * ROWS])
			cout << ".";
		else
			cout << "T";
		}
		cout << endl;
	}
	cout << endl;
}

/************************************************************************************************************/
/* Auxiliar function to display a path (vector)																*/
/************************************************************************************************************/
void printVector(vector<int> path){
	for(int j = 0; j < path.size(); j++){
		cout << " " << path[j] ;
	}
	cout << endl ;
}


/************************************************************************************************/
/* Preprocess the shortest paths and save in LUTs for all levels 								*/
/************************************************************************************************/
void preprocess_paths(){
	int upper_corner_index;
	int nbr_blocks = ROWSxCOLS;

	for (int i=1; i < vola_level; i++)
	{
		cout << "==================== Processing the shortest paths for level " << i+1 << "						====================" << endl;
		// the shortest paths are saved in LUT in goblal coords
		for(int index = 0; index < nbr_blocks; index++){
			for (int src=0; src < ROWSxCOLS; src++) {
				for (int dst = 0; dst <= src; dst++) {
					if ( g_grid[i][index].at(src) and g_grid[i][index].at(dst) ){ 							// Checking that start and goal are free cells
						upper_corner_index = (index/width[i-1] * width[i] * ROWS) + (index%width[i-1] * ROWS);		// offset to calculate global coords in the LUT
						callFindShortestPaths(g_grid[i][index], src, dst, G_LUT[i][index][src][dst], upper_corner_index, i, index);
						G_LUT_orig[i][index][src][dst] = G_LUT[i][index][src][dst];
					}
				}
			}
		}
		nbr_blocks *= ROWSxCOLS;
	}
}

/************************************************************************************************/
/* Preprocess the shortest paths and save in LUTs for all levels 								*/
/************************************************************************************************/
void restore_preprocess_paths(){
	int nbr_blocks = ROWSxCOLS;

	for (int i=1; i < vola_level; i++)
	{
		for(int index = 0; index < nbr_blocks; index++){
			for (int src=0; src < ROWSxCOLS; src++) {
				for (int dst = 0; dst <= src; dst++) {
						G_LUT[i][index][src][dst] = G_LUT_orig[i][index][src][dst];
				}
			}
		}
		nbr_blocks *= ROWSxCOLS;
	}
}
