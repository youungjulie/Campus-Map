# EE 538 Final Project Report

## TrojanMap

This project focuses on using data structures in C++ and implementing various graph algorithms to build a map application.

<p align="center"><img src="img/TrojanMap.png" alt="Trojan" width="500" /></p>

## Data Structure

To create the map, every location in the map is represented by Node. Every node is connected using graph data structure. 

### Classes

**Node**

Attributes:  

id: A unique id assign to each point  
lat: Latitude of this location  
lon: Longitude of this location  
name: Name of this location  
neighbors: List of the ids of all neighbor points  
attributes: List of the attributes of the location  

```cpp
class Node {
  public:
    Node(){};
    Node(const Node &n){id = n.id; lat = n.lat; lon = n.lon; name = n.name; neighbors = n.neighbors; attributes = n.attributes;};
    std::string id;    // A unique id assign to each point
    double lat;        // Latitude
    double lon;        // Longitude
    std::string name;  // Name of the location. E.g. "Bank of America".
    std::vector<std::string> neighbors;  // List of the ids of all neighbor points.
    std::unordered_set<std::string> attributes;  // List of the attributes of the location.
};
```

**Trojan Map**

Attributes:  

data: A map stores the id of location as key and its Node as value

*Function below are implementations of the map*


## Phase 1 report  

**GetLat()**  
Get the Latitude of a Node given its id  
- Time Complexity: O(1)  
- Space Complexity: O(1)  

**GetLon()**  
Get the Longitude of a Node given its id  
- Time Complexity: O(1)  
- Space Complexity: O(1)  

**GetName()**  
Get the name of a Node given its id  
- Time Complexity: O(n)  
- Space Complexity: O(1)  

**GetID**  
Get the id given its name  
- Time Complexity: O(1)  
- Space Complexity: O(1)  

**GetNeighborIDs**  
Get the neighbor ids of a Node  
- Time Complexity: O(1)  
- Space Complexity: O(1)  

**GetPosition**  
Given a location name, return its latitude and longitude. If id does not exist, return (-1, -1)    
- Time Complexity: O(1)  
- Space Complexity: O(1)  

**CalculateEditDistance**  
Calculate edit distance between two location names, that is how many substitution, deletion, and addition to modify one word to another.   
- Time Complexity: O(n*m)  
- Space Complexity: O(n*m)  
*n and m are the length of two words*  

**Autocomplete**  
Returns a vector of names given a partial name.  
- Time Complexity: O(n)  
- Space Complexity: O(n)  

Run time testing:  
![1](https://user-images.githubusercontent.com/63425702/166177570-711a82da-29ed-444a-9c20-22c09dd49183.png)

**FindClosestName**  
Iterate through the map and find the name with smallest edit distance.  
- Time Complexity: O(n * m * d)  
- Space Complexity: O(n * m * d)  
*n and m are the length of two words, d is the number of all data nodes*  

![image](https://user-images.githubusercontent.com/63425702/166177925-45bafe25-26ca-4e6c-b3aa-5d264865c557.png)


## Phase 2 report

**CalculateShortestPath_Dijkstra**

This function is using Dijkstra algorithm to calculate shortest distance from `location1_name` to `location2_name`.
```cpp
std::vector<std::string> TrojanMap::CalculateShortestPath_Dijkstra(std::string location1_name, std::string location2_name);
```
To realize this function, we use data structure below:
```cpp
std::vector<std::string> path;
std::priority_queue<std::pair<double, std::string>, std::vector<std::pair<double, std::string>> , std::greater<std::pair<double, std::string>> > pq;
std::unordered_map<std::string, double> shortest_map;
std::unordered_map<std::string, std::string> predecessor_map;
```
- Time Complexity: O(E + VlogV), where E is the number of edges and V is the number of vertices. 

Run time testing:  
- Input: `Ralphs`, `Target`
- Output: "2578244375","4380040154","4380040158","4380040167","6805802087","8410938469","6813416131","7645318201","6813416130","6813416129","123318563","452688940","6816193777","123408705","6816193774","452688933","452688931","123230412","6816193770","6787470576","4015442011","6816193692","6816193693","6816193694","4015377691","544693739","6816193696","6804883323","6807937309","6807937306","6816193698","4015377690","4015377689","122814447","6813416159","6813405266","4015372488","4015372487","6813405229","122719216","6813405232","4015372486","7071032399","4015372485","6813379479","6813379584","6814769289","5237417650",
The distance of the path is:0.927969 miles  
- Time taken by function: 118 ms  


**CalculateShortestPath_Bellman_Ford**

This function is using Bellman_Ford algorithm to calculate shortest distance from `

1_name` to `location2_name`.
```cpp
std::vector<std::string> TrojanMap::CalculateShortestPath_Bellman_Ford(std::string location1_name, std::string location2_name);
```
To realize this function, we use data structure below:
```cpp
std::vector<std::string> path;
std::unordered_map<std::string, std::vector<std::string>> neighbor_map;
std::unordered_map<std::string, double> shortest_map;
std::unordered_map<std::string, std::string> predecessor_map;
```
- Time Complexity: O(E*V), where E is the number of edges and V is the number of vertices. 

Run time testing:  
- Input: `Ralphs`, `Target`  
- Output: "2578244375","4380040154","4380040158","4380040167","6805802087","8410938469","6813416131","7645318201","6813416130","6813416129","123318563","452688940","6816193777","123408705","6816193774","452688933","452688931","123230412","6816193770","6787470576","4015442011","6816193692","6816193693","6816193694","4015377691","544693739","6816193696","6804883323","6807937309","6807937306","6816193698","4015377690","4015377689","122814447","6813416159","6813405266","4015372488","4015372487","6813405229","122719216","6813405232","4015372486","7071032399","4015372485","6813379479","6813379584","6814769289","5237417650",
The distance of the path is:0.927969 miles 
- Time taken by function: 7860 ms 
- Graph:
<img width="708" alt="trojanmap_shortestpath" src="https://user-images.githubusercontent.com/98196892/164145007-53d430d8-d53a-460b-8e27-ef9df43422b0.png">

**Summary: Run time based on input and methods**  
input | Dijkstra(/ms) | Bellman-Ford(/ms)
--- | --- | --- 
Target,  Ralphs | 84 | 1835
CAVA, USC Parking | 39 | 1610
Food 4 Less, USC Village Gym | 72 | 2632
7-Eleven, Ralphs | 149 | 1645
Chick-fil-A, Target | 63 | 1607

**CycleDetection**  
This function is to detect whether there exists a cycle in given area.  
```cpp
bool TrojanMap::CycleDetection(std::vector<std::string> &subgraph, std::vector<double> &square);
```
We used a helper function to detect whether there is a cycle:  
```cpp
bool TrojanMap::hasCycle(std::string current_id, std::map<std::string, bool> &visited, 
                std::string parent_id, std::vector<double> &square);
```
- Time Complexity: O(n)

Run time testing:  
- Input: `-118.299`,`-118.264`,`34.030`,`34.026` 
- Output: there exists a cycle in the subgraph 
- Time taken by function: 16 ms 
- Graph:
<img width="708" alt="cycle" src="https://user-images.githubusercontent.com/98196892/166170586-1769995f-b23d-416d-96f6-0985ed6afdc8.png">

**Topological Sort**  
This function uses Depth first search algorithm to sort the locations `input: location_names` based on their `dependencies`and return the sorted result.   

```cpp
std::vector<std::string> DeliveringTrojan(std::vector<std::string> &location_names,
                                            std::vector<std::vector<std::string>> &dependencies);
```

Run time testing:
```cpp
Input: 
location_names = {"Ralphs", "Chick-fil-A", "KFC", "Target", "CAVA"}
dependencies = {{ "Ralphs", "Chick-fil-A" }, { "Ralphs", "KFC" }, 
                                        { "Chick-fil-A", "KFC" }, { "Target", "CAVA" }, { "CAVA", "KFC" }}           
output: Target -> CAVA -> Ralphs -> Chick-fil-A -> KFC
Time taken by function: 0 ms
```

- Time Complexity: O(V + E) where V is number of vertexs and E is number of edges
- Graph:
![6_tpsort](https://user-images.githubusercontent.com/63425702/166176018-3b1d5a0e-96ca-4a16-96e3-6a090a9adfef.png)


## Phase 3 report


### Traveling Trojan Problem  
In this problem, given a list of `location_ids`, we need to find the shortest path that go through all locations and get back to the starting point.   
There are 3 methods used:   
1. Brute Force  
2. Early Backtracking
3. 2-opt


**TravellingTrojan_Brute_force**  

In this method, we generate all permutations, and returning the minimum using the given input `location_ids`.  
```c++
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_Brute_force(
      std::vector<std::string> location_ids);
```   
We used the helper function to develop the permutation:   
```c++
std::vector<std::vector<int>> TSP_aux(int start,
                                        std::vector<std::vector<double>> &weights,
                                        int cur_node, double cur_cost,
                                        std::vector<int> &cur_path, 
                                        std::pair<double, std::vector<std::vector<int>>> &records);
```

Data Structure:  
By using the tree data structure, and set each untravelled location as child nodes, 
we can find all permutations. The final path is stored at the leaf node.  
```cpp
std::vector<std::vector<double>> weights // record all location pairs distance
std::pair<double, std::vector<std::vector<std::string>>> records:
// records.first: The shortest distance
// records.second: Each path (permutation) we go through
```

- Time Complexity: O(n!) where n is the number of location_ids.

Run time testing:   
- input: "1855166322","9561828291","7863689394","4399914023","7875114139","122925411","122991336","6814820009",
- output: "1855166322","122925411","122991336","6814820009","7863689394","9561828291","7875114139","4399914023","1855166322",
The distance of the path is: 9.56062 miles         
Time taken by function: 32 ms
- Animation:
![output0](https://user-images.githubusercontent.com/98196892/166068709-70aaf4c7-80e9-4c57-927c-6ff69a28ea27.gif)

**TravellingTrojan_Backtracking**

In this method, given the input `location_ids`, we generate some permutation. This method is improved from the Brute Force that we abandon the path when the current distance is larger than current minimum path. Finally we return the minimum.   
```c++
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_Backtracking(
      std::vector<std::string> location_ids);
```  

We used the helper function to develop the permutation:   
```c++
std::vector<std::vector<int>> TSP_aux_Backtracking(int start,
                                        std::vector<std::vector<double>> &weights,
                                        int cur_node, double cur_cost,
                                        std::vector<int> &cur_path, 
                                        std::pair<double, std::vector<std::vector<int>>> &records);
```

Data Structure:   
By using the tree data structure, and set each untravelled location as child nodes, 
we can find all permutations. The final path is stored at the leaf node.
```cpp
std::vector<std::vector<double>> weights // record all location pairs distance
std::pair<double, std::vector<std::vector<std::string>>> records:
// records.first: The shortest distance
// records.second: Each path (permutation) we go through
```

- Time Complexity: O(n!) where n is the number of location_ids.

Run time testing:  
- input: "1855166322","9561828291","7863689394","4399914023","7875114139","122925411","122991336","6814820009",
- output: "1855166322","122925411","122991336","6814820009","7863689394","9561828291","7875114139","4399914023","1855166322",
The distance of the path is: 9.56062 miles
Time taken by function: 25 ms
- Animation:
![output0_backtracking](https://user-images.githubusercontent.com/98196892/166069076-b5061f29-faf5-4222-8351-4b75f0643e54.gif)


**TravellingTrojan_2opt**

In this method, given the input `location_ids`, we used 2 `for` loops to generate a sub part instead of using permutation. This method is improved that we reversed the sub part and checked if the distance becomes shorter. If yes, we updated the minmum distance we records. Finally we return the minimum.   
```c++
 std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_2opt(
      std::vector<std::string> location_ids);
```  

We used the helper function to realize the 2opt algorithm:   
```c++
std::vector<std::string> TwoOptSwap(std::vector<std::string> &route, int i, int k);
```

Data Structure:   
```cpp
  std::pair<double, std::vector<std::vector<std::string>>> records;
  std::vector<std::vector<std::string>> paths;
```

- Time Complexity: O(n^2) where n is the number of location_ids.

Run time testing:  
- input: "1855166322","9561828291","7863689394","4399914023","7875114139","122925411","122991336","6814820009",
- output: "1855166322","122925411","122991336","6814820009","7863689394","9561828291","7875114139","4399914023","1855166322",
The distance of the path is:9.56062 miles
Time taken by function: 2 ms
- Animation:
![output0_2opt](https://user-images.githubusercontent.com/98196892/166073112-58f09f1f-cef1-4df7-804c-8535004c4e91.gif)


**Summary Table**

Number of Points  | Distance | Brute Force  |  Backtracking | 2-opt 
------------- | ------------- | ------------- | ------------- |------------- 
8  | 7.01163 miles | 52 ms | 50 ms | 2 ms
9 | 8.51845 miles | 583 ms | 598 ms | 2 ms 
10  | 9.26329 miles  | 5197 ms | 5189 ms | 3 ms
11 | 8.27903 miles  | 81597 ms  | 81340 ms | 4 ms

### Find Nearby
In this method, we return the location names that satisfy the input of `attribute`, `location name`, `radius` and `k`.  
```c++
  std::vector<std::string> FindNearby(std::string, std::string, double, int);
```   
We construccted a new class to help:   
```c++
class CurNode{
  public:
  std::string name;
  double distance;
   CurNode(std::string name, double distance){
    this->name = name;
    this->distance = distance;
  }
 // Help to store the name and distance together
```

Data Structure:  
```cpp
  std::vector<std::string> res;
  std::string curid;
  std::string curname;
  std::priority_queue<CurNode, std::vector<CurNode>, keepmax> pq;
```

- Time Complexity: O(n).

Run time testing:   
- input: `supermarket`, `Ralphs`, `10`,`10`
- output: 

  1 Trader Joes
  
  2 Cal Mart Beer & Wine Food Store
  
  3 Food 4 Less
  
  Time taken by function: 75 ms
- Graph:
<img width="708" alt="截屏2022-04-29 下午3 38 38" src="https://user-images.githubusercontent.com/98196892/166077266-9a0ec246-e4e1-41be-b2e5-4d473cc2f72f.png">
