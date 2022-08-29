#include "trojanmap.h"
#include <algorithm>
//-----------------------------------------------------
// TODO: Student should implement the following:
//-----------------------------------------------------
/**
 * GetLat: Get the latitude of a Node given its id. If id does not exist, return -1.
 * 
 * @param  {std::string} id : location id
 * @return {double}         : latitude
 */
double TrojanMap::GetLat(const std::string& id) {
  auto search = data.find(id);
  if (search != data.end()) {
    return data[id].lat;
  }
  return -1;
}

/**
 * GetLon: Get the longitude of a Node given its id. If id does not exist, return -1.
 * 
 * @param  {std::string} id : location id
 * @return {double}         : longitude
 */
double TrojanMap::GetLon(const std::string& id) { 
  if (data.find(id) != data.end()) {
    return data[id].lon;
  }
  return -1;
}

/**
 * GetName: Get the name of a Node given its id. If id does not exist, return "NULL".
 * 
 * @param  {std::string} id : location id
 * @return {std::string}    : name
 */
std::string TrojanMap::GetName(const std::string& id) { 
  if (data.find(id) != data.end()) {
    return data[id].name;
  }
  return "NULL";
}

/**
 * GetNeighborIDs: Get the neighbor ids of a Node. If id does not exist, return an empty vector.
 * 
 * @param  {std::string} id            : location id
 * @return {std::vector<std::string>}  : neighbor ids
 */
std::vector<std::string> TrojanMap::GetNeighborIDs(const std::string& id) {
  if (data.find(id) != data.end()) {
    return data[id].neighbors;
  }
  return {};
}

/**
 * GetID: Given a location name, return the id. 
 * If the node does not exist, return an empty string. 
 *
 * @param  {std::string} name          : location name
 * @return {int}  : id
 */
std::string TrojanMap::GetID(const std::string& name) {
  std::string res = "";


  for (auto iter : data) {
    // iter.first: id
    // iter.second: Node
    auto search = iter.second.name;
    if (search == name) {
      return iter.first;
    }
  }

  return res;
}

/**
 * GetPosition: Given a location name, return the position. If id does not exist, return (-1, -1).
 *
 * @param  {std::string} name          : location name
 * @return {std::pair<double,double>}  : (lat, lon)
 */
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
  std::pair<double, double> results(-1, -1);

  auto find_id = GetID(name);
  if (find_id != "") {
    auto lat = GetLat(find_id);
    auto lon = GetLon(find_id); 
    std::pair<double, double> res(lat, lon);
    return res;
  }

  return results;
}

/**
 * CalculateEditDistance_naive: Calculate edit distance between two location names
 * Without using DP
 */
// int TrojanMap::CalculateEditDistance_naive(std::string a, std::string b, size_t a_length, size_t b_length){
//   // Time Complexity: O(n*m) 
//   // m and n is the length of two strings

//   if (a_length == 0) return b_length;
//   if (b_length == 0) return a_length;

//   if (a[a_length-1] == b[b_length-1]) {
//     return CalculateEditDistance_naive(a, b, a_length-1, b_length-1);
//   }

//   return 1 + std::min(std::min(CalculateEditDistance_naive(a, b, a_length, b_length-1), 
//                  CalculateEditDistance_naive(a, b, a_length-1, b_length)),
//                  CalculateEditDistance_naive(a, b, a_length-1, b_length-1));
// }


/**
 * CalculateEditDistance: Calculate edit distance between two location names
 * Using DP
 */
int TrojanMap::CalculateEditDistance(std::string a, std::string b){
  // Time Complexity: O(n*m) 
  // m and n is the length of two strings

  int m = a.length();
  int n = b.length();

  // Tabluation
  int dp[m+1][n+1];
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= n; j++) {
      if (i == 0) dp[i][j] = j;
      else if (j == 0) dp[i][j] = i;
      else {
        if (a[i-1] == b[j-1]) dp[i][j] = dp[i-1][j-1];
        else {
          dp[i][j] = std::min({dp[i-1][j], dp[i][j-1], dp[i-1][j-1]}) + 1;
        }
      }
    }
  }

  return dp[m][n];
}

/**
 * FindClosestName: Given a location name, return the name with smallest edit distance.
 *
 * @param  {std::string} name          : location name
 * @return {std::string} tmp           : similar name
 */
std::string TrojanMap::FindClosestName(std::string name) {
  // Iterate through the map and find the name with smallest edit distance 

  std::string tmp = "";
  int min_distance = name.length();

  for (auto iter : data) {
    auto temp_name = iter.second.name;
    int current_distance = CalculateEditDistance(name, temp_name);
    min_distance = std::min(min_distance, current_distance);
    if (min_distance == current_distance) tmp = temp_name;
  }
  return tmp;
}


/**
 * Autocomplete: Given a parital name return all the possible locations with
 * partial name as the prefix. The function should be case-insensitive.
 *
 * @param  {std::string} name          : partial name
 * @return {std::vector<std::string>}  : a vector of full names
 */
std::vector<std::string> TrojanMap::Autocomplete(std::string search_name){
  std::vector<std::string> results;

  
  for (auto iter : data) {
    std::string current_name = iter.second.name;
    if (search_name.size() > current_name.size()){
      continue;
    }

    std::string current_name_copy = iter.second.name;
    // Converting searching name and current name to lowercase
    std::transform(search_name.begin(), search_name.end(), search_name.begin(), ::tolower);
    std::transform(current_name.begin(), current_name.end(), current_name.begin(), ::tolower);
    int len = search_name.length();
    if (search_name == current_name.substr(0, len)) {
      results.push_back(current_name_copy);
    }
  }

  return results;
}

/**
 * CalculateDistance: Get the distance between 2 nodes. 
 * 
 * @param  {std::string} a  : a_id
 * @param  {std::string} b  : b_id
 * @return {double}  : distance in mile
 */
double TrojanMap::CalculateDistance(const std::string &a_id, const std::string &b_id) {
  // Do not change this function
  Node a = data[a_id];
  Node b = data[b_id];
  double dlon = (b.lon - a.lon) * M_PI / 180.0;
  double dlat = (b.lat - a.lat) * M_PI / 180.0;
  double p = pow(sin(dlat / 2),2.0) + cos(a.lat * M_PI / 180.0) * cos(b.lat * M_PI / 180.0) * pow(sin(dlon / 2),2.0);
  double c = 2 * asin(std::min(1.0,sqrt(p)));
  return c * 3961;
}

/**
 * CalculatePathLength: Calculates the total path length for the locations inside the vector.
 * 
 * @param  {std::vector<std::string>} path : path
 * @return {double}                        : path length
 */
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) {
  // Do not change this function
  double sum = 0;
  for (int i = 0;i < int(path.size())-1; i++) {
    sum += CalculateDistance(path[i], path[i+1]);
  }
  return sum;
}

/**
 * CalculateShortestPath_Dijkstra: Given 2 locations, return the shortest path which is a
 * list of id. Hint: Use priority queue.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Dijkstra(
    std::string location1_name, std::string location2_name) {
  std::vector<std::string> path;
  std::string start = GetID(location1_name);
  std::string goal = GetID(location2_name);

  std::priority_queue<std::pair<double, std::string>, std::vector<std::pair<double, std::string>> , std::greater<std::pair<double, std::string>> > pq;
  pq.push(make_pair(0, start));

  std::unordered_map<std::string, double> shortest_map;
  for(auto it = data.begin(); it != data.end(); it++){
    shortest_map[it->second.id] = INT_MAX;
  }
  shortest_map[start] = 0;

  std::unordered_map<std::string, std::string> predecessor_map;
  while(!pq.empty()) {
    double cur_dist = pq.top().first;
    std::string cur_node = pq.top().second;
    pq.pop();

    if(cur_node == goal) {
      break;
    } 

    if(cur_dist > shortest_map[cur_node]) {
      continue;
    } 

    for(auto neighbor : GetNeighborIDs(cur_node)){
      double new_dist = cur_dist + CalculateDistance(cur_node, neighbor);
      if(shortest_map[neighbor] > new_dist) {
        shortest_map[neighbor] = new_dist;
        predecessor_map[neighbor] = cur_node;
        pq.push(make_pair(new_dist, neighbor));
      }
    }
  }

  if(shortest_map[goal] != INT_MAX) {
    std::string cur_node = goal;
    while(cur_node != start) {
      path.push_back(cur_node);
      cur_node = predecessor_map[cur_node];
    }

    path.push_back(start);
    std::reverse(std::begin(path), std::end(path));
  }

  return path;
}

/**
 * CalculateShortestPath_Bellman_Ford: Given 2 locations, return the shortest path which is a
 * list of id. Hint: Do the early termination when there is no change on distance.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Bellman_Ford(
    std::string location1_name, std::string location2_name){
  std::vector<std::string> path;

  std::string start_node = GetID(location1_name);
  std::string end_node = GetID(location2_name);

  std::map<std::string, std::vector<std::string>> neighbor_map;
  neighbor_map[start_node] = GetNeighborIDs(start_node);

  std::map<std::string, double> shortest_map;
  for(auto it = data.begin(); it != data.end(); it++){
    shortest_map[it->second.id] = INT_MAX;
  }
  shortest_map[start_node] = 0;

  std::map<std::string, std::string> predecessor_map;
  bool stop = false;

  while(!stop){
    stop = true;
    for(auto it = neighbor_map.begin(); it != neighbor_map.end();) {
      for(auto neighbor : it->second){
        double new_dist = shortest_map[it->first] + CalculateDistance(it->first, neighbor);
        if(shortest_map[neighbor] > new_dist) {
          shortest_map[neighbor] = new_dist;
          predecessor_map[neighbor] = it->first;
          neighbor_map.insert(make_pair(neighbor, GetNeighborIDs(neighbor)));
          stop = false;
        }
      }
      neighbor_map.erase(it++);
    }
  }

  if(shortest_map[end_node] != INT_MAX) {
    std::string cur_node = end_node;
    while(cur_node != start_node) {
      path.push_back(cur_node);
      cur_node = predecessor_map[cur_node];
    }

    path.push_back(start_node);
    std::reverse(std::begin(path), std::end(path));
  }

  return path;
}


/**
 * Travelling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path
 */

std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_Brute_force(
                                    std::vector<std::string> location_ids) {
  // Brute Force method: find all permutations of the path and select the shortest distance method

  // Map: convert string id to int index
  std::unordered_map<int, std::string> map;
  for (size_t i = 0; i < location_ids.size(); i++) {
    std::pair<int, std::string> temp (i,  location_ids[i]);
    map.insert(temp);
  }
  
  // Initialize inputs
  std::pair<double, std::vector<std::vector<int>>> records;
  std::vector<std::vector<double>> weights = Distance_Map(location_ids);
  int start = 0;
  std::vector<int> cur_path = {start};

  // Generate permutations
  records.first = INT64_MAX;
  std::vector<std::vector<int>> temp_route = TSP_aux(start, weights, start, 0, cur_path, records);
  // print_matrix_int(temp_route);

  // Convert int index route to std::string route
  std::pair<double, std::vector<std::vector<std::string>>> records_final;
  print_matrix(records_final.second);
  records_final.first = records.first;
  for (size_t j = 0; j < temp_route.size(); j++) {
    std::vector<int> index = temp_route[j];
    records_final.second.push_back({});
    for (size_t m = 0; m < index.size(); m++) {
      records_final.second[j].push_back(map[index[m]]);
    }
    records_final.second[j].push_back(map[0]);
  }
  
  return records_final;
}

// Find all permutations of the input locations
std::vector<std::vector<int>> TrojanMap::TSP_aux(int start,
                                                std::vector<std::vector<double>> &weights,
                                                int cur_node, double cur_cost,
                                                std::vector<int> &cur_path, 
                                                std::pair<double, std::vector<std::vector<int>>> &records) {
  // If at the leaf
  if (cur_path.size() == weights.size()) {
    double final_cost = cur_cost + weights[cur_node][start];
    
    if (final_cost < records.first) {
      records.first = final_cost;
      records.second.push_back(cur_path);
    }
  }

  // Evaluvate all children
  for (size_t i = 0; i < weights.size(); i++) {
    if (std::find(cur_path.begin(), cur_path.end(), i) != cur_path.end()) {
      continue;
    }
    cur_path.push_back(i);
    std::vector<std::vector<int>> temp = TSP_aux(start, weights, i, 
                                                 cur_cost + weights[cur_node][i], cur_path, records);
    cur_path.pop_back();
  }

  return records.second;
}

// Distance Map
std::vector<std::vector<double>> TrojanMap::Distance_Map(
  const std::vector<std::string> location_ids) {
  int s = location_ids.size(); 
  std::vector<std::vector<double>> result(s, std::vector<double>(s, 0));
  if (location_ids.size() <= 1) {
    return result;
  }
  

  for (size_t i = 0; i < location_ids.size()-1 ; i++) {
    for (size_t j = i + 1; j < location_ids.size(); j++) {
      double dist = CalculateDistance(location_ids[i], location_ids[j]);
      result[i][j] = dist;
      result[j][i] = dist;
    }
  }
  return result;

}

void TrojanMap::print_matrix(const std::vector<std::vector<std::string>> input) {
  for (size_t i = 0; i < input.size(); i++) {
    for (size_t j = 0; j < input[0].size(); j++) {
      std::cout << input[i][j] << "->";
    }
    std::cout << std::endl;
  }
}

void TrojanMap::print_matrix_int(const std::vector<std::vector<int>> input) {
  for (size_t i = 0; i < input.size(); i++) {
    for (size_t j = 0; j < input[0].size(); j++) {
      std::cout << input[i][j] << "->";
    }
    std::cout << std::endl;
  }
}


std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_Backtracking(
                                    std::vector<std::string> location_ids) {
  // Map: convert string id to int index
  std::unordered_map<int, std::string> map;
  for (size_t i = 0; i < location_ids.size(); i++) {
    std::pair<int, std::string> temp (i,  location_ids[i]);
    map.insert(temp);
  }
  
  // Initialize inputs
  std::pair<double, std::vector<std::vector<int>>> records;
  std::vector<std::vector<double>> weights = Distance_Map(location_ids);
  int start = 0;
  std::vector<int> cur_path = {start};

  // Generate permutations
  records.first = INT64_MAX;
  std::vector<std::vector<int>> temp_route = TSP_aux_Backtracking(start, weights, start, 0, cur_path, records);
  // print_matrix_int(temp_route);

  // Convert int index route to std::string route
  std::pair<double, std::vector<std::vector<std::string>>> records_final;
  print_matrix(records_final.second);
  records_final.first = records.first;
  for (size_t j = 0; j < temp_route.size(); j++) {
    std::vector<int> index = temp_route[j];
    records_final.second.push_back({});
    for (size_t m = 0; m < index.size(); m++) {
      records_final.second[j].push_back(map[index[m]]);
    }
    records_final.second[j].push_back(map[0]);
  }
  
  return records_final;
}

// Find all permutations of the input locations
std::vector<std::vector<int>> TrojanMap::TSP_aux_Backtracking(int start,
                                                std::vector<std::vector<double>> &weights,
                                                int cur_node, double cur_cost,
                                                std::vector<int> &cur_path, 
                                                std::pair<double, std::vector<std::vector<int>>> &records) {
  // If at the leaf
  if (cur_path.size() == weights.size()) {
    double final_cost = cur_cost + weights[cur_node][start];
    
    if (final_cost < records.first) {
      records.first = final_cost;
      records.second.push_back(cur_path);
    }
  }

  // Early Backtracking
  if (cur_cost >= records.first) {return records.second;}

  // Evaluvate all children
  for (size_t i = 0; i < weights.size(); i++) {
    if (std::find(cur_path.begin(), cur_path.end(), i) == cur_path.end()) {
      cur_path.push_back(i);
      std::vector<std::vector<int>> temp = TSP_aux(start, weights, i, 
                                                 cur_cost + weights[cur_node][i], cur_path, records);
      cur_path.pop_back();
    }
  }

  return records.second;
}

// 2-opt
std::vector<std::string> TrojanMap::TwoOptSwap(std::vector<std::string> &route, int i, int k){
  std::vector<std::string> ans = route;
  reverse(ans.begin()+i,ans.begin()+k+1);
  return ans;
}

std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_2opt(
      std::vector<std::string> location_ids){
  std::pair<double, std::vector<std::vector<std::string>>> records;
  std::string start = location_ids[0];
  location_ids.push_back(start);
  std::vector<std::vector<std::string>> paths;
  paths.push_back(location_ids);
  double improvement;
  double best_distance;
  do{
    start_again:
    improvement = 0;
    best_distance = CalculatePathLength(location_ids);
    for(size_t i = 1; i <= location_ids.size()-3; i++){
        for(size_t k = i + 1; k <= location_ids.size()-2; k++){
        std::vector<std::string> newRoute = TwoOptSwap(location_ids, i, k);
        double newDistance = CalculatePathLength(newRoute);
        if(newDistance < best_distance){
          improvement = 1;
          best_distance = newDistance;
          location_ids = newRoute;
          paths.push_back(newRoute);
          goto start_again;
        }
      }
    }
  } while(improvement);
  records.first = best_distance;
  records.second = paths;
  return records;
}

/**
 * Given CSV filename, it read and parse locations data from CSV file,
 * and return locations vector for topological sort problem.
 *
 * @param  {std::string} locations_filename     : locations_filename
 * @return {std::vector<std::string>}           : locations 
 */
std::vector<std::string> TrojanMap::ReadLocationsFromCSVFile(std::string locations_filename){
  std::vector<std::string> location_names_from_csv;
  std::string line;
  std::fstream file(locations_filename, std::ios::in);
  if (file.is_open()) {
    // Skip the first line
    getline(file, line);

    // Read each line 
    while (getline(file, line)) {
      location_names_from_csv.push_back(line);
    }
  }
  return location_names_from_csv;
}

/**
 * Given CSV filenames, it read and parse dependencise data from CSV file,
 * and return dependencies vector for topological sort problem.
 *
 * @param  {std::string} dependencies_filename     : dependencies_filename
 * @return {std::vector<std::vector<std::string>>} : dependencies
 */
std::vector<std::vector<std::string>> TrojanMap::ReadDependenciesFromCSVFile(std::string dependencies_filename){
  std::vector<std::vector<std::string>> dependencies_from_csv;
  std::string line, word;
  std::vector<std::string> row;
  std::fstream file(dependencies_filename, std::ios::in);
  if (file.is_open()) {
    // Skip the first line
    getline(file, line);

    // Read each line 
    while (getline(file, line)) {
      row.clear();

      std::stringstream str(line);
      while (getline(str, word, ',')) {
        row.push_back(word);
      }

      dependencies_from_csv.push_back(row);
    }
  }

  return dependencies_from_csv;
}

/**
 * DeliveringTrojan: Given a vector of location names, it should return a sorting of nodes
 * that satisfies the given dependencies. If there is no way to do it, return a empty vector.
 *
 * @param  {std::vector<std::string>} locations                     : locations
 * @param  {std::vector<std::vector<std::string>>} dependencies     : prerequisites
 * @return {std::vector<std::string>} results                       : results
 */
std::vector<std::string> TrojanMap::DeliveringTrojan(std::vector<std::string> &locations,
                                                     std::vector<std::vector<std::string>> &dependencies){
  std::vector<std::string> result;
  std::map<std::string, int>result_position_map;
  std::map<std::string, std::vector<std::string>> edge_map;

  // Create the Directed Acyclic Graph
  for (auto pairs : dependencies) {
    edge_map[pairs[0]].push_back(pairs[1]);
  }
  
  // DFS
  // initialize the view map by 0 for all locations
  std::map<std::string, int> view;
  for (auto location : locations) {
    view[location] = 0;
  }

  for (auto location : locations) {
    auto location_position = std::find(result.begin(), result.end(), location);
    if (view[location]==0 && location_position == result.end()) {
      DeliveringTrojan_DFS_Helper(location, view, edge_map, result);
    }
  }
  std::reverse(result.begin(), result.end());

  // Detect Cycle
  // If a child node appears before its parent, then there is a cycle
  int index = 0;
  for (auto res:result) {
    result_position_map[res] = index;
    index++;
  }

  for (auto pair : dependencies) {
    if (result_position_map[pair[1]] < result_position_map[pair[0]]) {
      std::cout << "Detect Cycle, Empty result! " << std::endl;
      return {};
    }
  }

  return result;  
  
}

// DFS helper function
void TrojanMap::DeliveringTrojan_DFS_Helper(std::string root, std::map<std::string, int> &view, 
                                            std::map<std::string, std::vector<std::string>> edge_map, 
                                            std::vector<std::string> &result) {                                      
  view[root] = 1;
  for (const auto child : edge_map[root]) {
    if (view[child] != 1) {
      DeliveringTrojan_DFS_Helper(child, view, edge_map, result);
    }
  }
  result.push_back(root);
}


/**
 * inSquare: Give a id return whether it is in square or not.
 *
 * @param  {std::string} id            : location id
 * @param  {std::vector<double>} square: four vertexes of the square area
 * @return {bool}                      : in square or not
 */
bool TrojanMap::inSquare(std::string id, std::vector<double> &square) {
  auto lat = GetLat(id);
  auto lon = GetLon(id);
  if (lon >= square[0] && lon <= square[1] && lat <= square[2] && lat >= square[3]){
    return true;
  }
  return false;
}

/**
 * GetSubgraph: Give four vertexes of the square area, return a list of location ids in the squares
 *
 * @param  {std::vector<double>} square         : four vertexes of the square area
 * @return {std::vector<std::string>} subgraph  : list of location ids in the square
 */
std::vector<std::string> TrojanMap::GetSubgraph(std::vector<double> &square) {
  // include all the nodes in subgraph
  std::vector<std::string> subgraph;
  // example: square = {-118.299, -118.264, 34.032, 34.011}
  for (auto it = data.begin(); it != data.end();it++){
    if (inSquare(it->first, square)){
      subgraph.push_back(it->first);
    }
  }
  return subgraph;
}

/**
 * Cycle Detection: Given four points of the square-shape subgraph, return true if there
 * is a cycle path inside the square, false otherwise.
 * 
 * @param {std::vector<std::string>} subgraph: list of location ids in the square
 * @param {std::vector<double>} square: four vertexes of the square area
 * @return {bool}: whether there is a cycle or not
 */
bool TrojanMap::CycleDetection(std::vector<std::string> &subgraph, std::vector<double> &square) {
  std::map<std::string, bool> visited_map;
  for (auto ids: GetSubgraph(square)){
    visited_map[ids] = false;
  }
  
  for(auto it = visited_map.begin(); it != visited_map.end(); it++) {
    if(it->second == false) {
      auto parentid = it->first;
      if(hasCycle(it->first, visited_map, parentid, square) == true) {
        return true;
      }
    }
  }
  return false;
}
bool TrojanMap::hasCycle(std::string current_id, std::map<std::string, bool> &visited, 
                std::string parent_id, std::vector<double> &square){
  visited[current_id] = true;
  for(auto neighbor: GetNeighborIDs(current_id)) {
    if(inSquare(neighbor, square)) {

      if(visited[neighbor] == false) {
        auto parentid = current_id;
        if(hasCycle(neighbor, visited, parentid, square) == true) {
          return true;
        } 
      } else {
        if(neighbor != parent_id) {
          return true;
          }
      }
    }
  }             
  return false;              
}

/**
 * FindNearby: Given a class name C, a location name L and a number r, 
 * find all locations in class C on the map near L with the range of r and return a vector of string ids
 * 
 * @param {std::string} className: the name of the class
 * @param {std::string} locationName: the name of the location
 * @param {int} r: search radius
 * @param {int} k: search numbers
 * @return {std::vector<std::string>}: location name that meets the requirements
 */
std::vector<std::string> TrojanMap::FindNearby(std::string attributesName, std::string name, double r, int k) {
  std::vector<std::string> res;
  std::string curid;
  std::string curname;
  std::priority_queue<CurNode, std::vector<CurNode>, keepmax> pq;
  for(auto it = data.begin(); it != data.end(); it++){
    curid = it->first;
    curname = GetName(curid);
    // must be given attribute
    if (it->second.attributes.count(attributesName)){
      if (curname == name) continue;
      double dis = CalculateDistance(GetID(name),GetID(curname));
      CurNode *newnode = new CurNode(curid, dis);
      // within the radius
      if (dis <= r){
        pq.push(*newnode);
      }
      // smaller than k
      if (pq.size()>k){
        pq.pop();
      }
    }
  }
  while(!pq.empty()){
    CurNode cur = pq.top();
    pq.pop();
    std::string sub = cur.name;
    res.push_back(sub);
  }
  std::reverse(res.begin(), res.end());
  return res;
}

/**
 * CreateGraphFromCSVFile: Read the map data from the csv file
 * 
 */
void TrojanMap::CreateGraphFromCSVFile() {
  // Do not change this function
  std::fstream fin;
  fin.open("src/lib/data.csv", std::ios::in);
  std::string line, word;

  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);

    Node n;
    int count = 0;
    while (getline(s, word, ',')) {
      word.erase(std::remove(word.begin(), word.end(), '\''), word.end());
      word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '{'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '}'), word.end());
      if (count == 0)
        n.id = word;
      else if (count == 1)
        n.lat = stod(word);
      else if (count == 2)
        n.lon = stod(word);
      else if (count == 3)
        n.name = word;
      else {
        word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
        if (isalpha(word[0]))
          n.attributes.insert(word);
        if (isdigit(word[0]))
          n.neighbors.push_back(word);
      }
      count++;
    }
    data[n.id] = n;
  }
  fin.close();
}
