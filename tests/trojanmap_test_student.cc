#include "gtest/gtest.h"
#include "src/lib/trojanmap.h"

TEST(TrojanMapStudentTest, Test1) {
  EXPECT_EQ(true, true);
}

// Phase 1

// Test GetLat() function
TEST(TrojanMapTest, GetLat) {
  TrojanMap m;
  // Exist id
  EXPECT_EQ(m.GetLat("653725"), 34.0360852);
  // Non-exist id
  EXPECT_EQ(m.GetLat("100"), -1);
}

// Test GetLon() function
TEST(TrojanMapTest, GetLon) {
  TrojanMap m;
  // Exist id
  EXPECT_EQ(m.GetLon("653725"), -118.3212048);
  // Non-exist id
  EXPECT_EQ(m.GetLon("100"), -1);
}

// Test GetName() function
TEST(TrojanMapTest, GetName) {
  TrojanMap m;
  // Exist id, no name
  EXPECT_EQ(m.GetName("653725"), "");
  // Exist id, exist name
  EXPECT_EQ(m.GetName("368167117"), "Ahmanson Commons");
  // Non-exist id
  EXPECT_EQ(m.GetName("100"), "NULL");
}

// Test GetNeighborIDs() function
TEST(TrojanMapTest, GetNeighborIDs) {
  TrojanMap m;
  // Exist id
  std::vector<std::string> expected_result_1 = {"277327731", "1613324102"};
  EXPECT_EQ(m.GetNeighborIDs("653725"), expected_result_1);
  // Exist id
  std::vector<std::string> expected_result_2 = {"6816288727"};
  EXPECT_EQ(m.GetNeighborIDs("368167117"), expected_result_2);
  // Non-exist id
  std::vector<std::string> expected_result_3 = {};
  EXPECT_EQ(m.GetNeighborIDs("100"), expected_result_3);
}

// Test GetID Function
TEST(TrojanMapTest, GetID) {
  TrojanMap m;
  EXPECT_EQ(m.GetID("Ahmanson Commons"), "368167117");
  EXPECT_EQ(m.GetID("Chipotle"), "732641023");
  EXPECT_EQ(m.GetID("No place"), "");
}

// Test Autocomplete function
TEST(TrojanMapTest, Autocomplete) {
  TrojanMap m;
  // Test the simple case
  auto names = m.Autocomplete("Usc");
  std::unordered_set<std::string> gt = {"USC Village Gym", "USC Parking", "USC Fisher Museum of Art", "USC Roski Eye Institute", "USC Credit Union"}; // groundtruth for "USC"
  int success = 0;
  for (auto& n: names) {
    EXPECT_EQ(gt.count(n) > 0, true) << n + " is not in gt.";
    if (gt.count(n) > 0){
      success++;
    }
  }
  EXPECT_EQ(success, gt.size());
  // Test the lower case
  names = m.Autocomplete("usc");
  success = 0;
  for (auto& n: names) {
    EXPECT_EQ(gt.count(n) > 0, true) << n + " is not in gt.";
    if (gt.count(n) > 0){
      success++;
    }
  }
  EXPECT_EQ(success, gt.size());
  // Test the lower and upper case 
  names = m.Autocomplete("uSc"); 
  success = 0;
  for (auto& n: names) {
    EXPECT_EQ(gt.count(n) > 0, true) << n + " is not in gt.";
    if (gt.count(n) > 0){
      success++;
    }
  }
  EXPECT_EQ(success, gt.size());
  // Test the upper case 
  names = m.Autocomplete("USC"); 
  success = 0;
  for (auto& n: names) {
    EXPECT_EQ(gt.count(n) > 0, true) << n + " is not in gt.";
    if (gt.count(n) > 0){
      success++;
    }
  }
  EXPECT_EQ(success, gt.size());
}

// Test FindPosition function
TEST(TrojanMapTest, FindPosition) {
  TrojanMap m;
  
  // Test USC Village Gym
  auto position = m.GetPosition("USC Village Gym");
  std::pair<double, double> gt1(34.0252392, -118.2858186); // groundtruth for "USC Village Gym"
  EXPECT_EQ(position, gt1);
  // Test USC Parking
  position = m.GetPosition("USC Parking");
  std::pair<double, double> gt2(34.0238824, -118.2801114); // groundtruth for "USC Parking"
  EXPECT_EQ(position, gt2);
  // Test USC Fisher Museum of Art
  position = m.GetPosition("USC Fisher Museum of Art");
  std::pair<double, double> gt3(34.0186092, -118.2873476); // groundtruth for "USC Fisher Museum of Art"
  EXPECT_EQ(position, gt3);
  // Test Unknown
  position = m.GetPosition("XXX");
  std::pair<double, double> gt4(-1, -1);
  EXPECT_EQ(position, gt4);
}

// // Test CalculateEditDistance function
// TEST(TrojanMapTest, CalculateEditDistance_naive) {
//   TrojanMap m;
//   EXPECT_EQ(m.CalculateEditDistance_naive("abcdef", "bcf", 6, 3), 3);
//   EXPECT_EQ(m.CalculateEditDistance_naive("leaveylib", "dohenylib", 9, 9), 5);
// }

// Test CalculateEditDistance function
TEST(TrojanMapTest, CalculateEditDistance) {
  TrojanMap m;
  EXPECT_EQ(m.CalculateEditDistance("abcdef", "bcf"), 3);
  EXPECT_EQ(m.CalculateEditDistance("leaveylib", "dohenylib"), 5);
}

// Test FindClosestName function
TEST(TrojanMapTest, FindClosestName) {
  TrojanMap m;
  EXPECT_EQ(m.FindClosestName("Starbks"), "Starbucks");
  EXPECT_EQ(m.FindClosestName("Dulllce"), "Dulce");
}

// Test ReadLocationsFromCSVFile function
TEST(TrojanMapTest, ReadLocationsFromCSVFile) {
  TrojanMap m;

  // Input path must be changed depend on location on different computer
  std::string input_file_name = "/home/ee538/Documents/final-project-youungjulie/input/topologicalsort_locations.csv";
  std::vector<std::string> expected_output = { "Ralphs", "KFC", "Chick-fil-A", "Target", "CAVA" };
  EXPECT_EQ(m.ReadLocationsFromCSVFile(input_file_name), expected_output);
}

// Test ReadDependenciesFromCSVFile function
TEST(TrojanMapTest, ReadDependenciesFromCSVFile) {
  TrojanMap m;

  // Input path must be changed depend on location on different computer
  std::string input_file_name = "/home/ee538/Documents/final-project-youungjulie/input/topologicalsort_dependencies.csv";
  std::vector<std::vector<std::string>> expected_output = { { "Ralphs", "Chick-fil-A" }, { "Ralphs", "KFC" }, 
                                        { "Chick-fil-A", "KFC" }, { "Target", "CAVA" }, { "CAVA", "KFC" } };
  EXPECT_EQ(m.ReadDependenciesFromCSVFile(input_file_name), expected_output);
}
 
TEST(TrojanMapTest, DeliveringTrojan_1) {
  TrojanMap m;
  std::vector<std::string> locations = { "Ralphs", "KFC", "Chick-fil-A", "Target", "CAVA" };
  std::vector<std::vector<std::string>> dependencies =  {{ "Ralphs", "Chick-fil-A" }, { "Ralphs", "KFC" }, { "Chick-fil-A", "KFC" }, { "Target", "CAVA" }, { "CAVA", "KFC" } };
  std::vector<std::string> expected_result = {"Target", "CAVA", "Ralphs", "Chick-fil-A", "KFC"};
  EXPECT_EQ(m.DeliveringTrojan(locations, dependencies), expected_result);
}

TEST(TrojanMapTest, DeliveringTrojan_cycle) {
  TrojanMap m;
  std::vector<std::string> locations = {"Ralphs", "KFC", "Chick-fil-A"};
  std::vector<std::vector<std::string>> dependencies = {{"Ralphs","KFC"}, {"KFC", "Chick-fil-A"}, {"Chick-fil-A", "Ralphs"}};
  std::vector<std::string> expected_result = {};
  EXPECT_EQ(m.DeliveringTrojan(locations, dependencies), expected_result);
}

// Phase 2
// Test CalculateShortestPath_Dijkstra function
TEST(TrojanMapTest, CalculateShortestPath_Dijkstra) {
  TrojanMap m;
  
  // Test from Ralphs to Target
  auto path = m.CalculateShortestPath_Dijkstra("Ralphs", "Target");
  std::vector<std::string> gt{
      "2578244375","4380040154","4380040158","4380040167","6805802087","8410938469","6813416131",
      "7645318201","6813416130","6813416129","123318563","452688940","6816193777","123408705",
      "6816193774","452688933","452688931","123230412","6816193770","6787470576","4015442011",
      "6816193692","6816193693","6816193694","4015377691","544693739","6816193696","6804883323",
      "6807937309","6807937306","6816193698","4015377690","4015377689","122814447","6813416159",
      "6813405266","4015372488","4015372487","6813405229","122719216","6813405232","4015372486",
      "7071032399","4015372485","6813379479","6813379584","6814769289","5237417650"}; // Expected path
  // Print the path lengths
  std::cout << "My path length: "  << m.CalculatePathLength(path) << "miles" << std::endl;
  std::cout << "GT path length: " << m.CalculatePathLength(gt) << "miles" << std::endl;
  EXPECT_EQ(path, gt);
  
  // Reverse the input from Ralphs to Target
  path = m.CalculateShortestPath_Dijkstra("Target", "Ralphs");
  std::reverse(gt.begin(),gt.end()); // Reverse the path

  // Print the path lengths
  std::cout << "My path length: "  << m.CalculatePathLength(path) << "miles" << std::endl;
  std::cout << "GT path length: " << m.CalculatePathLength(gt) << "miles" << std::endl;
  EXPECT_EQ(path, gt);
}

// Test CalculateShortestPath_Bellman_Ford function
TEST(TrojanMapTest, CalculateShortestPath_Bellman_Ford) {
  TrojanMap m;
  
  // Test from Ralphs to Target
  auto path = m.CalculateShortestPath_Bellman_Ford("Ralphs", "Target");
  std::vector<std::string> gt{
    "2578244375","4380040154","4380040158","4380040167","6805802087","8410938469","6813416131",
    "7645318201","6813416130","6813416129","123318563","452688940","6816193777","123408705",
    "6816193774","452688933","452688931","123230412","6816193770","6787470576","4015442011",
    "6816193692","6816193693","6816193694","4015377691","544693739","6816193696","6804883323",
    "6807937309","6807937306","6816193698","4015377690","4015377689","122814447","6813416159",
    "6813405266","4015372488","4015372487","6813405229","122719216","6813405232","4015372486",
    "7071032399","4015372485","6813379479","6813379584","6814769289","5237417650"}; // Expected path
  // Print the path lengths
  std::cout << "My path length: "  << m.CalculatePathLength(path) << "miles" << std::endl;
  std::cout << "GT path length: " << m.CalculatePathLength(gt) << "miles" << std::endl;
  EXPECT_EQ(path, gt);
  
  // Reverse the input from Ralphs to Target
  path = m.CalculateShortestPath_Bellman_Ford("Target", "Ralphs");
  std::reverse(gt.begin(),gt.end()); // Reverse the path

  // Print the path lengths
  std::cout << "My path length: "  << m.CalculatePathLength(path) << "miles" << std::endl;
  std::cout << "GT path length: " << m.CalculatePathLength(gt) << "miles" << std::endl;
  EXPECT_EQ(path, gt);
}

// Test cycle detection function
TEST(TrojanMapTest, CycleDetection) {
  TrojanMap m;
  
  // Test case 1
  std::vector<double> square1 = {-118.260, -118.254, 34.032, 34.021};
  auto sub1 = m.GetSubgraph(square1);
  bool result1 = m.CycleDetection(sub1, square1);
  EXPECT_EQ(result1, true);

  // Test case 2
  std::vector<double> square2 = {-118.290, -118.289, 34.030, 34.021};
  auto sub2 = m.GetSubgraph(square2);
  bool result2 = m.CycleDetection(sub2, square2);
  EXPECT_EQ(result2, false);
}


// Test topologicalSort function
TEST(TrojanMapTest, TopologicalSort) {
  TrojanMap m;

  std::vector<std::string> location_names = {"Cardinal Gardens", "Coffee Bean1", "Aldewaniah", "CVS"};
  std::vector<std::vector<std::string>> dependencies = {{"Cardinal Gardens","Coffee Bean1"}, {"Cardinal Gardens","CVS"}, {"Coffee Bean1","CVS"}, {"CVS", "Aldewaniah"}};
  auto result = m.DeliveringTrojan(location_names, dependencies);
  std::vector<std::string> gt ={"Cardinal Gardens", "Coffee Bean1","CVS", "Aldewaniah"};
  EXPECT_EQ(result, gt);
}


// Test TravellingTrojan_Brute_force function 
TEST(TrojanMapTest, TravellingTrojan_Brute_force) {
  TrojanMap m;

  std::vector<std::string> location_ids = {"8201681442","6197156485","7786565237","6820972477", 
                                           "6807600525","1832234142","6819144993","1873055949"};
  std::pair<double, std::vector<std::vector<std::string>>> result = m.TravellingTrojan_Brute_force(location_ids);
  std::cout << "Total " << result.second.size() << " permutations" <<  std::endl;
  double actual = round(result.first * 100000.0) / 100000.0;
  double expect = 7.94756;
  
  EXPECT_EQ(actual, expect);
}


// Test TravellingTrojan_Backtracking function 
TEST(TrojanMapTest, TravellingTrojan_Backtracking) {
  TrojanMap m;

  std::vector<std::string> location_ids = {"8201681442","6197156485","7786565237","6820972477", 
                                           "6807600525","1832234142","6819144993","1873055949"};
  std::pair<double, std::vector<std::vector<std::string>>> result = m.TravellingTrojan_Backtracking(location_ids);
  std::cout << "Total " << result.second.size() << " permutations" <<  std::endl;
  double actual = round(result.first * 100000.0) / 100000.0;
  double expect = 7.94756;
  
  EXPECT_EQ(actual, expect);
}