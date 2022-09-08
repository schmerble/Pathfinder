// application.cpp <Starter Code>
// Name:Ryan Ly Joe
// UIN: 676047653
// NETID: RLYJOE2
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//
//
// HOW TO USE CREATIVE COMPONENT
// APP WILL ASK FOR 2 INPUTS, FOR BUILDING 1 AND BUILDING
// APP WILL THEN ASK FOR A INTEGER/DOUBLE INPUT, FOR MILES TO RUN
// INFO IS RELAYED, THEN PROCESS IS REPEATED UNTIL # IS INPUT FOR B1

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <stack>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <limits>

#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h"

using namespace std;
using namespace tinyxml2;

double INF = numeric_limits<double>::max();
typedef map<long long, Coordinates> nodesMap;
typedef vector<BuildingInfo> builVector;
typedef vector<FootwayInfo> footVector;
typedef graph<long long, double> Graph;
typedef long long ll;
typedef pair<ll, Coordinates> nodePair;
typedef pair<ll, double> llPair;

class prioritize {
  public:
  bool operator()(const llPair& p1, const llPair& p2) const {
    return p1.second > p2.second;
  }
};

typedef priority_queue<llPair , vector<llPair>, prioritize> pq;

//
// searchBuilding is a helper function --
// Uses string and vector of struct --
// Returns a building if found with string
// Otherwise return an N/A building
//
BuildingInfo searchBuilding(string query, vector<BuildingInfo>& Buildings) {
  // Checks if string is the abbrev of a building in vector
  for (auto& building : Buildings) {
    if (query == building.Abbrev) {
      return building;
    }
  }
  // Checks if string is atleast part of a building in vector
  for (auto& building : Buildings) {
    if (building.Fullname.find(query) != string::npos) {
      return building;
    }
  }
  // Return a null building
  BuildingInfo pogStop;
  pogStop.Abbrev = "N/A";
  pogStop.Fullname = "N/A";
  return pogStop;
}

//
// grabCenterCoords
// Takes into 2 buildings to calculate the midpoint
// Returns a coordinate struct of midpoint
//
Coordinates grabCenterCoords(BuildingInfo& b1, BuildingInfo& b2) {
  double c1Lat, c1Lon, c2Lat, c2Lon;
  c1Lat = b1.Coords.Lat;
  c1Lon = b1.Coords.Lon;
  c2Lat = b2.Coords.Lat;
  c2Lon = b2.Coords.Lon;
  // Save all the values to organize it better
  return centerBetween2Points(c1Lat, c1Lon, c2Lat, c2Lon);
}

//
// calcbuildDistnace
// Takes in a coordinate and a building
// Returns the distance between both points
//
double calcbuildDistance(Coordinates& midpoint, BuildingInfo& Building) {
  double c1Lat, c1Lon, c2Lat, c2Lon;
  c1Lat = midpoint.Lat;
  c1Lon = midpoint.Lon;
  c2Lat = Building.Coords.Lat;
  c2Lon = Building.Coords.Lon;
  // Save all the values to organize it better
  return distBetween2Points(c1Lat, c1Lon, c2Lat, c2Lon);
}

//
// nearestBuilding
// Takes in a coordinate, buildings, and N/A buildings
// Returns the nearest possible building
BuildingInfo nearestBuilding(Coordinates midpoint,
    vector<BuildingInfo>& Buildings,
    set<ll>& badB) {
  double min = INFINITY;
  BuildingInfo minBuilding;
  for (auto& building : Buildings) {
    // Go through every building and calculate distance
    double dist = calcbuildDistance(midpoint, building);
    ll ID = building.Coords.ID;
    // Check if building's distance is the lowest possible
    // ^ also checks if building is not in N/A set
    if (dist < min && (badB.count(ID) == 0)) {
      min = dist;
      minBuilding = building;
    }
  }
  // Returns nearest building
  return minBuilding;
}

//
// printPoint
// Takes in a building
// Prints out the building info
// Name, Lat, Lon
//
void printPoint(BuildingInfo& building) {
  cout << " " << building.Fullname << endl;
  cout << " (" << building.Coords.Lat << ", "
  << building.Coords.Lon << ")" << endl;
}

//
// printTripoints
// Takes in 3 buildings
// Prints out all building info
// Name, lat lon
//
void printTripoints(BuildingInfo& b1, BuildingInfo& b2,
    BuildingInfo& b3) {
  cout << "Person 1's point:"<< endl;
  printPoint(b1);
  cout << "Person 2's point:"<< endl;
  printPoint(b2);
  cout << "Destination Building:" << endl;
  printPoint(b3);
  cout << endl;
}

//
// nearestNode
// Grabs the nearest possible node from building
// Takes in a building, footways, and all possible nodes
// Returns ID of nearest node
//
ll nearestNode(BuildingInfo building, footVector& Footways, nodesMap allNodes) {
  double min = INFINITY;
  ll temp;
  for (auto& FootWay : Footways) {
    for (auto& FootNode : FootWay.Nodes) {
      // Gets to each possible node in a footway
      // Calculates distance from that node to the building
      double dist = calcbuildDistance(allNodes.at(FootNode), building);
      if (dist <= min) {
        min = dist;
        temp = FootNode;
      }
    }
  }
  // Returns nearest node;
  return temp;
}

//
// printNode
// Takes in an ID, and a node
// Prints out that Node ID's specific details
//
void printNode(ll& ID, nodesMap& Nodes) {
  cout << " " << ID << endl << " (" << Nodes.at(ID).Lat <<
  ", " << Nodes.at(ID).Lon << ")" << endl;
}

//
// printTriNodes
// Takes in 3 IDS, a building, all Nodes, and an integer
// iteNum determines if first time a building was calculated
// ^ if not, then this is not the original building
// Prints Nodes depending on iteNum
//
void printTriNodes(ll& p1, ll& p2, ll& pMid, BuildingInfo building,
    nodesMap& Nodes, int iteNum) {
  if (iteNum == 0) {
    // If first time, print all 3 nodes
    cout << "Nearest P1 node:" << endl;
    printNode(p1, Nodes);
    cout << "Nearest P2 node:" << endl;
    printNode(p2, Nodes);
    cout << "Nearest destination node:" << endl;
    printNode(pMid, Nodes);
  } else {
    // If any other time, print the NEW node
    cout << endl << "New destination building:" << endl;
    printPoint(building);
    cout << "Nearest destination node:" << endl;
    printNode(pMid, Nodes);
  }
}

//
// DijkstraShortestPath
// Takes in 2 IDS, a graph of all nodes, 2 maps, and a vector
// Changes predecessor and distances due to pass by reference
// Returns nothing
//
void DijkstraShortestPath(Graph& g, ll pID, ll destID,
    map<ll, ll>& predecessor, map<ll, double>& distances,
    vector<ll>& path) {
  predecessor.clear();  // <- Sanity check
  distances.clear();  // <- Sanity check
  pq unvisitedQueue;
  set<ll> visitedSet;
  // Starts by initiating every vertex for Dijkstra's
  // Adds them into an unvisitedQueue that is sorted.
  for (auto& vertice : g.getVertices()) {
    distances[vertice] = INFINITY;
    predecessor[vertice] = -1;
    unvisitedQueue.push(make_pair(vertice, INFINITY));
  }
  // Starts with the given start Vertex (pID)
  // ^ with distance = 0;
  // This should be first in the PQ;
  distances[pID] = 0;
  unvisitedQueue.push(make_pair(pID, 0));
  while (!unvisitedQueue.empty()) {
    // Go through each element in unvisited queue
    pair<ll, double> currentV = unvisitedQueue.top();
    unvisitedQueue.pop();
    // Grab top of queue and pop it.
    if (distances[currentV.first] == INFINITY) {
      // Don't run dijkstras if pair from queue is invalid
      break;
    } else if (visitedSet.count(currentV.first) == 1) {
      // If pair from queue is visited, skip past else
      continue;
    } else {
      // Add pair from queue to visited set
      visitedSet.insert(currentV.first);
      // Access neighbors of pair
      for (auto& adj : g.neighbors(currentV.first)) {
        double edgeLen = 0;
        g.getWeight(adj, currentV.first, edgeLen);
        double altPathDist = distances[currentV.first] + edgeLen;
        // Grabs weight, calculate distance from not to previous node
        if (altPathDist < distances[adj]) {
          // If a shorter path is found, change previous nodes.
          distances[adj] = altPathDist;
          predecessor[adj] = currentV.first;
          unvisitedQueue.push(make_pair(adj, altPathDist));
          // Add this pair BACK into the priority queue
        }
      }
    }
  }
}

//
// preparePath
// Takes in 2 maps, a path, and 2 IDS
// returns if destination ID is valid
// Changes the path of vector with pass by reference
//
bool preparePath(map<ll, ll>& predecessor,
    ll pID, ll destID,
    map<ll, double>& distances,
    vector<ll>& path) {
  path.clear();  // <-- Sanity Check
  // Checks if the distance of destID is invalid
  if (distances[destID] >= INFINITY) {
    return false;
  }
  stack<ll> predNodes;
  predNodes.push(destID);
  ll index = destID;
  // Creates a stack of predecessors
  while (predecessor.at(index) != -1) {
    index = predecessor.at(index);
    predNodes.push(index);
  }
  // Fills up stack by index of predecessor
  // modify index to get next predecessor
  while (!predNodes.empty()) {
    path.push_back(predNodes.top());
    predNodes.pop();
  }
  // While stack isn't empty
  // Adds a node to path, pops stack
  return true;
}

//
// printPath
// Takes in a path
// Prints out that path
//
void printPath(vector<ll>& path) {
  cout << "Path: ";
  // Iterates each ID node until second to last element
  for (int i = 0; i < path.size() - 1; i++) {
    cout << path[i] << "->";
  }
  // Print last element without arrow
  cout << path[path.size() - 1] << endl;
}

//
// calcDist
// Takes in a path, all nodes, and a map of distances
// returns the distance from a start node to a dest node
//
double calcDist(vector<ll>& path, nodesMap& allNodes,
    map<ll, double>& distances) {
  // If that path is 1 node or less
  // That means theres no distance, since it's one node
  // return 0;
  if (path.size() < 2) {
    return 0;
  }
  // returns distance from dest node to start
  return distances[path[path.size()-1]];
}

//
// parseInput
// Take in literally everything /s
// 3 IDS, 3 Buildings, Coordinates, All nodes, All Buildings
// cont: All Footways, all badBuildings, and an integer
// Chain of events to find ID of nodes.
void parseInput(ll& p1ID, ll& p2ID, ll& midID,
    BuildingInfo& p1, BuildingInfo& p2,
    BuildingInfo& pMid, Coordinates& midpoint,
    builVector& Buildings, nodesMap& Nodes, footVector& Footways,
    set<ll>& unreachableBuildings, int time) {
  // First, grab center point between two buildings
  midpoint = grabCenterCoords(p1, p2);
  // Next, grab the nearest building using midpoint
  // Account for unreachable buildings
  pMid = nearestBuilding(midpoint, Buildings, unreachableBuildings);
  // time = iteration amount
  // if first time (0) print tripoints
  if (time == 0) {
    printTripoints(p1, p2, pMid);
  }
  // Change the 3 given IDs based on nearestNodes
  // ^ of buildings.
  p1ID = nearestNode(p1, Footways, Nodes);
  p2ID = nearestNode(p2, Footways, Nodes);
  midID = nearestNode(pMid, Footways, Nodes);
}

//
// parseCreativeInput
// Takes in 2 buildings, a node, footways.
// Changes 2 ll ids according to respective building
// Slightly modified parse function from above
//
void parseCreativeInput(ll& p1ID, ll& p2ID,
    BuildingInfo& p1, BuildingInfo& p2,
    nodesMap& Nodes, footVector& Footways) {
    p1ID = nearestNode(p1, Footways, Nodes);
    p2ID = nearestNode(p2, Footways, Nodes);
    }
//
// CheckPath
// Takes in a map of node distances and an ID
// Returns if path exists
//
bool CheckPath(map<ll, double>& distances, ll p2ID) {
  // If distance of dest id is not valid return false
  // else return true!
  return (distances[p2ID] >= INFINITY) ? false : true;
}

//
// calcPaths
// Takes in a graph of all nodes, bad buildings,
// 3 node IDS, all Nodes, and building ID
// Returns if person1 and 2 can't reach mid building
bool calcPaths(Graph& g, ll& p1ID, ll& p2ID, ll& midID, nodesMap& Nodes,
    set<ll>& unreachableBuildings, ll& pMid) {
  map<ll, ll> predecessor; map<ll, double> distances; vector<ll> path;
  map<ll, ll> predecessor2; map<ll, double> distances2; vector<ll> path2;
  // ^--- Long ahh declarations for person 1 and person 2 paths
  DijkstraShortestPath(g, p1ID, midID, predecessor, distances, path);
  DijkstraShortestPath(g, p2ID, midID, predecessor2, distances2, path2);
  // Checks if both paths are valid to destination
  // Print if so
  if (CheckPath(distances, midID) && CheckPath(distances2, midID)) {
    preparePath(predecessor, p1ID, midID, distances, path);
    cout << endl << "Person 1's distance to dest: ";
    cout << calcDist(path, Nodes, distances) << " miles"<<  endl;
    printPath(path);
    preparePath(predecessor2, p2ID, midID, distances2, path2);
    cout << endl << "Person 2's distance to dest: ";
    cout << calcDist(path2, Nodes, distances2) << " miles"<<  endl;
    printPath(path2);
  // Checks if there's a path from person 1 to person 2
  // Tell em it's not reachable if so
  } else if (!CheckPath(distances, p2ID)) {
    cout << endl << "Sorry, destination unreachable. " << endl;
  // Checks if either person can't reach destination
  // return false if so
  } else if (!CheckPath(distances, midID) || !CheckPath(distances2, midID)) {
    cout << endl;
    cout << "At least one person was unable to reach the destination building. ";
    cout << "Finding next closest building..." << endl;
    unreachableBuildings.insert(pMid);
    // Put that building into unreachable status
    return false;
  }
  return true;
}

//
// creativecalcPath
// Takes in a graph, double for miles
// 2 buildingIDS, and a map of nodes;
// Calculates how many laps a user has to run to achieve
// input miles
// returns if a path is found;
//
bool creativecalcPaths(Graph& g, ll& p1ID, ll& p2ID, 
    nodesMap& Nodes, double miles) {
  map<ll, ll> predecessor; map<ll, double> distances; vector<ll> path;
  // ^--- Declarations
  DijkstraShortestPath(g, p1ID, p2ID, predecessor, distances, path);
  // Checks if path is valid
  // Print if so
  if (CheckPath(distances, p2ID)) {
    preparePath(predecessor, p1ID, p2ID, distances, path);
    cout << endl << "Distance to dest: ";
    cout << calcDist(path, Nodes, distances) << " miles"<<  endl;
    cout << "That's about " << ceil(miles/calcDist(path, Nodes, distances));
    cout << " laps! You'll be done in about ~~ ";
    cout << ceil(miles/8) << "hour(s) Go Go Go!" << endl;
    printPath(path);
  // Checks if there's a path from building 1 to 2;
  // Tell em it's not reachable if so
  } else if (!CheckPath(distances, p2ID)) {
    cout << endl << "Sorry, destination unreachable. " << endl;
  // Checks if either person can't reach destination
  // return false if so
  }
  return true;
}

//
// printTellInput
// Takes in an integer
// prints one of two statements
// I got lazy
//
void printTellInput(int i) {
  if (i == 0) {
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  } else {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
  }
}

//
// application
// Takes in all nodes, all footways, all buildings
// and a graph of nodes.
// Processes given input
//
void application(
    nodesMap& Nodes, footVector& Footways,
    builVector& Buildings,
    Graph& g) {
  string person1Building, person2Building; BuildingInfo pMid; Coordinates midpoint;
  ll p1ID, p2ID, midID;
  set<ll> badB;
  // ^-- Mad declarations
  cout << endl;
  printTellInput(0);
  // ^-- Prints to ask for input for p1
  getline(cin, person1Building);
  while (person1Building != "#") {
    bool pathFound = false;
    int time = 0;
    // Time represents amt iterated in while loop
    BuildingInfo p1 = searchBuilding(person1Building, Buildings);
    printTellInput(1);
    getline(cin, person2Building);
    // ^ Prints to ask for inpu for p2
    BuildingInfo p2 = searchBuilding(person2Building, Buildings);
    // ^ tries to get building info of given input
    while (p1.Abbrev == "N/A" || p2.Abbrev == "N/A") {
      if (p1.Abbrev == "N/A") {
        cout << "Person 1's building not found " << endl;
      } else if (p2.Abbrev == "N/A") {
        cout << "Person 2's building not found " << endl;
      }
      // ^-- Checks if buildings are valid
      printTellInput(0);
      getline(cin, person1Building);
      if (person1Building == "#") {
        break;
      }
      printTellInput(1);
      getline(cin, person2Building);
      p1 = searchBuilding(person1Building, Buildings);
      p2 = searchBuilding(person2Building, Buildings);
      // ^-- Asks for buildings against and reparses till valid
    }
    // Checks if string is # because I couldn't figure
    // out an intuitive way figure out exit command
    if (person1Building == "#") {
      continue;
    }
    cout << endl;
    badB.clear();
    // v-- While path isn't found // Parse input
    while (!pathFound) {
    parseInput(p1ID, p2ID, midID, p1, p2, pMid, midpoint, Buildings, Nodes, Footways, badB, time);
    printTriNodes(p1ID, p2ID, midID, pMid, Nodes, time);
    pathFound = calcPaths(g, p1ID, p2ID, midID, Nodes, badB, pMid.Coords.ID);
    if (pathFound) {
      break;
      }
    // ^-- restart if path is found
    time++;
    }
    // ^-- Iterates time if not found
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
    // ^-- Loop back to the start of the while loop
  }
}


//
// CREATIVE COMPONENT
// TITLE: Trackstar
// DESCRIPTION: Say you're practicing for the marathon, but you really
// love your campus, so much in fact, that you wanna practice on campus
// This app will ask you for 2 buildings, then ask for how many miles
// they would like to run that day, it will calculate the most optimal path
// then relay how many laps it would take, and how long if running at 8 miles/h
//
void creative(nodesMap& Nodes, footVector& Footways,
    builVector& Buildings,
    Graph& g) {
  string person1Building, person2Building;
  ll p1ID, p2ID;
  // ^-- Mad declarations
  cout << endl;
  cout << "Enter a building (partial name or abbreviation), or #> ";
  // ^-- Prints to ask for input for p1
  getline(cin, person1Building);
  while (person1Building != "#") {
    bool pathFound = false;
    // Time represents amt iterated in while loop
    BuildingInfo p1 = searchBuilding(person1Building, Buildings);
    cout << "Enter a building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    // ^ Prints to ask for inpu for p2
    BuildingInfo p2 = searchBuilding(person2Building, Buildings);
    // ^ tries to get building info of given input
    while (p1.Abbrev == "N/A" || p2.Abbrev == "N/A") {
      if (p1.Abbrev == "N/A") {
        cout << "Building 1 not found. " << endl;
      } else if (p2.Abbrev == "N/A") {
        cout << "Building 2 not found. " << endl;
      }
      // ^-- Checks if buildings are valid
      cout << "Enter a building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      if (person1Building == "#") {
        break;
      }
      cout << "Enter a building (partial name or abbreviation)> ";
      getline(cin, person2Building);
      p1 = searchBuilding(person1Building, Buildings);
      p2 = searchBuilding(person2Building, Buildings);
      // ^-- Asks for buildings against and reparses till valid
    }
    // Checks if string is # because I couldn't figure
    // out an intuitive way figure out exit command
    if (person1Building == "#") {
      continue;
    }
    cout << endl;
    cout << "How many miles do you want to run today?> " << endl;
    double miles;
    cin >> miles;
    // v-- While path isn't found // Parse input
    while (!pathFound) {
    parseCreativeInput(p1ID, p2ID, p1, p2, Nodes, Footways);
    pathFound = creativecalcPaths(g, p1ID, p2ID, Nodes, miles);
    if (pathFound) {
      break;
      }
    }
    cout << "Enter a building (partial name or abbreviation), or #> ";
    cin.ignore();
    // ^-- for some reason, cin has an extra \n breaking my loop input
    getline(cin, person1Building);
    // ^-- Loop back to the start of the while loop
  }
} 

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  nodesMap Nodes;
  // info about each footway, in no particular order
  footVector Footways;
  // info about each building, in no particular order
  builVector Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;


  //
  // TO DO: build the graph, output stats:
  //
  Graph G;
  for (auto& node : Nodes) {
    G.addVertex(node.first);
  }

  for (auto& index : Footways) {
    vector<long long> fwNodes = index.Nodes;
    for (int i = 0; i < fwNodes.size() - 1; i++) {
      double c1Lat = Nodes[index.Nodes[i]].Lat;
      double c1Lon = Nodes[index.Nodes[i]].Lon;
      double c2Lat = Nodes[index.Nodes[i+1]].Lat;
      double c2Lon = Nodes[index.Nodes[i+1]].Lon;
      double distance1 = distBetween2Points(c1Lat, c1Lon, c2Lat, c2Lon);
      G.addEdge(fwNodes[i], fwNodes[i+1], distance1);
      G.addEdge(fwNodes[i+1], fwNodes[i], distance1);
    }
  }

  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
        << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    // TO DO: add argument for the graph you make.
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    // TO DO: add arguments
    creative(Nodes, Footways, Buildings, G);
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
