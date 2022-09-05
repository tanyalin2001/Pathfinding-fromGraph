#include "ArgumentManager.h"
#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <stack>
#include <string>
#include <vector>
#include <algorithm>
#define INF 0x3f3f3f3f
#define NINF 0xFFFFFFFF
using namespace std;

class node {
  int v;
  double weight;

public:
  node(int _v, double _weight) {
    v = _v;
    weight = _weight;
  }
  int getV() { return v; }
  double getW() { return weight; }
};

struct CompareW {
  bool operator()(node &p1, node &p2) { return p1.getW() > p2.getW(); }
};

struct ComparenW {
  bool operator()(node &p1, node &p2) { return p1.getW() < p2.getW(); }
};

class Graph {
private:
  int V;
  list<node> *adj;
  list<node> *adjn;

public:
  Graph(int V);
  void addEdge(int u, int v, double weight);
  void addnEdge(int u, int v, double weight);
  void shortestPath(int s, int des, ofstream &out);
  void longestPath(int s, int des, ofstream &out);
  void longestPathUtil(int v, vector<bool> &visited, stack<int> &Stack);
  void longestPathUtil(int source, vector<bool> &visited, stack<int> &Stack,
                       vector<bool> &onstack, bool &isCyclic);
  int count_paths(int src, int dst);
  void path_counter(int src, int dst, int &path_count, vector<bool> &visited);
  void DFS(stack<int> &st, int source, vector<bool> &visted);
};

Graph::Graph(int V) {
  this->V = V;
  adj = new list<node>[V];
  adjn = new list<node>[V];
}

void Graph::addEdge(int u, int v, double weight) {
  node n(v, weight);
  adj[u].push_back(n);
}

void Graph::addnEdge(int u, int v, double weight) {
  node n(v, -1 * weight);
  adjn[u].push_back(n);
}

void Graph::shortestPath(int source, int des, ofstream &out) {
  priority_queue<node, vector<node>, CompareW> q;
  vector<double> dist(V, INF);
  node n(source, 0);
  q.push(n);
  dist[source] = 0;

  while (!q.empty()) {
    node p = q.top();
    int u = p.getV();
    q.pop();
    list<node>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i) {
      int v = (*i).getV();
      double weight = (*i).getW();
      if (dist[v] > dist[u] + weight) {
        dist[v] = dist[u] + weight;
        node no(v, dist[v]);
        q.push(no);
      }
    }
  }
  if (dist[des] == INF) {
    cout << "Infinite";
    out << "Infinite";
  } else {
    cout << (dist[des]) << endl;
    out << (dist[des]) << endl;
  }
}

void Graph::longestPathUtil(int source, vector<bool> &visited,
                            stack<int> &Stack) {
  visited[source] = true;
  list<node>::iterator i;
  for (i = adj[source].begin(); i != adj[source].end(); i++) {
    int v = (*i).getV();
    if (!visited[v])
      longestPathUtil(v, visited, Stack);
  }
  Stack.push(source);
}

void Graph::longestPath(int s, int des, ofstream &out) {
  vector<bool> visited(V, false);
  vector<bool> onstack(V, false);
  bool isCyclic = false;
  vector<double> dist(V, INF);
  dist[s] = 0;
  stack<int> st;
  for (int i = 0; i < V; i++)
    if (!visited[i])
      longestPathUtil(i, visited, st);
    cout << endl;
  while (!st.empty()) {
    int u = st.top();
    st.pop();
    if (dist[u] != INF) {
      for (node v : adj[u]) {
        if (dist[v.getV()] > dist[u] + v.getW() * -1){
          dist[v.getV()] = dist[u] + v.getW() * -1;
          }
      }
    }
  }
  if (dist[des] == INF) {
    cout << "Infinite" << endl;
    out << "Infinite" << endl;
  } else {
    cout << (dist[des] * -1) << endl;
    out << (dist[des] * -1) << endl;
  }
}

int Graph::count_paths(int src, int dst) {
  int count = 0;
  vector<bool> visited(V, false);
  path_counter(src, dst, count, visited);
  return count;
}

void Graph::path_counter(int src, int dst, int &count, vector<bool> &visited) {
  visited[src] = true;
  if (src == dst) {
    count++;
  } else {
    for (node u : adj[src]) {
      if (!visited[u.getV()])
        path_counter(u.getV(), dst, count, visited);
    }
  }
  visited[src] = false;
}

int main(int argc, char *argv[]) {
  ArgumentManager am(argc, argv);
  string infile = am.get("input");
  string outfile = am.get("output");
  string inpath = am.get("path");
  ifstream input(infile);
  ifstream path(inpath);
  ofstream output(outfile);

  string line;
  vector<double> l;
  vector<int> sta;
  vector<int> end;
  while (getline(input, line)) {
    if (line.length() == 0 || line.length() == 1) {
      continue;
    }
    int start = 0;
    int pos = 0;
    while (line[pos] != ' ') {
      pos++;
    }
    int i = stoi(line.substr(start, pos - start));
    start = pos + 1;
    pos++;
    while (line[pos] != ' ') {
      pos++;
    }
    int j = stoi(line.substr(start, pos - start));
    start = pos + 1;
    pos++;
    double length = stod(line.substr(start, line.size() - start));
    sta.push_back(i);
    end.push_back(j);
    l.push_back(length);
  }
  int s, des, max2;
  path >> s >> des;
  if(s>des)
    max2 = s;
  else
    max2 = des;
  int max = *max_element(sta.begin(), sta.end());
  int max1 = *max_element(end.begin(), end.end());
  
  if(max1 > max)
    max = max1;
  if(max2 > max)
    max = max2;
  Graph g(max+1);
  for(int i=0; i<sta.size(); i++){
    g.addEdge(sta[i], end[i], l[i]);
  }
  
  int count = g.count_paths(s, des);
  if (count == 0)
    output << "Infinite" << endl << "Infinite" << endl << "0";
  else if (s == des)
    output << '0' << endl << '0' << endl << '1' << endl;
  else {
    g.shortestPath(s, des, output);
    g.longestPath(s, des, output);
    output << count;
  }
}
