#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <cmath>

using namespace std;

// class Node {
// public:
//     Node(int value, bool set) {
//         m_value = value;
//         m_set = set;
//     }

//     int getValue() { return m_value; }

//     bool getSet() { return m_set; }
// private:
//     int m_value;
//     bool m_set;
// };

class Graph {
public:
    Graph(int n, int k, int b): m_n(n), m_k(k), m_b(b) {
        for (int i = 0; i < n; i++) {
            vector<pair<int, double>> v;
            m_graph.push_back(v);
        }
        m_lowerBound = 0;
        m_bestPrice = n * k / 2;
    }

    int getN() { return m_n; }

    int getK() { return m_k; }

    int getB() { return m_b; }

    void insertNode(int u, int v, double w) {
        m_graph[u].push_back(make_pair(v, w));
    }

    void insertExclusionPair(int u, int v) {
        m_exclusionPairs[u] = v;
    }

    void solveProblem() {
        vector<int> sol(m_n, -1);
        sol[0] = 0;
        bbDFS(0, -1, 0.0, 0, sol);
        

    }
private:
    int m_n, m_k, m_b;
    double m_lowerBound, m_bestPrice;
    vector<int> m_bestSolution;
    vector<vector<pair<int, double>>> m_graph;
    map<int, int> m_exclusionPairs;

    void bbDFS(int u, int p, double price, int cnt, vector<int> &sol) {
        if (cnt == m_n) {
            if (price < m_bestPrice) {
                m_bestSolution = sol;
                m_bestPrice = price;
            }
        }
        else {
            if (m_exclusionPairs.find(u) != m_exclusionPairs.end())
                sol[m_exclusionPairs[u]] = !sol[u];
            for (auto it: m_graph[u]) {
                if (it.first != p) {
                    if (sol[it.first] != -1) {
                        price = sol[it.first] == sol[u] ? price : price += it.second;
                        bbDFS(it.first, u, price, ++cnt, sol);
                    }
                    else {
                        sol[it.first] = 0;
                        price = sol[it.first] == sol[u] ? price : price += it.second;
                        if (price < m_bestPrice)
                            bbDFS(it.first, u, price, ++cnt, sol);

                        sol[it.first] = 1;
                        price = sol[it.first] == sol[u] ? price : price += it.second;
                        if (price < m_bestPrice)
                            bbDFS(it.first, u, price, ++cnt, sol);
                    }
                }
            }
        }
    }
};

template<typename T>
vector<T> split(const string& line) {
    istringstream is(line);
    return vector<T>(istream_iterator<T>(is), istream_iterator<T>());
}

Graph constructGraph() {
    string rawInput;
    getline(cin, rawInput);
    vector<int> inits = split<int>(rawInput);
    Graph g(inits[0], inits[1], inits[2]);

    int edgesNum = g.getN() * g.getK() / 2;
    for (int i = 0; i < edgesNum; i++) {
        getline(cin, rawInput);
        vector<double> nums = split<double>(rawInput);
        g.insertNode(int(nums[0]), int(nums[1]), nums[2]);
    }

    for (int i = 0; i < g.getB(); i++) {
        getline(cin, rawInput);
        vector<int> nums = split<int>(rawInput);
        g.insertExclusionPair(nums[0], nums[1]);
    }

    return g;
}

int main(int argc, char **argv) {
    Graph g = constructGraph();
    g.solveProblem();

    return 0;
}
