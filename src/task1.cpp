#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <cmath>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace chrono;

/**
 * Servers as a placeholder for a solution variables. Solution consists of best
 * price, number of recursive calls, execution time of B&B DFS algorithm, and
 * the vector with nodes that are separated into two disjoint sets.
 */
class Solution {
public:
    /**
     * Default constructor.
     */
    Solution(double price, vector<int> vec, double duration) {
        m_price = price;
        m_vec = vec;
        m_duration = duration;
    }

    /**
     * Returns best price.
     */
    double getPrice() { return m_price; }

    /**
     * Returns vector with nodes separated into two disjoint sets.
     */
    vector<int> getVec() { return m_vec; }

    /**
     * Returns execution time.
     */
    double getDuration() { return m_duration; }
private:
    double m_price, m_duration;
    vector<int> m_vec;
};

/**
 * Represents the graph and includes methods to construct the graph and to solve
 * the minimum exclusion cut.
 */
class Graph {
public:
    /**
     * Default constructor.
     */
    Graph(int n, int k, int b): m_n(n), m_k(k), m_b(b) {
        for (int i = 0; i < n; i++) {
            vector<pair<int, double>> v;
            m_graph.push_back(v);
        }
        m_bestPrice = n * k / 2;
        m_exclusionPairs = vector<int>(m_n, -1);
    }

    /**
     * Returns number of nodes in the graph.
     */
    int getN() { return m_n; }

    /**
     * Returns avarage degree of a node.
     */
    int getK() { return m_k; }

    /**
     * Returns number of exclusion pairs.
     */
    int getB() { return m_b; }

    /**
     * Returns graph representation. 
     */
    vector<vector<pair<int, double>>> getGraph() { return m_graph; }

    /**
     * Inserts edge between nodes u and v with the weight w to the graph.
     */
    void insertEdge(int u, int v, double w) {
        m_graph.at(u).push_back(make_pair(v, w));
        m_graph.at(v).push_back(make_pair(u, w));
    }

    /**
     * Inserts a new exclusion pair between u and v.
     */
    void insertExclusionPair(int u, int v) {
        m_exclusionPairs[v] = u;
    }

    /**
     * Solves the exclusion graph cut problem and returns solution.
     */
    Solution solveProblem() {
        vector<int> vec(m_n, -1);
        vec.at(0) = 0;

        auto start = high_resolution_clock::now(); 
        bbDFS(0, 0.0, vec);
        auto stop = high_resolution_clock::now();

        double duration =
            (double) duration_cast<milliseconds>(stop - start).count() / 1000;

        return Solution(m_bestPrice, m_bestVec, duration);
    }
private:
    int m_n, m_k, m_b;
    double m_bestPrice;
    vector<int> m_bestVec;
    vector<vector<pair<int, double>>> m_graph;
    vector<int> m_exclusionPairs;

    void bbDFS(int u, double price, vector<int>& vec) {
        int next = u + 1;

        if (next == m_n) {
            if (price < m_bestPrice) {
                m_bestPrice = price;
                m_bestVec = vec;
            }
            return;
        }

        double newPrice;

        if (m_exclusionPairs.at(next) != -1) {
            vec.at(next) = !vec.at(m_exclusionPairs.at(next));
            newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice)
                bbDFS(next, newPrice, vec);
        }
        else {
            vector<int> newVec;

            vec.at(next) = 0;
            newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice) {
                newVec = vec;
                bbDFS(next, newPrice, newVec);
            }
            
            vec.at(next) = 1;
            newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice) {
                newVec = vec;
                bbDFS(next, newPrice, newVec);
            }
        }
    }

    double recalculatePrice(int u, double price, const vector<int>& vec) {
        for (auto n: m_graph.at(u))
            if (vec.at(n.first) != -1 && vec.at(n.first) != vec.at(u))
                price += n.second;
        return price;
    }
};

/**
 * Helper function which splits a string through a whitespace and adds the
 * results to the vector.
 */
template<typename T>
vector<T> split(const string& line) {
    istringstream is(line);
    return vector<T>(istream_iterator<T>(is), istream_iterator<T>());
}

/**
 * Helper function which reads the input from the stdin and constructs a new
 * graph.
 */
Graph constructGraph() {
    string rawInput;
    getline(cin, rawInput);
    vector<int> inits = split<int>(rawInput);
    Graph g(inits.at(0), inits.at(1), inits.at(2));

    int edgesNum = g.getN() * g.getK() / 2;
    for (int i = 0; i < edgesNum; i++) {
        getline(cin, rawInput);
        vector<double> nums = split<double>(rawInput);
        g.insertEdge(int(nums.at(0)), int(nums.at(1)), nums.at(2));
    }

    for (int i = 0; i < g.getB(); i++) {
        getline(cin, rawInput);
        vector<int> nums = split<int>(rawInput);
        g.insertExclusionPair(nums.at(0), nums.at(1));
    }

    return g;
}

/**
 * Helper function which prints a solution in the human readable format.
 */
void printSolution(Solution& s, Graph &g) {
    auto graph = g.getGraph();
    auto vec = s.getVec();

    cout << "Price: " << s.getPrice() << endl;

    cout << "-------------------------" << endl;

    cout << "Execution time: " << s.getDuration() << endl;

    cout << "-------------------------" << endl;

    cout << "Set X: ";
    for (size_t i = 0; i < vec.size(); i++)
        if (vec.at(i) == 0)
            cout << i << " ";
    cout << endl;

    cout << "Set Y: ";
    for (size_t i = 0; i < vec.size(); i++)
        if (vec.at(i) == 1)
            cout << i << " ";
    cout << endl;

    cout << "-------------------------" << endl;

    cout << "Edges included in cut:" << endl;
    
    for (size_t i = 0; i < vec.size(); i++) {
        for (auto it: graph[i]) {
            if (vec[it.first] != vec[i] && (int) i < it.first)
                cout << "(" << i << ", " << it.first << ") ";
        }
    }
    cout << endl;
}

int main(int argc, char **argv) {
    Graph g = constructGraph();
    Solution s = g.solveProblem();
    printSolution(s, g);
    return 0;
}
