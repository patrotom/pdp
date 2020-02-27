#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <cmath>

using namespace std;

/**
 * Servers as a placeholder for a solution variables. Solution consists of best
 * price, edges included in the cut and the vector with nodes that are separated
 * into two disjoint sets.
 */
class Solution {
public:
    /**
     * Default constructor.
     */
    Solution(double price, vector<int> vec, size_t recCnt) {
        m_price = price;
        m_vec = vec;
        m_recCnt = recCnt;
    }

    /**
     * Returns best price.
     */
    double getPrice() { return m_price; }

    /**
     * Returns vector with nodes separated into two disjoint sets.
     */
    vector<int> getVec() { return m_vec; }
private:
    double m_price;
    size_t m_recCnt;
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
        m_recCnt = 0;
        bbDFS(0, 0.0, vec);

        return Solution(m_bestPrice, m_bestVec, m_recCnt);
    }
private:
    int m_n, m_k, m_b;
    size_t m_recCnt;
    double m_bestPrice;
    vector<int> m_bestVec;
    vector<vector<pair<int, double>>> m_graph;
    map<int, int> m_exclusionPairs;

    void bbDFS(int u, double price, vector<int> vec) {
        m_recCnt++;
        int next = u + 1;

        if (next == m_n) {
            if (price < m_bestPrice) {
                m_bestPrice = price;
                m_bestVec = vec;
            }
            return;
        }

        if (m_exclusionPairs.find(next) != m_exclusionPairs.end()) {
            vec.at(next) = !vec.at(m_exclusionPairs.at(next));
            double newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice)
                bbDFS(next, newPrice, vec);
        }
        else {
            vec.at(next) = 0;
            double newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice)
                bbDFS(next, newPrice, vec);
            
            vec.at(next) = 1;
            price = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice)
                bbDFS(next, price, vec);
        }
    }

    double recalculatePrice(int u, double price, const vector<int> &vec) {
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
void printSolution(Solution& s) {
    cout << "Best price: " << s.getPrice() << endl;

    cout << "-------------------------" << endl;

    cout << "Set X: ";
    for (size_t i = 0; i < s.getVec().size(); i++)
        if (s.getVec().at(i) == 0)
            cout << i << " ";
    cout << endl;

    cout << "Set Y: ";
    for (size_t i = 0; i < s.getVec().size(); i++)
        if (s.getVec().at(i) == 1)
            cout << i << " ";
    cout << endl;

    cout << "-------------------------" << endl;

    cout << "Edges included in cut:" << endl;
    // for (auto it: s.getEdges())
    //     cout << "(" << it.first << ", " << it.second << ")" << endl;
}

int main(int argc, char **argv) {
    Graph g = constructGraph();
    Solution s = g.solveProblem();
    printSolution(s);
    return 0;
}
