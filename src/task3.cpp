#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace chrono;

class Problem {
public:
    /**
     * Default constructor.
     */
    Problem(int n, int k, int b, int limit): m_n(n), m_k(k), m_b(b) {
        m_graph = vector<vector<pair<int, double>>>(
            m_n, vector<pair<int, double>>());
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
        m_exclusionPairs.at(v) = u;
    }
private:
    int m_n, m_k, m_b;
    vector<vector<pair<int, double>>> m_graph;
    vector<int> m_exclusionPairs;
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
 * Serves as a placeholder for a solution variables. Solution consists of best
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
 * Helper function which reads the input from the stdin and constructs a new
 * problem.
 */
Problem generateProblem(int limit) {
    string rawInput;
    getline(cin, rawInput);
    vector<int> inits = split<int>(rawInput);
    Problem p(inits.at(0), inits.at(1), inits.at(2), limit);

    int edgesNum = p.getN() * p.getK() / 2;
    for (int i = 0; i < edgesNum; i++) {
        getline(cin, rawInput);
        vector<double> nums = split<double>(rawInput);
        p.insertEdge(int(nums.at(0)), int(nums.at(1)), nums.at(2));
    }

    for (int i = 0; i < p.getB(); i++) {
        getline(cin, rawInput);
        vector<int> nums = split<int>(rawInput);
        p.insertExclusionPair(nums.at(0), nums.at(1));
    }

    return p;
}

/**
 * Helper function which prints a solution in the human readable format.
 */
void printSolution(Solution& s, Problem& p) {
    auto graph = p.getGraph();
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

int main() {
    return 0;
}
