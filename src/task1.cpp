#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>

using namespace std;

class Graph {
public:
    Graph(size_t n, size_t k, size_t b): m_n(n), m_k(k), m_b(b) {
        for (size_t i = 0; i < n; i++) {
            vector<pair<size_t, double>> v;
            m_graph.push_back(v);
        }
    }

    size_t getN() { return m_n; }

    size_t getK() { return m_k; }

    size_t getB() { return m_b; }

    void insertNode(size_t u, size_t v, double w) {
        m_graph[u].push_back(make_pair(v, w));
        m_graph[v].push_back(make_pair(u, w));
    }

    void insertExclusionPair(size_t u, size_t v) {
        m_exclusionPairs.emplace_back(make_pair(u, v));
    }
private:
    size_t m_n, m_k, m_b;
    vector<vector<pair<size_t, double>>> m_graph;
    vector<pair<size_t, size_t>> m_exclusionPairs;
};

template<typename T>
vector<T> split(const string& line) {
    istringstream is(line);
    return vector<T>(istream_iterator<T>(is), istream_iterator<T>());
}

Graph constructGraph() {
    string rawInput;
    getline(cin, rawInput);
    vector<size_t> inits = split<size_t>(rawInput);
    Graph g(inits[0], inits[1], inits[2]);

    size_t edgesNum = g.getN() * g.getK() / 2;
    for (size_t i = 0; i < edgesNum; i++) {
        getline(cin, rawInput);
        vector<double> nums = split<double>(rawInput);
        g.insertNode(size_t(nums[0]), size_t(nums[1]), nums[2]);
    }

    for (size_t i = 0; i < g.getB(); i++) {
        getline(cin, rawInput);
        vector<size_t> nums = split<size_t>(rawInput);
        g.insertExclusionPair(nums[0], nums[1]);
    }

    return g;
}

int main(int argc, char **argv) {
    Graph g = constructGraph();

    return 0;
}
