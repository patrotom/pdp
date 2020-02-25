#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <cmath>

using namespace std;

class Solution {
public:
    Solution(double price, vector<pair<int, int>> edges, vector<int> vec) {
        m_price = price;
        m_edges = edges;
        m_vec = vec;
    }

    double getPrice() { return m_price; }

    vector<pair<int, int>> getEdges() { return m_edges; }

    vector<int> getVec() { return m_vec; }
private:
    double m_price;
    vector<pair<int, int>> m_edges;
    vector<int> m_vec;
};

class Graph {
public:
    Graph(int n, int k, int b): m_n(n), m_k(k), m_b(b) {
        for (int i = 0; i < n; i++) {
            vector<pair<int, double>> v;
            m_graph.push_back(v);
        }
        m_bestPrice = n * k / 2;
    }

    int getN() { return m_n; }

    int getK() { return m_k; }

    int getB() { return m_b; }

    void insertNode(int u, int v, double w) {
        m_graph.at(u).push_back(make_pair(v, w));
        m_graph.at(v).push_back(make_pair(u, w));
    }

    void insertExclusionPair(int u, int v) {
        m_exclusionPairs[v] = u;
    }

    Solution solveProblem() {
        vector<int> vec(m_n, -1);
        vec.at(0) = 0;
        bbDFS(0, 0.0, vec, vector<pair<int, int>>());

        return Solution(m_bestPrice, m_bestEdges, m_bestVec);
    }
private:
    int m_n, m_k, m_b;
    double m_bestPrice;
    vector<int> m_bestVec;
    vector<pair<int, int>> m_bestEdges;
    vector<vector<pair<int, double>>> m_graph;
    map<int, int> m_exclusionPairs;

    void bbDFS(int u, double price, vector<int> vec, vector<pair<int, int>> edges) {
        int next = u + 1;

        if (next == m_n) {
            if (price < m_bestPrice) {
                m_bestPrice = price;
                m_bestVec = vec;
                m_bestEdges = edges;
            }
            return;
        }

        if (m_exclusionPairs.find(next) != m_exclusionPairs.end()) {
            vec.at(next) = !vec.at(m_exclusionPairs.at(next));
            double newPrice = recalculatePrice(next, price, vec, edges);
            if (newPrice < m_bestPrice)
                bbDFS(next, newPrice, vec, edges);
        }
        else {
            vec.at(next) = 0;
            vector<pair<int, int>> newEdges = edges;
            double newPrice = recalculatePrice(next, price, vec, newEdges);
            if (newPrice < m_bestPrice)
                bbDFS(next, newPrice, vec, newEdges);
            
            vec.at(next) = 1;
            price = recalculatePrice(next, price, vec, edges);
            if (newPrice < m_bestPrice)
                bbDFS(next, price, vec, edges);
        }
    }

    double recalculatePrice(int u, double price, const vector<int> &vec,
                            vector<pair<int, int>> &edges) {
        for (auto n: m_graph.at(u))
            if (vec.at(n.first) != -1 && vec.at(n.first) != vec.at(u)) {
                price += n.second;
                edges.push_back(make_pair(n.first, u));
            }
        return price;
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
    Graph g(inits.at(0), inits.at(1), inits.at(2));

    int edgesNum = g.getN() * g.getK() / 2;
    for (int i = 0; i < edgesNum; i++) {
        getline(cin, rawInput);
        vector<double> nums = split<double>(rawInput);
        g.insertNode(int(nums.at(0)), int(nums.at(1)), nums.at(2));
    }

    for (int i = 0; i < g.getB(); i++) {
        getline(cin, rawInput);
        vector<int> nums = split<int>(rawInput);
        g.insertExclusionPair(nums.at(0), nums.at(1));
    }

    return g;
}

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
    for (auto it: s.getEdges())
        cout << "(" << it.first << ", " << it.second << ")" << endl;
}

int main(int argc, char **argv) {
    Graph g = constructGraph();
    Solution s = g.solveProblem();
    printSolution(s);
    return 0;
}
