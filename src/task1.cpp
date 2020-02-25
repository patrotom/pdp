#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <cmath>

using namespace std;

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

    void solveProblem() {
        vector<int> sol(m_n, -1);
        sol.at(0) = 0;
        bbDFS(0, 0.0, sol);

        cout << m_bestPrice << endl;
    }
private:
    int m_n, m_k, m_b;
    double m_bestPrice;
    vector<int> m_bestSolution;
    vector<vector<pair<int, double>>> m_graph;
    map<int, int> m_exclusionPairs;

    void bbDFS(int u, double price, vector<int> sol) {
        int next = u + 1;

        if (next == m_n) {
            if (price < m_bestPrice) {
                m_bestPrice = price;
                m_bestSolution = sol;
            }
            return;
        }

        if (m_exclusionPairs.find(next) != m_exclusionPairs.end()) {
            sol.at(next) = !sol.at(m_exclusionPairs.at(next));
            double newPrice = recalculatePrice(next, price, sol);
            if (newPrice < m_bestPrice)
                bbDFS(next, newPrice, sol);
        }
        else {
            sol.at(next) = 0;
            double newPrice = recalculatePrice(next, price, sol);
            if (newPrice < m_bestPrice)
                bbDFS(next, newPrice, sol);
            
            sol.at(next) = 1;
            newPrice = recalculatePrice(next, price, sol);
            if (newPrice < m_bestPrice)
                bbDFS(next, newPrice, sol);
        }
    }

    double recalculatePrice(int u, double price, const vector<int> &sol) {
        for (auto n: m_graph.at(u))
            if (sol.at(n.first) != -1 && sol.at(n.first) != sol.at(u))
                price += n.second;
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

int main(int argc, char **argv) {
    Graph g = constructGraph();
    g.solveProblem();

    return 0;
}
