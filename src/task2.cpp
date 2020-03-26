#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <chrono>
#include <omp.h>

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
    Graph(int n, int k, int b, int limit): m_n(n), m_k(k), m_b(b), m_limit(limit) {
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
        #pragma omp parallel
        {
            #pragma omp single
                bbDFS(0, 0.0, vec);
        }
        auto stop = high_resolution_clock::now();

        double duration =
            (double) duration_cast<milliseconds>(stop - start).count() / 1000;

        return Solution(m_bestPrice, m_bestVec, duration);
    }
private:
    int m_n, m_k, m_b, m_limit;
    double m_bestPrice;
    vector<int> m_bestVec;
    vector<vector<pair<int, double>>> m_graph;
    vector<int> m_exclusionPairs;

    void bbDFS(int u, double price, vector<int>& vec) {
        int next = u + 1;

        if (next == m_n) {
            if (price < m_bestPrice) {
                #pragma omp critical
                {
                    if (price < m_bestPrice) {
                        m_bestPrice = price;
                        m_bestVec = vec;
                    }
                }
            }
            return;
        }

        double newPrice;

        if (m_exclusionPairs.at(next) != -1) {
            vec.at(next) = !vec.at(m_exclusionPairs.at(next));
            newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice) {
                #pragma omp task if (u < m_limit)
                    bbDFS(next, newPrice, vec);
            }
        }
        else {
            vector<int> newVec;

            vec.at(next) = 0;
            newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice) {
                newVec = vec;
                #pragma omp task if (u < m_limit)
                    bbDFS(next, newPrice, newVec);
            }
            
            vec.at(next) = 1;
            newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice) {
                newVec = vec;
                #pragma omp task if (u < m_limit)
                    bbDFS(next, newPrice, newVec);
            }
        }
        #pragma omp taskwait
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
Graph constructGraph(int limit) {
    string rawInput;
    getline(cin, rawInput);
    vector<int> inits = split<int>(rawInput);
    Graph g(inits.at(0), inits.at(1), inits.at(2), limit);

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
void printSolution(Solution& s, Graph& g) {
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

/**
 * Helper function to process input command line arguments
 */
int processArgs(int argc, char **argv) {
    if (argc != 2) {
        cerr << "Invalid number of arguments" << endl;
        return -1;
    }

    int x;
    string arg = argv[1];
    try {
        size_t pos;
        x = stoi(arg, &pos);
        if (pos < arg.size()) {
            cerr << "Trailing characters after number: " << arg << endl;
            return -1;
        }
    } catch (invalid_argument const &ex) {
        cerr << "Invalid number: " << arg << endl;
        return -1;
    } catch (out_of_range const &ex) {
        cerr << "Number out of range: " << arg << endl;
        return -1;
    }
    return x;
}

int main(int argc, char **argv) {
    int l = processArgs(argc, argv);
    if (l == -1)
        return 1;

    Graph g = constructGraph(l);
    Solution s = g.solveProblem();
    printSolution(s, g);

    return 0;
}
