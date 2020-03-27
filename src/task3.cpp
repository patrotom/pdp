#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <queue>
#include <iterator>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace chrono;

typedef vector<vector<pair<int, double>>> graph;

class State {
public:
    State(int size): price(0), depth(0) {
        vec = vector<int>(size, -1);
        vec.at(0) = 0;
    }

    State(State& prevState, int next) {
        vec = prevState.vec;
        depth = prevState.depth + 1;
        vec.at(depth) = next;
    }

    vector<int> vec;
    double price;
    int depth;
};

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

class Problem {
public:
    /**
     * Default constructor.
     */
    Problem(int n, int k, int b, int mlpCons):
        m_n(n), m_k(k), m_b(b), m_mlpCons(mlpCons),
        m_graph(graph(m_n, vector<pair<int, double>>())),
        m_exclusionPairs(vector<int>(m_n, -1)) {}

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
    graph getGraph() { return m_graph; }

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

    void solveProblem() {
        generateStates();

    }
private:
    int m_n, m_k, m_b, m_mlpCons;
    graph m_graph;
    vector<int> m_exclusionPairs;
    vector<State> m_states;

    void recalculatePrice(State& state) {
        for (auto n: m_graph.at(state.depth))
            if (state.vec.at(n.first) != -1 &&
                state.vec.at(n.first) != state.vec.at(state.depth))
                state.price += n.second;
    }

    void generateStates() {
        State s0(m_n);
        queue<State> q;
        q.push(s0);
        double bestPrice = m_n * m_k / 2;

        while(q.size() <= size_t(m_n * m_mlpCons) && !q.empty()) {
            State s = q.front();
            q.pop();
            if ((s.depth + 1) == m_n)
                break;

            State tmpState = State(s, 1);
            recalculatePrice(tmpState);
            if (tmpState.price < bestPrice)
                q.push(tmpState);
            
            tmpState = State(s, 0);
            recalculatePrice(tmpState);
            if (tmpState.price < bestPrice)
                q.push(tmpState);
        }

        while (!q.empty()) {
            m_states.push_back(q.front());
            q.pop();
        }


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
 * problem.
 */
Problem generateProblem(int mlpCons) {
    string rawInput;
    getline(cin, rawInput);
    vector<int> inits = split<int>(rawInput);
    Problem p(inits.at(0), inits.at(1), inits.at(2), mlpCons);

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
    int mlpCons = processArgs(argc, argv);
    if (mlpCons == -1)
        return 1;

    Problem p = generateProblem(mlpCons);
    p.solveProblem();

    return 0;
}
