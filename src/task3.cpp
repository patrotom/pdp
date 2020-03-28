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

    State() {}

    State(int size): price(0), depth(0), vec(vector<int>(size, -1)) {
        vec.at(0) = 0;
    }

    State(const State& prevState, int next, graph& graph): depth(prevState.depth + 1),
                                                           vec(prevState.vec) {
        vec.at(depth) = next;
        price = prevState.price;
        recalculatePrice(graph);
    }

    double price;
    int depth;
    vector<int> vec;
private:
    void recalculatePrice(graph& graph) {
        for (auto n: graph.at(depth))
            if (vec.at(n.first) != -1 &&
                vec.at(n.first) != vec.at(depth))
                price += n.second;
    }
};

/**
 * Serves as a placeholder for a solution variables. Solution consists of best
 * price, number of recursive calls, execution time of B&B DFS algorithm, and
 * the vector with nodes that are separated into two disjoint sets.
 */
class Solution {
public:

    Solution() {}

    Solution(int size, double initPrice): m_state(State(size)) {
        m_state.price = initPrice;
    }

    /**
     * Default constructor.
     */
    Solution(State state): m_state(state) {}

    State m_state;
    double m_duration;
};

class Problem {
public:
    /**
     * Default constructor.
     */
    Problem(int n, int k, int b, int mlpCons): m_n(n), m_k(k), m_b(b),
                                               m_mlpCons(mlpCons),
                                               m_bestPrice(n * k / 2),
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

    Solution solveProblem() {
        auto start = high_resolution_clock::now();
        Solution bestSol;
        double initPrice = m_bestPrice;
        generateStates();
        #pragma omp parallel for
            for (size_t i = 0; i < m_states.size(); i++) {
                Solution sol(m_n, initPrice);
                bbDFS(m_states[i], sol);
                #pragma omp critical
                    if (sol.m_state.price < m_bestPrice && (sol.m_state.depth + 1) == m_n) {
                        bestSol = sol;
                        m_bestPrice = sol.m_state.price;
                    }
            }
        auto stop = high_resolution_clock::now();

        bestSol.m_duration =
            (double) duration_cast<milliseconds>(stop - start).count() / 1000;
        return bestSol;
    }
private:
    int m_n, m_k, m_b, m_mlpCons;
    double m_bestPrice;
    graph m_graph;
    vector<int> m_exclusionPairs;
    vector<State> m_states;

    void generateStates() {
        State s0(m_n);
        queue<State> q;
        q.push(s0);

        while(q.size() <= size_t(m_n * m_mlpCons) && !q.empty()) {
            State s = q.front();
            q.pop();
            if ((s.depth + 1) == m_n)
                break;

            if (m_exclusionPairs.at(s.depth + 1) != -1) {
                int next = !s.vec.at(m_exclusionPairs.at(s.depth + 1));
                State tmpState = State(s, next, m_graph);
                q.push(tmpState);
            }
            else {
                State tmpState = State(s, 1, m_graph);
                q.push(tmpState);
                
                tmpState = State(s, 0, m_graph);
                q.push(tmpState);
            }
        }

        while (!q.empty()) {
            m_states.push_back(q.front());
            q.pop();
        }
    }

    void bbDFS(const State& state, Solution& sol) {
        if ((state.depth + 1) == m_n) {
            if (state.price < sol.m_state.price)
                sol.m_state = state;
            return;
        }

        if (m_exclusionPairs.at(state.depth + 1) != -1) {
            int next = !state.vec.at(m_exclusionPairs.at(state.depth + 1));
            State tmpState = State(state, next, m_graph);
            if (tmpState.price < sol.m_state.price)
                bbDFS(tmpState, sol);
        }
        else {
            State tmpState = State(state, 1, m_graph);
            if (tmpState.price < sol.m_state.price)
                bbDFS(tmpState, sol);
            
            tmpState = State(state, 0, m_graph);
            if (tmpState.price < sol.m_state.price)
                bbDFS(tmpState, sol);
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
    auto vec = s.m_state.vec;

    cout << "Price: " << s.m_state.price << endl;

    cout << "-------------------------" << endl;

    cout << "Execution time: " << s.m_duration << endl;

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
    Solution s = p.solveProblem();
    printSolution(s, p);

    return 0;
}
