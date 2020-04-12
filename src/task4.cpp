#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <queue>
#include <iterator>
#include <chrono>
#include <climits>
#include <omp.h>
#include <mpi.h>
#include <cstddef>

using namespace std;
using namespace chrono;

typedef vector<vector<pair<int, double>>> graph;

struct State {
    double price;
    int depth;
    vector<double> vec;
};

/**
 * Is responsible for creating graph represenation of a problem and solving this
 * problem.
 */
class MECSolver {
public:
    MECSolver() = default;

    /**
     * Constructor used to initialize all of the essential variables.
     */
    MECSolver(int n, int k, int b, int threadNum, int instNum): m_n(n), m_k(k), m_b(b),
                                                                m_threadNum(threadNum),
                                                                m_instNum(instNum),
                                                                m_graph(graph(m_n, vector<pair<int, double>>())),
                                                                m_exclusionPairs(vector<int>(m_n, -1)) {}

    /**
     * Returns number of nodes in the graph.
     */
    int getN() { return m_n; }

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

    void generateStates(queue<State>& q) {
        State s0;
        s0.price = 0;
        s0.depth = 0;
        s0.vec = vector<double>(m_n, -1);
        q.push(s0);

        while(q.size() <= size_t(m_instNum)) {
            State state = q.front();
            q.pop();
            if ((state.depth + 1) == m_n)
                break;

            if (m_exclusionPairs.at(state.depth + 1) != -1) {
                int next = !state.vec.at(m_exclusionPairs.at(state.depth + 1));
                updateState(state, next);
                q.push(state);
            }
            else {
                int depth = state.depth;
                double price = state.price;

                updateState(state, 0);
                q.push(state);

                state.depth = depth;
                state.price = price;
                updateState(state, 1);
                q.push(state);
            }
        }
    }

    /**
     * Solves the problem using OpenMP and data parallelism technique. Returns
     * the final solution.
     */
    State solve(State& state) {
        m_bestState.price = INT_MAX;

        #pragma omp parallel
        {
            #pragma omp single
                solveProblem(state);
        }

        return m_bestState;
    }
private:
    int m_n, m_k, m_b, m_threadNum, m_instNum;
    State m_bestState;
    graph m_graph;
    vector<int> m_exclusionPairs;

    void solveProblem(State& state) {
        int it = state.depth + 1;

        if (it == m_n) {
            if (state.price < m_bestState.price) {
                #pragma omp critical
                {
                    if (state.price < m_bestState.price)
                        m_bestState = state;
                }
            }
            return;
        }

        if (m_exclusionPairs.at(it) != -1) {
            updateState(state, !state.vec.at(m_exclusionPairs.at(it)));
            if (state.price < m_bestState.price) {
                #pragma omp task if ((it - 1) < m_threadNum)
                    solveProblem(state);
            }
        }
        else {
            state.vec.at(it) = 0;
            double price0 = recalculatePrice(it, state.price, state.vec);
            if (price0 < m_bestState.price) {
                State state0 = state;
                updateState(state0, price0);
                #pragma omp task if ((it - 1) < m_threadNum)
                    solveProblem(state0);
            }
            
            state.vec.at(it) = 1;
            double price1 = recalculatePrice(it, state.price, state.vec);
            if (price1 < m_bestState.price) {
                State state1 = state;
                updateState(state1, price1);
                #pragma omp task if ((it - 1) < m_threadNum)
                    solveProblem(state1);
            }
        }
        #pragma omp taskwait
    }

    void updateState(State& state, int next) {
        state.depth++;
        state.vec.at(state.depth) = next;
        state.price = recalculatePrice(state.depth, state.price, state.vec);
    }

    void updateState(State& state, double newPrice) {
        state.depth++;
        state.price = newPrice;
    }

    double recalculatePrice(int u, double price, const vector<double>& vec) {
        for (auto n: m_graph.at(u))
            if (vec.at(n.first) != -1 &&
                vec.at(n.first) != vec.at(u))
                price += n.second;
        return price;
    }
};

class ProcessHandler {
public:
    ProcessHandler(int argc, char **argv) {
        initMPI(argc, argv);
        initSolver(stoi(argv[1]), stoi(argv[2]), argv[3]);
        initStateType();
    }

    void solveProblem() {
        if (m_procNum == 0) {
            auto start = high_resolution_clock::now();
            doMasterWork();
            auto stop = high_resolution_clock::now();
            printSolution((double) duration_cast<milliseconds>(stop - start).count() / 1000);
        }
        else {
            doSlaveWork();
        }

        MPI_Finalize ();
    }
private:
    int m_provided, m_required, m_numProcs, m_procNum, m_threadNum;
    MECSolver m_solver;
    MPI_Datatype m_stateType;
    State m_bestState;

    static const int tag_work = 0;
    static const int tag_done = 1;
    static const int tag_finished = 2;

    void initMPI(int argc, char **argv) {
        m_provided = MPI_THREAD_FUNNELED;
        m_required = MPI_THREAD_FUNNELED;
        MPI_Init_thread(&argc, &argv, m_required, &m_provided);
        if (m_provided < m_required)
            throw runtime_error("MPI library does not provide required threading support");
        MPI_Comm_size(MPI_COMM_WORLD, &m_numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &m_procNum);
    }

    void initStateType() {
        const MPI_Aint displacements[3] = {offsetof(State, price), offsetof(State, depth), };
        const int lengths[3] = {1, 1, m_solver.getN()};
        MPI_Datatype types[3] = {MPI_DOUBLE, MPI_INT, MPI_INT};
        MPI_Type_create_struct(3, lengths, displacements, types, &m_stateType);
        MPI_Type_commit(&m_stateType);
    }

    void initSolver(const int& threadNum, const int& instNum, const string& fileName) {
        string rawInput;
        ifstream inFile(fileName);

        getline(inFile, rawInput);
        vector<int> inits = split<int>(rawInput);
        m_solver = MECSolver(inits.at(0), inits.at(1), inits.at(2), threadNum, instNum);

        int edgesNum = inits.at(0) * inits.at(1) / 2;
        for (int i = 0; i < edgesNum; i++) {
            getline(inFile, rawInput);
            vector<double> nums = split<double>(rawInput);
            m_solver.insertEdge(int(nums.at(0)), int(nums.at(1)), nums.at(2));
        }

        for (int i = 0; i < inits.at(2); i++) {
            getline(inFile, rawInput);
            vector<int> nums = split<int>(rawInput);
            m_solver.insertExclusionPair(nums.at(0), nums.at(1));
        }
    }

    vector<double> serializeState(State& state) {
        vector<double> v = state.vec;
        v.push_back(state.price);
        v.push_back(state.depth);
        return v;
    }

    State deserializeState(vector<double>& vec) {
        State s;
        s.depth = vec.at(vec.size() - 1);
        s.price = vec.at(vec.size() - 2);
        s.vec = vector<double>(vec.begin(), vec.begin() + vec.size() - 2);
        return s;
    }

    void doMasterWork() {
        MPI_Status status;
        queue<State> q;
        m_bestState.price = INT_MAX;
        m_solver.generateStates(q);
        int serSize = m_solver.getN() + 2;

        for (int dest = 1; dest < m_numProcs; dest++) {
            if (!q.empty()) {
                State state = q.front();
                q.pop();
                vector<double> ser = serializeState(state);
                MPI_Send(ser.data(), serSize, MPI_DOUBLE, dest, tag_work, MPI_COMM_WORLD);
            }
            else {
                break;
            }
        }

        int workingSlaves = m_numProcs - 1;

        while (workingSlaves > 0) {
            State state;
            vector<double> ser(serSize);
            MPI_Recv(ser.data(), serSize, MPI_DOUBLE, MPI_ANY_SOURCE, tag_done, MPI_COMM_WORLD, &status);
            state = deserializeState(ser);

            if (state.price < m_bestState.price)
                m_bestState = state;

            if (!q.empty()) {
                state = q.front();
                q.pop();
                ser = serializeState(state);
                MPI_Send(ser.data(), serSize, MPI_DOUBLE, status.MPI_SOURCE, tag_work, MPI_COMM_WORLD);
            }
            else {
                MPI_Send(ser.data(), serSize, MPI_DOUBLE, status.MPI_SOURCE, tag_finished, MPI_COMM_WORLD);
                workingSlaves--;
            }
        }
    }

    void doSlaveWork() {
        MPI_Status status;
        int serSize = m_solver.getN() + 2;
        while (true) {
            State state;
            vector<double> ser(serSize);
            MPI_Recv(ser.data(), serSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == tag_finished) {
                break;
            }
            else if (status.MPI_TAG == tag_work) {
                state = deserializeState(ser);
                State bestState = m_solver.solve(state);
                ser = serializeState(bestState);
                MPI_Send(ser.data(), serSize, MPI_DOUBLE, 0, tag_done, MPI_COMM_WORLD);
            }
            else {
                throw runtime_error("Unknown tag received");
            }
        }
    }

    void printSolution(double duration) {
        auto graph = m_solver.getGraph();
        auto price = m_bestState.price;
        auto vec = m_bestState.vec;

        cout << "Price: " << price << endl;

        cout << "-------------------------" << endl;

        cout << "Execution time: " << duration << endl;

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

    template<typename T>
    vector<T> split(const string& line) {
        istringstream is(line);
        return vector<T>(istream_iterator<T>(is), istream_iterator<T>());
    }
};

int main(int argc, char **argv) {
    ProcessHandler processHandler(argc, argv);
    processHandler.solveProblem();

    return 0;
}
