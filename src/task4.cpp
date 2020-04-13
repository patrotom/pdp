#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <queue>
#include <iterator>
#include <chrono>
#include <climits>
#include <omp.h>
#include <mpi.h>
#include <cstddef>
#include <stdexcept>

#define INST_SIZE 160

using namespace std;
using namespace chrono;

typedef vector<vector<pair<int, double>>> graph;

/**
 * Represents a state of problem we want to solve. State consists of price,
 * depth (current bit in vector) and a static array of bits.
 */
struct State {
    double price;
    int depth;
    array<int, INST_SIZE> vec;
};

/**
 * Is responsible for creating graph represenation of a problem and solving this
 * problem.
 */
class MECSolver {
public:
    /**
     * Default constructor
     */
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

    /**
     * Generates initial states(instances) of the problem using BFS algorithm.
     */
    void generateStates(queue<State>& q) {
        State s0;
        initState(s0);
        q.push(s0);

        while(q.size() <= size_t(m_instNum / m_n) && !q.empty()) {
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
     * Solves the problem using OpenMP and task parallelism technique. Returns
     * the best solution.
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
                #pragma omp task shared (state) if ((it - 1) < m_threadNum)
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

    double recalculatePrice(int u, double price, const array<int, INST_SIZE>& vec) {
        for (auto n: m_graph.at(u))
            if (vec.at(n.first) != -1 &&
                vec.at(n.first) != vec.at(u))
                price += n.second;
        return price;
    }

    void initState(State& state) {
        state.price = 0;
        state.depth = 0;

        for (int i = 0; i < m_n; i++)
            state.vec.at(i) = -1;

        state.vec.at(0) = 0;
    }
};

/**
 * Responsible for handling MPI communication. Divides work between master and
 * slave processes. Makes sure that passed messages are received and processed
 * by a correct process.
 */
class ProcessHandler {
public:
    /**
     * Constructor which initializes MPI, solver class, and custom structured
     * MPI type based on the State structure.
     */
    ProcessHandler(int argc, char **argv) {
        initMPI(argc, argv);
        initSolver(stoi(argv[1]), stoi(argv[2]), argv[3]);
        initStateType();
    }

    /**
     * Divides work between master and slave processes and solves the problem
     * by calling appropriate functions.
     */
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
    int m_provided, m_required, m_numProcs, m_procNum, m_threadNum, m_n;
    MECSolver m_solver;
    State m_bestState;
    MPI_Datatype m_stateType;

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

    void initSolver(const int& threadNum, const int& instNum, const string& fileName) {
        string rawInput;
        ifstream inFile(fileName);

        getline(inFile, rawInput);
        vector<int> inits = split<int>(rawInput);
        m_solver = MECSolver(inits.at(0), inits.at(1), inits.at(2), threadNum, instNum);
        m_n = inits.at(0);

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

    void initStateType() {
        const MPI_Aint displacements[3] = {offsetof(State, price),
                                           offsetof(State, depth),
                                           offsetof(State, vec)};
        const int lengths[3] = {1, 1, INST_SIZE};
        MPI_Datatype types[3] = {MPI_DOUBLE, MPI_INT, MPI_INT};
        MPI_Type_create_struct(3, lengths, displacements, types, &m_stateType);
        MPI_Type_commit(&m_stateType);
    }

    void doMasterWork() {
        MPI_Status status;
        queue<State> q;
        m_bestState.price = INT_MAX;
        m_solver.generateStates(q);

        for (int dest = 1; dest < m_numProcs; dest++) {
            if (!q.empty()) {
                State state = q.front();
                q.pop();
                MPI_Send(&state, 1, m_stateType, dest, tag_work, MPI_COMM_WORLD);
            }
            else {
                break;
            }
        }

        int workingSlaves = m_numProcs - 1;

        while (workingSlaves > 0) {
            State state;
            MPI_Recv(&state, 1, m_stateType, MPI_ANY_SOURCE, tag_done, MPI_COMM_WORLD, &status);

            if (state.price < m_bestState.price)
                m_bestState = state;

            if (!q.empty()) {
                state = q.front();
                q.pop();
                MPI_Send(&state, 1, m_stateType, status.MPI_SOURCE, tag_work, MPI_COMM_WORLD);
            }
            else {
                MPI_Send(&state, 1, m_stateType, status.MPI_SOURCE, tag_finished, MPI_COMM_WORLD);
                workingSlaves--;
            }
        }
    }

    void doSlaveWork() {
        MPI_Status status;
        while (true) {
            State state;
            MPI_Recv(&state, 1, m_stateType, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == tag_finished) {
                break;
            }
            else if (status.MPI_TAG == tag_work) {
                State bestState = m_solver.solve(state);
                MPI_Send(&bestState, 1, m_stateType, 0, tag_done, MPI_COMM_WORLD);
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
        for (size_t i = 0; i < m_n; i++)
            if (vec.at(i) == 0)
                cout << i << " ";
        cout << endl;

        cout << "Set Y: ";
        for (size_t i = 0; i < m_n; i++)
            if (vec.at(i) == 1)
                cout << i << " ";
        cout << endl;

        cout << "-------------------------" << endl;

        cout << "Edges included in cut:" << endl;
        
        for (size_t i = 0; i < m_n; i++) {
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
