#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <queue>
#include <iterator>
#include <chrono>
#include <omp.h>
#include <mpi.h>
#include <cstddef>

using namespace std;
using namespace chrono;

typedef vector<vector<pair<int, double>>> graph;

struct State {
    double price;
    vector<int> vec;
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
    MECSolver(int n, int k, int b, int limit): m_n(n), m_k(k), m_b(b),
                                               m_limit(limit),
                                               m_bestPrice(n * k / 2),
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

    /**
     * Solves the problem using OpenMP and data parallelism technique. Returns
     * the final solution.
     */
    State solve(vector<int>& vec) {
        #pragma omp parallel
        {
            #pragma omp single
                solveProblem(0, 0.0, vec);
        }
        State solution;
        solution.price = m_bestPrice;
        solution.vec = m_bestVec;

        return solution;
    }
private:
    int m_n, m_k, m_b, m_limit;
    double m_bestPrice;
    vector<int> m_bestVec;
    graph m_graph;
    vector<int> m_exclusionPairs;

    void solveProblem(int u, double price, vector<int>& vec) {
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
                    solveProblem(next, newPrice, vec);
            }
        }
        else {
            vector<int> newVec;

            vec.at(next) = 0;
            newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice) {
                newVec = vec;
                #pragma omp task if (u < m_limit)
                    solveProblem(next, newPrice, newVec);
            }
            
            vec.at(next) = 1;
            newPrice = recalculatePrice(next, price, vec);
            if (newPrice < m_bestPrice) {
                newVec = vec;
                #pragma omp task if (u < m_limit)
                    solveProblem(next, newPrice, newVec);
            }
        }
        #pragma omp taskwait
    }

    double recalculatePrice(int u, double price, const vector<int>& vec) {
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
        initSolver(stoi(argv[1]), argv[2]);
        initStateType();
    }

    void solveProblem() {
        MPI_Status status;

        if (m_procNum == 0)
            doMasterWork();
        else
            doSlaveWork();

        MPI_Finalize ();
    }
private:
    int m_provided, m_required, m_numProcs, m_procNum, m_threadNum;
    MECSolver m_solver;
    MPI_Datatype m_solutionType;

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
        const MPI_Aint displacements[2] = {offsetof(State, price), offsetof(State, vec)};
        const int lengths[2] = {1, m_solver.getN()};
        MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
        MPI_Type_create_struct(2, lengths, displacements, types, &m_solutionType);
        MPI_Type_commit(&m_solutionType);
    }

    void initSolver(const int& threadNum, const string& fileName) {
        string rawInput;
        ifstream inFile(fileName);

        getline(inFile, rawInput);
        vector<int> inits = split<int>(rawInput);
        m_solver = MECSolver(inits.at(0), inits.at(1), inits.at(2), threadNum);

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

    void doMasterWork() {
        cout << "I am master!" << endl;
    }

    void doSlaveWork() {
        cout << "I am slave!" << endl;
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

    // int mlpCons = processArgs(argc, argv);
    // if (mlpCons == -1)
    //     return 1;

    // Problem p = generateProblem(mlpCons);
    // Solution s = p.solveProblem();
    // printSolution(s, p);

    return 0;
}
