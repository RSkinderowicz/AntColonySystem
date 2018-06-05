/**
 * This is a single-file C++ 11 implementation of the Ant Colony System
 * algorithm for solving the TSP and ATSP, as described in:
 *
 * Dorigo, Marco, and Luca Maria Gambardella. "Ant colony system: a cooperative
 * learning approach to the traveling salesman problem." IEEE Transactions on
 * evolutionary computation 1.1 (1997): 53-66.
 *
 * It is intended mainly for educational puroposes and may not offer the best
 * possible performance.
 *
 * Licensed under terms of MIT license (see LICENSE)
 *
 * Copyright (c) 2018 Rafał Skinderowicz, rafal.skinderowicz@us.edu.pl
 */
#include <iostream>
#include <vector>
#include <cstdint>
#include <cassert>
#include <chrono>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

using namespace std;


std::default_random_engine & get_rng() {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static default_random_engine instance(seed);

    return instance;
}


uint32_t get_random_uint32(uint32_t min, uint32_t max_inclusive) {
    uniform_int_distribution<uint32_t> distribution(min, max_inclusive);
    return distribution( get_rng() );
}


double get_random_double(double from = 0.0, uint32_t to_exclusive = 1.0) {
    uniform_real_distribution<double> distribution(from, to_exclusive);
    return distribution( get_rng() );
}


struct ProblemInstance {
    uint32_t dimension_;
    bool is_symmetric_ = true;
    vector<double> distance_matrix_;
    vector<vector<uint32_t>> nearest_neighbor_lists_;


    ProblemInstance(uint32_t dimension,
                    const vector<double> &distance_matrix,
                    bool is_symmetric) :
        dimension_(dimension),
        is_symmetric_(is_symmetric),
        distance_matrix_(distance_matrix) {

        assert(dimension >= 2);
    }


    void initialize_nn_lists(uint32_t nn_list_size) {
        assert(dimension_ > 1);

        nn_list_size = min(dimension_ - 1, nn_list_size);

        nearest_neighbor_lists_.resize(dimension_);
        vector<uint32_t> neighbors(dimension_);
        for (uint32_t i = 0; i < dimension_; ++i) {
            neighbors[i] = i;
        }

        for (uint32_t node = 0; node < dimension_; ++node) {
            sort(neighbors.begin(), neighbors.end(),
                 [this, node] (uint32_t a, uint32_t b) {
                    return get_distance(node, a) < get_distance(node, b);
                 });
            assert( get_distance(node, neighbors.at(0))
                    <= get_distance(node, neighbors.at(1)) );

            auto &nn_list = nearest_neighbor_lists_.at(node);
            nn_list.clear();
            nn_list.reserve(nn_list_size);
            uint32_t count = 0;
            for (uint32_t i = 0; count < nn_list_size ; ++i) {
                if (neighbors[i] != node) {  // node is not its own neighbor
                    nn_list.push_back(neighbors[i]);
                    ++count;
                }
            }
        }
    }


    double get_distance(uint32_t from, uint32_t to) const {
        assert( (from < dimension_) && (to < dimension_) );

        return distance_matrix_[from * dimension_ + to];
    }

    const vector<uint32_t>& get_nearest_neighbors(uint32_t node) const {
        assert( node < nearest_neighbor_lists_.size() );

        return nearest_neighbor_lists_[node];
    }


    double calculate_route_length(const vector<uint32_t> &route) const {
        double distance = 0;
        if ( !route.empty() ) {
            auto prev_node = route.back();
            for (auto node : route) {
                distance += get_distance(prev_node, node);
                prev_node = node;
            }
        }
        return distance;
    }
};


/**
 * Tries to load a Traveling Salesman Problem (or ATSP) instance in TSPLIB
 * format from file at 'path'. Only the instances with 'EDGE_WEIGHT_TYPE:
 * EUC_2D' or 'EXPLICIT' are supported.
 *
 * Throws runtime_error if the file is in unsupported format or if an error was
 * encountered.
 *
 * Returns the loaded problem instance.
 */
ProblemInstance load_tsplib_instance(const char *path) {
    enum EdgeWeightType {
        EUC_2D, EXPLICIT
    };

    ifstream in(path);

    if (!in.is_open()) {
        throw runtime_error(string("Cannot open TSP instance file: ") + path);
    }

    string line;

    uint32_t dimension = 0;
    vector<double> distances;
    EdgeWeightType edge_weight_type { EUC_2D };
    bool is_symmetric = true;

    while (getline(in, line)) {
        cout << "Read line: " << line << endl;
        if (line.find("TYPE") == 0) {
            if (line.find(" TSP") != string::npos) {
                is_symmetric = true;
            } else if (line.find(" ATSP") != string::npos) {
                is_symmetric = false;
            } else {
                throw runtime_error("Unknown problem type");
            }
        } else if (line.find("DIMENSION") != string::npos) {
            istringstream line_in(line.substr(line.find(':') + 1));
            if ( !(line_in >> dimension) ) {
                throw runtime_error(string("Cannot read instance dimension"));
            }
        } else if (line.find("EDGE_WEIGHT_TYPE") != string::npos) {
            if (line.find(" EUC_2D") != string::npos) {
                edge_weight_type = EUC_2D;
            } else if (line.find(" EXPLICIT") != string::npos) {
                edge_weight_type = EXPLICIT;
            } else {
                throw runtime_error(string("Unsupported edge weight type"));
            }
        } else if (line.find("NODE_COORD_SECTION") != string::npos) {
            vector<pair<double, double>> coords;

            while (getline(in, line)) {
                if (line.find("EOF") == string::npos) {
                    istringstream line_in(line);
                    uint32_t id;
                    pair<double, double> point;
                    line_in >> id >> point.first >> point.second;
                    if (line_in.bad()) {
                        cerr << "Error while reading coordinates";
                    }
                    coords.push_back(point);
                } else {
                    break ;
                }
            }

            distances.resize(dimension * dimension, 0);

            for (uint32_t i = 0; i < dimension; ++i) {
                auto from = coords.at(i);

                for (uint32_t j = 0; j < dimension; ++j) {
                    if (i != j) {
                        auto to = coords.at(j);
                        auto dx = to.first - from.first;
                        auto dy = to.second - from.second;
                        double distance = int(sqrt(dx*dx + dy*dy) + 0.5);
                        distances.at(i * dimension + j) = distance;
                    }
                }
            }
        } else if (line.find("EDGE_WEIGHT_SECTION") != string::npos) {
            assert(dimension > 0);
            if (edge_weight_type != EXPLICIT) {
                throw runtime_error("Expected EXPLICIT edge weight type");
            }

            distances.reserve(dimension * dimension);
            while (getline(in, line)) {
                if (line.find("EOF") != string::npos) {
                    break ;
                }
                istringstream line_in(line);
                double distance;
                while (line_in >> distance) {
                    distances.push_back(distance);
                }
            }
            assert( distances.size() == dimension * dimension );
        }
    }
    in.close();

    assert(dimension > 2);

    return ProblemInstance(dimension, distances, is_symmetric);
}


struct Ant {
    vector<uint32_t> visited_;  // A list of visited nodes, i.e. a route
    vector<uint8_t> is_visited_;
    double cost_ = std::numeric_limits<double>::max();


    void initialize(uint32_t dimension) {
        visited_.clear();
        visited_.reserve(dimension);
        is_visited_.clear();
        is_visited_.resize(dimension, false);
    }


    void visit(uint32_t node) {
        assert( !is_visited_.at(node) );

        visited_.push_back(node);
        is_visited_.at(node) = true;
    }


    bool is_visited(uint32_t node) const {
        assert(node < is_visited_.size());

        return is_visited_[node];
    }


    bool all_visited() const {
        return find(is_visited_.begin(),
                    is_visited_.end(), false) == is_visited_.end();
    }
};


struct PheromoneMemory {
    uint32_t dimension_;
    vector<double> pheromone_values_;  // For every edge (a,b),
                                       // where 0 <= a, b < dimension_
    double min_pheromone_value_;


    PheromoneMemory(uint32_t dimension, double min_pheromone_value = 0) :
        dimension_(dimension),
        min_pheromone_value_(min_pheromone_value) {
        pheromone_values_.resize(dimension * dimension, min_pheromone_value);
    }


    double get(uint32_t from, uint32_t to) const {
        assert( (from < dimension_) && (to < dimension_) );

        return pheromone_values_[from * dimension_ + to];
    }


    void update(uint32_t from, uint32_t to,
                double evaporation_rate,
                double deposit,
                bool is_symmetric) {

        assert( (from < dimension_) && (to < dimension_) );

        auto &value = pheromone_values_[from * dimension_ + to];
        value = value * (1 - evaporation_rate) + evaporation_rate * deposit;

        if (is_symmetric) {
            pheromone_values_[to * dimension_ + from] = value;
        }
    }
};


/**
 * This are based on the article mentioned.
 */
struct ACSParameters {
    double q0_ = 0.9;
    double global_evaportaion_rate_ = 0.1;
    double local_evaporation_rate_ = 0.1;
    uint32_t ants_count_ = 10;
    double beta_ = 2.0;
    uint32_t cand_list_size_ = 15;
};


const uint32_t MaxCandListSize = 64;


/*
 * Moves 'ant' from its current node to a next one chosen according to the ACS
 * rules. Returns the selected node.
 */
uint32_t move_ant(const ProblemInstance &instance,
                  const PheromoneMemory &pheromone,
                  const vector<double> &heuristic,
                  Ant &ant, double q0) {
    assert( !ant.visited_.empty() );

    const auto dimension = instance.dimension_;
    const auto current_node = ant.visited_.back();
    const uint32_t offset = current_node * dimension;

    // A list of the nearest unvisited neighbors of 'current_node':
    uint32_t cand_list[MaxCandListSize];
    uint32_t cand_list_size = 0;
    double candidates_pheromone[MaxCandListSize];
    for (const auto node : instance.get_nearest_neighbors(current_node)) {
        const uint32_t valid = 1 - ant.is_visited(node);
        cand_list[cand_list_size] = node;
        candidates_pheromone[cand_list_size] = pheromone.get(current_node,
                                                             node);
        cand_list_size += valid;
    }

    uint32_t chosen_node = current_node;

    if (cand_list_size > 0) {  // Select from the closest nodes
        const auto q = get_random_double();
        if (q < q0) {   // Select node with the highest attractiveness
            chosen_node = cand_list[0];
            double max_product = candidates_pheromone[0]
                               * heuristic[offset + chosen_node];
            double max_product_trail = candidates_pheromone[0];

            for (uint32_t i = 1; i < cand_list_size; ++i) {
                // We exploit the fact that the nearest neighbors are sorted
                // according to the distance from 'current_node', thus we
                // read value of 'heuristic' only if necessary
                assert( i == 0
                        || ( heuristic.at(offset + cand_list[i])
                             <= heuristic.at(offset + cand_list[i-1]) ) );

                const auto trail = candidates_pheromone[i];
                if (trail > max_product_trail) {
                    const auto node = cand_list[i];
                    // Now its time to read heuristic for (current_node, node)
                    const auto product = trail * heuristic[offset + node];
                    if (product > max_product) {
                        max_product = product;
                        max_product_trail = trail;
                        chosen_node = node;
                    }
                }
            }
        } else {  // Pseudo-random proportional rule, i.e. a roulette wheel
            double products_prefix_sum[MaxCandListSize] = { 0 };
            double total = 0;
            for (uint32_t i = 0; i < cand_list_size; ++i) {
                const auto node = cand_list[i];
                const auto product = pheromone.get(current_node, node)
                                   * heuristic[offset + node];
                total += product;
                products_prefix_sum[i] = total;
            }

            chosen_node = cand_list[cand_list_size - 1];
            const auto r = get_random_double() * total;
            for (uint32_t i = 0; i < cand_list_size; ++i) {
                if (r < products_prefix_sum[i]) {
                    chosen_node = cand_list[i];
                    break ;
                }
            }
        }
    } else {  // Select from the rest of the unvisited nodes the one with the
              // maximum product of pheromone and heuristic
        double max_product = 0;

        for (uint32_t node = 0u; node < dimension; ++node) {
            if ( !ant.is_visited(node) ) {
                const auto product = pheromone.get(current_node, node)
                                   * heuristic[offset + node];
                if (product > max_product) {
                    max_product = product;
                    chosen_node = node;
                }
            }
        }
    }
    assert( chosen_node != current_node );

    ant.visit(chosen_node);
    return chosen_node;
}


/**
 * This creates a solution using nearest neighbor heuristic that always selects
 * a clostest of the (yet) unvisited nodes (cities).
 */
Ant create_solution_nn(const ProblemInstance &instance, uint32_t start_node = 0) {
    Ant ant;
    ant.initialize(instance.dimension_);

    uint32_t current_node = start_node;
    ant.visit(current_node);

    for (uint32_t i = 1; i < instance.dimension_; ++i) {
        uint32_t next_node = current_node;
        const auto &candidates = instance.get_nearest_neighbors(current_node);

        for (auto node : candidates) {
            if ( !ant.is_visited(node) ) {
                next_node = node;
                break ;
            }
        }

        if (next_node == current_node) {  // All closest nodes were visited,
                                          // we have to check the rest
            double min_distance = numeric_limits<double>::max();

            for (uint32_t node = 0; node < instance.dimension_; ++node) {
                if ( !ant.is_visited(node) ) {
                    auto distance = instance.get_distance(current_node, node);
                    if (distance < min_distance) {
                        min_distance = distance;
                        next_node = node;
                    }
                }
            }
        }

        assert( next_node != current_node );

        ant.visit(next_node);
        current_node = next_node;
    }
    return ant;
}


double calc_initial_pheromone_value(const ProblemInstance &instance) {
    const auto ant = create_solution_nn(instance);
    const auto cost = instance.calculate_route_length(ant.visited_);
    return 1. / (cost * instance.dimension_);
}


/**
 * Runs the ACS for the given number of iterations.
 * Returns the best solution (ant).
 */
Ant run_acs(const ProblemInstance &instance,
            const ACSParameters &params,
            uint32_t iterations) {
    const auto initial_pheromone = calc_initial_pheromone_value(instance);
    PheromoneMemory pheromone(instance.dimension_, initial_pheromone);

    vector<double> heuristic;
    heuristic.reserve(instance.dimension_ * instance.dimension_);

    for (auto distance : instance.distance_matrix_) {
        heuristic.push_back( 1 / pow(distance, params.beta_) );
    }

    vector<Ant> ants(params.ants_count_);
    Ant best_ant;

    for (uint32_t iteration = 0; iteration < iterations; ++iteration) {
        for (auto &ant : ants) {
            ant.initialize(instance.dimension_);
            auto start_node = get_random_uint32(0, instance.dimension_ - 1);
            ant.visit(start_node);
        }

        for (uint32_t j = 1; j < instance.dimension_; ++j) {
            for (auto &ant : ants) {
                const auto from_node = ant.visited_.back();
                const auto to_node = move_ant(instance, pheromone,
                                              heuristic, ant, params.q0_);
                pheromone.update(from_node,
                                 to_node,
                                 params.local_evaporation_rate_,
                                 initial_pheromone,
                                 instance.is_symmetric_);
            }
        }
        // Local evaporation on the route's closing edge
        for (const auto &ant : ants) {
            assert(ant.all_visited());

            uint32_t from_node = ant.visited_.back();
            uint32_t to_node = ant.visited_.front();

            pheromone.update(from_node, to_node,
                             params.local_evaporation_rate_,
                             initial_pheromone,
                             instance.is_symmetric_);
        }
        // Have we found an improved solution?
        for (auto &ant : ants) {
            ant.cost_ = instance.calculate_route_length(ant.visited_);
            if (ant.cost_ < best_ant.cost_) {
                best_ant = ant;

                cout << "New best solution found with the cost: "
                     << best_ant.cost_
                     << " at iteration " << iteration << endl;
            }
        }
        // Deposit pheromone on the edges belonging to the best ant's solution
        auto prev_node = best_ant.visited_.back();
        const double deposit = 1.0 / best_ant.cost_;
        for (auto node : best_ant.visited_) {
            // The global update of the pheromone trails
            pheromone.update(prev_node, node,
                             params.global_evaportaion_rate_,
                             deposit,
                             instance.is_symmetric_);
            prev_node = node;
        }
    }
    return best_ant;
}


int main(int argc, char *argv[]) {
    string path = "kroA100.tsp";
    if (argc >= 2) {
        path = argv[1];
    }
    try {
        ACSParameters params;
        auto instance = load_tsplib_instance(path.c_str());
        instance.initialize_nn_lists(params.cand_list_size_);

        run_acs(instance, params, 100000);
    } catch (runtime_error e) {
        cout << "An error has occurred: " << e.what() << endl;
    }
    return 0;
}
