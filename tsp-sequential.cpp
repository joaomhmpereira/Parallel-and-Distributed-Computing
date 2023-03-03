#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <queue>
#include <omp.h>
#include <iterator>
#include <algorithm>
#include "nqueue/queue.hpp"

using namespace std;

/* DATA STRUCTURES */

typedef struct city {
    int id;
    double min1;
    double min2;
} * City;

typedef struct node {
    struct node * parent_in_tour;
    double cost;
    double lower_bound;
    int    length;
    int    id;
} * Node;

struct cmp_op { bool operator()(const Node& n1, const Node& n2) {
    if (n1->lower_bound > n2->lower_bound) { 
        return true; 
    }
    else if (n1->lower_bound == n2->lower_bound) {
        return n1->id > n2->id;
    }
    else {
        return false;
    }
}};

/* GLOBAL VARIABLES */

double max_value, best_tour_cost, * matrix;
std::vector<Node> nodes_created;
std::vector<City> cities;
int n_cities, n_roads;
Node best_tour_node; 

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* AUXILIARIES */

double initialLowerBound() {
    double LB = 0;
    for (int i = 0; i < n_cities; i++) {
        LB += cities[i]->min1 + cities[i]->min2;
    }
    return LB / 2;
}

double newBound(City source, City target, double lower_bound) {
    double cost = matrix[source->id * n_cities + target->id];
    double cf = cost >= source->min2 ? source->min2 : source->min1;
    double ct = cost >= target->min2 ? target->min2 : target->min1;  
    
    double newLB = lower_bound + cost - (cf + ct) / 2;

    return newLB;
}

bool inTour(Node tour_node, int n) {
    Node parent = tour_node;

    while (parent) {
        if (parent->id == n) { return true; }
        else { parent = parent->parent_in_tour; }
    }
    return false;
}

void update_mins(int coord, int dist) {
    if (cities[coord]->min1 == -1) { cities[coord]->min1 = dist; } // if min1 is not set, set it
    else if (cities[coord]->min2 == -1) {  // if min2 is not set, set it
        if (dist < cities[coord]->min1) { // if dist is smaller than min1, set min2 to min1 and min1 to dist
            cities[coord]->min2 = cities[coord]->min1;
            cities[coord]->min1 = dist;
        }
        else { // if dist is bigger than min1, set min2 to dist
            cities[coord]->min2 = dist;
        }
    }
    else { // if min1 and min2 are set, check if dist is smaller than min1 or min2
        if (dist < cities[coord]->min1) { // if dist is smaller than min1, set min2 to min1 and min1 to dist
            cities[coord]->min2 = cities[coord]->min1;
            cities[coord]->min1 = dist;
        }
        else if (dist < cities[coord]->min2) { // if dist is smaller than min2, set min2 to dist
            cities[coord]->min2 = dist;
        }
    }
}

void free_nodes(std::vector<Node> nodes) {
    for (int i = 0; i < nodes.size(); i++) {
        free(nodes.at(i));
    }
}

/* MAINS */

void tsp(){
    Node root = (Node) calloc(1, sizeof(struct node));
    
    root->parent_in_tour = NULL;
    root->cost           = 0;
    root->lower_bound    = initialLowerBound();
    root->length         = 1;
    root->id             = 0;

    nodes_created.push_back(root);

    // initialize priority queue
    PriorityQueue<Node, cmp_op> queue;
    queue.push(root);

    best_tour_cost = max_value;

    while (!queue.empty()){
        Node node = queue.pop();
        int id = node->id;        

        // all remaining nodes worse than best
        if (node->lower_bound >= best_tour_cost) {
            return;
        }
        
        // tour complete, check if it is best
        if (node->length == n_cities) {
            if (node->cost + matrix[id * n_cities + 0] < best_tour_cost 
            && matrix[id * n_cities + 0] != 0) {
                best_tour_node = node;
                best_tour_cost = node->cost + matrix[id * n_cities + 0];
            }
        } else {
            for (int i = 0; i < n_cities; i++) {
                if (matrix[id * n_cities + i] != 0 && !inTour(node, i)) {
                    double new_bound_value = newBound(cities[id], cities[i], node->lower_bound);
                    
                    // cannot be better than best so far
                    if(new_bound_value > best_tour_cost) {
                        continue;
                    }

                    Node newNode = (Node) calloc(1, sizeof(struct node));
                    newNode->parent_in_tour = node;
                    newNode->cost           = node->cost + matrix[id * n_cities + i];
                    newNode->lower_bound    = new_bound_value;
                    newNode->length         = node->length + 1;
                    newNode->id             = i;

                    queue.push(newNode);
                    nodes_created.push_back(newNode);
                }
            }
        }
    }
    return;
}

void parse_inputs(string inputFile) {
    ifstream myFile;
    myFile.open(inputFile, ios::in);
    if(!myFile) {
        cout << "error opening file";
        exit(1);
    }
    int coord1, coord2;
    double dist;

    myFile >> n_cities >> n_roads;

    // initialize cities
    for (int i = 0; i < n_cities; i++) {
        City city = (City) calloc(1, sizeof(struct city));
        city->id = i;
        city->min1 = -1;
        city->min2 = -1;
        cities.push_back(city);
    }
    matrix = (double *) calloc(n_cities * n_cities, sizeof(double));

    // initialize matrix and every city with min1 and min2 values
    while (myFile >> coord1 >> coord2 >> dist) {
        matrix[coord1 * n_cities + coord2] = dist;
        matrix[coord2 * n_cities + coord1] = dist;

        update_mins(coord1, dist);
        update_mins(coord2, dist);
    }
    
    myFile.close();
}

void print_result() {
    if (best_tour_cost > max_value || (best_tour_cost == max_value && !best_tour_node) || (best_tour_cost == max_value && best_tour_node->length != n_cities)){
        cout << "NO SOLUTION" << endl;
    } else {
        fprintf(stdout, "%.1f\n", best_tour_cost);
        Node parent = best_tour_node;

        std::vector<int> path;
        path.push_back(0);
        while (parent != NULL) {
            path.push_back(parent->id);
            parent = parent->parent_in_tour;
        }
        
        reverse(path.begin(), path.end());
        for (int i = 0; i < path.size(); i++) {
            if (i == 0) { cout << path[i]; }
            else {cout << " " << path[i]; }
        }
        cout << endl;
    }
}

int main(int argc, char *argv[]) {
    double exec_time;

    if (argc < 3) {
        std::cerr << "Usage: tsp <input_file> <max-value>\n";
        return 1;
    }
    parse_inputs(argv[1]);
    max_value = std::stod(argv[2]);

    exec_time = -omp_get_wtime();

    tsp();
    
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result();

    free(matrix);
    for (int i = 0; i < cities.size(); i++) {
        free(cities.at(i));
    }
    free_nodes(nodes_created);
    return 0;
}
