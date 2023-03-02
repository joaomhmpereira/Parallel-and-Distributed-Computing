#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <queue>
#include <omp.h>
#include "nqueue/queue.hpp"

using namespace std;

typedef struct city {
    int id;
    double min1;
    double min2;
} * City;

typedef struct node {
    std::vector<int> tour;
    double cost;
    double lower_bound;
    int length;
} * Node;

struct cmp_op { bool operator()(const Node& n1, const Node& n2) {
    if (n1->lower_bound > n2->lower_bound) { 
        return true; 
    }
    else if (n1->lower_bound == n2->lower_bound) {
        return n1->tour.back() > n2->tour.back();
    }
    else {
        return false;
    }
}};

double max_value, best_tour_cost;
std::vector<int> best_tour; 

int n_nodes, n_edges;

double * matrix;
std::vector<City> cities;

double initialLowerBound(){
    double LB = 0;
    for (int i = 0; i < n_nodes; i++) {
        LB += cities[i]->min1 + cities[i]->min2;
    }
    return LB / 2;
}

double newBound(City source, City target, double lower_bound){
    double cost = matrix[source->id * n_nodes + target->id];
    double cf = cost >= source->min2 ? source->min2 : source->min1;
    double ct = cost >= target->min2 ? target->min2 : target->min1;  
    
    return lower_bound + cost - (cf + ct) / 2;
}

bool in_tour(std::vector<int> tour, int v) {
    for (int i = 0; i < tour.size(); i++) {
        if (tour.at(i) == v) { return true; }
    }
    return false;
}

void tsp(){
    Node root = (Node) calloc(1, sizeof(struct node));
    std::vector<Node> nodes_created;
    
    root->tour.push_back(0);
    root->cost = 0;
    root->lower_bound = initialLowerBound();
    root->length = 1;

    //nodes_created.push_back(root);

    // initialize priority queue
    PriorityQueue<Node, cmp_op> queue;
    queue.push(root);

    best_tour_cost = max_value;

    while (!queue.empty()){
        Node node = queue.pop();
        int id = node->tour.back();        
        
        //cout << "Queue size: " << queue.size() << endl;

        // All remaining nodes worse than best
        if (node->lower_bound >= best_tour_cost) {
            //for (int i = 0; i < nodes_created.size(); i++) {
            //    free(nodes_created[i]);
            //}
            free(node);
            return;
        }
        
        // Tour complete, check if it is best
        if (node->length == n_nodes) {
            if (node->cost + matrix[id * n_nodes + 0] < best_tour_cost) { 
                best_tour = node->tour;
                best_tour.push_back(0);
                best_tour_cost = node->cost + matrix[id * n_nodes + 0];
            }
        } else {
            for (int i = 0; i < n_nodes; i++) {
                if (matrix[id * n_nodes + i] != 0 && !in_tour(node->tour, i)) {
                    
                    double new_bound_value = newBound(cities[id], cities[i], node->lower_bound);
                    if(new_bound_value > best_tour_cost) {
                        continue;
                    }
                    Node newNode = (Node) calloc(1, sizeof(struct node));
                    newNode->tour = node->tour;
                    newNode->tour.push_back(i);
                    newNode->cost = node->cost + matrix[id * n_nodes + i];
                    newNode->lower_bound = new_bound_value;
                    newNode->length = node->length + 1;
                    queue.push(newNode);
                    //nodes_created.push_back(newNode);
                }
            }
        }
        free(node);
    }
    //for (int i = 0; i < nodes_created.size(); i++) {
    //    free(nodes_created[i]);
    //}
    
    return;
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

void readInputFile(string inputFile) {
    //cout << "Reading input file" << endl;
    ifstream myFile;
    myFile.open(inputFile, ios::in);
    if(!myFile) {
        cout << "error opening file";
        exit(1);
    }
    int coord1, coord2;
    double dist;

    myFile >> n_nodes >> n_edges;

    for (int i = 0; i < n_nodes; i++) {
        City city = (City) calloc(1, sizeof(struct city));
        city->id = i;
        city->min1 = -1;
        city->min2 = -1;
        cities.push_back(city);
    }
    matrix = (double *) calloc(n_nodes * n_nodes, sizeof(double));
    // initialize matrix and every city with min1 and min2 values
    while (myFile >> coord1 >> coord2 >> dist) {
        matrix[coord1 * n_nodes + coord2] = dist;
        matrix[coord2 * n_nodes + coord1] = dist;

        update_mins(coord1, dist);
        update_mins(coord2, dist);
    }
    
    myFile.close();
}

void print_result() {
    if (best_tour_cost > max_value || best_tour.size() == 0){
        cout << "NO SOLUTION" << endl;
    } else {
        fprintf(stdout, "%.1f\n", best_tour_cost);
        for (int i = 0; i < best_tour.size(); i++) {
            cout << best_tour.at(i) << " ";
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
    readInputFile(argv[1]);
    max_value = std::stod(argv[2]);

    exec_time = -omp_get_wtime();

    tsp();
    
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result();

    free(matrix);
    for (int i = 0; i < cities.size(); i++) {
        free(cities[i]);
    }
    return 0;
}