#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
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
    int current_city;
    friend bool operator>(struct node left, struct node right) { return left.lower_bound > right.lower_bound 
                        || (left.lower_bound == right.lower_bound && left.current_city > right.current_city); }
} * Node;

double max_value, best_tour_cost;
std::vector<int> best_tour; 

int n_nodes, n_edges;

double * matrix;
std::vector<City> cities;

double initialLowerBound(){
    double LB;
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
    Node root = (Node) malloc(sizeof(struct node));
    std::vector<Node> nodes_created;
    root->tour.push_back(0);
    root->cost = 0;
    root->lower_bound = initialLowerBound();
    root->length = 1;
    root->current_city = 0;
    nodes_created.push_back(root);

    cout << "Initial lower bound: " << root->lower_bound << endl;

    // initialize priority queue
    PriorityQueue<Node> queue;
    queue.push(root);

    best_tour_cost = max_value;

    cout << "CHECK 1" << endl;

    while (!queue.empty()){
        Node node = queue.pop();
        int id = node->current_city;
        
        cout << "CHECK 2" << endl;

        // All remaining nodes worse than best
        if (node->lower_bound >= best_tour_cost) {
            return;
        }
        
        // Tour complete, check if it is best
        if (node->length == n_nodes) {
            cout << "Tour complete" << endl;
            if (node->cost + matrix[id * n_nodes + 0] < best_tour_cost) { 
                cout << "New best tour found" << endl; 
                best_tour = node->tour;
                best_tour.push_back(0);
                best_tour_cost = node->cost + matrix[id * n_nodes + 0];
            }
        } else {
            cout << "CHECK 3" << endl;
            for (int i = 0; i < n_nodes; i++) {
                if (matrix[id * n_nodes + i] != 0 && !in_tour(node->tour, matrix[id * n_nodes + i])) {
                    cout << "NODE " << i << endl;
                    
                    double new_bound_value = newBound(cities[id], cities[i], node->lower_bound);
                    cout << "New bound value: " << new_bound_value << endl;
                    if(new_bound_value > best_tour_cost) {
                        continue;
                    }
                    cout << "CHECK 4" << endl;
                    
                    Node newNode = (Node) malloc(sizeof(struct node));
                    newNode->tour = node->tour;
                    newNode->tour.push_back(i);
                    newNode->cost = node->cost + matrix[id * n_nodes + i];
                    newNode->lower_bound = new_bound_value;
                    newNode->length = node->length + 1;
                    newNode->current_city = i;
                    queue.push(newNode);
                    nodes_created.push_back(newNode);
                    //print node tour
                    for (int j = 0; j < newNode->tour.size(); j++) {
                        cout << newNode->tour.at(j) << " ";
                    }
                    cout << endl;
                    cout << "CHECK 5" << endl;
                }
            }
        }
    }
    
    for (int i = 0; i < nodes_created.size(); i++) {
        free(nodes_created.at(i));
    }
    return;
}

void readInputFile(string inputFile) {
    cout << "Reading input file" << endl;
    ifstream myFile;
    myFile.open(inputFile, ios::in);
    if(!myFile) {
        cout << "NO SOLUTION";
        exit(1);
    }
    int coord1, coord2;
    double dist;

    myFile >> n_nodes >> n_edges;

    for (int i = 0; i < n_nodes; i++) {
        City city = (City) malloc(sizeof(struct city));
        city->id = i;
        city->min1 = -1;
        city->min2 = -1;
        cities.push_back(city);
    }
    cout << "For" << endl;
    matrix = (double *) calloc(n_nodes * n_nodes, sizeof(double));
    cout << "Matrix initialized" << endl;
    // initialize matrix and every city with min1 and min2 values
    while (myFile >> coord1 >> coord2 >> dist) {
        matrix[coord1 * n_nodes + coord2] = dist;
        matrix[coord2 * n_nodes + coord1] = dist;

        if (cities[coord1]->min1 == -1) { cities[coord1]->min1 = dist; } // if min1 is not set, set it
        else if (cities[coord1]->min2 == -1) {  // if min2 is not set, set it
            if (dist < cities[coord1]->min1) { // if dist is smaller than min1, set min2 to min1 and min1 to dist
                cities[coord1]->min2 = cities[coord1]->min1;
                cities[coord1]->min1 = dist;
            }
            else { // if dist is bigger than min1, set min2 to dist
                cities[coord1]->min2 = dist;
            }
        }
        else { // if min1 and min2 are set, check if dist is smaller than min1 or min2
            if (dist < cities[coord1]->min1) { // if dist is smaller than min1, set min2 to min1 and min1 to dist
                cities[coord1]->min2 = cities[coord1]->min1;
                cities[coord1]->min1 = dist;
            }
            else if (dist < cities[coord1]->min2) { // if dist is smaller than min2, set min2 to dist
                cities[coord1]->min2 = dist;
            }
        }

        if (cities[coord2]->min1 == -1) { cities[coord2]->min1 = dist; }
        else if (cities[coord2]->min2 == -1) { 
            if (dist < cities[coord2]->min1) {
                cities[coord2]->min2 = cities[coord2]->min1;
                cities[coord2]->min1 = dist;
            }
            else {
                cities[coord2]->min2 = dist;
            }
        }
        else {
            if (dist < cities[coord2]->min1) {
                cities[coord2]->min2 = cities[coord2]->min1;
                cities[coord2]->min1 = dist;
            }
            else if (dist < cities[coord2]->min2) {
                cities[coord2]->min2 = dist;
            }
        }
    }
    
    myFile.close();

    //print matrix
    for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
            cout << matrix[j * n_nodes + i] << " ";
        }
        cout << endl;
    }

    //print cities
    for (int i = 0; i < n_nodes; i++) {
        cout << "City " << i << ": " << cities[i]->min1 << " " << cities[i]->min2 << endl;
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
    //print_result();       // to the stdout! TODO

    free(matrix);
    for (int i = 0; i < cities.size(); i++) {
        free(cities.at(i));
    }
    return 0;
}