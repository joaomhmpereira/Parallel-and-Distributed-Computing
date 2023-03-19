/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                   TRAVELING SALESMAN PROBLEM                  */
/*                    CPD 2022/2023 - Group 15                   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <omp.h>
#include "nqueue/queue.hpp"

using namespace std;

/* STRUCTURES */

typedef struct city { /* contain's a city's information */
    int id;                /* city number */
    double min1;           /* 1st smallest edge */
    double min2;           /* 2nd smallest edge */
} * City;

typedef struct node { /* a node in the search tree */
    int * tour;            /* tour so far (its ancestors) */
    double cost;           /* tour cost so far */
    double lower_bound;    /* node's lower bound */
    int length;            /* tour size */
} * Node;

/* determines whether a node is "smaller" than another node */
struct cmp_op { bool operator()(const Node& n1, const Node& n2) {
    return (n1->lower_bound > n2->lower_bound) 
    || (n1->lower_bound == n2->lower_bound && n1->tour[n1->length - 1] > n2->tour[n2->length - 1]);
}};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* calculates the initial lower bound (node 0) */
double initialLowerBound(int n_cities, City * cities) {
    double LB = 0;
    for (int i = 0; i < n_cities; i++) {
        LB += cities[i]->min1 + cities[i]->min2;
    }
    return LB / 2;
}

/* calculates the lower bound of a node (coming from a source) */
double newBound(City source, City target, double lower_bound, int n_cities, double * matrix){
    double cost = matrix[source->id * n_cities + target->id];
    double cf = cost >= source->min2 ? source->min2 : source->min1;
    double ct = cost >= target->min2 ? target->min2 : target->min1;

    return lower_bound + cost - (cf + ct) / 2;
}

/* updates the min1 and min2 of the city given by coord */
void update_mins(int coord, double dist, City ** cities) {
    if ((*cities)[coord]->min1 == -1) { (*cities)[coord]->min1 = dist; } // if min1 is not set, set it
    else if ((*cities)[coord]->min2 == -1) {  // if min2 is not set, set it
        if (dist < (*cities)[coord]->min1) { // if dist is smaller than min1, set min2 to min1 and min1 to dist
            (*cities)[coord]->min2 = (*cities)[coord]->min1;
            (*cities)[coord]->min1 = dist;
        }
        else { // if dist is bigger than min1, set min2 to dist
            (*cities)[coord]->min2 = dist;
        }
    }
    else { // if min1 and min2 are set, check if dist is smaller than min1 or min2
        if (dist < (*cities)[coord]->min1) { // if dist is smaller than min1, set min2 to min1 and min1 to dist
            (*cities)[coord]->min2 = (*cities)[coord]->min1;
            (*cities)[coord]->min1 = dist;
        }
        else if (dist < (*cities)[coord]->min2) { // if dist is smaller than min2, set min2 to dist
            (*cities)[coord]->min2 = dist;
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void tsp(double * best_tour_cost, int max_value, int n_cities, int ** best_tour, double * matrix, City * cities){
    Node root = (Node) calloc(1, sizeof(struct node));

    root->tour = (int *) calloc(1, sizeof(int));
    root->tour[0] = 0;

    root->cost = 0;
    root->lower_bound = initialLowerBound(n_cities, cities);
    root->length = 1;

    //printf("Alloc node id=%d\n", root->tour[root->length - 1]);

    // initialize priority queue
    // PriorityQueue<Node, cmp_op> queue;
    // queue.push(root);

    (*best_tour_cost) = max_value;
    (*best_tour) = (int *) calloc(n_cities + 1, sizeof(int));
    
    int n_threads;
    int n_emptied = 0;
    
    #pragma omp parallel 
    {
        //printf("Alloc tour_nodes for thread %d\n", omp_get_thread_num());
        bool * tour_nodes = (bool *) calloc(n_cities, sizeof(bool));
        PriorityQueue<Node, cmp_op> private_queue;
        bool emptied = false;
        
        #pragma omp master
        {
            n_threads = omp_get_num_threads();
            private_queue.push(root);
        }
            
        while (n_emptied != n_threads) {
            Node node = NULL;
            
            if (!private_queue.empty()){
                node = private_queue.pop();
                if (emptied) {
                    #pragma omp atomic
                    n_emptied -= 1;
                }
                emptied = false;
                printf("!!! Thread %d got the node !!!\n", omp_get_thread_num());
            }
            else {
                printf("!!! Thread %d found the queue empty !!!\n", omp_get_thread_num());
                if (!emptied) {
                    emptied = true;
                    #pragma omp atomic
                    n_emptied += 1;
                }
            }
            
            int id = 0;
            if (!emptied) 
                id = node->tour[node->length - 1];

            bool aux = false;
            #pragma omp critical(best_tour_lock)
            if (!emptied)
                aux = (node->lower_bound >= (*best_tour_cost));

            // All remaining nodes worse than best
            if (!emptied && aux) {
                //printf("Free node id=%d\n", node->tour[node->length - 1]);
                free(node->tour);
                free(node);

                while (!private_queue.empty()) {
                    Node n = private_queue.pop();
                    //printf("Free node id=%d\n", n->tour[n->length - 1]);
                    free(n->tour);
                    free(n);
                }

                emptied = true;
                #pragma omp atomic
                n_emptied += 1;
            }

            // Tour complete, check if it is best
            if (!emptied && node->length == n_cities) {
                #pragma omp critical(best_tour_lock)
                if (node->cost + matrix[id * n_cities + 0] < (*best_tour_cost) && matrix[id * n_cities + 0] >= 0.0001) {
                    int i;
                    for (i = 0; i < n_cities; i++) {
                        (*best_tour)[i] = node->tour[i];
                    }
                    (*best_tour)[i] = 0;

                    (*best_tour_cost) = node->cost + matrix[id * n_cities + 0];
                }
            } 
            else if (!emptied) {
                printf("Thread %d is going to check neighbours\n", omp_get_thread_num());
                for (int i = 0; i < node->length; i++) {
                    tour_nodes[node->tour[i]] = true;
                }
                
                #pragma omp for
                for (int i = 0; i < n_cities; i++) {
                    printf("Thread %d is checking node %d\n", omp_get_thread_num(), i);
                    if (matrix[id * n_cities + i] >= 0.0001 && !tour_nodes[i]) {
                        double new_bound_value = newBound(cities[id], cities[i], node->lower_bound, n_cities, matrix);

                        bool aux_2 = false;
                        #pragma omp critical(best_tour_lock)
                        aux_2 = (new_bound_value > (*best_tour_cost));

                        if (!aux_2) {
                            Node newNode = (Node) calloc(1, sizeof(struct node));
                            newNode->tour = (int *) calloc(node->length + 1, sizeof(int));
                            for (int j = 0; j < node->length; j++) {
                                newNode->tour[j] = node->tour[j];
                            }
                            newNode->tour[node->length] = i;
                            newNode->cost = node->cost + matrix[id * n_cities + i];
                            newNode->lower_bound = new_bound_value;
                            newNode->length = node->length + 1;

                            private_queue.push(newNode);
                            printf("Thread %d pushed node %d onto its stack\n", omp_get_thread_num(), i);
                        }
                    }
                }
                memset(tour_nodes, false, n_cities * sizeof(bool));
            }

            if (!emptied) {
                free(node->tour);
                free(node);
            }
        }

        free(tour_nodes);
    }
    
    return;
}

void readInputFile(string inputFile, int * n_cities, int * n_roads, double ** matrix, City ** cities) {
    ifstream myFile;
    myFile.open(inputFile, ios::in);
    if(!myFile) {
        cerr << "error opening file";
        exit(1);
    }
    int coord1, coord2;
    double dist;

    myFile >> *n_cities >> *n_roads;

    *cities = (City *) calloc((*n_cities), sizeof(City));

    for (int i = 0; i < (*n_cities); i++) {
        City city = (City) calloc(1, sizeof(struct city));
        city->id = i;
        city->min1 = -1;
        city->min2 = -1;
        (*cities)[i] = city;
    }

    *matrix = (double *) calloc((*n_cities) * (*n_cities), sizeof(double));

    // initialize matrix and every city with min1 and min2 values
    while (myFile >> coord1 >> coord2 >> dist) {
        (*matrix)[coord1 * (*n_cities) + coord2] = dist;
        (*matrix)[coord2 * (*n_cities) + coord1] = dist;

        update_mins(coord1, dist, cities);
        update_mins(coord2, dist, cities);
    }

    myFile.close();
}

void print_result(double best_tour_cost, double max_value, int n_cities, int * best_tour) {
    if (best_tour_cost > max_value || best_tour[1] == 0) {
        cout << "NO SOLUTION" << endl;
    } else {
        fprintf(stdout, "%.1f\n", best_tour_cost);
        for (int i = 0; i < n_cities + 1; i++) {
            if (i == n_cities) 
                cout << best_tour[i];
            else
                cout << best_tour[i] << " ";
        }
        cout << endl;
    }
    free(best_tour);
}

int main(int argc, char *argv[]) {
    double max_value;       /* max. solution cost accepted */
    double best_tour_cost;  /* best solution cost */
    double * matrix = NULL; /* adjacencies matrix */
    double exec_time;       /* execution time */

    int * best_tour;        /* best solution path */
    int n_cities;           /* number of cities */
    int n_roads;            /* number of roads - links between cities */

    City * cities = NULL;   /* cities in the map */

    /* * * * * * * * * * * * * * * * * * */

    if (argc < 3) {
        cerr << "Usage: tsp <input_file> <max-value>\n";
        return 1;
    }

    readInputFile(argv[1], &n_cities, &n_roads, &matrix, &cities);

    max_value = stod(argv[2]);

    exec_time = -omp_get_wtime();

    tsp(&best_tour_cost, max_value, n_cities, &best_tour, matrix, cities);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);

    print_result(best_tour_cost, max_value, n_cities, best_tour);

    free(matrix);
    for (int i = 0; i < n_cities; i++) {
        free(cities[i]);
    }
    free(cities);
    return 0;
}