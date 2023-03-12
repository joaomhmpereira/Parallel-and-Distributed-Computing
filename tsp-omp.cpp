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

/* determines if city v is in the tour */
bool in_tour(int * tour, int length, int v) {
    for (int i = 0; i < length; i++) {
        if (tour[i] == v) { return true; }
    }
    return false;
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

bool toReturn(bool to_return[], int size) {
    for(int i = 0; i < size; i++){
        if (to_return[i]) {
            printf("< thread %d is still running >\n", i);
            return true;
        }
    }
    printf("All threads are done\n");
    return false;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void tsp(double * best_tour_cost, int max_value, int n_cities, int ** best_tour, double * matrix, City * cities, int n_roads){
    Node root = (Node) calloc(1, sizeof(struct node));

    root->tour = (int *) calloc(n_cities, sizeof(int));
    root->tour[0] = 0;

    root->cost = 0;
    root->lower_bound = initialLowerBound(n_cities, cities);
    root->length = 1;

    // initialize priority queue
    PriorityQueue<Node, cmp_op> queue;
    queue.push(root);

    (*best_tour_cost) = max_value;
    (*best_tour) = (int *) calloc(n_cities + 1, sizeof(int));

    /* TODO: o que e melhor: verificar esta flag em todo o lado ou meter nested ifs? */

    bool * to_return;
    int n_threads;
    #pragma omp parallel shared(to_return, n_threads)
    {
    n_threads = omp_get_num_threads();
    to_return = (bool*) calloc(n_threads, sizeof(int));

    const int thread_id = omp_get_thread_num();
    while (!toReturn(to_return, n_threads)) {
        Node node;
        #pragma omp critical(queue_lock)
        if (!queue.empty()){
            printf("--- Thread %d got the node ---\n", thread_id);
            to_return[thread_id] = false;
            node = queue.pop();
        }
        else {
            printf("--- Thread %d is setting flag to false ---\n", thread_id);
            to_return[thread_id] = true;
        }

        int id;
        if (!to_return[thread_id])
            id = node->tour[node->length - 1];
        
        // All remaining nodes worse than best
        bool aux = false;
        #pragma omp critical(best_tour_lock)
        if (!to_return[thread_id]) {
            aux = (node->lower_bound >= (*best_tour_cost));
            //printf("AUX INSIDE: %d Thread id: %d\n", aux, thread_id);
        }
        //printf("AUX OUTSIDE: %d Thread id: %d\n", aux, thread_id);
        
        if (!(to_return[thread_id]) && aux) {
            free(node->tour);
            free(node);

            #pragma omp critical(queue_lock)
            while (!queue.empty()) {
                Node n = queue.pop();
                free(n->tour);
                free(n);
            }
            
            printf("--- Setting all thread flags to false --- (made by thread: %d)\n", thread_id);
            for (int i = 0; i < n_threads; i++)
                to_return[i] = true;
        }
        else if (!to_return[thread_id]) {
            // Tour complete, check if it is best
            if (node->length == n_cities) {
                /* CRITICAL: need to check and update only one at each time */
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
            else {
                /* PARALLEL FOR */
                /* testar condicoes:
                    - n_roads
                    - n_cities
                    - max_value
                */
                //#pragma omp parallel for //schedule(dynamic, 2)
                for (int i = 0; i < n_cities; i++) {
                    if (matrix[id * n_cities + i] != 0 && !in_tour(node->tour, node->length, i)) {
                        double new_bound_value = newBound(cities[id], cities[i], node->lower_bound, n_cities, matrix);
                        
                        /* CRITICAL: cannot check while someone is writing here */
                        bool aux_2 = false;
                        #pragma omp critical(best_tour_lock)
                        aux_2 = (new_bound_value <= (*best_tour_cost));
                        
                        if (aux_2) {
                            Node newNode = (Node) calloc(1, sizeof(struct node));
                            newNode->tour = (int *) calloc(n_cities, sizeof(int));
                            int j;
                            for (j = 0; j < node->length; j++) {
                                newNode->tour[j] = node->tour[j];
                            }
                            newNode->tour[j] = i;
                            newNode->cost = node->cost + matrix[id * n_cities + i];
                            newNode->lower_bound = new_bound_value;
                            newNode->length = node->length + 1;

                            /* CRITICAL: access to the queue */
                            #pragma omp critical(queue_lock)
                            queue.push(newNode);
                        }
                    }
                }
            }

            free(node->tour);
            free(node);
            
        }
    }        
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

    tsp(&best_tour_cost, max_value, n_cities, &best_tour, matrix, cities, n_roads);

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