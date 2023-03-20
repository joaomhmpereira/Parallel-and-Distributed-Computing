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

bool there_are_nodes(Node * nodes_in_processing, int n_threads) {
    for (int i = 0; i < n_threads; i++) {
        if (nodes_in_processing[i]) return true;
    }
    return false;
}

void merge_best_tour(double * array_best_tour_cost, int ** array_best_tour, int n_threads, int n_cities, PriorityQueue<Node, cmp_op> * queue_array) {
    // discover best tour cost
    double best_tour_cost = array_best_tour_cost[0];
    int best_tour_index = 0;
    for (int i = 1; i < n_threads; i++) {
        if (array_best_tour_cost[i] < best_tour_cost) {
            best_tour_cost = array_best_tour_cost[i];
            best_tour_index = i;
        }
    }

    // for all other threads, update best 
    for (int i = 0; i < n_threads; i++) {
        if (i != best_tour_index) {
            for (int j = 0; j < n_cities + 1; j++) {
                array_best_tour[i][j] = array_best_tour[best_tour_index][j];
            }
        }
        array_best_tour_cost[i] = best_tour_cost;

        Node node = queue_array[i].top();
        bool worse_than_best = (node->lower_bound >= best_tour_cost);
                
        if (worse_than_best) {
            while (!queue_array[i].empty()) {
                Node n = queue_array[i].pop();
                free(n->tour);
                free(n);
            }
        }
    }
}


void tsp(double * best_tour_cost, int max_value, int n_cities, int ** best_tour, double * matrix, City * cities){
    Node root = (Node) calloc(1, sizeof(struct node));

    root->tour = (int *) calloc(1, sizeof(int));
    root->tour[0] = 0;

    root->cost = 0;
    root->lower_bound = initialLowerBound(n_cities, cities);
    root->length = 1;

    (*best_tour_cost) = max_value;
    (*best_tour) = (int *) calloc(n_cities + 1, sizeof(int));
    
    omp_lock_t best_tour_lock;
    omp_lock_t * queue_locks;
    
    omp_init_lock(&best_tour_lock);
    PriorityQueue<Node, cmp_op> * queue_array;
    
    int n_threads;    
    Node * nodes_in_processing =  (Node *) calloc(4, sizeof(Node));
    bool thereAreNodes = true;

    double * array_best_tour_cost;
    int ** array_best_tour;
    
    //comecam aqui 4 threads
    #pragma omp parallel shared(nodes_in_processing, thereAreNodes, queue_locks, queue_array, array_best_tour_cost, array_best_tour, best_tour_cost)
    {
        bool * tour_nodes = (bool *) calloc(n_cities, sizeof(bool));
        const int thread_id = omp_get_thread_num();
        n_threads = omp_get_num_threads(); 
        bool freed = false;
        int iterations = 0;

        #pragma omp single
        {
            array_best_tour_cost = (double *) calloc(n_threads, sizeof(double));
            array_best_tour = (int **) calloc(n_threads, sizeof(int *));
            for (int i = 0; i < n_threads; i++) {
                array_best_tour_cost[i] = max_value;
                array_best_tour[i] = (int *) calloc(n_cities + 1, sizeof(int));
            }

            queue_array = (PriorityQueue<Node, cmp_op> *) calloc(n_threads, sizeof(PriorityQueue<Node, cmp_op>));
            queue_locks = (omp_lock_t *) calloc(n_threads, sizeof(omp_lock_t));
            for (int l = 0; l < n_threads; l++) {
                omp_init_lock(&queue_locks[l]);
            }
        }
        
        nodes_in_processing[0] = root;
        for (int i = 1; i < n_threads; i++) {
            nodes_in_processing[i] = NULL;
        }
        
        #pragma omp barrier

        while (thereAreNodes) {
            //printf("Thread %d in\n", thread_id);
            #pragma omp single
            {
                iterations++;
                if (iterations % 100 == 0) {
                    merge_best_tour(array_best_tour_cost, array_best_tour, n_threads, n_cities, queue_array);
                }   
            }
            Node node = nodes_in_processing[thread_id];
            if (node) {
                //printf("Thread %d is processing node %d\n", thread_id, node->tour[node->length - 1]);
			
                int id = node->tour[node->length - 1];

                bool worse_than_best = (node->lower_bound >= array_best_tour_cost[thread_id]);
                
                // all remaining nodes worse than best
                if (worse_than_best) {
                    //printf("Thread %d is freeing node %d level %d\n", thread_id, node->tour[node->length - 1], node->length);
                    free(node->tour);
                    free(node);

                    nodes_in_processing[thread_id] = NULL;
                    omp_set_lock(&queue_locks[thread_id]);
                    while (!queue_array[thread_id].empty()) {
                        Node n = queue_array[thread_id].pop();
                        //printf("Thread %d is freeing node %d level %d\n", thread_id, n->tour[n->length - 1], n->length);
                        free(n->tour);
                        free(n);
                    }
                    omp_unset_lock(&queue_locks[thread_id]);
                    freed = true;
                }
                else {
                    // Tour complete, check if it is best
                    if (node->length == n_cities) {
                        if (node->cost + matrix[id * n_cities + 0] < array_best_tour_cost[thread_id] && matrix[id * n_cities + 0] >= 0.0001) {
                            int i;
                            for (i = 0; i < n_cities; i++) {
                                array_best_tour[thread_id][i] = node->tour[i];
                            }
                            array_best_tour[thread_id][i] = 0;

                            array_best_tour_cost[thread_id] = node->cost + matrix[id * n_cities + 0];
                        }
                    } 
                }
            }

            #pragma omp barrier 

            for (int j = 0; j < n_threads; j++) {
                Node shared_node = nodes_in_processing[j];

                if (shared_node) {
                    //printf("Thread %d will now process its portion of node %d\n", thread_id, shared_node->tour[shared_node->length - 1]);

                    for (int i = 0; i < shared_node->length; i++) {
                        tour_nodes[shared_node->tour[i]] = true;
                    }

                    int id = shared_node->tour[shared_node->length - 1];
                    #pragma omp for
                    for (int i = 0; i < n_cities; i++) {
                        if (matrix[id * n_cities + i] != 0 && !tour_nodes[i]) {
                            double new_bound_value = newBound(cities[id], cities[i], shared_node->lower_bound, n_cities, matrix);
                            
                            bool better_than_best = (new_bound_value <= array_best_tour_cost[thread_id]);

                            if (better_than_best) {
                                Node newNode = (Node) calloc(1, sizeof(struct node));
                                newNode->tour = (int *) calloc(shared_node->length + 1, sizeof(int));
                                for (int j = 0; j < shared_node->length; j++) {
                                    newNode->tour[j] = shared_node->tour[j];
                                }
                                newNode->tour[shared_node->length] = i;
                                newNode->cost = shared_node->cost + matrix[id * n_cities + i];
                                newNode->lower_bound = new_bound_value;
                                newNode->length = shared_node->length + 1;
                                //printf("Thread %d is storing node %d level %d\n", thread_id, newNode->tour[newNode->length - 1], newNode->length);
                                int random_index = rand() % n_threads;
                                omp_set_lock(&queue_locks[random_index]);
                                //printf("storing in index %d\n", random_index);
                                queue_array[random_index].push(newNode);
                                omp_unset_lock(&queue_locks[random_index]);
                            }
                        }
                    }
                    
                    memset(tour_nodes, false, n_cities * sizeof(bool));
                }

            }
            
            #pragma omp barrier

            if (node && !freed) {
                //printf("Thread %d is freeing node %d level %d\n", thread_id, node->tour[node->length - 1], node->length);
                free(node->tour);
                free(node);
            }

            omp_set_lock(&queue_locks[thread_id]);
            if (!queue_array[thread_id].empty())
                nodes_in_processing[thread_id] = queue_array[thread_id].pop();
            else
                nodes_in_processing[thread_id] = NULL;
            omp_unset_lock(&queue_locks[thread_id]);

            freed = false;
            
            //printf("Thread %d out\n", thread_id);
            #pragma omp barrier
            #pragma omp single
            {
                thereAreNodes = there_are_nodes(nodes_in_processing, n_threads);
            }
            #pragma omp barrier
        }
        //printf("Thread %d finished, no more nodes\n", thread_id);

        free(tour_nodes);
    }

    for (int i = 0; i < n_threads; i++) {
        if (array_best_tour_cost[i] < (*best_tour_cost)) {
            (*best_tour_cost) = array_best_tour_cost[i];
            for (int j = 0; j < n_cities + 1; j++) {
                (*best_tour)[j] = array_best_tour[i][j];
            }
        }
        free(array_best_tour[i]);
        omp_destroy_lock(&queue_locks[i]);
    }

    free(queue_locks);
    free(array_best_tour);
    free(array_best_tour_cost);
    free(queue_array);
    free(nodes_in_processing);
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

    srand (time(NULL)); // initialize random seed

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
