/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                   TRAVELING SALESMAN PROBLEM                  */
/*                    CPD 2022/2023 - Group 15                   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstddef>
#include <mpi.h>
#include "nqueue/queue.hpp"

using namespace std;


#define COST_TAG 1
#define TOUR_TAG 2
#define WHITE 3
#define BLACK 4
#define TOKEN_TAG 5
#define IDLE_TAG 6
#define WORK_TAG 7
#define INFO_TAG 8
#define SOLUTION_TAG 9

/* STRUCTURES */

typedef struct city { /* contain's a city's information */
    int id;                /* city number */
    double min1;           /* 1st smallest edge */
    double min2;           /* 2nd smallest edge */
} * City;

typedef struct node { /* a node in the search tree */
    int length;            /* tour size */
    int * tour;            /* tour so far (its ancestors) */
    double cost;           /* tour cost so far */
    double lower_bound;    /* node's lower bound */
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void balance_LB(int me, int other, int my_load, int other_load, int delta_LB, PriorityQueue<Node, cmp_op> queue, MPI_Comm comm, int * color, int n_cities) {
    if (other_load - my_load > delta_LB) {
        for (int i = 0; i < 1; i++) {
            if (queue.empty()) return; 
            if (other < me) *(color) = BLACK;
            Node node, aux;
            aux = queue.pop();
            node = queue.pop();
            queue.push(aux);

            /* serialize and send node */
            int length = node->length;
            char * buffer = (char *) calloc(1, sizeof(double)*2 + sizeof(int) + sizeof(int)*length);
            memcpy(buffer, &length, sizeof(int));
            memcpy(buffer + sizeof(int), node->tour, sizeof(int) * length);
            memcpy(buffer + sizeof(int) + sizeof(int) * length, &node->cost, sizeof(double));
            memcpy(buffer + sizeof(int) + sizeof(int) * length + sizeof(double), &node->lower_bound, sizeof(double));

            MPI_Send(buffer, sizeof(double)*2 + sizeof(int) + sizeof(int)*length, MPI_BYTE, other, WORK_TAG, comm);
            
            free(node->tour);
            free(node);
            free(buffer);
        }
    }
}

void balance_amount(int me, int other, int my_load, int other_load, int delta_amount, double send_amount_rate, PriorityQueue<Node, cmp_op> queue, MPI_Comm comm, int * color, int n_cities) {
    if (my_load - other_load > delta_amount) {
        int n = (int) ((my_load - other_load) * send_amount_rate);
        for (int i = 0; i < 1; i++) {
            if (queue.empty()) return; 
            if (other < me) *(color) = BLACK;

            Node node = queue.get_buffer().back();

            /* serialize and send node */
            int length = node->length;
            char * buffer = (char *) calloc(1, sizeof(double)*2 + sizeof(int) + sizeof(int)*length);
            memcpy(buffer, &length, sizeof(int));
            memcpy(buffer + sizeof(int), node->tour, sizeof(int) * length);
            memcpy(buffer + sizeof(int) + sizeof(int) * length, &node->cost, sizeof(double));
            memcpy(buffer + sizeof(int) + sizeof(int) * length + sizeof(double), &node->lower_bound, sizeof(double));

            MPI_Send(buffer, sizeof(double)*2 + sizeof(int) + sizeof(int)*length, MPI_BYTE, other, WORK_TAG, comm);
            
            queue.get_buffer().pop_back();
            free(buffer);
        }
    }
}

void inform_neighbors(int me, int n_tasks, int my_load_LB, int my_load_amount, MPI_Comm comm, int * color) {
    for (int i = 0; i < n_tasks; i++) {
        if (i != me) {
            if (i < me) *(color) = BLACK;
            int info[2] = {my_load_LB, my_load_amount};
            MPI_Send(info, 2, MPI_INT, i, INFO_TAG, comm);
        }
    }
}

void tsp(double * best_tour_cost, int max_value, int n_cities, int ** best_tour, double * matrix, City * cities, int n_tasks, int id){
    /* INITIALIZATION */
    bool finished = false;
    bool * tour_nodes = (bool *) calloc(n_cities, sizeof(bool));
    
    /* RING TERMINATION */
    int color = BLACK;
    int token = BLACK;
    int next_rank = (id + 1) % n_tasks;
    int prev_rank = (id + n_tasks - 1) % n_tasks;
    MPI_Request request_broadcast, request_receive, request_allreduce;
    int terminate = -1;     
    int flag, flag2;
    int empty_queue = 2;

    //printf("[TASK %d] Executing sequential part...\n", id);
    /* init root */
    Node root = (Node) calloc(1, sizeof(struct node));
    root->tour = (int *) calloc(1, sizeof(int));
    root->tour[0] = 0;
    root->cost = 0;
    root->lower_bound = initialLowerBound(n_cities, cities);
    root->length = 1;

    /* init best */
    (*best_tour_cost) = max_value;
    (*best_tour) = (int *) calloc(n_cities + 1, sizeof(int));

    /* SEQUENTIAL PART */
    bool * tour_nodes_init = (bool *) calloc(n_cities, sizeof(bool));
    PriorityQueue<Node, cmp_op> initial_queue;
    initial_queue.push(root);
    while (initial_queue.size() < 100000 || initial_queue.empty()) {
        Node node = initial_queue.pop();
        int node_id = node->tour[node->length - 1];

        // All remaining nodes worse than best
        if (node->lower_bound >= (*best_tour_cost)) {
            free(node->tour);
            free(node);
            while (!initial_queue.empty()) {
                Node n = initial_queue.pop();
                free(n->tour);
                free(n);
            }
            free(tour_nodes_init);
            break;
        }

        // Tour complete, check if it is best
        if (node->length == n_cities) {
            if (node->cost + matrix[node_id * n_cities + 0] < (*best_tour_cost) && matrix[node_id * n_cities + 0] >= 0.0001) {
                int i;
                for (i = 0; i < n_cities; i++) {
                    (*best_tour)[i] = node->tour[i];
                }
                (*best_tour)[i] = 0;
                (*best_tour_cost) = node->cost + matrix[node_id * n_cities + 0];
            }
        } 
        else {
            for (int i = 0; i < node->length; i++) {
                tour_nodes_init[node->tour[i]] = true;
            }
            for (int i = 0; i < n_cities; i++) {
                if (matrix[node_id * n_cities + i] != 0 && !tour_nodes_init[i]) {
                    double new_bound_value = newBound(cities[node_id], cities[i], node->lower_bound, n_cities, matrix);
                    if(new_bound_value > (*best_tour_cost)) {
                        continue;
                    }
                    Node newNode = (Node) calloc(1, sizeof(struct node));
                    newNode->tour = (int *) calloc(n_cities + 1, sizeof(int));
                    for (int j = 0; j < node->length; j++) {
                        newNode->tour[j] = node->tour[j];
                    }
                    newNode->tour[node->length] = i;
                    newNode->cost = node->cost + matrix[node_id * n_cities + i];
                    newNode->lower_bound = new_bound_value;
                    newNode->length = node->length + 1;
                    initial_queue.push(newNode);
                }
            }
            memset(tour_nodes_init, false, n_cities * sizeof(bool));
        }
        free(node->tour);
        free(node);
    }  

    if (initial_queue.empty()) return;

    PriorityQueue<Node, cmp_op> queue;
    int i = 0;
    while (!initial_queue.empty()) {
        if (i % n_tasks != id) {
            Node n = initial_queue.pop();
            free(n->tour);
            free(n);
        }
        else {            
            queue.push(initial_queue.pop());
        }
        i++;
    }

    /* PARALLEL PART - EACH PROCESS WORKS ON ITS SUBPROBLEM */
    bool broadcasted = false, reduced = false;
    bool * idle_processes;
    MPI_Status inform_neighbors_status;
    int num_idle_processes = 0;
    MPI_Status solution_status;
    MPI_Status balance_status;
    MPI_Status bcast_status;
    int iterations = 0;
    double min_cost;
    
    /* minimum delay since last communication: approx. time to process a node */
    double delay = 1;
    /* deltas and rates */
    int delta_LB = 1;
    int delta_amount = 3;
    double send_amount_rate = 0.65;
    int send_LB_rate = 3;

    /* last time processes interchanged each type of message */
    double last_comm = MPI_Wtime();
    /* amount of load each process has (depending on type) */
    int * load_LB;
    int * load_amount;

    idle_processes = (bool *) calloc(n_tasks, sizeof(bool));
    load_LB = (int *) calloc(n_tasks, sizeof(int));
    load_amount = (int *) calloc(n_tasks, sizeof(int));

    while (1) {
        if (MPI_Wtime() - last_comm > delay) {
            if (!queue.empty())
                inform_neighbors(id, n_tasks, queue.top()->lower_bound, queue.size(), MPI_COMM_WORLD, &color);
            else
                inform_neighbors(id, n_tasks, *(best_tour_cost), 0, MPI_COMM_WORLD, &color);
            last_comm = MPI_Wtime();
        } 
        
        MPI_Iprobe(MPI_ANY_SOURCE, INFO_TAG, MPI_COMM_WORLD, &flag, &inform_neighbors_status);
        if (flag) {
            int source = inform_neighbors_status.MPI_SOURCE;
            int info[2];
            MPI_Recv(&info, 2, MPI_INT, source, INFO_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            load_LB[source] = info[0];
            load_amount[source] = info[1];
            if (!queue.empty())
                balance_LB(id, source, queue.top()->lower_bound, info[0], delta_LB, queue, MPI_COMM_WORLD, &color, n_cities);
            else
                balance_LB(id, source, *(best_tour_cost), info[0], delta_LB, queue, MPI_COMM_WORLD, &color, n_cities);
            balance_amount(id, source, queue.size(), info[1], delta_amount, send_amount_rate, queue, MPI_COMM_WORLD, &color, n_cities);
        }

        flag = false;
        MPI_Iprobe(MPI_ANY_SOURCE, WORK_TAG, MPI_COMM_WORLD, &flag, &balance_status);
        if (flag) {
            fprintf(stderr, "Process %d: received work request\n", id);
            int source = balance_status.MPI_SOURCE;
            Node subproblem = (Node) calloc(1, sizeof(struct node));
            size_t size = sizeof(double) * 2 + sizeof(int) + sizeof(int) * (n_cities + 1);
            char * placeholder_buffer = (char *) calloc(size, sizeof(char));

            /* receive and deserialize node */
            MPI_Recv(placeholder_buffer, size, MPI_BYTE, source, WORK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            memcpy(&subproblem->length, placeholder_buffer, sizeof(int));

            subproblem->tour = (int *) calloc(subproblem->length, sizeof(int));
            memcpy(subproblem->tour, placeholder_buffer + sizeof(int), sizeof(int) * subproblem->length);
            memcpy(&subproblem->cost, placeholder_buffer + sizeof(int) + sizeof(int) * subproblem->length, sizeof(double));
            memcpy(&subproblem->lower_bound, placeholder_buffer + sizeof(int) + sizeof(int) * subproblem->length + sizeof(double), sizeof(double));
            
            free(placeholder_buffer);
            
            //fprintf(stderr, "111Process %d: received work from %d\n", id, source);
            
            //fprintf(stderr, "receiving tour\n");
            //MPI_Recv(subproblem->tour, n_cities + 1, MPI_INT, source, WORK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            queue.push(subproblem);
            // fprintf(stderr, "222Process %d: received work from %d\n", id, source);
            // fprintf(stderr, "Node: cost = %f, lower_bound = %f, length = %d, tour = %d\n", subproblem->cost, subproblem->lower_bound, subproblem->length, subproblem->tour);

            // for(int i = 0; i < subproblem->length; i++) {
            //     fprintf(stderr, "%d ", subproblem->tour[i]);
            // }
            // fprintf(stderr, "\n");
        }

        flag = false;
        MPI_Iprobe(MPI_ANY_SOURCE, SOLUTION_TAG, MPI_COMM_WORLD, &flag, &solution_status);
        if (flag) {
            int source = solution_status.MPI_SOURCE;
            double cost;
            MPI_Recv(&cost, 1, MPI_DOUBLE, source, SOLUTION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (cost < *(best_tour_cost)) {
                *(best_tour_cost) = cost;
            }
        }

        if (!queue.empty()) {
            Node node = queue.pop();
            int node_id = node->tour[node->length - 1];

            // All remaining nodes worse than best
            if (node->lower_bound >= (*best_tour_cost)) {
                free(node->tour);
                free(node);
                while (!queue.empty()) {
                    Node n = queue.pop();
                    free(n->tour);
                    free(n);
                }
                continue;
            }

            // Tour complete, check if it is best
            if (node->length == n_cities) {
                if (node->cost + matrix[node_id * n_cities + 0] < (*best_tour_cost) && matrix[node_id * n_cities + 0] >= 0.0001) {
                    int i;
                    for (i = 0; i < n_cities; i++) {
                        (*best_tour)[i] = node->tour[i];
                    }
                    (*best_tour)[i] = 0;

                    (*best_tour_cost) = node->cost + matrix[node_id * n_cities + 0];

                    for (int i = 0; i < n_tasks; i++) {
                        if (i != id) {
                            if (i < id) color = BLACK;
                            MPI_Send(best_tour_cost, 1, MPI_DOUBLE, i, SOLUTION_TAG, MPI_COMM_WORLD);
                        }
                    }
                }
            } 
            else {
                
                for (int i = 0; i < node->length; i++) {
                    tour_nodes[node->tour[i]] = true;
                }
                
                for (int i = 0; i < n_cities; i++) {
                    if (matrix[node_id * n_cities + i] != 0 && !tour_nodes[i]) {
                        double new_bound_value = newBound(cities[node_id], cities[i], node->lower_bound, n_cities, matrix);
                        if(new_bound_value > (*best_tour_cost)) {
                            continue;
                        }
                        Node newNode = (Node) calloc(1, sizeof(struct node));
                        newNode->tour = (int *) calloc(node->length + 1, sizeof(int));
                        for (int j = 0; j < node->length; j++) {
                            newNode->tour[j] = node->tour[j];
                        }
                        newNode->tour[node->length] = i;
                        newNode->cost = node->cost + matrix[node_id * n_cities + i];
                        newNode->lower_bound = new_bound_value;
                        newNode->length = node->length + 1;
                        queue.push(newNode);
                    }
                }

                memset(tour_nodes, false, n_cities * sizeof(bool));

            }
            //fprintf(stderr, "Process %d: FREEING node %d\n", id, node);
            free(node->tour);
            free(node);
        }         
        else {
            if (!broadcasted) {
                //MPI_Ibcast(&empty_queue, 1, MPI_INT, id, MPI_COMM_WORLD, &request_broadcast);
                for (int i = 0; i < n_tasks; i++) {
                    if (i != id) {
                        MPI_Isend(&empty_queue, 1, MPI_INT, i, IDLE_TAG, MPI_COMM_WORLD, &request_broadcast);
                    }
                }
                broadcasted = true;
                num_idle_processes++;
                idle_processes[id] = true;
                //fprintf(stderr, "[TASK %d] ::: broadcasted idle tag! ::: \n", id);
            }
            
            //MPI_Test(&request_broadcast, &flag, &bcast_status);
            // MPI_Iprobe checks if there are any messages with IDLE_TAG ready to be received
            MPI_Iprobe(MPI_ANY_SOURCE, IDLE_TAG, MPI_COMM_WORLD, &flag, &bcast_status);
            if (flag) {
                int source = bcast_status.MPI_SOURCE;
                MPI_Recv(&empty_queue, 1, MPI_INT, source, IDLE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                idle_processes[source] = true;
                num_idle_processes++;
                //fprintf(stderr, "[TASK %d] Received broadcast from %d, Idle: %d\n", id, source, num_idle_processes);
            }
            //flag = 0;
        }
        
        //fprintf(stderr, "[TASK %d] Num idle processes before checking: %d\n", id, num_idle_processes);

        if (num_idle_processes == n_tasks) {
            //fprintf(stderr, "[TASK %d] Initiating ring termination!\n", id);
            color = WHITE;
            if (!id) {
                /* P0 sends a white token to P1 */
                token = WHITE; 
                //fprintf(stderr, "[TASK %d] Sending token to %d : token color -> %d \n", id, next_rank, token);
                MPI_Send(&token, 1, MPI_INT, next_rank, TOKEN_TAG, MPI_COMM_WORLD);
                /* P0 is waiting for a token from P(n_tasks - 1) */
                MPI_Recv(&token, 1, MPI_INT, prev_rank, TOKEN_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //MPI_Test(&request_receive, &flag, MPI_STATUS_IGNORE);
                flag = 1;
                if (flag) {
                    //fprintf(stderr, "[TASK %d] Received token from %d : token color -> %d \n", id, prev_rank, token);
                    /* if P0 receives a black token, it will pass a white token */
                    if (token == BLACK) token = WHITE;
                    /* if P0 receives a white token, computation can terminate */
                    else {
                        terminate = 0;
                    }
                }
            }
            else {
                /* when a process finishes, it waits to receive the token */
                //fprintf(stderr, "[TASK %d] Waiting for token...\n", id);
                MPI_Recv(&token, 1, MPI_INT, prev_rank, TOKEN_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //MPI_Test(&request_receive, &flag, MPI_STATUS_IGNORE);
                flag = 1;
                if (flag) {
                    //fprintf(stderr, "[TASK %d] Received token from %d : token color -> %d \n", id, prev_rank, token);
                    /* if the color of the process is black, it will pass a black token */
                    if (color == BLACK) token = BLACK;
                    //fprintf(stderr, "[TASK %d] Sending token to %d : token color -> %d \n", id, next_rank, token);
                    MPI_Send(&token, 1, MPI_INT, next_rank, TOKEN_TAG, MPI_COMM_WORLD);
                    /* a black process becomes white when it passes the token */
                    color = WHITE;
                }
            }

            if (terminate == 0 || id) {
                MPI_Bcast(&terminate, 1, MPI_INT, 0, MPI_COMM_WORLD);
                break;
            }
        }
    }

    free(idle_processes);
    free(load_LB);
    free(load_amount);
    free(tour_nodes);
    // ring termination (so qd fizermos o load balancing)

    /* in the end, determine the best solution */
    double best_tour_cost_neighbor;
    int best_tour_neighbor[n_cities + 1];

    if (id) {
        MPI_Send(best_tour_cost, 1, MPI_DOUBLE, 0, COST_TAG, MPI_COMM_WORLD);
        MPI_Send((*best_tour), n_cities + 1, MPI_INT, 0, TOUR_TAG, MPI_COMM_WORLD);
    }
    else {
        MPI_Status status;
        for (int i = 1; i < n_tasks; i++) {
            MPI_Recv(&best_tour_cost_neighbor, 1, MPI_DOUBLE, i, COST_TAG, MPI_COMM_WORLD, &status);
            
            MPI_Recv(best_tour_neighbor, n_cities + 1, MPI_INT, i, TOUR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (best_tour_cost_neighbor < (*best_tour_cost)) {
                (*best_tour_cost) = best_tour_cost_neighbor;
                memcpy((*best_tour), best_tour_neighbor, (n_cities + 1) * sizeof(int));
            }
        }
    }
        
    // gather e dar a melhor tour (na variavel global)
    // mpi finish
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
    int n_tasks, id;
    
    City * cities = NULL;   /* cities in the map */

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    
    exec_time = - MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    /* * * * * * * * * * * * * * * * * * */

    if (argc < 3) {
        cerr << "Usage: tsp <input_file> <max-value>\n";
        return 1;
    }

    readInputFile(argv[1], &n_cities, &n_roads, &matrix, &cities);

    
    max_value = stod(argv[2]);


    tsp(&best_tour_cost, max_value, n_cities, &best_tour, matrix, cities, n_tasks, id);
    
    
    exec_time += MPI_Wtime();
    
    if(!id){
        fprintf(stderr, "%.1fs\n", exec_time);
        print_result(best_tour_cost, max_value, n_cities, best_tour);
    }

    free(matrix);
    for (int i = 0; i < n_cities; i++) {
        free(cities[i]);
    }
    free(cities);

    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

//if(!finished){
    //    char ** buffer = (char **) calloc(n_tasks*n_tasks, sizeof(char *));
    //    for (int i = 0; i < n_tasks*n_tasks; i++){
    //        Node node = initial_queue.pop();
    //        int length = node->length;
    //        buffer[i] = (char *) calloc(1, sizeof(double)*2 + sizeof(int) + sizeof(int)*length);
    //        memcpy(buffer[i], &node->cost, sizeof(double));
    //        memcpy(buffer[i] + sizeof(double), &node->lower_bound, sizeof(double));
    //        memcpy(buffer[i] + sizeof(double)*2, &length, sizeof(int));
    //        memcpy(buffer[i] + sizeof(double)*2 + sizeof(int), node->tour, sizeof(int)*length);
    //        MPI_Send(buffer[i], sizeof(double)*2 + sizeof(int) + sizeof(int)*length, MPI_CHAR, i % n_tasks, 0, MPI_COMM_WORLD);
    //    }
    //}