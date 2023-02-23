#ifndef _TSP_QUEUE_H
#define _TSP_QUEUE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// A queue where the elements are stored in an increasing order.
// This implementation uses a binary heap.
typedef struct
{
	void** buffer;
	size_t size;
	size_t max_size;
        char (*cmpfn)(void *, void *);
} priority_queue_t;

// Create a new priority queue
priority_queue_t *queue_create(char (*)(void *, void *));

// Delete an existing priority
void queue_delete(priority_queue_t *queue);

// Insert a new element in the queue and then sort its contents.
void queue_push(priority_queue_t *queue, void* new_element);

// Return the element with the lowest value in the queue, after removing it.
void* queue_pop(priority_queue_t *queue);

// Print the contents of the priority queue
void queue_print(priority_queue_t *queue, FILE *, void (*)(FILE *, void*));

#endif //_TSP_QUEUE_H
