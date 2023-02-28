/*************************************************************************
	Copyright (C) 2023 Nicolas L. Guidotti. All rights reserved.

	This file is part of the tsp library, which is licensed under 
	the terms contained in the LICENSE file.
 **************************************************************************/

#ifndef _TSP_QUEUE_HPP
#define _TSP_QUEUE_HPP

// A queue where the elements are stored in an increasing order.
// This implementation uses a binary heap.
//
// Parameters:
// - "T": Type of each element in the queue
// - "Compare": Type of the callable object with the compare function
//
// Usage Example:
// PriorityQueue<int> queue;       // Creates a new queue of ints
// queue.push(10);                 // Add the element "10" to the queue.
// queue.push(5);                  // Add the element "5" to the queue. The program will automatically sort the queue.
// int val = queue.pop();          // Remove the lowest element in the queue. In this case, "val = 5"
//
// Note: this structure assumes that the operator ">" is defined for "T", i.e.,
// bool operator>(const T& left, const T& right) { return left > right; }			// Change here how the elements compare to each other
//
// Alternatively, you can also define a struct that hold the compare function:
// struct cmp_op { bool operator()(const T& left, const T& right) { return left > right; } }
// and then pass it to the template parameters of the priority queue:
// PriorityQueue<int, cmp_op> queue;
//
template<typename T, typename Compare = std::greater<T>>
class PriorityQueue
{
    private:
		using index_t = int;

		// Array for storing the elements in the binary heap
		std::vector<T> _buffer;
		
		// Callable object with the compare function
		Compare _cmp;
		
		// Return the index of the parent node
		static index_t parent_of(index_t i)
		{
			return (i - 1) / 2;
		}
		
		// Bubble-down the element to the correct position
		// (i.e., compare it to its child and then swap them if necessary).
		// Assume that all the elements in the subtree is already sorted.
		void bubble_down(index_t node)
		{
			index_t left_child = 2 * node + 1;
			index_t right_child = 2 * node + 2;
			index_t i = node;
			
			// Compare with the left node
			if (left_child < _buffer.size() and _cmp(_buffer[node], _buffer[left_child]))
			{
				i = left_child;
			}
			
			// Compare with the right node
			if (right_child < _buffer.size() and _cmp(_buffer[i], _buffer[right_child]))
			{
				i = right_child;
			}
			
			// If node is not in the correct position, swap and then sort the subtree
			if (i != node)
			{
				std::swap(_buffer[i], _buffer[node]);
				bubble_down(i);
			}
		}
		
    public:
        PriorityQueue() = default;
        virtual ~PriorityQueue() = default;
		PriorityQueue& operator=(const PriorityQueue& other) = default;
		
		// Check if the priority is empty
		bool empty()
		{
			return _buffer.empty();
		}
		
		// Return the number of elements in the queue
		index_t size()
		{
			return _buffer.size();
		}
		
		// Return a raw pointer to contents of the queue
		T* data()
		{
			return _buffer.data();
		}
		
		// Insert a new element using the copy semantics. This routine automatically sorts the queue
		void push(const T& new_element)
		{
			// Insert the new_element at the end of the buffer
			index_t node = _buffer.size();
			_buffer.push_back(new_element);
			
			// Bubble-up the new element to the correct position
			// (i.e., compare it to the parent and then swap them if necessary)
			while (node > 0 and _cmp(_buffer[parent_of(node)], _buffer[node]))
			{
				index_t parent = parent_of(node);
				std::swap(_buffer[node], _buffer[parent]);
				node = parent;
			}
		}
		
		// Insert a new element using the move semantics. This routine automatically sorts the queue.
		void push(T&& new_element)
		{
			// Insert the new_element at the end of the buffer
			index_t node = _buffer.size();
			_buffer.push_back(std::forward<T>(new_element));
			
			// Bubble-up the new element to the correct position
			// (i.e., compare it to the parent and then swap them if necessary)
			while (node > 0 and _cmp(_buffer[parent_of(node)], _buffer[node]))
			{
				index_t parent = parent_of(node);
				std::swap(_buffer[node], _buffer[parent]);
				node = parent;
			}
		}
		
		// Return the element with the lowest value, after removing it from the queue
		T pop()
		{
			// If the queue is empty, calls the default constructor
			if (empty()) return T();
			
			// Stores the lowest element in a temporary
			T top_val = std::move(_buffer.front());
			
			// Put the last element in the queue in the front.
			_buffer[0] = std::move(_buffer.back());
			
			// Remove the duplicated element in the back.
			_buffer.pop_back();
			
			// Sort the queue based on the value of the nodes.
			bubble_down(0);
			
			return top_val;
		}
		
		// Print the contents of the queue. Note that you can use any callable object as parameter,
		// including lambdas, std::function, structs with the operator(), etc.
		//
		// Usage:
		// PriorityQueue<std::pair<int, float>> queue;
		// ...
		// queue.print([](const std::pair<int, float>& elem){ fprintf(stderr, "%d %f\n", elem.first, elem.second); });
		//
		template<typename Print>
		void print(Print print_per_element)
		{
			PriorityQueue queue_copy = *this;
			
			while (not queue_copy.empty())
			{
				T elem = queue_copy.pop();
				print_per_element(elem);
			}
		}
};

#endif //_TSP_QUEUE_HPP
