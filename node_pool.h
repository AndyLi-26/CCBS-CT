#pragma once

#include "cpool.h"
#include "structs.h"
#include <stdint.h>
namespace node_pool_ns
{
	static const uint64_t NBS = 8; // node block size; set this >= 8
	static const uint64_t LOG2_NBS = 3;
	static const uint64_t NBS_MASK = 7;
}

class node_pool
{
	public:
		node_pool();
        node_pool(size_t num_nodes):blocks_(0)
		{init(num_nodes);};
		~node_pool();

		// return a warthog::search_node object corresponding to the given id.
		// if the node has already been generated, return a pointer to the
		// previous instance; otherwise allocate memory for a new object.
		Node* generate(int node_id);

        // return a pre-allocated pointer. if the corresponding node has not
        // been allocated yet, return null
        Node* get_ptr(int node_id);
		size_t mem();

	private:
        void init(size_t nblocks);

		size_t num_blocks_;
		Node** blocks_;
		cpool* blockspool_;
};
