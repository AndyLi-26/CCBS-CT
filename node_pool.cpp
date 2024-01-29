#include "node_pool.h"
#include "structs.h"



void node_pool::init(size_t num_nodes)
{
	num_blocks_ = ((num_nodes) >> node_pool_ns::LOG2_NBS)+1;
	blocks_ = new Node*[num_blocks_];
	for(size_t i=0; i < num_blocks_; i++)
	{
		blocks_[i] = 0;
	}

    // by default:
    // allocate one chunk of memory of size
    // warthog::mem::DEFAULT_CHUNK_SIZE and assign addresses
    // from that pool in order to generate blocks of nodes. when the pool is
    // full, cpool pre-allocates more, one chunk at a time.
    size_t block_sz =
        node_pool_ns::NBS * sizeof(Node);
        node_pool_ns::NBS * sizeof(Node);
    blockspool_ = new cpool(block_sz, 1);
}

node_pool::~node_pool()
{
	blockspool_->reclaim();
	delete blockspool_;

	for(size_t i=0; i < num_blocks_; i++)
	{
		if(blocks_[i] != 0)
		{
			//std::cerr << "deleting block: "<<i<<std::endl;
			blocks_[i] = 0;
		}
	}
    delete [] blocks_;
}
Node* node_pool::generate(int node_id)
{
	int block_id = node_id >> node_pool_ns::LOG2_NBS;
	int list_id = node_id &  node_pool_ns::NBS_MASK;

    // id outside the pool address range
	if(block_id > num_blocks_) { return 0; }

    // add a new block of nodes if necessary
	if(!blocks_[block_id])
	{
		//std::cerr << "generating block: "<<block_id<<std::endl;
        blocks_[block_id] = new (blockspool_->allocate()) Node[node_pool_ns::NBS];

        // initialise memory
        int current_id = node_id - list_id;
		for( uint32_t i  = 0; i < node_pool_ns::NBS; i+=8)
		{
            new (&blocks_[block_id][i])   Node(current_id++);
            new (&blocks_[block_id][i+1]) Node(current_id++);
            new (&blocks_[block_id][i+2]) Node(current_id++);
            new (&blocks_[block_id][i+3]) Node(current_id++);
            new (&blocks_[block_id][i+4]) Node(current_id++);
            new (&blocks_[block_id][i+5]) Node(current_id++);
            new (&blocks_[block_id][i+6]) Node(current_id++);
            new (&blocks_[block_id][i+7]) Node(current_id++);
		}
	}

	// return the node from its position in the assocated block
    return &(blocks_[block_id][list_id]);
}

Node* node_pool::get_ptr(int node_id)
{
	int block_id = node_id >> node_pool_ns::LOG2_NBS;
	int list_id = node_id &  node_pool_ns::NBS_MASK;

    // id outside the pool address range
	if(block_id > num_blocks_) { return 0; }

	if(!blocks_[block_id])
    {
        return 0;
    }
    return &(blocks_[block_id][list_id]);
}

size_t node_pool::mem()
{
	size_t bytes =
        sizeof(*this) +
        blockspool_->mem() +
		num_blocks_*sizeof(void*);

	return bytes;
}
