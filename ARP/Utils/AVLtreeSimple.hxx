/*
	This file contains AVL tree definitions.

	Kalev Kask, March 2002.
*/

#ifndef AVLSimple_HXX_INCLUDED
#define AVLSimple_HXX_INCLUDED

#include <inttypes.h>

#define AVL_DEFAULT_SIZE	1024
#define AVL_SIZE_LIMIT		1048576

/*
	This structure is a simple node in the AVL tree. Simple node has only id, no data.
*/

typedef struct _MauiAVLNodeSimple {

	// the id of this entity (or the key)
	int32_t m_id ;

	// a linked list of elements in use, or not in use
	// for used nodes, 
	//		this is a doubly linked list, sorted (from left to right) in increasing order of keys
	// for unsed nodes, 
	//		this if a singly linked list (only m_next_entity is used) of unused nodes.
	uint32_t m_prev_entity, m_next_entity ;

	// left child and right child
	uint32_t m_LC, m_RC ;
	// height is the length of the int32_test path from this node to a leaf.
	// we want to use AVLtree[0].m_height = -1, so that we can do 'AVLtree[parent] = ++AVLtree[child].m_height' and it is correct when child=0.
	signed short m_height ;
	// balance = height of right subtree - height of left subtree
	signed char m_balance ;

	// 0 if not in use, 1 if in use
	unsigned char m_in_use ;

public :

	void operator<<(const _MauiAVLNodeSimple & Node) ;

} MauiAVLNodeSimple ;

/*
	This class is used to represent an AVL tree.
*/

class CMauiAVLTreeSimple
{

protected :

	// pointer to the actual AVL tree
	MauiAVLNodeSimple *m_db_tree ;

	// this is the amount of space allocated for the AVL tree
	uint32_t m_allocated_space ;

	// When reallocating the tree, system will reserve that much more memory.
	// If this variable is 0, system will never reallocate.
	// That means that when the tree is full and the user calls Insert() on the tree, it will fail.
	// Initially equal to AVL_DEFAULT_SIZE.
	uint32_t m_re_allocation_size ;

	// this is the actual (current) size of the AVL tree (ie. number of currently used blocks)
	uint32_t m_size ;

	// number of free blocks
	uint32_t m_free_size ;

	// points to the first free block
	uint32_t m_first_free_block ;

	// points to the node that has the smallest key.
	// this key is the first node in the ordered list of keys.
	// this list if build using m_next_entity and m_prev_entity pointers stored in each node.
	uint32_t m_first_used_block ;

	// a pointer to the root of the AVL tree.
	// if root pointer is less than 1, (eg. 0 or -1) the root node is undefined.
	// this happens when the tree is empty, or the memory for the tree is not allocated.
	uint32_t m_root ;

	// these variables are used to store the stack (usually) when traversing the tree
	// since the height of a AVL tree is always no more than 1.44 log(n) we don't need very much space.
	// assuming we don't allow more than 1G objects, the depth is bounded by 1.44*30.
	int32_t Left[64], Right[64], Middle[64] ;

	// this variable is used during balancing of the AVL tree.
	uint32_t path_len ;

private :

	bool operator<<(const CMauiAVLTreeSimple & AVLtree) ;

private :

	// Management functions are private.

	void R_rotation(uint32_t k) ;
	void L_rotation(uint32_t k) ;
	void BalanceInsertTree(void) ;
	void BalanceDeleteTree(void) ;
	// initial memory allocation
	int AllocateMemory(uint32_t size) ;
	// when the tree is full, then system tries to reallocate the entire tree
	int ReAllocateMemory(void) ;
	void PrintTree(void) ;

// *****************************************************************************************************
// Timer.
// *****************************************************************************************************

public :

	inline uint32_t GetSize(void) const { return m_size ; }
	inline uint32_t GetSpace(void) const { return m_allocated_space ; }
	inline uint32_t GetFirstUsedBlockIndex(void) const { return m_first_used_block ; }

	void Empty(void) ;
	void EmptyQuick(void) ;

	// returns TRUE iff the key in in the tree
	int32_t Find(int32_t key) ;
	int32_t FindNext(int32_t key, int32_t *next_key) ;

	// returns FALSE iff it fails (no memory). if key is in use, it will replace.
	int32_t Insert(int32_t key) ;

	// Remove a key from the tree.
	// Note that it is not a good idea to delete keys while the tree is being traversed using the
	// GetNext() function, because the selection/current pointer might be invalid.
	int32_t Remove(int32_t key) ; // returns TRUE iff everything is fine
	int32_t RemoveFirst(int32_t *key /* output */) ; // returns TRUE iff everything is fine

	void SetReallocationSize(uint32_t reallocation_size) ; // has to be non-negative. default AVL_DEFAULT_SIZE.

	// to find the smallest positive key not in the tree.
	// return value is > 0.
	int32_t FindSmallestPositiveKeyNotUsed(void) ;

/*
	Member functions to traverse the tree.
	Returns FALSE iff the tree is empty or no next element.
	Note that the keys are returned in an increasing order of keys.
*/
	// Note that this function does set the current pointer to the first element in the tree
	int32_t GetFirst(int32_t *key /* can be NULL */, int32_t & current) const ;
	// if current = -1 is passed in, this function returns the first key.
	// otherwise this function returns the key next to the current in the increasing order of keys.
	int32_t GetNext(int32_t *key /* can be NULL */, int32_t & current) const ;
	// this function returns the key/data pointed to by 'current' 
	// (if -1==current, then this means the first, if 0==current, then this means NONE).
	// it also forwars the 'current' pointer to the next element in the set.
	int32_t GetCurrentAndAdvance(int32_t *key /* can be NULL */, int32_t & current) const ;

	void PrintString(char *) { }
	int32_t CheckTree(char **pErrorStr = NULL) ; // returns TRUE iff the tree is OK.
	int32_t TestTree(int numIterations, char **pErrorStr = NULL) ; // returns TRUE iff the tree is OK.
	bool TestConsistency(int32_t checkbitvector = 255) const ;

// *****************************************************************************************************
// Construction and destruction.
// *****************************************************************************************************

public :

	// Default constructor creates an empty tree of size AVL_DEFAULT_SIZE.
	CMauiAVLTreeSimple(void) ;

	// This constructor creates a tree from a set of keys.
	CMauiAVLTreeSimple(uint32_t size /* number of keys (has to be >= 0) */, 
			int32_t *key /* an array of keys */, 
			uint32_t space /* how much space to allocate for the AVL tree (> 0) */, 
			uint32_t re_allocation_size = AVL_DEFAULT_SIZE) ;

	// Destructor has to release the memory that the AVL tree uses.
	~CMauiAVLTreeSimple(void) ;

} ;

#endif // AVLSimple_HXX_INCLUDED
