/* MIST copyright 2016 by Yujin Chung and Jody Hey */

#include <iostream>
#include <algorithm>
#include "coaltree.hpp"
#include "misc.hpp"


/// YC 6/11/2014
/// Mutiply branch lengths by treeHeightScaler (ranges from 2/3 to 3/2)
/// The tree will be extended or shrinked
void node::rescale_treeHeight(double treeHeightScaler)
{
  age *= treeHeightScaler;
  if(isTip == 0)
    {
      desc[0]->rescale_treeHeight(treeHeightScaler);
      desc[1]->rescale_treeHeight(treeHeightScaler);
    }
  return;
}


node* node::search_node2detach(int nodeID)
{
	// NOTE: A node on a tree has its unique ID,
	// although node ID is not a member of class 'node'.
	// The roles of node ID assignment are
	//		1) A node has a smaller ID than its descendants. (e.g., the ID of the root is zero)
	//		2) If node A is lighter (has fewer descendants) than its sibling node B,
	//			then all the descendants (nodes) of A have smaller ID than those of B.
	//      3) If sibling nodes A and B have the same number of descendants,
	//         but A is older than B, then all the descendants (nodes) of A
	//		   have smaller ID than those of B.
	// In other words,	we can list gene trees in a unique way. For every pair of sibling subtrees
	// on a tree, we put lighter subtree (subtree with fewer descendants) on the left side.
	// If the sizes are the same, then we put older subtree on the left side.
	// Then we label nodes from top to bottom and from left to right.

	// node to detach
	// node* node_detach;

	// REMOVE
	//std::cout << "in node::search_node2detach()\n";
	//std::cout <<"Current node is ";
	//print_coaltree();
	//std::cout << "nodeID = " << nodeID <<"\n";

	if(nodeID == 0)
		// the current node is the node to detach
	{
		// REMOVE
		// print_nodeMembers();

		if(isTip == 0)
		{
			return desc[0]->par; //get_node();
			//node_detach = desc[0]->par; //get_node();
		}
		else
		{
			// REMOVE
			//std::cout << "The sibling order of the following subtree is " << siblingOrder <<"\n";
			//print_coaltree();

			return par->desc[siblingOrder];
			//node_detach = par->desc[siblingOrder];
		}

	}
	else
		// the node to detach is a descendant of the current node.
	{
		node* left_tree;
		node* right_tree;
		int size1 = desc[0]->size_tree();
		int size2 = desc[1]->size_tree();

		// Find which child node is the left tree.
		if(size1 == size2)
		{
			if(desc[0]->age > desc[1]->age)
			//	then desc[0] is the left subtree
			{
				left_tree = desc[0];
				right_tree = desc[1];
			}
			else if(desc[0]->isTip && desc[1]->isTip && desc[0]->tipID < desc[1]->tipID)
			{
				left_tree = desc[0];
				right_tree = desc[1];
			}
			else
			{
				left_tree = desc[1];
				right_tree = desc[0];
			}
		}
		else
		{
			if(size1 < size2)
			{
				left_tree = desc[0];
				right_tree = desc[1];
			}
			else
			{
				left_tree = desc[1];
				right_tree = desc[0];
			}
		}
		//REMOVE
		// std::cout <<"left_tree->size_tree() = " <<left_tree->size_tree()<<"\n";
		// std::cout <<"right_tree->size_tree() = " <<right_tree->size_tree()<<"\n";

		// Find the node with 'nodeID'
		if(2*(left_tree->size_tree())-1 >= nodeID)
		// the node to detach in 'left_tree' or its descendants
			return left_tree->search_node2detach(nodeID-1);
		else
		{// the node to detach in 'right_tree' or its descendants
			int newNodeID = nodeID-2*(left_tree->size_tree());

			//REMOVE testing
			// cout << "newNodeID = "<< newNodeID <<"\n";

			return right_tree->search_node2detach(newNodeID);
					//node_detach = right_tree->search_node2detach(newNodeID);
			// that is, the new node ID is the old node ID - the number of nodes in the left tree (2*nTips -1) - 1 (for the current node)
		}
	}

	// REMOVE
	//		std::cout << "Node picked\n";
	//		node_detach->print_coaltree();
	//		std::cout << "Parent is ";
	//		node_detach->par->print_coaltree();

	//return node_detach;
}



node* node::search_node2detach_fromRoot(int nodeID)
{
	node* node_detach;
	node* left_tree;
	node* right_tree;
	int size1 = desc[0]->size_tree();
	int size2 = desc[1]->size_tree();


	// NOTE: A node on a tree has its unique ID,
	// although node ID is not a member of class 'node'.
	// The roles of node ID assignment are
	//		1) A node has a smaller ID than its descendants. (e.g., the ID of the root is zero)
	//		2) If node A is lighter (has fewer descendants) than its sibling node B,
	//			then all the descendants (nodes) of A have smaller ID than those of B.
	//      3) If sibling nodes A and B have the same number of descendants,
	//         but A is older than B, then all the descendants (nodes) of A
	//		   have smaller ID than those of B.
	// In other words,	we can list gene trees in a unique way. For every pair of sibling subtrees
	// on a tree, we put lighter subtree (subtree with fewer descendants) on the left side.
	// If the sizes are the same, then we put older subtree on the left side.
	// Then we label nodes from top to bottom and from left to right.

	// Find which child node is the left tree.
	if(size1 == size2)
	{
		if(desc[0]->age > desc[1]->age)
			//	then desc[0] is the left subtree
		{
			left_tree = desc[0];
			right_tree = desc[1];
		}
		else
		{
			left_tree = desc[1];
			right_tree = desc[0];
		}
	}
	else
	{
		if(size1 < size2)
		{
			left_tree = desc[0];
			right_tree = desc[1];
		}
		else
		{
			left_tree = desc[1];
			right_tree = desc[0];
		}
	}

	// Find the node with 'nodeID'
	if(2*left_tree->size_tree()-1 >= nodeID) // the node to detach in 'left_tree' or its descendants
		node_detach = left_tree->search_node2detach(nodeID-1);
	else // the node to detach in 'right_tree' or its descendants
		node_detach = right_tree->search_node2detach(nodeID-2*(left_tree->size_tree()));


	return node_detach;
}

/***
 *  Randomly pick a node to detach
 */
node* node::randomPick_node()
{
	// REMOVE
	//std::cout << "in node::randomPick_node()\n";
	//std::cout << "give tree is ";
	//print_coaltree();

	// The number of candidate nodes on a tree.
	// The root node is excluded.
	int n_nodes2pick = 2*size_tree()-2;

	// Randomly pick a node ID to detach.
	// Candidate node IDs are 1, 2, ..., n_nodes2pick.
	// Note. node IDs are 0 to 2*n-2, where n is the number of tips.
	int nodeID = runiform_discrete(n_nodes2pick)+1;
	// REMOVE
	// std::cout << "nodeID = " <<nodeID <<"\n";

	// Find the node whose ID is 'nodeID'.
	node* node_detach = search_node2detach_fromRoot(nodeID);

	// REMOVE
	//std::cout << "Node picked\n";
	//node_detach->print_coaltree();
	//std::cout << "Parent is ";
	//node_detach->par->print_coaltree();

	return node_detach;
}


/**
 *
 */
void node::slider(double dist)
{
	// REMOVE
	//std::cout << "in node::slider()\n";
	//print_coaltree();
	//std::cout << "siblingOrder is "<<siblingOrder <<"\n\n";

	//int sisID = find_sisterID();
	//int nodeID = 1-sisID;

  double brlen = par->age - age;
  double brlen_sister = par->age - par->desc[1-siblingOrder]->age;
  
  // REMOVE
  //	std::cout << "in node::slider()\n";
	//std::cout << "dist = " << dist << "\n";
	//	std::cout << "brlen = " << brlen <<" and brlen_sister = " << brlen_sister << "\n";
	//	std::cout << "sister node tree: ";
	//	par->desc[sisID]->print_coaltree();
  
  if(dist < 0.0 )
    {
      // std::cout << "Distance is negative.\n";
      
      // The branch from the current node to its parent has a tendency to be shorten.
      // That is, the parent node moves toward the current node and the sister.
      
      dist *= -1.0;
      if(dist < std::min(brlen, brlen_sister))
	// if dist < both of branch lengths,
	// then branch lengths are shorten. No tree topology change.
	// that is, lengthen the branch from the parent node to the grand parent node
	{
	  par-> age -= dist;
	  par->makeNULL_ancestorsLik();
	  
	  if(go2root()->validTree() == 0)
			{
			  std::cout << "in node::slider(). Not valid tree is proposed (Error1)\n";
			}
	}
      else
	{
	  if(brlen <= brlen_sister)
	    {
	      // Since parent node can't move behind the current node,
	      // move down the parent node to the current node its sister node and move it back (reflect)
	      
	      par->age = age; // Move down the parent node to its sister node.
	      // Then	the branch length from the current node to its parent = 0
	      par->makeNULL_ancestorsLik();
	      
	      if(go2root()->validTree() == 0)
		{
		  std::cout << "in node::slider(). Not valid tree is proposed (Error2)\n";
		}
	      
	      dist -= brlen; // The remaining distance after the move.
	      slider(dist); // So far tree topology doesn't change, but it can be changed in the recursion
	    }
	  else
			  // if brlen_sister < brlen (and brlen_sister < dist),
			  // then the current node can attach to a descendant edge of its sister edge.
			  {
			    // pick a child node of the sister node
				int childID = 0;
				if( runiform() < 0.5)
					childID = 1;

				//REMOVE
				//std::cout << "picked childID =" <<childID <<"\n";
				//std::cout << "and the node tree is ";
				//par->desc[1-siblingOrder]->desc[childID]->print_coaltree();

				// Disregard the node (parent node) connecting the grand parent node and the sister node
				if(par->isRoot == 0)
				{
					unsigned int newSiblingOrder_par =par->get_siblingOrder();
					// REMOVE
					//std::cout << "newSiblingOrder4sister is" << newSiblingOrder_par <<"\n";
					try{
						(par->par)->desc[newSiblingOrder_par] = par->desc[1-siblingOrder];
					}catch (std::exception &e) {
						std::cout << "In node::slider()\n Can't access (par->par)->desc[newSiblingOrder4sister]"
							" or par->desc[1-siblingOrder] - array index out of bounds\n"
							"newSiblingOrder4sister = " << newSiblingOrder_par
							<<" and siblingOrder = "<< siblingOrder <<"\n";
					}
					(par->desc[1-siblingOrder])->par = par->par;
					(par->desc[1-siblingOrder])->assign_siblingOrder(newSiblingOrder_par);

					// Move the parent node to the edge connecting the sister node and its child node.
					(par->desc[1-siblingOrder]->desc[childID])->par = par;
					par->desc[1-siblingOrder] = par->desc[1-siblingOrder]->desc[childID];
					par->desc[1-siblingOrder]->assign_siblingOrder(1-siblingOrder);
					(par->par->desc[newSiblingOrder_par])->desc[childID] = par;
					par->par =par->par->desc[newSiblingOrder_par];
					par->assign_siblingOrder(childID);
					par->age = par->par->age;
				}
				else // The sister node is the new root.
				{
					// REMOVE
					//std::cout << "sister node is \n";
					//par->desc[1-siblingOrder]->print_coaltree();

					par->par = par->desc[1-siblingOrder]; // Sister node is the new parent
					par->isRoot = 0;

					// Move the parent node to the edge connecting the sister node and its child node.
					(par->desc[1-siblingOrder]->desc[childID])->par = par;
					par->desc[1-siblingOrder] = par->desc[1-siblingOrder]->desc[childID];
					par->desc[1-siblingOrder]->assign_siblingOrder(1-siblingOrder);
					// new root
					(par->par)->desc[childID] = par;
					par->assign_siblingOrder(childID);
					(par->par)->isRoot = 1;
					(par->par)->par = 0;
					(par->par)->assign_siblingOrder(2);
					par->age = par->par->age;
				}


				// Add a node on the edge between sister and sister's child node with 'childID'
				// The new node has the sister node as parent
				// and has the current node and sister's child as children.
				/*
				node* new_node = new node;
				new_node->initialization();
				//node tmp;
				//node* new_node = &tmp;
				new_node -> isTip = 0;
				new_node->isRoot = 0;
				new_node->isLikelihoodNULL = 1;
				new_node -> desc[0] = par->desc[1-siblingOrder]->desc[childID];
				new_node -> par = par->desc[1-siblingOrder];
				new_node -> desc[1] = par->desc[siblingOrder];
				new_node -> age = par->desc[1-siblingOrder]->age;


				// sister's node has a new parent
				par->desc[1-siblingOrder]->desc[childID]->par = new_node;
				// sister has a new child
				par->desc[1-siblingOrder]->desc[childID] = new_node;

				// std::cout <<par->isRoot << " here\n";
				// The root changes
				if(par->isRoot)
				{
					new_node->par->isRoot =1;
					//par->desc[sisID]->par = NULL;
					par->desc[1-siblingOrder]->par->isRoot = 1;
				}
				else
				{
					par->desc[1-siblingOrder]->par = par->par; // Connect the sister node to its grand parent (sister has new parent)
					// Remove the original parent node
					int parentID = 1- par->find_sisterID();
					par->par->desc[parentID] = par->desc[1-siblingOrder];
				}
				new_node -> desc[1]->par = new_node; // update the current node again

				par = node_node; // current node has a new parent.
				// delete new_node;
				 */
				par->makeNULL_ancestorsLik();
				if(go2root()->validTree() == 0)
				{
					std::cout << "in node::slider(). Not valid tree is proposed (Error3)\n";
				}

				// REMOVE
				// std::cout << "Newly add node tree is ";
				// new_node->print_coaltree();
				// std::cout << "par->age = " << par->age <<"\n";
				// std::cout << "parent node tree is ";
				// par->print_coaltree();
				// if(!(par->isRoot))
				// {
				// 	std::cout << "grand parent node tree is ";
				//	par->par->print_coaltree();
				// }
				// std::cout << "updated current node: ";
				// print_coaltree();
				// std::cout << "The whole tree is \n";
				// go2root()->print_coaltree();

				// The remaining distance after the current node is connected to the child edge of its original sister
				dist -= brlen_sister;
				dist *= -1;
				slider(dist);
			}
		}
	}
	else
	{
		// std::cout << "Distance is not negative.\n";

		// if dist >= 0, the move toward the root.

		if((par->isRoot) == 0)
		{
			brlen = par->par->age - par->age;
		}

		if(dist < brlen || par->isRoot==1)
		{
			// If the grand parent node is root or the amount of move ('dist') is smaller than the parent edge,
			// then the parent branch is lengthen.

			// REMOVE
			//std::cout << "parent node is ";
			//par->print_coaltree();

			par->age += dist;
			par->makeNULL_ancestorsLik();

			//std::cout << "parent node after update is ";
			//par->print_coaltree();

			if(go2root()->validTree() == 0)
			{
				std::cout << "in node::slider(). Not valid tree is proposed (Error4)\n";
			}

		}
		else
		{
			// YC
			// If the grand parent node (not root) is reached,
			// then randomly pick one between the other child edge or the parent edge of the grand parent node,
			// and keep moving along the edge.

			// Remaining distance to move
			dist -= brlen;

			//REMOVE testing
			//std::cout << "Remaining dist: "<<dist <<"\n";

			// Randomly pick one between the sister edge
			// or the parent edge of the grand parent node.
			int nodeID2go = 0;
			if(runiform() <0.5)
				nodeID2go = 1; //the parent edge of the grand parent node is picked

			//REMOVE testing
			//std::cout << "nodeID2go = " << nodeID2go <<"\n";

			// int parentID = 1 - par->find_sisterID();

			par -> age = (par->par)->age;

			if(nodeID2go) // the parent edge of the grand parent node is picked
			{
				//REMOVE testing
				// std::cout << "Move along the parent edge of the grand parent\n";



				if(par->par->isRoot == 1) // new root
				{
					par->isRoot = 1;
					(par->par)->isRoot = 0;

					// Remove the parent node between the grandparent node and the sister node
					unsigned int originalSiblingOrder_par =par->get_siblingOrder();
					(par->par)->par = par;
					(par->par)->assign_siblingOrder(1-siblingOrder);
					(par->par)->desc[originalSiblingOrder_par] = par->desc[1-siblingOrder];
					(par->desc[1-siblingOrder])->par = par->par;
					(par->desc[1-siblingOrder])->assign_siblingOrder(originalSiblingOrder_par);
					// The original parent node is the new root
					par->desc[1-siblingOrder] = par->par;
					//delete par->par;
					par->par = 0;
					par->assign_siblingOrder(2);

					if(go2root()->validTree() == 0)
					{
						std::cout << "in node::slider(). Not valid tree is proposed (Error5)\n";
					}
				}
				else
				{
					// REMOVE
					//std::cout << "grand-grand parent ";
					//(par->par->par)->print_coaltree();
					//std::cout << "grand-grand parent's child is ";
					//(par->par->par)->desc[(par->par)->get_siblingOrder()]->print_coaltree();

					// grand-grand-parent has the parent node as a child
					(par->par->par)->desc[(par->par)->get_siblingOrder()] = par;

					// REMOVE
					//std::cout << "has a new child\n";
					//(par->par->par)->print_coaltree();


					// Disregard the parent node between the grandparent node and the sister node
					unsigned int originalSiblingOrder_par =par->get_siblingOrder();
					unsigned int originalSiblingOrder_grandPar =(par->par)->get_siblingOrder();
					(par->desc[1-siblingOrder])->par = par->par;
					(par->desc[1-siblingOrder])->assign_siblingOrder(originalSiblingOrder_par);
					(par->par)->desc[originalSiblingOrder_par] = par->desc[1-siblingOrder];

					// REMOVE
					//std::cout << "Updated sister is ";
					//(par->desc[1-siblingOrder])->print_coaltree();
					//std::cout << "Updated grand parent is ";
					//(par->par)->print_coaltree();


					par->par = (par->par->par);
					// REMOVE
					//std::cout << "Updated parent is";
					//(par)->print_coaltree();
					//std::cout << "Updated parent's par is";
					//(par->par)->print_coaltree();
					//
					//(par->par)->par = par;
					//(par->par)->assign_siblingOrder(1-siblingOrder);
					(par->desc[1-siblingOrder]->par)->par = par; // the original grandparent's new parent is the original parent
					(par->desc[1-siblingOrder]->par)->assign_siblingOrder(1-siblingOrder);
					par->desc[1-siblingOrder] = (par->desc[1-siblingOrder])->par;
					par->assign_siblingOrder(originalSiblingOrder_grandPar);

					// REMOVE
					//std::cout << "Updated sister is ";
					//(par->desc[1-siblingOrder])->print_coaltree();
					//std::cout << "The final updated grand parent is ";
					//(par->par)->print_coaltree();

					if(go2root()->validTree() == 0)
					{
						std::cout << "in node::slider(). Not valid tree is proposed (Error6)\n";
					}

				}



				/*
				node* grandPar = par->par;


				// Create a new node on the parent edge of the grand parent node
				// The new node's parent is the parent of the grand parent node,
				// and children are the current node and the grand parent node.

				//node tmpNode;
				//node* new_par = &tmpNode;
				node* new_par = new node;
				new_par->initialization();
				new_par->isLikelihoodNULL = 1;
				if(grandPar->isRoot)
				{
					new_par->isRoot =1;
				}
				else
				{
					new_par->isRoot = 0;
					new_par -> par = grandPar->par; // parent is the grand-grand parent
				}
				new_par -> desc[0] = grandPar; // child is the grand parent
				new_par -> desc[1] = par->desc[siblingOrder]; //child is the current node
				new_par->age = grandPar->age;
				new_par->isTip = 0;

				// Connect the grand-grand parent node to the new node
				if(!(grandPar->isRoot))
				{
					int grandparentID = 1 - par->par->find_sisterID();
					par->par->par->desc[grandparentID] = new_par;
				}

				// The grandparent node's parent is the new node.
				grandPar->par = new_par;
				if(grandPar->isRoot)
					grandPar->isRoot = 0;

				// Connect the sister node to the grand parent node
				// That is, remove the parent node between them.
				par->desc[1-siblingOrder]->par = grandPar;
				grandPar->desc[parentID] = par->desc[1-siblingOrder];

				// The current node's new parent is the new_node
				// Remove the parent node foreaver.
				par = new_par;
				new_par -> desc[1] ->par = new_par;
				*/

				par->desc[1-siblingOrder]->makeNULL_ancestorsLik();

				// par->initialize_lik();

				// REMOVE
				// std::cout << "The whole tree is \n";
				// go2root()->print_coaltree();


				// Move toward the root more.
				slider(dist);

			}
			else
			{ // the the other child edge of the grand parent node is picked.

				unsigned int originalSiblingOrder_par =par->get_siblingOrder();

				// Assign new parent node to parent's sister
				(par->par->desc[1-originalSiblingOrder_par])->par = par;
				(par->par->desc[1-originalSiblingOrder_par])->assign_siblingOrder(1-siblingOrder);

				// REMOVE
				//std::cout << "new grand parent's descendants:\n";
				//(par->par)->print_coaltree();

				// REMOVE
				//	std::cout << "Test:\n";
				//	(par->par)->desc[originalSiblingOrder_par]->print_coaltree();
				//	(par->desc[1-siblingOrder])->print_coaltree();
				//	(par->desc[siblingOrder])->print_coaltree();
				//	std::cout << "siblingOrder = " << siblingOrder <<"\n";
				//	std::cout << "sister's siblingOrder = " << (par->desc[0])->get_siblingOrder()<< "\n";
				//	std::cout << "sister's siblingOrder = " << (par->desc[1])->get_siblingOrder()<< "\n";


				// Remove the parent node between the grand parent node and sister node
				(par->desc[1-siblingOrder])->par = par->par;
				(par->desc[1-siblingOrder])->assign_siblingOrder(originalSiblingOrder_par);
				(par->par)->desc[originalSiblingOrder_par] =(par->desc[1-siblingOrder]);


				// REMOVE
				//std::cout << "new grand parent's descendants:\n";
				//(par->par)->print_coaltree();

				// Assign new sister node
				par->desc[1-siblingOrder] =	(par->par->desc[1-originalSiblingOrder_par]);

				// Assign new child node (the original parent node) to the original grand parent node
				(par->par)->desc[1-originalSiblingOrder_par] = par;
				par->assign_siblingOrder(1-originalSiblingOrder_par);



				if(go2root()->validTree() == 0)
				{
					std::cout << "in node::slider(). Not valid tree is proposed (Error7)\n";
				}

				// REMOVE remove
				// std::cout <<"Move to parent sister's node\n";
				// std::cout << "parentID = "<< parentID << "\n";

				/*

				node* grandPar = par->par;
				node* parSis = grandPar->desc[1-parentID];
				node* sis =par->desc[1-siblingOrder];

				//REMOVE testing
				//std::cout << "grand parent ";
				//grandPar->print_coaltree();
				//std::cout << "sister of parent ";
				//parSis->print_coaltree();
				//std::cout << "sister ";
				//sis->print_coaltree();

				// Create a new node on the edge between the grand parent node and the sister of the parent node.
				// The new node's parent is the grand parent node, and
				// children are the sister of the parent node and the current node.


				//node tmp2;
				//node* new_par = &tmp2;
				 node* new_par = new node;
				 new_par->initialization();
				new_par->isLikelihoodNULL=1;
				new_par->isTip = 0;
				new_par->isRoot = 0;
				new_par->par = grandPar; // par->par->par;
				new_par->desc[0] = parSis; //par->par->desc[1-parentID];
				new_par->desc[1] = par->desc[siblingOrder];
				new_par->age = grandPar->age; // par->par->age;

				// the sister of the parent node has the new node as parent
				parSis->par = new_par;

				// The grand parent's children are the new node and the sister node
				grandPar->desc[1-parentID] = new_par;
				grandPar->desc[parentID] = sis;

				// Sister's parent should the grand parent node;
				sis->par = grandPar;

				// the current node's parent is the new node
				// Remove the original parent node
				par = new_par;
				new_par -> desc[1] ->par = new_par;
				*/
				par->desc[1-siblingOrder]->makeNULL_ancestorsLik();


				// par->initialize_lik();

				//REMOVE testing
				//std::cout << "Newly add node tree is ";
				//new_par->print_coaltree();
				//std::cout << "par->age = " << par->age <<"\n";
				//std::cout << "parent node tree is ";
				//par->print_coaltree();
				//std::cout << "grand parent node tree is ";
				//par->par->print_coaltree();
				//std::cout << "updated current node: ";
				//print_coaltree();
				//std::cout << "The whole tree is \n";
				//go2root()->print_coaltree();

				// Move down along the child edge of the grand parent node
				dist *= -1;
				slider(dist);
			}
		}
	}
}

node* node::go2root()
{
	//REMOVE
	//cout << "in go2root(): ";
	//print_coaltree();

	if(isRoot)
		return desc[0]->par;
	else
		return par->go2root();
}



/**
 * Propose a new coalescent tree from the current tree.
 * Sliding window approach is used.
 * The code is logically same as the tree proposal part of function 'updategenealogy' in IMa2.
 */
node* node::propose_coaltree(double slidedist)
{

	// REMOVE
	//std::cout <<"In node::propose_coaltree()\n";
	//print_coaltree();

	// pick a node to detach
  node* node_detach = randomPick_node();


  // sliding
  node_detach->slider(slidedist);
  
  // Move to the root node
  node* new_tree = node_detach->go2root();

	// REMOVE
	//std::cout << "Therefore the newly propose tree is \n";
	//new_tree->print_coaltree();
	// std::cout << "Proposed tree's right child is";
	// new_tree->desc[1]->print_coaltree();
	// std::cout << "isLikelihoodNULL: " << new_tree->desc[1]->isLikelihoodNULL <<"\n";

  // for likelihood calculation later
  new_tree->computeTotalLengthsOfLineages();

  return new_tree;
}




/**
 * YC 3/4/2014
 * Compute the total coalescent rate of a tree
 * The total coalescent rate is similar to Eq.(16) in Hey and Nielsen (2007).
 * For example, for a time interval t_i, there are n_i lineages.
 * Then the corresponding coalescent rate used in this function is
 *     t_i * n_i*(n_i-1)/2.
 * Note the rate used in Hey and Nielsen (2007) is
 *     t_i * n_i*(n_i-1)
 * This function returns the TOTAL coalescent rate.
 */
// YC 2/23/2015
// changed the total coalescent rate as that used in Hey and Nielsen (2007)
// that is, t_i * n_i*(n_i-1)
void node::compute_totalCoalescentRate(std::list<double> coalTimes)
{
	totalCoalRate = 0.0;

	// The coalescent times on a tree are obtained.
	// Here the coalescent times mean the ages of internal nodes (including the root)
	// and it is sorted in ascending order (i.e., the age of the root comes at last).
	coalTimes.sort(); // sort the elements in ascending order
	list<double>::iterator iter_set, iter_prev;

	unsigned int n_lineages = size_tree();
	iter_set = coalTimes.begin();
	totalCoalRate = *iter_set * n_lineages *(n_lineages-1);
	iter_set++;
	for( ;iter_set != coalTimes.end(); ++iter_set)
	{
		n_lineages--;
		iter_prev = iter_set;
		--iter_prev;
		totalCoalRate += (*iter_set- *(iter_prev)) * n_lineages *(n_lineages-1);
	}

	// YC 2/23/2015
	// totalCoalRate /= 2.0;

	return;
}


/**
 * YC 3/4/2014
 * Compute the total coalescent rate of a tree
 * The total coalescent rate is similar to Eq.(16) in Hey and Nielsen (2007).
 * For example, for a time interval t_i, there are n_i lineages.
 * Then the corresponding coalescent rate used in this function is
 *     t_i * n_i*(n_i-1)/2.
 * Note the rate used in Hey and Nielsen (2007) is
 *     t_i * n_i*(n_i-1)
 * This function returns the TOTAL coalescent rate.
 */
// YC 2/23/2015
// changed the total coalescent rate as that used in Hey and Nielsen (2007)
// that is, t_i * n_i*(n_i-1)
void node::compute_totalCoalescentRate()
{
	totalCoalRate = 0.0;

	// The coalescent times on a tree are obtained.
	// Here the coalescent times mean the ages of internal nodes (including the root)
	// and it is sorted in ascending order (i.e., the age of the root comes at last).
	list<double> coalTimes; //(0.0,size_tree());
	list<double>::iterator iter_set, iter_prev;
	coalTimes = get_coalescentTimes(coalTimes); // extract the coalescent times from the tree
	coalTimes.sort(); // sort the elements in ascending order

	unsigned int n_lineages = size_tree();
	iter_set = coalTimes.begin();
	totalCoalRate = *iter_set * n_lineages *(n_lineages-1);
	iter_set++;
	for( ;iter_set != coalTimes.end(); ++iter_set)
	{
		n_lineages--;
		iter_prev = iter_set;
		--iter_prev;
		totalCoalRate += (*iter_set- *(iter_prev)) * n_lineages *(n_lineages-1);
	}

	// YC 2/23/2015	
	// totalCoalRate /= 2.0;

	return;
}

