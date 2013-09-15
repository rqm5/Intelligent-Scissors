/* iScissor.cpp */
/* Main file for implementing project 1.  See TODO statments below
 * (see also correlation.cpp and iScissor.h for additional TODOs) */

#include <assert.h>

#include "correlation.h"
#include "iScissor.h"
#include "PriorityQueue.h"

const double linkLengths[8] = { 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2 };

// two inlined routines that may help;

inline Node& NODE(Node* n, int i, int j, int width)
{
    return *(n + j * width + i);
}

inline unsigned char PIXEL(const unsigned char* p, int i, int j, int c, int width)
{
    return *(p + 3 * (j * width + i) + c);
}

void assignCoordinates (Node* nodes, int imgWidth, int imgHeight);
void computeDerivatives (Node* nodes, const unsigned char* img, int imgWidth, int imgHeight);
double getMaxD (Node* nodes, int imgWidth, int imgHeight);
void computeCost (Node* nodes, int imgWidth, int imgHeight, double maxD);
void initNodeState(Node* nodes, int width, int height);

/************************ TODO 1 ***************************
 *InitNodeBuf
 *	INPUT:
 *		img:	a RGB image of size imgWidth by imgHeight;
 *		nodes:	a allocated buffer of Nodes of the same size, one node corresponds to a pixel in img;
 *  OUPUT:
 *      initializes the column, row, and linkCost fields of each node in the node buffer.
 */

void InitNodeBuf(Node* nodes, const unsigned char* img, int imgWidth, int imgHeight)
{
	assignCoordinates(nodes, imgWidth, imgHeight);
	computeDerivatives(nodes, img, imgWidth, imgHeight);
	double maxD = getMaxD(nodes, imgWidth, imgHeight);
	computeCost(nodes, imgWidth, imgHeight, maxD);
}


/*
 *assignCoordinates
 *	INPUT:
 *		nodes:	a allocated buffer of Nodes of the same size, one node corresponds to a pixel in img. the linkCosts in
 *		        each of the nodes is the magnitude of the derivative in each of the 8 directions at that pixel.
 *		img:	a RGB image of size imgWidth by imgHeight;
 *  OUPUT:
 *      returns the buffer of Nodes with row and column values assigned for each Node.
 */

void assignCoordinates (Node* nodes, int imgWidth, int imgHeight) {
	for (int x = 0; x < imgWidth; ++x) {
		for (int y = 0; y < imgHeight; ++y) {
			int nodeIndex = y*imgWidth+x;
			nodes[nodeIndex].column = x;
			nodes[nodeIndex].row = y;
		}
	}
}


/*
 *computeDerivatives
 *	INPUT:
 *		nodes:	a allocated buffer of Nodes of the same size, one node corresponds to a pixel in img. the linkCosts in
 *		        each of the nodes is the magnitude of the derivative in each of the 8 directions at that pixel.
 *		img:	a RGB image of size imgWidth by imgHeight;
 *  OUPUT:
 *      returns the buffer of Nodes with the derivatives in each of the 8 directions stored in linkCost for each Node.
 */

void computeDerivatives (Node* nodes, const unsigned char* img, int imgWidth, int imgHeight) {
	int kernelWidth = 3;
	int kernelHeight = 3;
	double scale = 1.0000;
	double offset = 0.0000;
	const unsigned char* selection = NULL;
	for (int linkIndex = 0; linkIndex < 8; ++linkIndex) {
		double* rsltImg = new double[3*imgHeight*imgWidth];
		image_filter(rsltImg, img, selection, imgWidth, imgHeight, kernels[linkIndex], kernelWidth, kernelHeight, scale, offset);
		for (int x = 0; x < imgWidth; ++x) {
			for (int y = 0; y < imgHeight; ++y) {
				int nodeIndex = y*imgWidth+x;
				nodes[nodeIndex].linkCost[linkIndex] = (pow(rsltImg[3*(y*imgWidth+x)], 2.0000) + pow(rsltImg[3*(y*imgWidth+x)+1], 2.0000) + 
					                                   pow(rsltImg[3*(y*imgWidth+x)+2], 2.0000))/3;
			}
		}
		delete[] rsltImg;
	}
}


/*
 *getMaxD
 *	INPUT:
 *		nodes:	a allocated buffer of Nodes of the same size, one node corresponds to a pixel in img. the linkCosts in
 *		        each of the nodes is the magnitude of the derivative in each of the 8 directions at that pixel.
 *		img:	a RGB image of size imgWidth by imgHeight;
 *  OUPUT:
 *      returns the maximum derivative magnitude.
 */

double getMaxD (Node* nodes, int imgWidth, int imgHeight) {
	double maxD = 0.0000;
	for (int x = 0; x < imgWidth; ++x) {
		for (int y = 0; y < imgHeight; ++y) {
			int nodeIndex = y*imgWidth+x;
			for (int linkIndex = 0; linkIndex < 8; ++linkIndex) {
				if (nodes[nodeIndex].linkCost[linkIndex] > maxD) {
					maxD = nodes[nodeIndex].linkCost[linkIndex];
				}
			}
		}
	}
	return maxD;
}


/*
 *computeCost
 *	INPUT:
 *		nodes:	a allocated buffer of Nodes of the same size, one node corresponds to a pixel in img. the linkCosts in
 *		        each of the nodes is the magnitude of the derivative in each of the 8 directions at that pixel.
 *		img:	a RGB image of size imgWidth by imgHeight;
 *		maxD:	the maximum derivative magnitude among all of the node links.
 *  OUPUT:
 *      returns the buffer of Nodes with the correct linkCost values.
 */

void computeCost (Node* nodes, int imgWidth, int imgHeight, double maxD) {
	for (int x = 0; x < imgWidth; ++x) {
		for (int y = 0; y < imgHeight; ++y) {
			int nodeIndex = y*imgWidth+x;
			for (int linkIndex = 0; linkIndex < 8; ++linkIndex) {
				double length = 1.0000;
				if (linkIndex%2 == 1) {
					length = SQRT2;
				}
				nodes[nodeIndex].linkCost[linkIndex] = (maxD - nodes[nodeIndex].linkCost[linkIndex]) * length;
			}
		}
	}
}

/************************ END OF TODO 1 ***************************/

static int offsetToLinkIndex(int dx, int dy)
{
    int indices[9] = { 3, 2, 1, 4, -1, 0, 5, 6, 7 };
    int tmp_idx = (dy + 1) * 3 + (dx + 1);
    assert(tmp_idx >= 0 && tmp_idx < 9 && tmp_idx != 4);
    return indices[tmp_idx];
}

/************************ TODO 4 ***************************
 *LiveWireDP:
 *	INPUT:
 *		seedX, seedY:	seed position in nodes
 *		nodes:			node buffer of size width by height;
 *      width, height:  dimensions of the node buffer;
 *		selection:		if selection != NULL, search path only in the subset of nodes[j*width+i] if selection[j*width+i] = 1;
 *						otherwise, search in the whole set of nodes.
 *		numExpanded:		compute the only the first numExpanded number of nodes; (for debugging)
 *	OUTPUT:
 *		computes the minimum path tree from the seed node, by assigning
 *		the prevNode field of each node to its predecessor along the minimum
 *		cost path from the seed to that node.
 */

void LiveWireDP(int seedX, int seedY, Node* nodes, int width, int height, const unsigned char* selection, int numExpanded)
{
	Node* seed = &NODE(nodes, seedX, seedY, width);
	seed->prevNode = NULL;
	seed->totalCost = 0;

	if (selection != NULL && selection[seedY * width + seedX] == 0)
		return; //Exit function if the seed isn't in the selection.

	initNodeState(nodes, width, height);

	CTypedPtrHeap<Node> pq;
	pq.Insert(seed);

	while (!pq.IsEmpty()) {
		Node* minNode = pq.ExtractMin();
		minNode->state = EXPANDED;

		for (int linkIndex = 0; linkIndex < 8; ++linkIndex) {

			int offsetX, offsetY;
			minNode->nbrNodeOffset(offsetX, offsetY, linkIndex);

			int neighborX = minNode->column + offsetX;
			int neighborY = minNode->row + offsetY;

			if ((0 <= neighborX && neighborX < width) && (0 <= neighborY && neighborY < height)) {
				Node* neighborNode = &(NODE(nodes, neighborX, neighborY, width));

				//Mildly convuoluted if statement to avoid NullPointerExceptions: statement will break 
				//before trying to dereference selection if selection == NULL.
				if (neighborNode->state != EXPANDED && !(selection != NULL && selection[neighborY * width + neighborX] != 1)) {
					if (neighborNode->state == INITIAL) {
						neighborNode->state = ACTIVE;

						neighborNode->totalCost = minNode->totalCost + minNode->linkCost[linkIndex];
						neighborNode->prevNode = minNode;
						pq.Insert(neighborNode);
					}
					else {
						double currentCost = minNode->totalCost + minNode->linkCost[linkIndex];

						if (currentCost < neighborNode->totalCost) {
							neighborNode->totalCost = currentCost;
							neighborNode->prevNode = minNode;
							pq.Update(neighborNode);
						}
					}
				}
			}
		}
	}
	
}
/************************ END OF TODO 4 ***************************/

/************************ TODO 5 ***************************
 *MinimumPath:
 *	INPUT:
 *		nodes:				a node buffer of size width by height;
 *		width, height:		dimensions of the node buffer;
 *		freePtX, freePtY:	an input node position;
 *	OUTPUT:
 *		insert a list of nodes along the minimum cost path from the seed node to the input node.
 *		Notice that the seed node in the buffer has a NULL predecessor.
 *		And you want to insert a *pointer* to the Node into path, e.g.,
 *		insert nodes+j*width+i (or &(nodes[j*width+i])) if you want to insert node at (i,j), instead of nodes[nodes+j*width+i]!!!
 *		after the procedure, the seed should be the head of path and the input code should be the tail
 *  NB PROF. SNAVELY:
 *		If the seed is invalid (outside the selection), this will return a one-point path consisting of the input node.
 */

void MinimumPath(CTypedPtrDblList <Node>* path, int freePtX, int freePtY, Node* nodes, int width, int height)
{
	int inputNodeIndex = freePtY * width + freePtX;
	Node* currentNode = &(NODE(nodes, freePtX, freePtY, width));
	while (currentNode != NULL) {
		path->AddHead(currentNode);
		currentNode = currentNode->prevNode;
	}
}


/*
 *initNodeState
 *	INPUT:
 *		nodes:	a allocated buffer of Nodes of the same size, one node corresponds to a pixel in img. the linkCosts in
 *		        each of the nodes is the magnitude of the derivative in each of the 8 directions at that pixel.
 *  OUPUT:
 *      returns the buffer of Nodes with the states of the Nodes set to INITIAL.
 */

void initNodeState(Node* nodes, int width, int height) {
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			NODE(nodes, x, y, width).state = INITIAL;
		}
	}
}


/*
 *getSeedIndex
 *	INPUT:
 *		nodes:	a allocated buffer of Nodes of the same size, one node corresponds to a pixel in img. the linkCosts in
 *		        each of the nodes is the magnitude of the derivative in each of the 8 directions at that pixel.
 *  OUPUT:
 *      returns the index of the seed Node. if none are found, returns -1.
 */

int getSeedIndex(Node* nodes, int width, int height) {
	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			int nodeIndex = y*width+x;
			if (nodes[nodeIndex].prevNode == NULL) {
				return nodeIndex;
			}
		}
	}
	return -1;
}

/************************ END OF TODO 5 ***************************/

/************************ An Extra Credit Item ***************************
 *SeedSnap:
 *	INPUT:
 *		img:				a RGB image buffer of size width by height;
 *		width, height:		dimensions of the image buffer;
 *		x,y:				an input seed position;
 *	OUTPUT:
 *		update the value of x,y to the closest edge based on local image information.
 */

void SeedSnap(int& x, int& y, unsigned char* img, int width, int height)
{
    printf("SeedSnap in iScissor.cpp: to be implemented for extra credit!\n");
}

//generate a cost graph from original image and node buffer with all the link costs;
void MakeCostGraph(unsigned char* costGraph, const Node* nodes, const unsigned char* img, int imgWidth, int imgHeight)
{
    int graphWidth = imgWidth * 3;
    int graphHeight = imgHeight * 3;
    int dgX = 3;
    int dgY = 3 * graphWidth;

    int i, j;
    for (j = 0; j < imgHeight; j++) {
        for (i = 0; i < imgWidth; i++) {
            int nodeIndex = j * imgWidth + i;
            int imgIndex = 3 * nodeIndex;
            int costIndex = 3 * ((3 * j + 1) * graphWidth + (3 * i + 1));

            const Node* node = nodes + nodeIndex;
            const unsigned char* pxl = img + imgIndex;
            unsigned char* cst = costGraph + costIndex;

            cst[0] = pxl[0];
            cst[1] = pxl[1];
            cst[2] = pxl[2];

            //r,g,b channels are grad info in seperate channels;
            DigitizeCost(cst	   + dgX, node->linkCost[0]);
            DigitizeCost(cst - dgY + dgX, node->linkCost[1]);
            DigitizeCost(cst - dgY      , node->linkCost[2]);
            DigitizeCost(cst - dgY - dgX, node->linkCost[3]);
            DigitizeCost(cst	   - dgX, node->linkCost[4]);
            DigitizeCost(cst + dgY - dgX, node->linkCost[5]);
            DigitizeCost(cst + dgY	   ,  node->linkCost[6]);
            DigitizeCost(cst + dgY + dgX, node->linkCost[7]);
        }
    }
}

