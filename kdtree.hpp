/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */

    if (curDim >= Dim) {
      return false;
    }

    if (first[curDim] < second[curDim]) {
      return true;
    }
    if (first[curDim] > second[curDim]) {
      return false;
    }

    if (first[curDim] == second[curDim]) {
      return (first < second);
    }

    return false;
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     * 
     */
    unsigned potentialToTarget = 0.0;
    unsigned currentToTarget = 0.0;

    unsigned i = 0;
    while (i < Dim) { 
      potentialToTarget += (target[i] - potential[i]) * (target[i] - potential[i]);
      currentToTarget += (currentBest[i] - target[i]) * (currentBest[i] - target[i]);
      i++;
    }

    potentialToTarget = sqrt(potentialToTarget);
    currentToTarget = sqrt(currentToTarget);

    if (potentialToTarget < currentToTarget) { 
      return true;
    }
    if (potentialToTarget > currentToTarget) {
      return false;
    } else {
      return (potential < currentBest);
    }
    
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
    if(!(newPoints.empty())){
      vector<Point<Dim>> points = newPoints;
       size = 0;
       constructHelper(points, 0, points.size() - 1, 0, root);
    }
    else {
      root = NULL;
      size = 0;
    }
}

template <int Dim>
void KDTree<Dim>::constructHelper(vector<Point<Dim>>& newPoints, int left, int right, int dimension, KDTreeNode*& subroot){
  if (left > right) {
    return;
  }
  size_t median = (left + right) / 2;
  subroot = new KDTreeNode(select(newPoints, left, right, median, dimension));
  size++;
  constructHelper(newPoints, left, median - 1, (dimension + 1) % Dim, subroot->left);
  constructHelper(newPoints, median + 1, right, (dimension + 1) % Dim, subroot->right);
}

template <int Dim>
int KDTree<Dim>::partition(vector<Point<Dim>>& list, int left, int right, size_t pivotIndex, int dimension){
  Point<Dim> pivotValue = list[pivotIndex];
    Point<Dim> pivotTemp = list[pivotIndex];
  	list[pivotIndex] = list[right];
  	list[right] = pivotTemp;
  size_t storeIndex = left;
  for (int i = left; i < right; i++) {
    if (smallerDimVal(list[i], pivotValue, dimension)) {
      pivotTemp = list[storeIndex];
      list[storeIndex] = list[i];
      list[i] = pivotTemp;
      storeIndex++;
    }
  }
  pivotTemp = list[storeIndex];
  list[storeIndex] = list[right];
  list[right] = pivotTemp;
  return storeIndex;
}

template <int Dim>
Point<Dim> KDTree<Dim>::select(vector<Point<Dim>>& list, int left, int right, size_t k, int dimension){
  if (left == right) {  
    return list[left]; 
  }
  size_t pivotIndex = (left + right) / 2;
  pivotIndex = partition(list, left, right, pivotIndex, dimension);
  if (k == pivotIndex){
    return list[k];
  }
  else if (k < pivotIndex) {
    return select(list, left, pivotIndex - 1, k, dimension);
  }
  else {
    return select(list, pivotIndex + 1, right, k, dimension);
  }
}


template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
  copy(this->root, other.root);
  size = other.size;

}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */
  remove(root);
  copy(root, rhs->root);
  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
  remove(root);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    Point<Dim> best;
    findNN(root, query, 0, best);
    return best;
}

template <int Dim>
void KDTree<Dim>::findNN(KDTreeNode* child, const Point<Dim>& query, int dimension, Point<Dim>& best) const {
  best = child->point;
  bool check;
  if (child->left == NULL && child->right == NULL) {
    best = child->point;
    return;
  }

  if (smallerDimVal(query, child->point, dimension)) {
    if (child->left == NULL) { 
      findNN(child->right, query, (dimension + 1) % Dim, best);
    } else {
      findNN(child->left, query, (dimension + 1) % Dim, best);
    }
    check = true;
  } else {
    if (child->right == NULL) { 
      findNN(child->left, query, (dimension + 1) % Dim, best);
    } else {
      findNN(child->right, query, (dimension + 1) % Dim, best);
    }
    check = false;
  }
  if (shouldReplace(query, best, child->point)) {
    best = child->point;
  }

  double distance = getDistance(query, child->point, dimension);
  double radius = 0;
  int index = 0;
  while (index < Dim) {
    radius += pow(query[index] - best[index], 2);
    index++;
  }
  
  
  if (distance <= radius) {
    KDTreeNode * c = NULL;
    if (check) {
      c = child->right;
    } else {
      c = child->left;
    }
		if (c != NULL) {
			Point<Dim> other;
      findNN(c, query, (dimension + 1) % Dim, other);
			if (shouldReplace(query, best, other)) {
        best = other;
      } 
		}
	}
  return;
}


template <int Dim>
double KDTree<Dim>::getDistance(const Point<Dim> &q, const Point<Dim> &s, int dim) const {
  double distance = 0.0;
      distance = pow((q[dim] - s[dim]), 2);
      return distance;
}

template <int Dim>
void KDTree<Dim>::copy(KDTreeNode * subroot, KDTreeNode * other) {
	subroot = new KDTreeNode();
	subroot->point = other->point;
	copy(subroot->left, other->left);
	copy(subroot->right, other->right);
}

template <int Dim>
void KDTree<Dim>::remove(KDTreeNode * subroot) {
	if (subroot == NULL) return;
	if (subroot->left != NULL) {
    remove(subroot->left);
  } 
	if (subroot->right != NULL) {
    remove(subroot->right);
  }
	delete subroot;
	subroot = NULL;
}

