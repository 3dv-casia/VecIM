#include "base.h"
#include "kdtree/kdTree.h"


void find_closest_K_points(
	kdtree::KdTree* tree, const Point_3 p, unsigned int k, vector<unsigned int>& neighbors, vector<double>& squared_distances
	){
	kdtree::Vector3D v3d( p.x(), p.y(), 0 );
	tree->setNOfNeighbours( k );
	tree->queryPosition( v3d );

	unsigned int num = tree->getNOfFoundNeighbours();
	if (num == k) {
	  neighbors.resize(k);
	  squared_distances.resize(k);
	  for (unsigned int i=0; i<k; ++i) {
		neighbors[i] = tree->getNeighbourPositionIndex(i);
		squared_distances[i] = tree->getSquaredDistance(i);
	  }		
	} else
		cout << "less than " << k << " points found" << endl;
}

float density(int size, vector<Point_3>& pointsets){

   int points_num_ = pointsets.size();
   kdtree::Vector3D* points = new kdtree::Vector3D[points_num_];
   for (size_t i = 0; i < pointsets.size(); ++i){
	points[i].x = pointsets[i].x();
	points[i].y = pointsets[i].y();
	points[i].z = 0;
   }

   unsigned int maxBucketSize = 16 ;	// number of points per bucket
   kdtree::KdTree* tree = new kdtree::KdTree(points, points_num_, maxBucketSize );
   delete [] points;
	
   double total = 0;
   vector<unsigned int> neighbors;
   vector<double> sqr_distances;
   for (int i = 0; i < points_num_; i++){
	find_closest_K_points(tree, pointsets[i], size, neighbors, sqr_distances);
	double avg = 0;
	for (unsigned int k = 0; k < neighbors.size(); ++k)
	   avg += sqrt(sqr_distances[k]);
	total += (avg / neighbors.size());
   }
 
   return static_cast<float>(total / points_num_);
}




