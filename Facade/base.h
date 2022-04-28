#ifndef _BASE_H_
#define _BASE_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <time.h>
#include <cmath>
#include "../basic/Config.h"
#include "approximate_squared_distance.h"

#include <boost/shared_ptr.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>
#include <CGAL/algorithm.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Bbox_2.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel		K;
typedef CGAL::Exact_predicates_exact_constructions_kernel		K_epec;
typedef CGAL::Cartesian_converter<K_epec, K>		E2I;
typedef CGAL::Cartesian_converter<K, K_epec>		I2E;
typedef K::Point_2		Point_2;
typedef K::Vector_2		Vector_2;
typedef K::Segment_2		Segment_2;
typedef K::Line_2		Line_2;
typedef K::Intersect_2		Intersect_2;
typedef K::Ray_2		Ray_2;
typedef CGAL::Bbox_2		Bbox_2;
typedef pair<Point_2, Vector_2>		Point_with_normal2 ;
typedef vector<Point_with_normal2>		Input_range;
typedef K::Triangle_3		Triangle_3;
typedef K::Point_3		Point_3;
typedef K::Vector_3		Vector_3;
typedef K::Line_3 		Line_3;
typedef K::Plane_3 		Plane_3;
typedef pair<Point_3, Vector_3>		Point_with_normal;
typedef vector<Point_with_normal>		Pwn_vector;


/*!
 * \description: Cylinder class.
 */
class Circle
{
public:
	Point_2 center; 		// circle center
	double radius; 		// circle radius
	int number; 		// supporting points number
	int flag = 1; 		// valid (==1) or not (==0)
	double hmax; 		// cylinder max height
	double hmin; 		// cylinder min height

	Circle (Point_2 c, double r, double hx, double hn, int n=0):center(c), radius(r), number(n), hmax(hx), hmin(hn){};
	// descending order
	bool operator < (const Circle& s) { 
		return number > s.number;
	}
};

/*!
 * \description: 2D Line class.
 */
class Line
{
public:
	Line_2 line2; 		// line passing through the segment {s, t}
	vector<Point_3> pointset; 		// epsilon-nerghbor supporting points
	int ID = 0; 		// line id
	Point_2 s; 		// segment source
	Point_2 t; 		// segment target
	Vector_2 normal; 		// line oritation 
	double distance; 		// average distance between points and line
	double height0 = numeric_limits<double>::min(); 		// highest point z value
	double height1 = numeric_limits<double>::max(); 		// lowest point z value
	 
	// descending order
	bool operator < (const Line& line) { 
	   return pointset.size() > line.pointset.size();
	}
	bool operator == (const Line& line) const {
           return line2 == line.line2;   
	}
	
};

/*!
 * \description: 2D segment intersection class.
 */
class IntersectPoint
{
public:
	Point_2 p;
	vector<int> IDs; 		// segment ids
	vector<int> IDl; 		// line ids
	int flag = 0; 		// boundary (==1) or not (==0)
	Point_2 heights; 		// the height of point {high, low}
	IntersectPoint(){}
	IntersectPoint(Point_2 p_, Point_2 h, vector<int>&& ids):p(p_), heights(h), IDs(ids){};
};

/*!
 * \description: 2D segment class.
 */
class Segment
{
public:
	Segment_2 segment2;
	int ID = 0; 		// id of line which the segment belongs to 
	int flags = 0; 		// source point is not on the boundary
	int flagt = 0; 		// target point is not on the boundary

	Segment(Segment_2 seg, int id, int s=0, int t=0){
	  ID = id;
	  segment2 = seg;
	  flags = s;
	  flagt = t;
	}
	Segment(Segment_2 seg){
	  segment2 = seg;
	}
	bool operator < (const Segment& s) { 
	   return segment2.source().x() < s.segment2.source().x() || segment2.source().x() == s.segment2.source().x() && segment2.source().y() < s.segment2.source().y();
	}
};

/*!
 * \description: Facade class for potential facades.
 */
class Facade
{
public:
   K_epec::Segment_2 seg; 		// exact 2D projection segment {s, t}
   vector<double> heights; 		// facade height {s-high, s-low, t-high, t-low}
   Vector_2 normal; 		// segment normal
   int flag = 0; 		// selected potential facade (==1) or not (==0)
   Facade(K_epec::Segment_2 s, vector<double> h):seg(s), heights(h){};
   Facade(){}
};

/*!
 * \description: Mesh class for storage.
 */
template <class T>
class Mesh
{
public:
	vector<T> vertices;
	vector<T> faces;
};


/*!
 * \description: Judge number range
 */
template <class T>
inline bool number_range(T x){
	return x < numeric_limits<T>::max() && x > numeric_limits<T>::lowest();
}

/*!
 * \description: Judge the validation of 2D point
 */
inline bool point_2_range(Point_2 p){
  return number_range<double>(p.x()) && number_range<double>(p.y());
}


/******************************FileIO.cpp******************************/
/*!
 * \description: Load 3D pointcloud from ply file
 * \param fileName pointcloud file name to be loaded
 * \param points 3D pointcloud with normal
 */
bool PLYPointLoad3(const string fileName, Pwn_vector& points);

/*!
 * \description: Save 3D pointcloud into ply file
 * \param fileName pointcloud file name to be saved
 * \param points 3D pointcloud with normal
 * \param type storage format (ascii-1, binary_big_endian-2, binary_little_endian-3)
 * \return 
 */
bool PLYPointSave(const string fileName, Pwn_vector& points, int type = 3);

/*!
 * \description: Save 3D mesh into ply file
 * \param fileName mesh file name to be saved
 * \param mesh 3D mesh
 * \param type storage format (ascii-1, binary_big_endian-2, binary_little_endian-3)
 * \return 
 */
bool PLYMeshSave(const string fileName, Mesh<Point_3>& mesh, int type = 3);

/*!
 * \description: Save 3D mesh into off file
 * \param fileName mesh file name to be saved
 * \param mesh 3D mesh
 * \return 
 */
bool OFFMeshSave(const string fileName, Mesh<K_epec::Point_3>& mesh);

/*!
 * \description: Load 3D mesh from ply file
 * \param fileName file name to be loaded
 * \param mesh 3D mesh
 * \return 
 */
bool PLYMeshLoad(string fileName, Mesh<Point_3>& mesh);

/*!
 * \description: Save 3D mesh into ply file with adaptive triangle facets
 * \param mesh 3D triangle mesh
 * \param config configuration file
 * \param filename mesh file name to be saved
 * \return
 */
bool PLYTriMeshSave(Mesh<Point_3>& mesh, const cm::Config& config, string& filename);


/******************************ransac_detection.cpp******************************/
/*!
 * \description: detect 3D facade planes from facade pointcloud by RANSAC
 * \param points facade pointcloud
 * \param probability RANSAC parameter
 * \param min_points RANSAC parameter
 * \param epsilon RANSAC parameter
 * \param cluster_epsilon RANSAC parameter
 * \param normal_threshold RANSAC parameter
 * \param cos_angle RANSAC parameter
 * \param lines facade projection lines
 * \param wdir file save directory
 * \return 
 */
bool ransac_detection_p(Pwn_vector& points, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, vector<Line>& lines, string wdir);

/*!
 * \description: detect 3D cylinder from cylinder pointcloud by RANSAC
 * \param points cylinder pointcloud
 * \param radius cylinder radius threshold
 * \param probability RANSAC parameter
 * \param min_points RANSAC parameter
 * \param epsilon RANSAC parameter
 * \param cluster_epsilon RANSAC parameter
 * \param normal_threshold RANSAC parameter
 * \param cos_angle RANSAC parameter
 * \param circles facade projection lines
 * \param wdir file save directory
 * \return 
 */
bool ransac_detection_c(Pwn_vector& points, float radius, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, vector<Circle>& circles, string wdir);


bool ransac_detection_fc(Pwn_vector& points, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, string wdir, string object);

bool ransac_detection_floor(Pwn_vector& points,  vector<Plane_3>& floors, vector<Line_2>& intersects, const cm::Config& config, string wdir);

bool ransac_detection_ceiling(Pwn_vector& points,  vector<Plane_3>& ceilings, vector<Line_2>& intersects, const cm::Config& config, string wdir);

// hypothesis_generator.cpp
bool generation(vector<Line>& lines, vector<Segment>& bsegments, vector<Segment>& segments, unordered_map<Point_2, int>& points_idx, vector<IntersectPoint>& points, float angle_thetad, float num_threshold, float extend_length, string wdir);

// segment_selection.cpp
void optimize(vector<Segment_2>& selected, unordered_map<Point_2, int>& points_idx, vector<Point_2>& heights, vector<Vector_2>& oritations, vector<Line>& lines, vector<Segment>& segments, vector<IntersectPoint>& points, vector<Segment>& bsegments, double lambda_data_fitting, double lambda_model_coverage, double lambda_model_complexity, double epsilon, double alpha, string wdir);
void optimize2(vector<Facade>& potential_facade, vector<IntersectPoint>& points, vector<int>& suppoints_num, vector<double>& area_ratio, double lambda_data_fitting, double lambda_model_coverage, double lambda_model_complexity, string wdir);
void normal_adjust(vector<Segment_2> segments, vector<Vector_2>& oritations);
void producemesh(vector<Segment_2>& selected, Mesh<K_epec::Point_3>& planes, vector<Point_2>& heights, vector<Vector_2>& oritations, double top=0, double bottom=1);
bool segment_wall(vector<Segment_2> selected, const cm::Config& config, string wdir);

// kdtree_search.cpp
float density(int size, vector<Point_3>& pointsets);

//generate_fc_map.cpp
bool generate_fc(vector<Segment_2> & segments, const Point_3 corner, const Point_2 length, const cm::Config& config, string wdir);
bool generate_fc(vector<Segment_2> segments, double top, double bottom, string wdir);
bool generate_floor(vector<Plane_3>& floors, cv::Mat& hfmap, cv::Mat& nfmap, vector<Line_2>& intersectf, const cm::Config& config, Point_3 corner, double step, string wdir);
bool generate_ceiling(vector<Plane_3>& ceilings, cv::Mat& hcmap, cv::Mat& ncmap, vector<Line_2>& intersectc, const cm::Config& config, Point_3 corner, double step, string wdir);
bool generate_map(string wdir, double size);
bool generate_modeling(string wdir);

/*!
 * \description: Approximate segment to segment squared distance.
 * \param {s1} first segment source point 
 * \param {t1} first segment target point
 * \param {s2} second segment source point
 * \param {t2} second segment target point
 * \return {*} squared distance
 */
inline double robust_squared_distance(const Point_2 &s1, const Point_2 &t1, const Point_2 &s2, const Point_2 &t2) {
  return segment_segment_squared_distance(
    s1.x(), s1.y(), t1.x(), t1.y(),
    s2.x(), s2.y(), t2.x(), t2.y());
}

#endif // _BASE_H_







