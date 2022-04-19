/****** CGAL ******/
#include <boost/shared_ptr.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>
#include <CGAL/algorithm.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Bbox_2.h>

/****** OPENCV ******/
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

/****** SRANDARD ******/
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

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel K_epec;
typedef K::Triangle_3 Triangle_3;
typedef K::Point_2 		Point_2;
typedef K::Vector_2		Vector_2;
typedef K::Segment_2 		Segment_2;
typedef K::Line_2 		Line_2;
typedef K::Intersect_2 		Intersect_2;
typedef CGAL::Bbox_2		Bbox_2;
using Point_with_normal2 = pair<Point_2, Vector_2>;
using Input_range       = vector<Point_with_normal2>;

typedef K::Point_3		Point_3;
typedef K::Vector_3		Vector_3;
typedef K::Line_3 		Line_3;
typedef K::Plane_3 		Plane_3;
typedef pair<Point_3, Vector_3>         Point_with_normal;
typedef vector<Point_with_normal>                       Pwn_vector;



/****** key class ******/

class Circle
{
public:
  Point_2 center;
  double radious;
  int number;
  int flag = 1;
  double hmax;
  double hmin;
  Circle (Point_2 c, double r, double hx, double hn, int n=0):center(c), radious(r), number(n), hmax(hx), hmin(hn){};
  bool operator < (const Circle& s) { 
	   return number > s.number;
	}
};


class Line
{
public:
	Line_2 line2;
	vector<Point_3> pointset; // epsilon-nerghbor points (or store points index)
	int ID; // line id number
	Point_2 s; // segment source
	Point_2 t; // target
	Vector_2 normal; // oritation 
	float distance; // ave distance between points and line
	float height0 = -1e9; // highest point z axis
	float height1 = 1e9; // lowest point z axis
	 
	Line(){ 
	 ID = 0;
	}
	// decreasing order
	bool operator < (const Line& line) { 
	   return pointset.size()>line.pointset.size();
	}
	bool operator == (const Line& line) const {
           return line2 == line.line2;   
	}
	
};

class IntersectPoint
{
public:
	Point_2 p;
	vector<int> IDs; // segments id
	vector<int> IDl; // line id
	int flag = 0; // boundary(1) or not
	Point_2 heights; // point height {high, low}
	IntersectPoint(){}
	IntersectPoint(Point_2 p_, Point_2 h, vector<int> ids):p(p_), heights(h), IDs(ids){};
};

class Segment
{
public:
	Segment_2 segment2;
	int ID = 0; // id number of line which it belong to 
	int flags = 0; // source is not on boundary
	int flagt = 0; // target is not on boundary

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

class Facade
{
public:
   K_epec::Segment_2 seg; // {s,t}
   vector<double> heights; // {s-high,s-low,t-high,t-low}
   Vector_2 normal;
   int flag = 0; // selected potential_facade -> flag = 1
   Facade(K_epec::Segment_2 s, vector<double> h):seg(s),heights(h){};
   Facade(){}
};

class Mesh
{
public:
	vector<Point_3> vertices;
	vector<Point_3> faces;
};


/****** function defination ******/

// judge number range
template <class T>
inline bool number_range(T x){
	return x < numeric_limits<T>::max() && x > numeric_limits<T>::lowest();
}

inline bool point_2_range(Point_2 p){
  return number_range<double>(p.x()) && number_range<double>(p.y());
}

// FileIO.cpp
bool PLYPointLoad3(const string fileName, Pwn_vector& points);
bool PLYPointSave(const string fileName, Pwn_vector& points, int type=3);
bool PLYMeshSave(const string fileName, Mesh& mesh, int type=3);
bool OFFMeshSave(const string fileName, Mesh& mesh);
bool PLYMeshLoad(string fileName, Mesh& mesh);
bool PLYTriMeshSave(Mesh& mesh, const cm::Config& config, string& name);

// ransac_detection.cpp
bool ransac_detection_p(Pwn_vector& points, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, vector<Line>& lines, string wdir);
bool ransac_detection_c(Pwn_vector& points, float radious, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, vector<Circle>& circles, string wdir);
bool ransac_detection_fc(Pwn_vector& points, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, string wdir, string object);
bool ransac_detection_floor(Pwn_vector& points,  vector<Plane_3>& floors, vector<Line_2>& intersects, const cm::Config& config, string wdir);
bool ransac_detection_ceiling(Pwn_vector& points,  vector<Plane_3>& ceilings, vector<Line_2>& intersects, const cm::Config& config, string wdir);

// hypothesis_generator.cpp
bool generation(vector<Line>& lines, vector<Segment>& bsegments, vector<Segment>& segments, unordered_map<Point_2, int>& points_idx, vector<IntersectPoint>& points, float angle_thetad, float num_threshold, float extend_length, string wdir);

// segment_selection.cpp
void optimize(vector<Segment_2>& selected, unordered_map<Point_2, int>& points_idx, vector<Point_2>& heights, vector<Vector_2>& oritations, vector<Line>& lines, vector<Segment>& segments, vector<IntersectPoint>& points, vector<Segment>& bsegments, double lambda_data_fitting, double lambda_model_coverage, double lambda_model_complexity, double epsilon, double alpha, string wdir);
void optimize2(vector<Facade>& potential_facade, vector<IntersectPoint>& points, vector<int>& isolate_seg, vector<int>& suppoints_num, vector<double>& area_ratio, double lambda_data_fitting, double lambda_model_coverage, double lambda_model_complexity, string wdir);
void normal_adjust(vector<Segment_2> segments, vector<Vector_2>& oritations);

// segment_selection.cpp: plane mesh
void producemesh(vector<Segment_2>& selected, Mesh& planes, vector<Point_2>& heights, vector<Vector_2>& oritations, double top=0, double bottom=1);

// kdtree_search.cpp
float density(int size, vector<Point_3>& pointsets);

//generate_fc_map.cpp
bool generate_fc(vector<Segment_2> & segments, const Point_3 corner, const Point_2 length, const cm::Config& config, string wdir);
bool generate_fc(vector<Segment_2> segments, double top, double bottom, string wdir);
bool generate_floor(vector<Plane_3>& floors, cv::Mat& hfmap, cv::Mat& nfmap, vector<Line_2>& intersectf, const cm::Config& config, Point_3 corner, double step, string wdir);
bool generate_ceiling(vector<Plane_3>& ceilings, cv::Mat& hcmap, cv::Mat& ncmap, vector<Line_2>& intersectc, const cm::Config& config, Point_3 corner, double step, string wdir);
bool generate_map(string wdir, double size);
bool generate_modeling(string wdir);

// segment_selection.cpp
bool segment_wall(vector<Segment_2> selected, const cm::Config& config, string wdir);

/*!
 * \brief Approximate segment to segment squared distance.
 * \note Possible underflow and overflow of float number type.
 * \param s1 first segment source point 
 * \param t1 first segment target point
 * \param s2 second segment source point
 * \param t2 second segment target point
 * \return squared distance
 */
inline double robust_squared_distance(const Point_2 &s1, const Point_2 &t1, const Point_2 &s2, const Point_2 &t2) {
  return segment_segment_squared_distance(
    s1.x(), s1.y(), t1.x(), t1.y(),
    s2.x(), s2.y(), t2.x(), t2.y());
}








