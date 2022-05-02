#include "base.h"
#include <glog/logging.h>

#define pi 3.1415926 

/*!
 * \description: Calculate line parameters by least square method.
 * \param k slope
 * \param b intercept
 * \return
 */
void least_square_line(double& k, double& b, vector<Point_2>& pointset)
{
   double A=0; double B=0; double C=0; double D=0;
   int N=0;
   double dis=0;
   for (int i = 0; i < pointset.size(); i++) 
   {
		Point_2 p = pointset[i];
		A+=p.x();
		B+=p.y();
		C+=p.x()*p.x();
		D+=p.x()*p.y();
		N++;
    }
    // compute line parameter
    k = (A*B-N*D)/(A*A-N*C);
    b = (A*D-B*C)/(A*A-N*C);
}

/*!
 * \description: Merge line1 and line2.
 * \param l1 line1
 * \param l2 line2
 * \param lines total lines
 * \param extend segment extension scale
 * \return
 */
bool merge(Line& l1, Line& l2, vector<Line>& lines, float extend)
{
	Line m;
	Segment_2 s1(l1.s, l1.t);
	Segment_2 s2(l2.s, l2.t);
	Point_2 mid1 = CGAL::midpoint(l1.s, l1.t);
	Point_2 mid2 = CGAL::midpoint(l2.s, l2.t);
	double len1 = CGAL::approximate_sqrt(s1.squared_length());
	double len2 = CGAL::approximate_sqrt(s2.squared_length());
	if(len1 < 1e-4 || len2 < 1e-4)
		LOG(INFO) << "Too small vector.";
	Point_2 midn((len1*mid1.x() + len2*mid2.x())/(len1+len2), (len1*mid1.y() + len2*mid2.y())/(len1+len2));
	Vector_2 dir = s1.to_vector() * s2.to_vector() > 0 ? (s1.to_vector() + s2.to_vector()) : (s1.to_vector() - s2.to_vector());
	m.line2 = Line_2(midn, dir);
	double dis = 0;
	double xmin = 1e9;
	double xmax = -1e9;
	double x1, y1, x2, y2;
	for (int i = 0; i < l1.pointset.size(); i++) 
	{
		Point_2 p = Point_2(l1.pointset[i].x(), l1.pointset[i].y());
		m.pointset.push_back(l1.pointset[i]);
		Point_2 pp = m.line2.projection(p);
		dis += CGAL::sqrt((p-pp).squared_length());
		if (xmin > pp.x() || xmin == pp.x() && y1 > pp.y()) { xmin = pp.x(); x1 = p.x(); y1 = p.y(); }
		if (xmax < pp.x() || xmax == pp.x() && y2 < pp.y()) { xmax = pp.x(); x2 = p.x(); y2 = p.y(); }
	}
	for (int i = 0; i < l2.pointset.size(); i++) 
	{
		Point_2 p = Point_2(l2.pointset[i].x(), l2.pointset[i].y());
		m.pointset.push_back(l2.pointset[i]);
		Point_2 pp = m.line2.projection(p);
		dis += CGAL::sqrt((p-pp).squared_length());
		if (xmin > pp.x() || xmin == pp.x() && y1 > pp.y()) { xmin = pp.x(); x1 = p.x(); y1 = p.y(); }
		if (xmax < pp.x() || xmax == pp.x() && y2 < pp.y()) { xmax = pp.x(); x2 = p.x(); y2 = p.y(); }
	}
	double d1 = l1.line2.perpendicular(Point_2(0, 0)).to_vector() * l1.normal;
	double d2 = l1.line2.perpendicular(Point_2(0, 0)).to_vector() * l2.normal;
	if(d1*d2 >0 || d1==0 && d2>0 || d2==0 && d1>0) // the normal of the supporting point is in the same direction
	{
	//if(true) {
		if(m.pointset.size() == 0){
			LOG(INFO) << "empty pointset.";
			return false;
		}
		m.distance = dis / m.pointset.size();
		Point_2 s = m.line2.projection(Point_2(x1, y1));
		Point_2 t = m.line2.projection(Point_2(x2, y2));
		vector<Point_2> end;
		double mid_length = sqrt(point_point_squared_distance(s.x(),s.y(),t.x(),t.y()));
		Point_2 midpoint = CGAL::midpoint(s,t);
		Vector_2 dir = m.line2.to_vector();
		if(dir.squared_length() < 1e-8)
			LOG(INFO) << "Too small vector.";
		dir = dir/sqrt(dir.squared_length());
		m.s = midpoint + (extend*mid_length) * dir;
		m.t = midpoint - (extend*mid_length) * dir;
		m.height0 = max(l1.height0, l2.height0);
		m.height1 = min(l1.height1, l2.height1);
		// oritation
		m.normal = l1.normal;
		lines.push_back(m);
		vector<Line>::iterator pos = find(lines.begin(), lines.end(), l1);
		if (pos != lines.end()) {
			lines.erase(pos);
		}
		else {
			LOG(INFO) << "Fatal error: line doesn't exist" ;
		}
		pos = find(lines.begin(), lines.end(), l2);
		if (pos != lines.end()) {
			lines.erase(pos);
		}
		else {
			LOG(INFO) << "Fatal error: line doesn't exist" ;
		}
		return true;
	}
	else
		return false;
}

/*!
 * \description: The number of points whose distance to the line is within the threshold.
 * \param pointset pointset
 * \param l line
 * \param dist distance threshold
 * \return number of points that meet the conditions
 */
int points_on_line(vector<Point_3>& pointset, Line& l, float dist)
{
	int num = 0;
	for(int i = 0; i < pointset.size(); i++)
	{
		Point_2 p = Point_2(pointset[i].x(), pointset[i].y());
		float sdist = (l.line2.projection(p) - p).squared_length();
		if (sqrt(sdist) <= dist)
		num++;	
	}
	return num;
}

/*!
 * \description: Refine lines according to the above conditions:
	1.have enough common points or have enough intersection length
	2.have sigma-small angle and beta-close distance
	3.the direction of the normal is close
 * \param lines line set
 * \param angle_thred angle threshold at merge
 * \param num_threshold number threshold of common points at merge
 * \param extend segment extension scale
 * \return
 */
void refine_lines(vector<Line>& lines, float angle_thred, float num_threshold, float extend)
{
	float angle_theta = static_cast<float>(pi * angle_thred / 180.0f); // in radian
	float dist_theshold = 0; // average distance
	int num0 = lines.size();
	for(int i = 0; i < num0; i++)
		dist_theshold += lines[i].distance;
	dist_theshold /= num0; // ave distance between points and lines
	LOG(INFO) << "Ave distance between points and lines: " << dist_theshold;
	// extend segments
	for (int i = 0; i < lines.size(); i++)
	{
		vector<Point_2> end;
		double mid_length = sqrt(point_point_squared_distance(lines[i].s.x(),lines[i].s.y(), lines[i].t.x(),lines[i].t.y()));
		Point_2 midpoint = CGAL::midpoint(lines[i].s,lines[i].t);
		Vector_2 dir = lines[i].line2.to_vector();
		if(dir.squared_length() < 1e-8)
			LOG(INFO) << "Too small vector.";
		dir = dir/sqrt(dir.squared_length());
		lines[i].s = midpoint + (extend*mid_length) * dir;
		lines[i].t = midpoint - (extend*mid_length) * dir;
	}
	int con1 = 0;
	int con2 = 0;
	bool merged = false;
	do
	{
		merged = false;
		sort(lines.begin(), lines.end()); // decreasing order
		for(int i = 0; i < lines.size(); i++)
		{
			Line l1 = lines[i];
			auto v1 = l1.line2.to_vector();
			float len1 = CGAL::approximate_sqrt(v1.squared_length());
			float num1 = l1.pointset.size() / num_threshold; 
			for(int j = i+1; j < lines.size(); j++)
			{
				Line l2 = lines[j];
				auto v2 = l2.line2.to_vector();
				float len2 = CGAL::approximate_sqrt(v2.squared_length());
				float num2 = l2.pointset.size() / num_threshold;
				if(len1 < 1e-4 || len2 < 1e-4)
					LOG(INFO) << "Too small vector.";
				float cos2 = abs((v1*v2)/(len1*len2));
				if (cos2 > cos(angle_theta) && sqrt(robust_squared_distance(l1.s, l1.t, l2.s, l2.t)) < dist_theshold*2)
				{ // small angle and close distance
					con1++;
					// compute overlaps
					int num1on2 = points_on_line(l1.pointset, l2, l2.distance*2); 
					int num2on1 = points_on_line(l2.pointset, l1, l1.distance*2);
					if (num1on2 > num1 && num2on1 > num2) 
					{ // enough common points
						con2++;
						merged = merge(l1, l2, lines, extend);
						if (merged == true)
							break;
					}
				}
			}
			if (merged)
				break;
		}
	} while (merged);
	LOG(INFO) << con1 << " times satisify close distance and angle. " << con2 << " times satisify enough points.";

	sort(lines.begin(), lines.end());
	for (int i = 0; i < lines.size(); i++) // update ID number
		lines[i].ID = i;
	LOG(INFO) << num0 - lines.size() << " lines are merged.";
}

/*!
 * \description: Clip segments according to bbox.
 * \param bsegments bbox
 * \param lines line set
 * \param segments initial segments
 * \param wdir file save/load directory
 * \return 
 */
bool initial_seg(vector<Segment>& bsegments, vector<Line>& lines, vector<Segment>& segments, string wdir)
{ 
    // initial segments
    for (int i = 0; i < lines.size(); i++)
	{
		vector<Point_2> end;
		Segment_2 ex_seg(lines[i].s, lines[i].t); 
		for (int  j = 0; j < bsegments.size(); j++)
		{
			auto result = CGAL::intersection(bsegments[j].segment2, ex_seg);    
			if (result)
				if(const Point_2* pt = boost::get<Point_2>(&*result))
				{		
					vector<Point_2>::iterator pos = find(end.begin(), end.end(), *pt);
					if (pos == end.end()) // no exist
						end.push_back(*pt);
				}
		}
		int flag = 0;
		for(auto& p: end)
		{
			if(!point_2_range(p))
			{
				LOG(INFO) << "invalid intersection line: "<< lines[i].s << " " << lines[i].t;
				flag = 1;
				break;
			}
		}
		if(flag == 1)
		{
			segments.push_back(Segment(ex_seg,i,1,1));
			continue;
		}
		if (end.size() > 2)
		{
			LOG(INFO) << "Initial segments fails!" ;
			return false;
		}
		else if ( end.size() == 0 )
			segments.push_back(Segment(ex_seg,i,0,0)); // i -  line id number 
		else if ( end.size() == 2 )
			segments.push_back(Segment(Segment_2(end[0],end[1]),i,1,1)); // i -  line id number
		else
		{ // one intersection point, the other is an end point.
			double xmin = bsegments[0].segment2.source().x();
			double ymin = bsegments[0].segment2.source().y();
			double xmax = bsegments[3].segment2.source().x();
			double ymax = bsegments[3].segment2.source().y();
			if(ex_seg.source().x() >= xmin && ex_seg.source().x() <= xmax && ex_seg.source().y() >= ymin && ex_seg.source().y() <= ymax)
				segments.push_back(Segment(Segment_2(end[0],ex_seg.source()),i,1,0)); 
			else
				segments.push_back(Segment(Segment_2(end[0],ex_seg.target()),i,1,0)); 
		}
    }
    // initial segments save
    ofstream file(wdir + "segments_initial", ios::out);
	for(int i = 0; i < segments.size(); i++)
	{
	  file << setiosflags(ios::fixed) << setprecision(8) << segments[i].segment2.source().x() << " " << segments[i].segment2.target().x() << " " << segments[i].segment2.source().y() << " "<< segments[i].segment2.target().y() << endl;
	}
	file.close();
    return true;
}

/*!
 * \description: Generate hypothesis segments
 * \param segments hypothesis segments
 * \param points_idx point index in pointcloud
 * \param points intersection between segments
 * \param initial_seg initial segments
 * \param bsegments bbox  
 * \return {*}
 */
bool hypothesis_seg(vector<Segment>& segments, unordered_map<Point_2, int>& points_idx, vector<IntersectPoint>& points, vector<Segment>& initial_seg, vector<Segment>& bsegments)
{
   int initial_num = initial_seg.size();
   if (initial_num == 0) 
   {
		LOG(INFO) << "Initial segments no exist!";
        return false;
   } 
   int *vertag = new int[initial_num](); // end points intersection or not: 0-none 1-s 2-t 3-both
   for (int i = 0; i < initial_num; i++)
   {
		for (int j = i+1; j < initial_num; j++)
		{
			auto result = CGAL::intersection(initial_seg[j].segment2, initial_seg[i].segment2);
			Point_2 pt;
			if (result && CGAL::assign(pt, result))
			{
				IntersectPoint ip;
				ip.p = pt;
				ip.IDl.push_back(i);
				ip.IDl.push_back(j);
				points.push_back(ip);	
				if (bsegments[0].segment2.has_on(pt) || bsegments[1].segment2.has_on(pt) || bsegments[2].segment2.has_on(pt) || bsegments[3].segment2.has_on(pt))
				{
					if(initial_seg[i].segment2.source() != pt && initial_seg[i].segment2.target() != pt && initial_seg[j].segment2.source() != pt && initial_seg[j].segment2.target() != pt)
					{
						LOG(INFO) << "Error: Should on the segment!";
						return false;
					}
					if (initial_seg[i].segment2.source() == pt)
					{
						if(vertag[i]==0) vertag[i]=1;
						else if (vertag[i]!=1) vertag[i]=3;
					}
					else{
						if(vertag[i]==0) vertag[i]=2;
						else if (vertag[i]!=2) vertag[i]=3;			 
					}
					if (initial_seg[j].segment2.source() == pt)
					{
						if(vertag[j]==0) vertag[j]=1;
						else if (vertag[j]!=1) vertag[j]=3;
					}
					else{
						if(vertag[j]==0) vertag[j]=2;
						else if (vertag[j]!=2) vertag[j]=3;			 
					}
				}
			}		
		}	
    }
    for (int i = 0; i < initial_num; i++)
	{ // add end point
		if (vertag[i]==0)
		{ // add s, t
			IntersectPoint ip, ip2;
			ip.p = initial_seg[i].segment2.source();
			ip.IDl.push_back(i);
			if (initial_seg[i].flags ==1) ip.flag = 1;
			points.push_back(ip);

			ip2.p = initial_seg[i].segment2.target();
			ip2.IDl.push_back(i);
			if (initial_seg[i].flagt ==1) ip2.flag = 1;
			points.push_back(ip2);
		}
		else if (vertag[i]==1)
		{ // add t
			IntersectPoint ip;
			ip.p = initial_seg[i].segment2.target();
			ip.IDl.push_back(i);
			if (initial_seg[i].flagt ==1) ip.flag = 1;
			points.push_back(ip);
		}
		else if (vertag[i]==2)
		{ // add s
			IntersectPoint ip;
			ip.p = initial_seg[i].segment2.source();
			ip.IDl.push_back(i);
			if (initial_seg[i].flags ==1) ip.flag = 1;
			points.push_back(ip);
		}	
    }
    delete [] vertag;
    // points stores all intersection points and line ID
    
    // sort points by x,y 
    sort(points.begin(), points.end(),
    [](const IntersectPoint &a, const IntersectPoint &b) {
      return a.p.x() != b.p.x()?a.p.x() > b.p.x():a.p.y() > b.p.y();
    });

    Point_2* pp = new Point_2[initial_num];
    int* tag = new int[initial_num]();
    int* idxl = new int[initial_num]();

    for(int i = 0; i < points.size(); i++)
	{
		points_idx[points[i].p] = i;
		int nums = points[i].IDl.size();
		if (nums == 0) 
		{ 
			LOG(INFO) <<"Error: should have lines!";
			return false;
		}
		for(int j = 0; j < nums; j++)
		{
			int idx = points[i].IDl[j];
			if (tag[idx] == 0) 
			{
				pp[idx] = points[i].p;
				tag[idx] = 1;
				idxl[idx] = i;
			}
			else
			{
				segments.push_back(Segment(Segment_2(pp[idx],points[i].p),idx));
				pp[idx] = points[i].p;
				points[i].IDs.push_back(segments.size()-1);
				points[idxl[idx]].IDs.push_back(segments.size()-1);
				idxl[idx] = i;
			}
		}
    }
    delete[] tag;
    delete[] pp;
	delete[] idxl;
 
    LOG(INFO) << "Generate " << segments.size() << " hypothesis segments, " <<  points.size() << " intersecting points.";
    return true;
}

// main function
bool generation(vector<Line>& lines, vector<Segment>& bsegments, vector<Segment>& segments, unordered_map<Point_2, int>& points_idx, vector<IntersectPoint>& points, float angle_thred, float num_threshold, float extend, string wdir)
{
	// timer
	clock_t start,end;
	start = clock();

	refine_lines(lines, angle_thred, num_threshold, extend);
	vector<Segment> initial_segs;
	bool j1 = initial_seg(bsegments, lines, initial_segs, wdir);
	bool j2 = hypothesis_seg(segments, points_idx, points, initial_segs, bsegments);

	end = clock();
	cout << "Hypothesis generation time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" << endl;
	return j1 && j2;
}


