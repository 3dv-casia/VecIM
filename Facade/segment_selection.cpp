#include "base.h"
#include <glog/logging.h>
#include "method_common.h"
#include "../math/linear_program.h"
#include "../math/linear_program_solver.h"
#include <malloc.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>

typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CGAL::Partition_traits_2<K>                         Traits;
typedef Traits::Point_2                                     Pointp2;
typedef Traits::Polygon_2                                   Polygon2;
typedef std::list<Polygon2>     Polygon_list;

void optimize(vector<Segment_2>& selected, unordered_map<Point_2, int>& points_idx, vector<Point_2>& heights, vector<Vector_2>& oritations, vector<Line>& lines, vector<Segment>& segments, vector<IntersectPoint>& points, vector<Segment>& bsegments, double lambda_data_fitting, double lambda_model_coverage, double lambda_model_complexity, double alpha, string wdir)
{
	// timer
	clock_t start, end;
	start = clock();

	// linear solve
	LinearProgramSolver::SolverName solver_name = LinearProgramSolver::SCIP; // scip 
	LinearProgram program_;

	double total_points = 0;
	double dist_threshold = 0; // max distance between points and lines
	vector<Point_3> pointsets;
	for (int i = 0; i < lines.size(); i++)
	{
		total_points += lines[i].pointset.size();
		dist_threshold = max(lines[i].distance, dist_threshold);
		for (int j = 0; j < lines[i].pointset.size(); j++)
			pointsets.push_back(lines[i].pointset[j]);
	}
	LOG(INFO) << "Max distance between points and lines after merging: " << dist_threshold;

	// compute point density 
	/* nerghbor size : 20 */
	float den = density(20, pointsets);
	LOG(INFO) << "Points density: " << den;
	float radius = alpha*den; // for point coverage
	vector<Point_3>().swap(pointsets); // clear vector

	// binary variables:
	// x[0] ... x[num_segments - 1] : binary labels of all the input segments
	// x[num_segments] ... x[num_segments + num_points] : binary labels of all the intersecting points (remain or not)
	// x[num_segments + num_points] ... x[num_segments + num_points + num_points] : binary labels of corner points (sharp point or not)
	int num_segments = segments.size();
	int num_points = points.size(); 
	int num_sharp_points = points.size();
	int total_variables = num_segments + num_points + num_sharp_points;
   
	// better scale
	double bbox_length = 2*(sqrt(bsegments[0].segment2.squared_length()) + sqrt(bsegments[1].segment2.squared_length())); 
	double coeff_data_fitting = lambda_data_fitting;
	double coeff_coverage = total_points * lambda_model_coverage / bbox_length ;
	double coeff_complexity = total_points * lambda_model_complexity / double(points.size());
	LOG(INFO) << "bbox_length: " << bbox_length;
	LOG(INFO) << "total_points: " << total_points;
	LOG(INFO) << "coeff_data_fitting: " << lambda_data_fitting;
	LOG(INFO) << "coeff_coverage: " << coeff_coverage;
	LOG(INFO) << "coeff_complexity: " << coeff_complexity;

	// point supporting
	vector<double> support_num(segments.size(), 0);
	vector<vector<Point_2>> support_pp(segments.size()); // 2D supporting points
	E2I to_inexact;
	I2E to_exact;
	for(int i = 0; i < segments.size(); i++)
	{
		int id = segments[i].ID; // line index
		int n = 0;
		for (int j = 0; j < lines[id].pointset.size(); j++)
		{
			K_epec::Point_2 p = K_epec::Point_2(lines[id].pointset[j].x(), lines[id].pointset[j].y());
			K_epec::Line_2 l(to_exact(segments[i].segment2.source()), to_exact(segments[i].segment2.target()));
			K_epec::Point_2 pp = l.projection(p);
			double minx = segments[i].segment2.min().x();
			double maxx = segments[i].segment2.max().x();
			double miny = segments[i].segment2.min().y();
			double maxy = segments[i].segment2.max().y();
			K_epec::Point_2 pmin = to_exact(segments[i].segment2.min());
			K_epec::Point_2 pmax = to_exact(segments[i].segment2.max());
			if(pmin <= pp && pp <= pmax)
			{
				n++;
				support_num[i] += 1;
				support_pp[i].push_back(to_inexact(pp));
			}
		}
	}    
	LOG(INFO) << "Point supporting down.";

	//test
	// record support point
   ofstream file_sp(wdir+"ppoints");
   for(auto sp:support_pp){
       if(sp.size()==0)
     continue;
     for(auto p:sp)
     file_sp << p <<"\n";
   }
   file_sp.close();

	// point coverage
	vector<double> uncovered_length(segments.size(), 0);
	for(int i = 0; i < segments.size(); i++)
	{
		double covered_length = 0;
		double segment_length = sqrt(segments[i].segment2.squared_length());
		if(support_pp[i].size()==0)
		{
			uncovered_length[i] = segment_length;
			continue;
		} 
		sort(support_pp[i].begin(), support_pp[i].end());
		for(int j = 1; j < support_pp[i].size(); j++)
		{
			if(sqrt((support_pp[i][j]-support_pp[i][j-1]).squared_length()) < radius)
				covered_length += sqrt((support_pp[i][j]-support_pp[i][j-1]).squared_length());
		}
		uncovered_length[i] = (segment_length - covered_length)  < 0 ? 0 : (segment_length - covered_length);
	}
	LOG(INFO) << "Point coverage down";

	// linear
	program_.clear();
	LinearObjective* objective = program_.create_objective(LinearObjective::MINIMIZE);

	const vector<Variable*>& variables = program_.create_n_variables(total_variables);
	for (size_t i = 0; i < total_variables; ++i) 
	{
		Variable* v = variables[i];
		v->set_variable_type(Variable::BINARY);
	}
	LOG(INFO) << "#Total variables: " << program_.num_variables() ;

	for(int i = 0; i < num_segments; i++)
	{
		if(i > total_variables)
		{
			LOG(INFO) << "Error: variables undefine." ;
			return ;
		}
		objective->add_coefficient(i, -coeff_data_fitting * support_num[i]); // accumulate data fitting term
		objective->add_coefficient(i, coeff_coverage * uncovered_length[i]); // accumulate model coverage term
	} 
	for(int i = num_segments + num_points; i < total_variables; i++)
	{
		if(i > total_variables)
		{
			LOG(INFO) << "Error: variables undefine.";
			return;
		}
		objective->add_coefficient(i, coeff_complexity); // accumulate model complexity term
	}  

	// 1.Add constraints: the number of segments associated with a point must be either 2,3,4 or 0
	// -(e1+e2+e3+e4) >= -4xi && e1+e2+e3+e4 >= 2xi (0,2,3,4)
	// (e1+e2+e3+e4) = xi (0,1)
	// (e1+e2+e3+e4) = 2xi (0,2)
	for(int i = 0; i < points.size(); i++)
	{
		LinearConstraint* c1 = program_.create_constraint();
		LinearConstraint* c2 = program_.create_constraint();
		for(int j = 0; j < points[i].IDs.size(); j++){
			if(points[i].IDs[j] > total_variables){
				LOG(INFO) << "Error: variables undefine.";
				return;
			}
			c1->add_coefficient(points[i].IDs[j], 1.0);
			c2->add_coefficient(points[i].IDs[j], -1.0);
		}
		if(num_segments+i > total_variables){
			LOG(INFO) << "Error: variables undefine.";
			return;
		}
		c1->add_coefficient(num_segments+i, -2.0);
		c1->set_bound(LinearConstraint::LOWER, 0.0);
		c2->add_coefficient(num_segments+i, 4.0);
		c2->set_bound(LinearConstraint::LOWER, 0.0);
	}
	/*
	// boundary point must be either 1 or 0, non-boundary point must be either 2 or 0
	for(int i = 0; i < points.size(); i++)
	{
		if (points[i].flag == 1){
			if(points[i].IDs.size()!=1){
				cout << "Error: boundary points!" << endl;
				return;
			}	
			LinearConstraint* c = program_.create_constraint(LinearConstraint::FIXED,0.0,0.0);
			for(int j = 0; j<points[i].IDs.size(); j++){
				if(points[i].IDs[j] > total_variables){
					cout << "Error: variables undefine." << endl;
					return ;
				}
				c->add_coefficient(points[i].IDs[j], 1.0);
			}
			c->add_coefficient(num_segments+i, -1.0);
			continue;
		}

		LinearConstraint* c1 = program_.create_constraint(LinearConstraint::FIXED,0.0,0.0);	
		for(int j = 0; j<points[i].IDs.size(); j++){
			if(points[i].IDs[j] > total_variables){
				LOG(INFO) <<  "Error: variables undefine." ;
				return;
			}
			c1->add_coefficient(points[i].IDs[j], 1.0);	  
		}
		if(num_segments+i > total_variables){
			LOG(INFO) << "Error: variables undefine." ;
			return;
		}
		c1->add_coefficient(num_segments+i, -2.0);
	}
	*/

	// 2.Add constraints: for the sharp points. 
	double M = 1.0;
	for (int i = 0; i < points.size(); i++)
	{
		// if a point is sharp, the point must be selected first:
		// X[var_point_usage_idx] >= X[var_point_sharp_idx]
		LinearConstraint* c = program_.create_constraint();
		int var_point_usage_idx = num_segments+i;
		int var_point_sharp_idx = num_segments+num_points+i;
		if(var_point_usage_idx > total_variables || var_point_sharp_idx > total_variables){
			LOG(INFO) << "Error: variables undefine.";
			return ;
		}

		c->add_coefficient(var_point_usage_idx, 1.0);
		c->add_coefficient(var_point_sharp_idx, -1.0);
		c->set_bound(LinearConstraint::LOWER, 0.0);

		for (int j = 0; j<points[i].IDs.size(); j++){
			int s1 = points[i].IDs[j];
			Line_2 l1 = lines[segments[s1].ID].line2;
			for (int k = j+1; k<points[i].IDs.size(); k++){
				int s2 = points[i].IDs[k];
				Line_2 l2 = lines[segments[s2].ID].line2;
				if(l1!=l2 && l1!=l2.opposite()){ // non-colinear
					// the constraint is:
					//X[var_point_sharp_idx] + M * (3 - (X[s1] + X[s2] + X[var_point_usage_idx])) >= 1
					c = program_.create_constraint();
					if(s1 > total_variables || s2 > total_variables){
						LOG(INFO) << "Error: variables undefine.";
						return ;
					}
					c->add_coefficient(var_point_sharp_idx, 1.0);
					c->add_coefficient(s1, -M);
					c->add_coefficient(s2, -M);
					c->add_coefficient(var_point_usage_idx, -M);
					c->set_bound(LinearConstraint::LOWER, 1.0 - 3.0 * M);
				}
			}
		} 
	}

	LOG(INFO) << "#Total constraints: " << program_.constraints().size();

	// optimize model
	LOG(INFO) << "Solving the binary program. Please wait..." ;

	LinearProgramSolver solver;
	if (solver.solve(&program_, solver_name)) 
	{
		LOG(INFO) << "Solving the binary program done. ";
		vector<int> index;
		vector<Point_2> l;
		vector<Point_2> r;
		vector<Segment> segs;
		// mark result
		const vector<double>& X = solver.solution();
		for(int i = 0; i < segments.size(); i++)
		{
			if (static_cast<int>(std::round(X[i])) == 0) 
				continue;
			else
			{
				Point_2 s,t;
				s = segments[i].segment2.source();
				t = segments[i].segment2.target();
				if(s.x() > t.x() || s.x() == t.x() && s.y() > t.y())
					segs.push_back(Segment(Segment_2(t,s),segments[i].ID));
				else 
					segs.push_back(Segment(segments[i].segment2,segments[i].ID));				
			}
		}
		LOG(INFO)<< segs.size() << " segments are selected." ;

		// merge collinear segs
		sort(segs.begin(),segs.end());
		for(int i = 0; i < segs.size(); i++)
		{
			Point_2 ss = segs[i].segment2.source();
			Point_2 tt = segs[i].segment2.target();  
			vector<Point_2>::iterator re = find(r.begin(), r.end(), ss);
			if(re != r.end())
			{
				if(!points_idx.count(ss))
				{
					LOG(INFO) << "Intersect point is not found. " ;
					return;
				}
				int idx = points_idx[ss];
				int sum = 0;
				for(auto s: points[idx].IDs)
					if(static_cast<int>(std::round(X[s])) == 1) sum++;
				if(sum == 3 || sum == 4)
				{
					l.push_back(ss);
					r.push_back(tt);
					index.push_back(segs[i].ID);
					continue;			
				}
				vector<Point_2>::iterator re2 = find(re+1,r.end(),ss);
				if(segs[i].ID == index[distance(r.begin(),re)])
					*re = tt;
				else if (re2 != r.end() && segs[i].ID == index[distance(r.begin(),re2)])
					*re2 = tt;
				else
				{
					l.push_back(ss);
					r.push_back(tt);
					index.push_back(segs[i].ID);
				}
			}
			else
			{
				l.push_back(ss);
				r.push_back(tt);
				index.push_back(segs[i].ID);
			}
		}
		
		// record corresponding 3D points index for each segments
		// formal: segments.size() num_points points
		ofstream file(wdir + "3dpointsindex", ios::out);
		Pwn_vector facade_points;
		file << l.size() << " ";
		for(int i=0; i<l.size(); i++)
		{
			selected.push_back(Segment_2(l[i],r[i]));
			heights.push_back(Point_2(lines[index[i]].height0, lines[index[i]].height1)); // highest and lowest
			file << lines[index[i]].pointset.size() << " ";
			file << lines[index[i]].normal.x() << " " << lines[index[i]].normal.y() <<  " ";
			for(auto &p: lines[index[i]].pointset)
			{
				file << p << " ";
				facade_points.push_back(Point_with_normal({p,Vector_3(0,1,0)}));
			}
			oritations.push_back(lines[index[i]].normal);
		}
		file.close();
		PLYPointSave(wdir + "facade_points_selected.ply", facade_points);
	}
	else
	{ 
		LOG(INFO) << "solving the binary program failed. " ;
	}
	end = clock();
	LOG(INFO) << " Optimization Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" ;
}	

void optimize2(vector<Facade>& potential_facades, vector<IntersectPoint>& points, vector<int>& suppoints_num, vector<double>& area_ratio, double lambda_data_fitting, double lambda_model_coverage, double lambda_model_complexity, string wdir)
{

	if(potential_facades.size() != suppoints_num.size())
	{
		LOG(INFO) << "Variables Error. ";
		return;
	}
	// timer
	clock_t start, end;
	start = clock();
	// linear solve
	LinearProgramSolver::SolverName solver_name = LinearProgramSolver::SCIP; // scip 
	LinearProgram program_;

	double total_points = 0;
	for (int i = 0; i< suppoints_num.size(); i++)
		total_points += suppoints_num[i];

	// binary variables:
	// x[0] ... x[num_facades - 1] : binary labels of all the potential facades
	// x[num_facades] ... x[num_facades + num_points] : binary labels of all the intersecting points (remain or not)
	// x[num_facades + num_points] ... x[num_facades + num_points + num_points] : binary labels of corner points (sharp point of not)
	int num_facades = potential_facades.size();
	int num_points = points.size(); 
	int num_sharp_points = points.size();
	int total_variables = num_facades + num_points + num_sharp_points;

	// better scale
	double coeff_data_fitting = lambda_data_fitting;
	double coeff_coverage = total_points * lambda_model_coverage;
	double coeff_complexity = total_points * lambda_model_complexity / double(points.size());
	LOG(INFO) << "coeff_data_fitting: " << lambda_data_fitting;
	LOG(INFO) << "coeff_coverage: " << coeff_coverage;
	LOG(INFO) << "coeff_complexity: " << coeff_complexity;

	//test
    ifstream file(wdir + "energy");
	E2I to;
	for(int i = 0; i < suppoints_num.size(); i++){
		LOG(INFO) << potential_facades[i].seg << " " << suppoints_num[i] << " " << -coeff_data_fitting * suppoints_num[i] << " " << area_ratio[i] << " " << coeff_coverage * area_ratio[i] << "\n";
	}
	file.close();

	// linear
	program_.clear();
	LinearObjective* objective =  program_.create_objective(LinearObjective::MINIMIZE);

	const vector<Variable*>& variables = program_.create_n_variables(total_variables);
	for (size_t i = 0; i < total_variables; ++i) 
	{
		Variable* v = variables[i];
		v->set_variable_type(Variable::BINARY);
	}
	LOG(INFO) << "#Total variables: " << program_.num_variables() << ", " << num_facades << " facades, " << num_points << " points.";

	for(int i = 0; i < num_facades; i++)
	{
		if(i > total_variables)
		{
			LOG(INFO) << "Error: variables undefine." ;
			return ;
		}
		objective->add_coefficient(i, -coeff_data_fitting * suppoints_num[i]);// accumulate data fitting term
		objective->add_coefficient(i, coeff_coverage * area_ratio[i]);// accumulate model coverage term
	} 
	for(int i = num_facades + num_points; i < total_variables; i++)
	{
		if(i > total_variables)
		{
			LOG(INFO) << "Error: variables undefine.";
			return ;
		}
		objective->add_coefficient(i, coeff_complexity);// accumulate model complexity term
	}  

	// 1.Add constraints: the number of segments associated with a point must be either 2,3,4 or 0
	// -(e1+e2+e3+e4) >= -4xi && e1+e2+e3+e4 >= 2xi (0,2,3,4)
	// (e1+e2+e3+e4) = xi (0,1)
	// (e1+e2+e3+e4) = 2xi (0,2)
	for(int i = 0; i<points.size(); i++)
	{
		// boundary point dont have constraints 
		if (points[i].flag == 1)	
			continue;
		// Add constraints: points 0 2 3 4
		LinearConstraint* c1 = program_.create_constraint();
		LinearConstraint* c2 = program_.create_constraint();
		for(int j = 0; j < points[i].IDs.size(); j++)
		{
			if(points[i].IDs[j] > num_facades)
			{
				LOG(INFO) << "Error: variables undefine.";
				return;
			}
			c1->add_coefficient(points[i].IDs[j], 1.0);
			c2->add_coefficient(points[i].IDs[j], -1.0);
		}
		if(num_facades+i > total_variables)
		{
			LOG(INFO) << "Error: variables undefine.";
			return ;
		}
		c1->add_coefficient(num_facades+i, -2.0);
		c1->set_bound(LinearConstraint::LOWER, 0.0);
		c2->add_coefficient(num_facades+i, 4.0);
		c2->set_bound(LinearConstraint::LOWER, 0.0);
	}

	LOG(INFO) << "#Total constraints: " << program_.constraints().size();

	// optimize model
	LOG(INFO) << "Solving the binary program. Please wait..." ;

	LinearProgramSolver solver;
	if (solver.solve(&program_, solver_name)) 
	{
		LOG(INFO) << "Solving the binary program done. ";
		// mark result
		int nums = 0;
		const vector<double>& X = solver.solution();
		for(int i = 0; i < num_facades; i++)
		{
			if(static_cast<int>(std::round(X[i])) == 0)
				continue;
			else
			{
				potential_facades[i].flag = 1;
				nums++;
			}
		}
		LOG(INFO)<< nums << " facades are selected." ;
	}
	else
	{ 
		LOG(INFO) << "solving the binary program failed. " ;
	}
	end = clock();
	LOG(INFO) << "Optimization Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" ;
}	


bool segment_wall(vector<Segment_2> selected, const cm::Config& config, string wdir)
{

    ifstream file(wdir + "3dpointsindex", ios::in);
	int segs_num; int points_num;
    file >> segs_num;
    double x,y,z;
    if(segs_num != selected.size()) 
	{
      LOG(INFO) << "facade record error.";
      return false;
    } 
    vector<vector<Point_3>> points(segs_num);
    for(auto &p_set: points)
	{
        file >> points_num;
        file >> x >> y;    
        for(int i = 0; i < points_num; i++)
		{
            file >> x >> y >> z;
            p_set.push_back(Point_3(x,y,z));
        }
    } 
	auto step = config.get_optional<double>("door.step");
	auto thre_num = config.get_optional<double>("door.thre_num");
	auto weight = config.get_optional<double>("door.weight");
	auto ratio = config.get_optional<double>("door.ratio");
	auto h_ratio = config.get_optional<double>("door.h_ratio");
	vector<Segment_2> seg_new;
	for(int k = 0; k < selected.size(); k++)
	{
		auto seg = selected[k];
		if(seg.squared_length() < (*weight)*(*weight))
		{
		   seg_new.push_back(seg);
		   continue;
		}
		Point_2 s = seg.source() < seg.target() ? seg.source() : seg.target();
		Point_2 t = seg.source() > seg.target() ? seg.source() : seg.target();  
		int bin_num = ceil(sqrt(seg.squared_length()) / (*step));
		struct p3{
			int x=0;
			double y=-1e9,z=1e9;
		};
		vector<p3> bins(bin_num); // num, zmax, zmin
		double zmax = -1e9;
		double zmin = 1e9; 
        Vector_2 dir = t - s;
		dir = dir / sqrt(dir.squared_length());
        for(auto p: points[k])
		{
			Point_2 pp = seg.supporting_line().projection(Point_2(p.x(), p.y()));
			if(pp >=s && pp <=t)
			{
				int idx = sqrt((pp - s).squared_length()) / (*step);
				bins[idx].x = bins[idx].x + 1;
				bins[idx].y = max(bins[idx].y, p.z());
				bins[idx].z = min(bins[idx].z, p.z());
				zmax = max(zmax, bins[idx].y);
				zmin = min(zmin, bins[idx].z);
			}
		}
		int black = 0;
		double thre_zmax = zmax - (zmax - zmin) *(*h_ratio);
		double thre_zmin = zmin + (zmax - zmin) * (*h_ratio);
		vector<int> tag(bin_num, 0); // 1-black 2-top 3-bottom
		for(int i = 0; i < bin_num; i++)
		{
			if(bins[i].x < *thre_num)
			{
				black++;
				tag[i] = 1;
			}
			else if(bins[i].y <= thre_zmax)
			{
				tag[i] = 3;
			}
			else if(bins[i].z >= thre_zmin)
			{
				tag[i] = 2;
			}

		}
        if(black >= bin_num * (*ratio))
		{ // no operation
			seg_new.push_back(seg);
		    continue;
		}
		int line_bin = (*weight)/(*step);
		map<pair<int,int>,int> record; // increasing order
		vector<int> start(4, -1);
		for(int i = 0; i < bin_num; i++)
		{
           if(i == 0 && tag[i] != 0)
		   {
              start[tag[i]] = i;
		   }
		   else
		   {
			   if(tag[i] != tag[i-1])
			   {
				   if(tag[i-1] != 0 && i - start[tag[i-1]] >= line_bin){
					   record[{start[tag[i-1]], i-1}] = tag[i-1];
				   }
				   if(tag[i] != 0){
                       start[tag[i]] = i;
				   }
			   }
		   }
		}
		if(tag[bin_num-1] != 0 && bin_num - start[tag[bin_num-1]] >= line_bin)
				record[{start[tag[bin_num-1]], bin_num-1}] = tag[bin_num-1];
		// divide segment
		int prev_idx = -1;
		Point_2 prev = s;
		Point_2 end = s;
		int num = 0;
		for(auto iter: record)
		{
			num++;
			if(iter.first.first != prev_idx + 1)
			{
                  end = prev + (*step) * (iter.first.first - prev_idx - 1) * dir;
				  seg_new.push_back({prev, end});
				  if(end >= t){
					  LOG(INFO) << "Division error.";
					  LOG(INFO) << dir.squared_length() << " " << iter.first.first << " " << prev_idx << "  " << prev << " " << end;
					  LOG(INFO) << s << " " << t;
				  }
				  prev = end;
				  prev_idx = iter.first.first - 1;
			}
			if(iter.second != 1)
			{
                  end = prev + (*step) * (1 + iter.first.second - iter.first.first) * dir;
				  end = end < t ? end : t; // end
				  if(end == t && num != record.size()){
				    LOG(INFO) << "segment division error.";
					return false;
					}
				  seg_new.push_back({prev, end});
				  prev = end;
				  prev_idx = iter.first.second;
			}
			else{
				end = prev + (*step) * (1 + iter.first.second - iter.first.first) * dir;
				end = end < t ? end : t; // end
				prev = end;
				prev_idx = iter.first.second;
			}
		}
	    // end process!
		if(end != t)
		   seg_new.push_back({end,t});
	}

	ofstream file2(wdir + "segments_division", ios::out);
    for(int i = 0; i<seg_new.size(); i++)
	{
      file2 << double(seg_new[i].source().x())<<" "<< double(seg_new[i].target().x())<<" "<< double(seg_new[i].source().y())<<" "<< double(seg_new[i].target().y())<< endl;
    }
	  file2.close();
    return true;
}

void producemesh(vector<Segment_2>& selected, Mesh<K_epec::Point_3>& planes, vector<Point_2>& heights, vector<Vector_2>& oritations, double top, double bottom)
{
    K_epec::Point_3 p1, p2, p3, p4;
	for(int i = 0; i < selected.size(); i++)
	{
		if (top < bottom)
		{
			p1 = K_epec::Point_3(selected[i].source().x(),selected[i].source().y(),heights[i].y());
			p2 = K_epec::Point_3(selected[i].source().x(),selected[i].source().y(),heights[i].x());
			p3 = K_epec::Point_3(selected[i].target().x(),selected[i].target().y(),heights[i].y());
			p4 = K_epec::Point_3(selected[i].target().x(),selected[i].target().y(),heights[i].x());
		}
		else
		{
			p1 = K_epec::Point_3(selected[i].source().x(),selected[i].source().y(),bottom);
			p2 = K_epec::Point_3(selected[i].source().x(),selected[i].source().y(),top);
			p3 = K_epec::Point_3(selected[i].target().x(),selected[i].target().y(),bottom);
			p4 = K_epec::Point_3(selected[i].target().x(),selected[i].target().y(),top);
		}
		int idx = planes.vertices.size();
		planes.vertices.push_back(p1); planes.vertices.push_back(p2);
		planes.vertices.push_back(p3); planes.vertices.push_back(p4);
		planes.vertices.push_back(p2); planes.vertices.push_back(p3);
		K_epec::Vector_3 v;
		if(i < oritations.size())
			v = K_epec::Vector_3(oritations[i].x(),oritations[i].y(),0);
		else 
			v = K_epec::Vector_3(0,0,0);

		K_epec::Plane_3 plane(p1 ,p2, p3);
		if (plane.orthogonal_vector()*v > 0 )
		{
			planes.faces.push_back(K_epec::Point_3(idx, idx+1, idx+2));
			planes.faces.push_back(K_epec::Point_3(idx+5, idx+4, idx+3));
		}
		else
		{
			planes.faces.push_back(K_epec::Point_3(idx+1, idx, idx+2));
			planes.faces.push_back(K_epec::Point_3(idx+4, idx+5, idx+3));
		}
	}
}


