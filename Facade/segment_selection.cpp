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

void optimize(vector<Segment_2>& selected, unordered_map<Point_2, int>& points_idx, vector<Point_2>& heights, vector<Vector_2>& oritations, vector<Line>& lines, vector<Segment>& segments, vector<IntersectPoint>& points, vector<Segment>& bsegments, double lambda_data_fitting, double lambda_model_coverage, double lambda_model_complexity, double epsilon, double alpha, string wdir){

   // timer
   clock_t start,end;
   start = clock();

   // linear solve
   LinearProgramSolver::SolverName solver_name = LinearProgramSolver::SCIP; // scip 
   LinearProgram	program_;


   double total_points = 0;
   float dist_thethold = 0; // max distance between points and lines
   vector<Point_3> pointsets;
   for (int i = 0; i<lines.size(); i++){
	total_points += lines[i].pointset.size();
	dist_thethold = max(lines[i].distance,dist_thethold);
	for (int j = 0; j<lines[i].pointset.size(); j++)
	   pointsets.push_back(lines[i].pointset[j]);
   }
   LOG(INFO) << "Max distance between points and lines after merging: " << dist_thethold;
   
    ofstream file2(wdir + "2dpoints", ios::out);
	for(int i=0; i<pointsets.size(); i++){
	  file2 << pointsets[i].x() << " " << pointsets[i].y() << endl;
	}
	file2.close();

   // compute point density 
   /* nerghbor size : 20 */
   float den = density(20, pointsets);
   LOG(INFO) << "Points density: " << den;
   float radius = alpha*den; // for point coverage
   vector<Point_3>().swap(pointsets); // clear vector

   // binary variables:
	// x[0] ... x[num_segments - 1] : binary labels of all the input segments
	// x[num_segments] ... x[num_segments + num_points] : binary labels of all the intersecting points (remain or not)
	// x[num_segments + num_points] ... x[num_segments + num_points + num_points] : binary labels of corner points (sharp point of not)
   int num_segments = segments.size();
   int num_points = points.size(); // *****consider 4 segments
   int num_sharp_points = points.size();
   int total_variables = num_segments + num_points + num_sharp_points;
   
   // better scale
   double bbox_length =  2*(sqrt(bsegments[0].segment2.squared_length()) + sqrt(bsegments[1].segment2.squared_length())); // *1 one side,*2 two side
   double coeff_data_fitting = lambda_data_fitting;
   double coeff_coverage = total_points * lambda_model_coverage / bbox_length ;
   double coeff_complexity = total_points * lambda_model_complexity / double(points.size());
   LOG(INFO) << "bbox_length: " << bbox_length;
   LOG(INFO) << "total_points: " << total_points;
   LOG(INFO) << "coeff_data_fitting: " << lambda_data_fitting;
   LOG(INFO) << "coeff_coverage: " << coeff_coverage ;
   LOG(INFO) << "coeff_complexity: " << coeff_complexity;

   //TODO:TEST
   double lens = 0;
   for(auto s: segments){
      lens += sqrt(s.segment2.squared_length());
	  //LOG(INFO) << s.segment2.source() << " " << s.segment2.target() << " " << s.segment2.squared_length(); 
   }
    LOG(INFO) << "total_candidate_length: " << lens;

   // data fitting term: support points of segments
   //double *support_num = new double[segments.size()]();
   vector<double> support_num(segments.size(), 0);
   vector<vector<Point_2>> support_pp(segments.size());// support points projection
   // fixed: accurate computation
   	typedef CGAL::Cartesian_converter<K_epec, K> To_geom;
    To_geom to_geom;
   	typedef CGAL::Cartesian_converter<K, K_epec> To_geom2;
    To_geom2 to_geom2;
   for(int i = 0; i < segments.size(); i++){
	int id = segments[i].ID; // line
	int n = 0;
	for (int j = 0; j < lines[id].pointset.size(); j++){
	  K_epec::Point_2 p = K_epec::Point_2(lines[id].pointset[j].x(), lines[id].pointset[j].y());
	  K_epec::Line_2 l(to_geom2(segments[i].segment2.source()), to_geom2(segments[i].segment2.target()));
	  K_epec::Point_2 pp = l.projection(p);
	  double s = segments[i].segment2.source().x();
	  double t = segments[i].segment2.target().x();
	  double minx = s > t ? t:s;
	  double maxx = s < t ? t:s;// unsure 1e-4 1e-8
	  double s2 = segments[i].segment2.source().y();
	  double t2 = segments[i].segment2.target().y();
	  double miny = s2 > t2 ? t2:s2;
	  double maxy = s2 < t2 ? t2:s2;// unsure 1e-4 1e-8
	  //TODO: x,y constraints
	  bool con =  CGAL::certainly(CGAL::is_zero(l.b())) ? (pp.y()>=(miny-1e-4) && pp.y()<=(maxy+1e-4)) 
	  : ( CGAL::certainly(CGAL::is_zero(l.a())) ? (pp.x()>=(minx-1e-4) && pp.x()<=(maxx+1e-4)) :
	   (l.x_at_y(pp.y())>=(minx-1e-4) && l.x_at_y(pp.y())<=(maxx+1e-4) && 
	    l.y_at_x(pp.x())>=(miny-1e-4) && l.y_at_x(pp.x())<=(maxy+1e-4)));
	  if(con){ //&& sqrt((p-pp).squared_length()) < dist_thethold){ // support points
		n++;
		//support_num[i] += 1-(sqrt((p-pp).squared_length()))/(dist_thethold);
		support_num[i] += 1;
		support_pp[i].push_back(to_geom(pp));
	  }
	}
   }    
    LOG(INFO) << " Point support down.";
   // point coverage
   double *uncovered_length = new double[segments.size()]();
   for(int i = 0; i < segments.size(); i++){
	double covered_length = 0;
	double segment_length = sqrt(segments[i].segment2.squared_length());
	if(support_pp[i].size()==0){
	 uncovered_length[i] = segment_length;
	   //uncovered_length[i] = 1.0;
	 continue;
	} 
	sort(support_pp[i].begin(),support_pp[i].end());
	for(int j=1; j<support_pp[i].size(); j++){
	  if (sqrt((support_pp[i][j]-support_pp[i][j-1]).squared_length()) < radius)
		covered_length += sqrt((support_pp[i][j]-support_pp[i][j-1]).squared_length());
	}
	uncovered_length[i] = (segment_length - covered_length)  < 0 ? 0 : (segment_length - covered_length);
   }
   LOG(INFO) << "Point coverage down";

   // record support point
   ofstream file_sp(wdir+"ppoints");
   for(auto sp:support_pp){
       if(sp.size()==0)
	   continue;
	   for(auto p:sp)
	   file_sp << p <<"\n";
   }
   file_sp.close();

   // linear
   program_.clear();
   LinearObjective* objective =      program_.create_objective(LinearObjective::MINIMIZE);

 // set
   const vector<Variable*>& variables =    program_.create_n_variables(total_variables);
   for (size_t i = 0; i < total_variables; ++i) {
		Variable* v = variables[i];
		v->set_variable_type(Variable::BINARY);
   }
  LOG(INFO) << "#Total variables: " << program_.num_variables() ;

   for(int i = 0; i < num_segments; i++){
	if(i > total_variables){
	  LOG(INFO) << "Error: variables undefine." ;
	  return ;
	}
	objective->add_coefficient(i, -coeff_data_fitting * support_num[i]);// accumulate data fitting term
	objective->add_coefficient(i, coeff_coverage * uncovered_length[i]);// accumulate model coverage term
   } 
   for(int i = num_segments + num_points; i < total_variables; i++){
	if(i > total_variables){
	  LOG(INFO) << "Error: variables undefine.";
	  return ;
	}
	objective->add_coefficient(i, coeff_complexity);// accumulate model complexity term

   }  

   // 1.Add constraints: the number of segments associated with a point must be either 2,3,4 or 0
   // -(e1+e2+e3+e4) >= -4xi && e1+e2+e3+e4 >= 2xi (0,2,3,4)
   // (e1+e2+e3+e4) = xi (0,1)
   // (e1+e2+e3+e4) = 2xi (0,2)
   for(int i = 0; i<points.size(); i++){

	// boundary point must be either 1 or 0
	/*
		// boundary points 0 1
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
	*/
	/*
		// points 0 2
		LinearConstraint* c1 = program_.create_constraint(LinearConstraint::FIXED,0.0,0.0);	
		for(int j = 0; j<points[i].IDs.size(); j++){
		if(points[i].IDs[j] > total_variables){
		LOG(INFO) <<  "Error: variables undefine." ;
		return ;
		}
		c1->add_coefficient(points[i].IDs[j], 1.0);	  
		}
		if(num_segments+i > total_variables){
		LOG(INFO) << "Error: variables undefine." ;
		return ;
		}
		c1->add_coefficient(num_segments+i, -2.0);
	*/

		// points 0 2 3 4
		LinearConstraint* c1 = program_.create_constraint();
		LinearConstraint* c2 = program_.create_constraint();
		for(int j = 0; j<points[i].IDs.size(); j++){
		if(points[i].IDs[j] > total_variables){
		LOG(INFO) << "Error: variables undefine.";
		return ;
		}
		c1->add_coefficient(points[i].IDs[j], 1.0);
		c2->add_coefficient(points[i].IDs[j], -1.0);
		}
		if(num_segments+i > total_variables){
		LOG(INFO) << "Error: variables undefine.";
		return ;
		}
		c1->add_coefficient(num_segments+i, -2.0);
		c1->set_bound(LinearConstraint::LOWER, 0.0);

		c2->add_coefficient(num_segments+i, 4.0);
		c2->set_bound(LinearConstraint::LOWER, 0.0);

   }


  // 2.Add constraints: for the sharp points. 
   double M = 1.0;
   for (int i = 0; i<points.size(); i++){
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
	     if( l1!=l2 && l1!=l2.opposite() ){ // non-colinear
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
   if (solver.solve(&program_, solver_name)) {
	LOG(INFO) << "Solving the binary program done. ";
	vector<int> index;
	vector<Point_2> l;
	vector<Point_2> r;
	vector<Segment> segs;
	//TODO:test
	ofstream test(wdir + "results");
	// mark result
	const vector<double>& X = solver.solution();
	for(int i = 0; i<segments.size(); i++){
		test << segments[i].segment2.source() << " " << segments[i].segment2.target() << " " << coeff_coverage * uncovered_length[i] << " " << uncovered_length[i] << " " << -coeff_data_fitting * support_num[i] << " " << support_num[i] << "\n";
	   if (static_cast<int>(std::round(X[i])) == 0) 
		continue;
	   else{
		Point_2 s,t;
		s = segments[i].segment2.source();
		t = segments[i].segment2.target();
		if(s.x() > t.x() || s.x() == t.x() && s.y() > t.y())
		  segs.push_back(Segment(Segment_2(t,s),segments[i].ID));
		else 
		  segs.push_back(Segment(segments[i].segment2,segments[i].ID));
		//index.push_back(i);
	    }
    }
	LOG(INFO)<< segs.size() << " segments are selected." ;

	//test
	ofstream fp(wdir + "points"); 
	for(int i = segments.size(); i<segments.size()+points.size(); i++){
            if (static_cast<int>(std::round(X[i])) == 0) continue;
			fp << i-segments.size() << " " <<points[i-segments.size()].p << "\n";
			//LOG(INFO) << points[i-segments.size()].p;
			//for(auto s: points[i-segments.size()].IDs)
			 //LOG(INFO) << s << " " <<static_cast<int>(std::round(X[s])) << " " <<  segments[s].segment2;
	} 
    fp.close();

	// merge collinear segs
	sort(segs.begin(),segs.end());
	for(int i = 0; i < segs.size(); i++){
	  Point_2 ss = segs[i].segment2.source();
	  Point_2 tt = segs[i].segment2.target();  
	  vector<Point_2>::iterator re = find(r.begin(),r.end(),ss);
	  if(re != r.end()){
		//unsure: reserve corner point
		if(!points_idx.count(ss)){
		     LOG(INFO) << "Intersect point is not found. " ;
			 return;
		}
		int idx = points_idx[ss];
		int sum = 0;
		for(auto s: points[idx].IDs)
		   if(static_cast<int>(std::round(X[s])) == 1) sum++;
		//LOG(INFO) << idx << " " << ss <<" " << tt << " " <<points[idx].p << " " <<  sum;   
		if(sum == 3 || sum == 4){
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
	    else{
	     l.push_back(ss);
	     r.push_back(tt);
	     index.push_back(segs[i].ID);
	    }
	  }
	  else{
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
	for(int i=0; i<l.size(); i++){
	  selected.push_back(Segment_2(l[i],r[i]));
	  heights.push_back(Point_2(lines[index[i]].height0, lines[index[i]].height1));// highest and lowest
	  file << lines[index[i]].pointset.size() << " ";
	  file << lines[index[i]].normal.x() << " " << lines[index[i]].normal.y() <<  " ";
	  for(auto &p: lines[index[i]].pointset){
	      file << p << " ";
		  facade_points.push_back(Point_with_normal({p,Vector_3(0,1,0)}));
	  }
	  oritations.push_back(lines[index[i]].normal);
	}
	file.close();
	PLYPointSave(wdir + "facade_points_selected.ply", facade_points);
	
	// adjust face normal 
	//normal_adjust(selected, oritations);


   }
   else{ 
        LOG(INFO) << "solving the binary program failed. " ;
   }

   end = clock();
   LOG(INFO) << " Optimization Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" ;
  
}	

void optimize2(vector<Facade>& potential_facades, vector<IntersectPoint>& points, vector<int>& isolate_seg, vector<int>& suppoints_num, vector<double>& area_ratio, double lambda_data_fitting, double lambda_model_coverage, double lambda_model_complexity, string wdir){
   
   if(potential_facades.size() != suppoints_num.size()){
	   LOG(INFO) << "Variables Error. ";
	   return;
   }
   // timer
   clock_t start,end;
   start = clock();

   // linear solve
   LinearProgramSolver::SolverName solver_name = LinearProgramSolver::SCIP; // scip 
   LinearProgram program_;

  
   double total_points = 0;
   for (int i = 0; i< suppoints_num.size(); i++)
	total_points += suppoints_num[i];

/*
   // test
   for (int i = 0; i< suppoints_num.size(); i++){
	   LOG(INFO) << potential_facades[i].seg << " " << suppoints_num[i]/total_points << " " << area_ratio[i];
   }
*/

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
   LOG(INFO) << "coeff_coverage: " << coeff_coverage ;
   LOG(INFO) << "coeff_complexity: " << coeff_complexity;

   // linear
   program_.clear();
   LinearObjective* objective =  program_.create_objective(LinearObjective::MINIMIZE);

   const vector<Variable*>& variables = program_.create_n_variables(total_variables);
   for (size_t i = 0; i < total_variables; ++i) {
		Variable* v = variables[i];
		v->set_variable_type(Variable::BINARY);
   }
  LOG(INFO) << "#Total variables: " << program_.num_variables() << ", " << num_facades << " facades, " << num_points << " points.";

   for(int i = 0; i < num_facades; i++){
	if(i > total_variables){
	  LOG(INFO) << "Error: variables undefine." ;
	  return ;
	}
	objective->add_coefficient(i, -coeff_data_fitting * suppoints_num[i]);// accumulate data fitting term
	objective->add_coefficient(i, coeff_coverage * area_ratio[i]);// accumulate model coverage term
	//LOG(INFO) << -coeff_data_fitting * suppoints_num[i]  << " " << coeff_coverage * area_ratio[i]; //TEST
   } 
   for(int i = num_facades + num_points; i < total_variables; i++){
	if(i > total_variables){
	  LOG(INFO) << "Error: variables undefine.";
	  return ;
	}
	objective->add_coefficient(i, coeff_complexity);// accumulate model complexity term

   }  

   // 1.Add constraints: the number of segments associated with a point must be either 2,3,4 or 0
   // -(e1+e2+e3+e4) >= -4xi && e1+e2+e3+e4 >= 2xi (0,2,3,4)
   // (e1+e2+e3+e4) = xi (0,1)
   // (e1+e2+e3+e4) = 2xi (0,2)
   //test
   int test = 0;
   vector<int> test2;
   for(int i = 0; i<points.size(); i++){


	// boundary point dont have constraints 
	if (points[i].flag == 1)	
 	    continue;
    
	

   	// Add constraints: points 0 2 3 4
	LinearConstraint* c1 = program_.create_constraint();
	LinearConstraint* c2 = program_.create_constraint();
	//test
	if(points[i].p.x() < 282 && points[i].p.x() > 274.26 && points[i].p.y() < 271 && points[i].p.y() > 264.34){
		LOG(INFO) << "TEST4: " << points[i].p;
		test = num_facades+i;
		for(auto k:points[i].IDs)
		   LOG(INFO) << k << " " << isolate_seg[k] << " " << CGAL::to_double(potential_facades[k].seg.squared_length());
	}
	for(int j = 0; j < points[i].IDs.size(); j++){
	  if(points[i].IDs[j] > num_facades){
	  LOG(INFO) << "Error: variables undefine.";
	  return ;
	  }
    if(isolate_seg[points[i].IDs[j]] == 1 && CGAL::to_double(potential_facades[points[i].IDs[j]].seg.squared_length()) < 25e-6 || isolate_seg[points[i].IDs[j]] == -1) // the source and target of segments overlap
        continue;
    //test
	 if(points[i].p.x() < 282 && points[i].p.x() > 274.26 && points[i].p.y() < 271 && points[i].p.y() > 264.34){
	  LOG(INFO) << i << " " <<  points[i].IDs[j];
	  test2.push_back( points[i].IDs[j]);
	}
	  c1->add_coefficient(points[i].IDs[j], 1.0);
	  c2->add_coefficient(points[i].IDs[j], -1.0);
	}
	if(num_facades+i > total_variables){
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
   if (solver.solve(&program_, solver_name)) {
	LOG(INFO) << "Solving the binary program done. ";

	// mark result
    int nums = 0;
	const vector<double>& X = solver.solution();
	for(int i = 0; i < num_facades; i++){
		//LOG(INFO) << X[i]; // test
	   if (isolate_seg[i] == 1 && CGAL::to_double(potential_facades[i].seg.squared_length()) < 25e-6 || isolate_seg[i] == -1 || static_cast<int>(std::round(X[i])) == 0) 
		continue;
	   else{
		potential_facades[i].flag = 1;
		nums++;
	   }
	}
	LOG(INFO)<< nums << " facades are selected." ;


   }
   else{ 
        LOG(INFO) << "solving the binary program failed. " ;
   }

   end = clock();
   LOG(INFO) << "Optimization Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" ;
  
}	


bool segment_wall(vector<Segment_2> selected, const cm::Config& config, string wdir){

    ifstream file(wdir + "3dpointsindex", ios::in);
	int segs_num; int points_num;
    file >> segs_num;
    double x,y,z;
    if(segs_num != selected.size()) {
      LOG(INFO) << "facade record error.";
      return false;
    } 
    vector<vector<Point_3>> points(segs_num);
    for(auto &p_set: points){
        file >> points_num;
        file >> x >> y;    
        for(int i = 0; i < points_num; i++){
            file >> x >> y >> z;
            p_set.push_back(Point_3(x,y,z));
        }
    } 
    //
	 auto step = config.get_optional<double>("door.step");
	 auto thre_num = config.get_optional<double>("door.thre_num");
	 auto weight = config.get_optional<double>("door.weight");
	 auto ratio = config.get_optional<double>("door.ratio");
	 auto h_ratio = config.get_optional<double>("door.h_ratio");
	vector<Segment_2> seg_new;
	for(int k = 0; k < selected.size(); k++){
		auto seg = selected[k];
		if(seg.squared_length() < (*weight)*(*weight)){
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
        for(auto p: points[k]){
			Point_2 pp = seg.supporting_line().projection(Point_2(p.x(), p.y()));
			if(pp >=s && pp <=t){
				int idx = sqrt((pp - s).squared_length()) / (*step);
				bins[idx].x = bins[idx].x + 1;
				bins[idx].y = max(bins[idx].y, p.z());
				bins[idx].z = min(bins[idx].z, p.z());
				zmax = max(zmax, bins[idx].y);
				zmin = min(zmin, bins[idx].z);
			}
		}
		//LOG(INFO) << "Bins finished.";
		int black = 0;
		double thre_zmax = zmax - (zmax - zmin) *(*h_ratio);
		double thre_zmin = zmin + (zmax - zmin) * (*h_ratio);
		vector<int> tag(bin_num, 0); // 1-black 2-top 3-bottom
		for(int i = 0; i < bin_num; i++){
		  if(bins[i].x < *thre_num){
             black++;
			 tag[i] = 1;
		  }
		  else if(bins[i].y <= thre_zmax){
             tag[i] = 3;
		  }
		  else if(bins[i].z >= thre_zmin){
             tag[i] = 2;
		  }

		}
        if(black >= bin_num * (*ratio)){ // no operation
			seg_new.push_back(seg);
		    continue;
		}
				//test
		if(abs(s.x() - 62.3785)< 0.001 && abs(s.y() - 68.3534)< 0.001)
		   for(auto t:tag) cout << t;
		int line_bin = (*weight)/(*step);
		map<pair<int,int>,int> record; // increasing order
		vector<int> start(4, -1);
		for(int i = 0; i < bin_num; i++){
           if(i == 0 && tag[i] != 0){
              start[tag[i]] = i;
		   }
		   else{
			   if(tag[i] != tag[i-1]){
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
		for(auto iter: record){
			num++;
			if(iter.first.first != prev_idx + 1){
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
			if(iter.second != 1){
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
    for(int i = 0; i<seg_new.size(); i++){
      file2 << double(seg_new[i].source().x())<<" "<< double(seg_new[i].target().x())<<" "<< double(seg_new[i].source().y())<<" "<< double(seg_new[i].target().y())<< endl;
    }
	  file2.close();

    return true;
}

void producemesh(vector<Segment_2>& selected, Mesh& planes, vector<Point_2>& heights, vector<Vector_2>& oritations, double top, double bottom){
  Point_3 p1, p2, p3, p4;
	for(int i = 0; i < selected.size(); i++){
		if (top < bottom){
			p1 = Point_3(selected[i].source().x(),selected[i].source().y(),heights[i].y());
			p2 = Point_3(selected[i].source().x(),selected[i].source().y(),heights[i].x());
			p3 = Point_3(selected[i].target().x(),selected[i].target().y(),heights[i].y());
			p4 = Point_3(selected[i].target().x(),selected[i].target().y(),heights[i].x());
		}
		else{
			p1 = Point_3(selected[i].source().x(),selected[i].source().y(),bottom);
			p2 = Point_3(selected[i].source().x(),selected[i].source().y(),top);
			p3 = Point_3(selected[i].target().x(),selected[i].target().y(),bottom);
			p4 = Point_3(selected[i].target().x(),selected[i].target().y(),top);
		}
		/*
		int idx1,idx2,idx3,idx4;
		vector<Point_3>::iterator it = find(planes.vertices.begin(),planes.vertices.end(),p1);
		if(it!=planes.vertices.end()) idx1 = distance(planes.vertices.begin(), it);
		else{
			idx1 = planes.vertices.size();
			planes.vertices.push_back(p1);
		}

		it = find(planes.vertices.begin(),planes.vertices.end(),p2);
		if(it!=planes.vertices.end()) idx2 = distance(planes.vertices.begin(), it);
		else{
		idx2 = planes.vertices.size();
		planes.vertices.push_back(p2);
		}

		it = find(planes.vertices.begin(),planes.vertices.end(),p3);
		if(it!=planes.vertices.end()) idx3 = distance(planes.vertices.begin(), it);
		else{
		idx3 = planes.vertices.size();
		planes.vertices.push_back(p3);
		}

		it = find(planes.vertices.begin(),planes.vertices.end(),p4);
		if(it!=planes.vertices.end()) idx4 = distance(planes.vertices.begin(), it);
		else{
		idx4 = planes.vertices.size();
		planes.vertices.push_back(p4);
		}
		*/
		int idx = planes.vertices.size();
		planes.vertices.push_back(p1); planes.vertices.push_back(p2);
		planes.vertices.push_back(p3); planes.vertices.push_back(p4);
		planes.vertices.push_back(p2); planes.vertices.push_back(p3);
		Vector_3 v;
		if(i < oritations.size())
		v = Vector_3(oritations[i].x(),oritations[i].y(),0);
		else 
		Vector_3 Vector_3(0,0,0);
		Plane_3 plane(p1,p2,p3);

		if ( plane.orthogonal_vector()*v > 0 ){
		planes.faces.push_back(Point_3(idx,idx+1,idx+2));
		planes.faces.push_back(Point_3(idx+5,idx+4,idx+3));
		}
		else{
		planes.faces.push_back(Point_3(idx+1,idx,idx+2));
		planes.faces.push_back(Point_3(idx+4,idx+5,idx+3));
		}
	}
}


struct TriIndex
{
    unsigned value[3];
    TriIndex () {
		for (unsigned i = 0; i < 3; ++i)
				value[i] = -1;}
    TriIndex (unsigned i1, unsigned i2, unsigned i3){
			value[0] = i1;
			value[1] = i2;
			value[2] = i3;}
		
	unsigned &operator[](unsigned i){return value[i];}	
};
typedef std::vector<TriIndex> TriIndices;

unsigned findVertexIndex(Point_2 pt,std::vector<Point_2> points)
{
    std::vector<Point_2>::iterator item = std::find(points.begin(), points.end(), pt);
	if (item != points.end())
	{
		unsigned pos = std::distance(points.begin(), item);
		return pos;
    }
	return -1;
}
 
TriIndices  GetTriIndices (vector<Point_2> points)
{
    Polygon_list polys;
	CGAL::approx_convex_partition_2(points.begin(),
    points.end(),std::back_inserter(polys));
	LOG(INFO) << polys.size() << " convex polygons.";

    CDT cdt;
	Polygon2 poly;
	unsigned vertIndex;
	TriIndex triIndex;
	TriIndices triIndices;

	for (Polygon_list::iterator it = polys.begin(); it != polys.end();++it)
	{
		
		//约束三角剖分
		cdt.clear();
		cdt.insert(it->vertices_begin(), it->vertices_end());
		//cdt.insert(points.begin(), points.begin() + 5);

	for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it)
		{
			for (unsigned i = 0; i < 3; ++i)
			{
				Point_2 point = it->vertex(i)->point();
				vertIndex = findVertexIndex(point, points);
				triIndex[i] = vertIndex;
			}
	    triIndices.push_back(triIndex);
	    } 
	}
    return triIndices;
	//}
}

unsigned findthirdvertex(vector<Point_2>& points, Segment_2& seg, TriIndices& tris){
	int flag[3]={0,0,0};
	for(auto tri: tris){
		for(int k = 0; k <3; k++)
		if(seg.source() == points[tri[k]] || seg.target() == points[tri[k]]) flag[k] = 1;
		if(flag[0] + flag[1] + flag[2] == 2)
			return flag[0] == 0? tri[0]: (flag[1] == 0? tri[1]: tri[2]);
		flag[0] = flag[1] = flag[2] = 0;
	}
    return -1;
}

void normal_adjust(vector<Segment_2> segments, vector<Vector_2>& oritations){

    map<Point_2,int> points;
    vector<Point_2> single_points;
	int num_seg = segments.size();

    double xmin = 1e9; double ymin = 1e9;
    double xmax = -1e9; double ymax = -1e9;
 
   for(auto s:segments){
       if (points.count(s.source()))
           points[s.source()]++;
       else
           points.insert(make_pair(s.source(),1));
       if (points.count(s.target()))
           points[s.target()]++;
       else
           points.insert(make_pair(s.target(),1));
       
   }

   map<Point_2,int>::iterator iter;
    for (iter=points.begin(); iter!=points.end(); iter++){
       if(iter->second == 1)
         single_points.push_back(iter->first);
       if (iter->first.x() < xmin) xmin = iter->first.x();
       if (iter->first.x() > xmax) xmax = iter->first.x();
       if (iter->first.y() < ymin) ymin = iter->first.y();
       if (iter->first.y() > ymax) ymax = iter->first.y();  
	}

    // clear points
    map<Point_2,int>().swap(points);
    malloc_trim(0);

    if(single_points.size()%2==1){
        LOG(INFO) << "Facade error.";
        return;
    }

	vector<Point_2> xmins;	vector<Point_2> xmaxs;
	vector<Point_2> ymins; 	vector<Point_2> ymaxs;
	for(auto p: single_points){
		if (abs(p.x()-xmin) < 1e-5) {xmins.push_back(p);continue;}
		if (abs(p.x()-xmax) < 1e-5) {xmaxs.push_back(p);continue;}
		if (abs(p.y()-ymin) < 1e-5) {ymins.push_back(p);continue;}
		if (abs(p.y()-ymax) < 1e-5) {ymaxs.push_back(p);continue;}
	}

	if(xmins.size()%2 != 0 || ymins.size()%2 != 0 || xmaxs.size()%2 != 0 || ymaxs.size()%2 != 0 ){
		LOG(INFO) << "Point error.";
		LOG(INFO) << "Boundary size: " << xmins.size() << " " << ymins.size() << " " << xmaxs.size() << " "<< ymaxs.size();
        return;
	}
	for(int i = 0; i < xmins.size(); i=i+2)
       segments.push_back(Segment_2(xmins[i],xmins[i+1]));
	for(int i = 0; i < ymins.size(); i=i+2)
       segments.push_back(Segment_2(ymins[i],ymins[i+1]));
	for(int i = 0; i < xmaxs.size(); i=i+2)
       segments.push_back(Segment_2(xmaxs[i],xmaxs[i+1]));
	for(int i = 0; i < ymaxs.size(); i=i+2)
       segments.push_back(Segment_2(ymaxs[i],ymaxs[i+1]));
	// closed segments

	vector<Segment_2> segs;
	int * flag = new int[segments.size()];
	flag[0] = 1;
	Point_2 tt = segments[0].target();
	segs.push_back(segments[0]);
 	for(int i = 1; i < segments.size(); i++){
		for(int j = 0; j < segments.size(); j++){
            if(flag[j]==1) continue;
			if(segments[j].source() == tt){
				flag[j] = 1; 
				segs.push_back(segments[j]);
				tt = segments[j].target();
				break;
			}
			if(segments[j].target() == tt){
				flag[j] = 1; 
				segs.push_back(segments[j].opposite());
				tt = segments[j].source();
				break;
			}
		}	
	}
	delete [] flag;
	// assume one closed facades
	/*
	if(segs.size() != segments.size()){
		LOG(INFO) << "Facade error.";
		return;
	}*/
	
	int counter;//0-顺 1-逆
	Point_2 last,cur,next;
	vector<Point_2> counterwise_p;
	for(int s = 0; s < segs.size(); s++)
	   if(segs[s].target().x() == xmin){
		   last = segs[s].source();
		   cur = segs[s].target();
		   next = segs[(s+1)%segs.size()].target();
		   if(next.x() > xmin)
		      break;
	   }
	if ((cur - last).x() * (next - cur).y() - (next - cur).x() * (cur - last).y() > 0) // counterwise
	   counter = 1;
	else
	   counter = 0;
	for(int i = 0; i <segs.size(); i++){
		if(counter == 1)
		counterwise_p.push_back(segs[i].target());
		else
		counterwise_p.push_back(segs[segs.size()-1-i].target());
	}

	TriIndices tris = GetTriIndices(counterwise_p);

	for(int i = 0; i < num_seg; i++){
		int index = findthirdvertex(counterwise_p, segments[i], tris);
		if (index == -1){
			LOG(INFO) << "Normal maybe inaccurate.";
			//return;
			oritations.push_back(Vector_2(0, 0));
		}
		else
		oritations.push_back(Vector_2(CGAL::midpoint(segments[i].source(),segments[i].target()), counterwise_p[index]));
	}	
}