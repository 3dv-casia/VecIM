#include "base.h"
#include <glog/logging.h>
#include <limits>

#include <CGAL/property_map.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/Regularization/regularize_planes.h>

// Type declarations
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<K,
  Pwn_vector, Point_map, Normal_map>            Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits> Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<Traits>            Plane;
typedef CGAL::Shape_detection_3::Cylinder<Traits>         Cylinder;

// intersection projection
void PtoS(Plane* plane,Pwn_vector& points,int num,vector<Line>& lines, vector<Segment_2>& seg)
{
   Line line;
   line.ID = num;
   double sum_distances=0;
  
   int N=0;
   double xmin=1e9;double xmax=-1e9;double x1,y1,x2,y2;
   // compute line parameter
   double k = -plane->plane_normal().x()/plane->plane_normal().y();
   double b = 0;
   Vector_2 normal(plane->plane_normal().x(), plane->plane_normal().y());
   Line_2 dir(Point_2(0,0), normal.perpendicular(CGAL::POSITIVE));
  

   vector<size_t>::const_iterator
      index_it = plane->indices_of_assigned_points().begin();
    double vx=0, vy=0;
    double disx = 0; // normal: 1, 0, 0
    while (index_it != plane->indices_of_assigned_points().end()) {
      
      // Retrieves point
      const Point_with_normal &p = *(points.begin() + (*index_it));  
      sum_distances += CGAL::sqrt(plane->squared_distance(p.first));
      Point_2 p2(Point_2(p.first.x(),p.first.y()));
      line.pointset.push_back(p.first);
      if(line.height0 < p.first.z()) // max
	      line.height0 = p.first.z();
      if(line.height1 > p.first.z()) // min
      	line.height1 = p.first.z();
      b += -k*p.first.x()+p.first.y();
      vx += p.second.x();
      vy += p.second.y();
      if(plane->plane_normal().y() == 0) disx += p.first.x();
      // Proceeds with next point.
      index_it++;
      N++;

      // seg record     
      if (xmin>(dir.projection(p2)).x() || xmin==(dir.projection(p2)).x() && p2.y() < y1) {xmin=dir.projection(p2).x(); x1 = p2.x(); y1 = p2.y();}
      if (xmax<(dir.projection(p2)).x() || xmax==(dir.projection(p2)).x() && p2.y() > y2) {xmax=dir.projection(p2).x(); x2 = p2.x(); y2 = p2.y();}
    }
    //cout <<"Euclidean distance: "<<sum_distances<<endl;
    if(plane->plane_normal().y() == 0) {
      disx /= N;
      line.line2 = Line_2(Point_2(disx, 1), Point_2(disx, -1));
    }
    else
      line.line2 = Line_2(k,-1,b/N);
    line.distance = sum_distances/N;
    line.s = line.line2.projection(Point_2(x1,y1));
    line.t = line.line2.projection(Point_2(x2,y2));
    line.normal = Vector_2(vx/N,vy/N);
    if(!point_2_range(line.s) || !point_2_range(line.t)){
      LOG(INFO) << "invalid line: "<< line.s << " " << line.t;
      return;
    }
    lines.push_back(line);
    seg.push_back(Segment_2(line.s,line.t));
    //LOG(INFO) << "normal test: " << Segment_2(line.s,line.t) << " " << line.normal;
}

// clean circles intersection
void Clean(vector<Circle>& circles, float r){
     sort(circles.begin(),circles.end());// 降序
     for(int i = 0; i < circles.size();i++){
       if(circles[i].radius > r){
         circles[i].flag = 0;
         continue;
       }
       if(circles[i].flag == 0) continue;
       for(int j = i+1; j < circles.size();j++){
          if(circles[j].radius > r){
           circles[j].flag = 0;
           continue;
          }
          if(circles[j].flag == 0) continue;        
          if (sqrt((circles[i].center-circles[j].center).squared_length()) < circles[i].radius + circles[j].radius)
             circles[j].flag = 0;
       }
     }
}


// a point is in a polygon or not
bool isin(vector<Point_3> &points, vector<int>& poly, Point_2 point, double& h){
    int tag = 0;
    if(poly.size() < 3)
    return false;
    for(int i = 0; i < poly.size(); i++){
        Segment_2 seg(Point_2(points[poly[i]].x(), points[poly[i]].y()), Point_2(points[poly[(i+1)%poly.size()]].x(), points[poly[(i+1)%poly.size()]].y()));
        if(seg.has_on(point)){
          tag = 1;
          break;
        }
        auto result = CGAL::intersection(seg, K::Ray_2(point,Vector_2(0,1)));
        if(result)
        if(auto p = boost::get<Point_2>(&*result))
          tag++;
    }
    if(tag%2 == 1){
      Point_3 p1 = points[poly[0]], p2 = points[poly[0]], p3 = points[poly[0]];
      for(auto p:poly)
       if((points[p]-p1).squared_length() > 1e-6){
         p2 = points[p];
         break;
       }
      for(auto p:poly)
       if((points[p]-p1).squared_length() > 1e-6 && (points[p]-p2).squared_length() > 1e-6){
         if((Line_3(p1, p2).projection(points[p]) - points[p]).squared_length() < 1e-6) continue;
         p3 = points[p];
         break;
       }
       if(p1 == p2 && p1 == p3)
          return false;
        Plane_3 pl(p1, p2, p3);
        auto result = CGAL::intersection(Line_3(Point_3(point.x(), point.y(), 0), Vector_3(0, 0, 1)), pl);
        if(result)
          if(auto *p = boost::get<Point_3>(&*result)){
               h = p->z();
               return true;
          }
        return false;
    }
    else
          return false;

}

bool isin(vector<Point_3> &points, Point_3& poly, Point_2 point, double& h){
    int tag = 0;
    for(int i = 0; i < 3; i++){
        Segment_2 seg(Point_2(points[poly[i]].x(), points[poly[i]].y()), Point_2(points[poly[(i+1)%3]].x(), points[poly[(i+1)%3]].y()));
        if(seg.has_on(point)){
          tag = 1;
          break;
        }
        auto result = CGAL::intersection(seg, K::Ray_2(point,Vector_2(0,1)));
        if(result)
        if(auto p = boost::get<Point_2>(&*result))
          tag++;
    }
    if(tag%2 == 1){
      Point_3 p1 = points[poly[0]], p2 = points[poly[0]], p3 = points[poly[0]];
      for(int i = 0; i < 3; i++)
       if((points[poly[i]]-p1).squared_length() > 1e-6){
         p2 = points[poly[i]];
         break;
       }
      for(int i = 0; i < 3; i++)
       if((points[poly[i]]-p1).squared_length() > 1e-6 && (points[poly[i]]-p2).squared_length() > 1e-6){
         if((Line_3(p1, p2).projection(points[poly[i]]) - points[poly[i]]).squared_length() < 1e-6) continue;
         p3 = points[poly[i]];
         break;
       }
       if(p1 == p2 && p1 == p3)
          return false;
        Plane_3 pl(p1, p2, p3);
        auto result = CGAL::intersection(Line_3(Point_3(point.x(), point.y(), 0), Vector_3(0, 0, 1)), pl);
        if(result)
          if(auto *p = boost::get<Point_3>(&*result)){
               h = p->z();
               return true;
          }
        return false;
    }
    else
          return false;

}

// circles to planes (8)
bool Toplane(vector<Circle>& circles, string wdir){

  int  num = 0;
  vector<Point_2> vertexes;
  // generate cylinder file
  ofstream ofss(wdir + "cylinder_final");
  for(int i = 0; i < circles.size(); i ++){
    if (circles[i].flag == 0) continue;
    num++;
    // 8 points
    Vector_2 add[8] ={Vector_2(0, circles[i].radius), Vector_2(circles[i].radius/1.414,circles[i].radius/1.414), 
                      Vector_2(circles[i].radius, 0), Vector_2(circles[i].radius/1.414, -circles[i].radius/1.414), 
                      Vector_2(0, -circles[i].radius), Vector_2(-circles[i].radius/1.414, -circles[i].radius/1.414), 
                      Vector_2(-circles[i].radius, 0), Vector_2(-circles[i].radius/1.414, circles[i].radius/1.414)};

    for(int j = 0; j < 8; j++)
      vertexes.push_back(circles[i].center + add[j]);
    for(int j = 8*i + 0; j < 8*i + 8 -1 ; j++)
     ofss << vertexes[j].x() << " " << vertexes[(j+1)%(8*(i+1))].x() << " " << vertexes[j].y() << " " << vertexes[(j+1)%(8*(i+1))].y() << "\n";
    ofss << vertexes[8*i + 7].x() << " " << vertexes[8*i + 0].x() << " " << vertexes[8*i + 7].y() << " " << vertexes[8*i + 0].y() << "\n";
  }
  LOG(INFO) << num << " cylinders remained."; 
  ofss.close();
/*
  // generate height according to floor.off ceiling.off
  vector<Point_2> height(vertexes.size()); // {max, min}
  ifstream floor(wdir + "floor_vec.ply");
  ifstream ceil(wdir + "ceiling_vec.ply");
  vector<Point_3> pointsf;
  vector<vector<int>> polysf;
  vector<Point_3> pointsc;
  vector<vector<int>> polysc;
  int vn, fn, o;
  string coff;
  double x,y,z;
  // floor.off
  floor >> coff;
  floor >> vn >> fn >> o;
  while(vn-->0){
     floor >> x >> y >> z;
     pointsf.push_back(Point_3(x, y, z));
  }
  while(fn-->0){
    floor >> o;
    vector<int> poly;
    while(o-->0){
      floor >> vn;
      poly.push_back(vn);
    }
    polysf.push_back(poly);
  }
  // ceiling.off
  ceil >> coff;
  ceil >> vn >> fn >> o;
  while(vn-->0){
     ceil >> x >> y >> z;
     pointsc.push_back(Point_3(x, y, z));
  }
  while(fn-->0){
    ceil >> o;
    vector<int> poly;
    while(o-->0){
      ceil >> vn;
      poly.push_back(vn);
    }
    polysc.push_back(poly);
  }
  */

   // generate height according to floor.ply ceiling.ply
  vector<Point_2> height(vertexes.size()); // {max, min}
  string floor(wdir + "floor_vec.ply");
  string ceiling(wdir + "ceiling_vec.ply");
  vector<Point_3> pointsf;
  vector<vector<int>> polysf;
  vector<Point_3> pointsc;
  vector<vector<int>> polysc;
  Mesh<Point_3> meshf,meshc;
  PLYMeshLoad(floor, meshf);
  PLYMeshLoad(ceiling, meshc);

  for(int i = 0; i < vertexes.size(); i++){
       double hc=1e9, hf=1e9;
       for(int j = 0; j < meshf.faces.size(); j++){
         if(isin(meshf.vertices, meshf.faces[j], vertexes[i], hf))
            break;
       }
       for(int j = 0; j < meshc.faces.size(); j++){
         if(isin(meshc.vertices, meshc.faces[j], vertexes[i], hc)){
            break;
         }
       }
       if(hc==1e9 || hf==1e9){
         LOG(INFO) << vertexes[i] << " Cylinder failed.";
         return false;
       }
       height[i] = Point_2(hc, hf);
  }

  // obtain height for each vertexes, 8 vertexes -> 1 cylinder
  vector<Point_3> points;
  vector<vector<int>> polygons;
  for(int i = 0; i < vertexes.size(); i=i+8){
    for(int j = 0; j < 8; j++){
      Point_3 p1 = Point_3(vertexes[i+j].x(), vertexes[i+j].y(), height[i+j].x());
      Point_3 p2 = Point_3(vertexes[i+j].x(), vertexes[i+j].y(), height[i+j].y());
      Point_3 p3 = Point_3(vertexes[i+(j+1)%8].x(), vertexes[i+(j+1)%8].y(), height[i+(j+1)%8].y());
      Point_3 p4 = Point_3(vertexes[i+(j+1)%8].x(), vertexes[i+(j+1)%8].y(), height[i+(j+1)%8].x());
      int id1,id2,id3,id4;
      auto it = find(points.begin(), points.end(), p1);
      id1 = std::distance(points.begin(), it);
      if(it == points.end())
       points.push_back(p1);

      it = find(points.begin(), points.end(), p2);
      id2 = std::distance(points.begin(), it);
      if(it == points.end())
       points.push_back(p2);

      it = find(points.begin(), points.end(), p3);
      id3 = std::distance(points.begin(), it);
      if(it == points.end())
       points.push_back(p3);

      it = find(points.begin(), points.end(), p4);
      id4 = std::distance(points.begin(), it);
      if(it == points.end())
       points.push_back(p4);

      polygons.push_back({id1, id4, id3, id2});
    }
  }
/*
  // generate cylinder.off
  ofstream ofs(wdir + "cylinder.off");
  ofs << "COFF\n" << points.size() << ' ' << polygons.size() << ' ' << "0\n";
  for (const auto &p : points)
    ofs << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
  for (const auto &plg : polygons) {
    ofs << plg.size();
    for (const auto &p : plg)
      ofs << ' ' << p;
    ofs << '\n';
  }
  ofs.close();
*/
  string ofs = wdir + "cylinder_vec.ply";
  Mesh<Point_3> mesh;
  for(auto& p: points)
    mesh.vertices.push_back(p);
  for(auto& poly: polygons){
    for(int i = 2; i < poly.size(); i++){
        mesh.faces.push_back({poly[0], poly[i-1], poly[i]});
    }
  }
  PLYMeshSave(ofs, mesh, 1);

  return true;
}

bool ransac_detection_p(Pwn_vector& points, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, vector<Line>& lines, string wdir)
{
  // timer
  clock_t start,mid,end;
  start = clock(); 
  
  //test random
  //std::random_shuffle(points.begin(), points.end());
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>();// Register shapes for detection 
  Efficient_ransac::Parameters parameters;// Sets parameters for shape detection.
  parameters.probability = probability; 
  parameters.min_points = min_points;
  parameters.epsilon = epsilon;  
  parameters.normal_threshold = normal_threshold;   
  if (cluster_epsilon)
    parameters.cluster_epsilon = cluster_epsilon; 
  
  double re_angle = cm::get_config().get<double>("facade.ransac.angle");

  // Detects shapes
  for(int i = 0; i < 20; i++){
    ransac.detect(parameters);
    if(ransac.shapes().end() - ransac.shapes().begin() > 0) continue;
  }
  // Prints number of detected shapes and unassigned points.
  LOG(INFO) << ransac.shapes().end() - ransac.shapes().begin() << " detected planes, "
     << ransac.number_of_unassigned_points()
     << " unassigned points." ;

  if (ransac.shapes().end() - ransac.shapes().begin() ==0)
	return 0;

  Efficient_ransac::Plane_range planes = ransac.planes();

  LOG(INFO) << "ransac regularization angle: " << re_angle;
                   
  // Regularize detected planes.
  CGAL::regularize_planes(points,
                          Point_map(),
                          planes,
                          CGAL::Shape_detection::Plane_map<Traits>(),
                          CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                          true,  // regularize parallelism
                          true,  // regularize orthogonality
                          true, //  regularize coplanarity
                          true,  // regularize Z-symmetry (default)
                          re_angle,
                          0.02);   // 10 degrees of tolerance for parallelism / orthogonality
                      

  vector<Segment_2> seg;
  vector<Point_2> ps;
  vector<Point_2> support_pp;

  int num=0;
  Efficient_ransac::Plane_range::iterator it = planes.begin();

  // record possible ceilings
  ofstream fce(wdir + "ceilings_facade", ios::out);
  Pwn_vector cpoints;
  while (it != planes.end()) {

    if (Plane* plane = dynamic_cast<Plane*>(it->get()))
      {
        Vector_3 normal = plane->plane_normal(); 
        // remain horizantal ceiling (thethold=80) 
        //if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z()))> 1 - cos_angle) {
        // remain non-facade ceiling (thethold=80) 
        if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z()))> cos_angle) {
          // find ave z height and record ceiling points
          double z = 0;
          for(auto p = plane->indices_of_assigned_points().begin(); p!= plane->indices_of_assigned_points().end();p++){
            z += (points.begin() + *p)->first.z();
            cpoints.push_back(points[*p]);
          }
          z /= plane->indices_of_assigned_points().size();
          if(plane->plane_normal() == Vector_3(0,0,1) || plane->plane_normal() == Vector_3(0,0,-1))
              fce << plane->plane_normal() <<" "<< plane->d() << "\n";
          else{
              if(plane->plane_normal().z() > 0)
                fce << plane->plane_normal() <<" "<< plane->d() << "\n";//fce << 0 << " " << 0 << " " << 1 <<" "<< z << "\n";
              else
                fce << plane->plane_normal() <<" "<< plane->d() << "\n";//fce << 0 << " " << 0 << " " << -1 <<" "<< z << "\n";  
          }
         
        } 
 
        if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z()))>cos_angle) {it++; continue;} // remove non-facade (thethold=80)

        // store point sets and 2D lines ->facade plane
       	PtoS(plane, points, num, lines, seg);
        num++;
      }  
    // Proceeds with next detected shape.
    it++;
  }
  fce.close();

  //TODO: the second detection
  Efficient_ransac ransac2;
  Pwn_vector points2;
  for(auto &p:ransac.indices_of_unassigned_points())
      points2.push_back(points[p]);
  ransac2.set_input(points2);
  ransac2.add_shape_factory<Plane>();// Register shapes for detection 
  // Detects shapes
  parameters.min_points = 50;
  ransac2.detect(parameters);
  Efficient_ransac::Plane_range planes2 = ransac2.planes();
   // Regularize detected planes.
  CGAL::regularize_planes(points2,
                          Point_map(),
                          planes2,
                          CGAL::Shape_detection::Plane_map<Traits>(),
                          CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points2, planes2),
                          true,  // regularize parallelism
                          true,  // regularize orthogonality
                          false, // do not regularize coplanarity
                          true,  // regularize Z-symmetry (default)
                          re_angle,
                           0.02);   // 10 degrees of tolerance for parallelism / orthogonality          

    it = planes2.begin();
    while (it != planes2.end()) {
      if (Plane* plane = dynamic_cast<Plane*>(it->get()))
        {
          Vector_3 normal = plane->plane_normal();
          if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z()))>cos_angle) {it++; continue;} // remove non-facade (thethold=80)

          // store point sets and 2D lines ->facade plane
       	  PtoS(plane, points2, num, lines, seg);
          num++;
        }  
      // Proceeds with next detected shape.
      it++;
    }
    LOG(INFO) << ransac2.shapes().end() - ransac2.shapes().begin() << " detected planes";
  

  PLYPointSave(wdir + "ceilings_facade.ply", cpoints);
  LOG(INFO) << num << " lines are generated. ";
  end = clock();
  LOG(INFO) << "Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s";
  if(num != lines.size())
    LOG(INFO) << num << " " << lines.size() << " not match.";
  // record 2d segments and facade points
  Pwn_vector facade_points;
  ofstream file(wdir + "2dsegments", ios::out);
	for(int i = 0; i < lines.size(); i++){
	  file << lines[i].s.x() << " " << lines[i].t.x() << " "<< lines[i].s.y() << " " << lines[i].t.y() << endl;
	  for(auto p: lines[i].pointset)
      facade_points.push_back(Point_with_normal({p,Vector_3(0,1,0)}));
	  
  }
	file.close();
  PLYPointSave(wdir + "facade_points_ransac.ply", facade_points);
	

  return 1;
}

bool ransac_detection_c(Pwn_vector& points, float radius, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, vector<Circle>& circles, string wdir){
// timer
  clock_t start, mid, end;
  start = clock(); 
  
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Cylinder>();// Register shapes for detection 
  Efficient_ransac::Parameters parameters;// Sets parameters for shape detection.
  parameters.probability = probability; 
  parameters.min_points = min_points;
  parameters.epsilon = epsilon;  
  parameters.normal_threshold = normal_threshold;   
  if (cluster_epsilon)
    parameters.cluster_epsilon = cluster_epsilon; 

  // Detects shapes
  ransac.detect(parameters);
  // Prints number of detected shapes and unassigned points.
  LOG(INFO) << ransac.shapes().end() - ransac.shapes().begin() << " detected cylinders, "
     << ransac.number_of_unassigned_points()
     << " unassigned points." ;

  if (ransac.shapes().end() - ransac.shapes().begin() == 0){
    LOG(INFO) << "No cylinder is detected.";
    return 0;
  }

  Efficient_ransac::Shape_range cylinders = ransac.shapes();
  Pwn_vector cylinder_points;
  int num=0;
  Efficient_ransac::Shape_range::iterator it = cylinders.begin();
  while (it != cylinders.end()) {

    if (Cylinder* cylinder = dynamic_cast<Cylinder*>(it->get()))
      {
        Line_3 axis = cylinder->axis();   
        if (abs(axis.to_vector().z()/sqrt(axis.to_vector().x()*axis.to_vector().x()+axis.to_vector().y()*axis.to_vector().y()+axis.to_vector().z()*axis.to_vector().z()))<cos_angle) {it++; continue;} // remove non-facade (thethold=80)
        double mid = 0;
        int num_c = 0;// find the middle height !
        for(auto p:cylinder->indices_of_assigned_points()){
           cylinder_points.push_back(points[p]);
           mid += points[p].first.z();
           num_c++;
        }    
        Plane_3 mid_p(Point_3(0,0,mid/num_c), Vector_3(0,0,1)); 
        auto result = CGAL::intersection(mid_p, axis);    
	      if (result)
	        if(const Point_3* pt = boost::get<Point_3>(&*result))	
	            circles.push_back(Circle(Point_2(pt->x(), pt->y()),cylinder->radius(),cylinder->indices_of_assigned_points().end()-cylinder->indices_of_assigned_points().begin(), 0, 0));
           else
             LOG(INFO) << "Cylinder detection is faled! ";       
      }  
    // Proceeds with next detected shape.
    it++;
    num++;
  }
  PLYPointSave(wdir + "cylinders.ply", cylinder_points);
  LOG(INFO) << num << " cylinders are generated. ";
  
  // exclude intersection circles
  Clean(circles, radius);

  // valid circles produce eight planes
  if(!Toplane(circles, wdir))
    return 0;

  // save cylinder file
  string fname = wdir +"cylinder";
  ofstream ofs(fname);
  for(int i = 0; i < circles.size(); i++){
    if (circles[i].flag == 0) continue;
    ofs << circles[i].center.x() << " " << circles[i].center.y() <<" "<< circles[i].radius << "\n";
  }
  ofs.close();

  end = clock();
  LOG(INFO) << "Cylinder RANSAC Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s";
  return 1;
}

// ransac one floor ceiling
bool ransac_detection_fc(Pwn_vector& points, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, string wdir, string object){
  // timer
  clock_t start,mid,end;
  start = clock(); 
  
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>();// Register shapes for detection 
  Efficient_ransac::Parameters parameters;// Sets parameters for shape detection.
  parameters.probability = probability; 
  parameters.min_points = min_points;
  parameters.epsilon = epsilon;  
  parameters.normal_threshold = normal_threshold;   
  if (cluster_epsilon)
    parameters.cluster_epsilon = cluster_epsilon; 

  // Detects shapes
  ransac.detect(parameters);
  // Prints number of detected shapes and unassigned points.
  LOG(INFO) << ransac.number_of_unassigned_points() << " unassigned points." ;

  if (ransac.shapes().end() - ransac.shapes().begin() == 0)
	return 0;

  double height = 0;
  int num = 0;
  Efficient_ransac::Shape_range planes = ransac.shapes();
  Plane* p = dynamic_cast<Plane*>(ransac.shapes().begin()->get());
  vector<size_t>::const_iterator
      index_it = p->indices_of_assigned_points().begin();
  while (index_it != p->indices_of_assigned_points().end()) {
     height += (*(points.begin() + (*index_it))).first.z();
     num ++;
     index_it++;
  }    
  height /= num;
  ofstream ofs(wdir + object);
  ofs << height;
  ofs.close();
  return 1;    
}

// ransac maybe multi floors
bool ransac_detection_floor(Pwn_vector& points,  vector<Plane_3>& floors, vector<Line_2>& intersects, const cm::Config& config, string wdir){
  // timer
  clock_t start,mid,end;
  start = clock(); 
  
  vector<Plane_3> verticals;// record vertical floors
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>();// Register shapes for detection 
  Efficient_ransac::Parameters parameters;// Sets parameters for shape detection.
  parameters.probability = *(config.get_optional<double>("floor.ransac.probability")); 
  parameters.min_points = *(config.get_optional<double>("floor.ransac.min_points"));
  parameters.epsilon = *(config.get_optional<double>("floor.ransac.epsilon"));  
  parameters.normal_threshold = *(config.get_optional<double>("floor.ransac.normal_threshold"));   
  if (*(config.get_optional<double>("floor.ransac.cluster_epsilon")))
    parameters.cluster_epsilon = *(config.get_optional<double>("floor.ransac.cluster_epsilon")); 

  // Detects shapes
  ransac.detect(parameters);

  
  if (ransac.shapes().end() - ransac.shapes().begin() ==0)
	return 0;

  Efficient_ransac::Plane_range planes = ransac.planes();


   // Regularize detected planes.
  CGAL::regularize_planes(points,
                          Point_map(),
                          planes,
                          CGAL::Shape_detection::Plane_map<Traits>(),
                          CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                          true,  // regularize parallelism
                          true,  // regularize orthogonality
                          false, // do not regularize coplanarity
                          true,  // regularize Z-symmetry (default)
                          15,
                           0.15);   // 10 degrees of tolerance for parallelism / orthogonality
                      

  int num=0;
  Efficient_ransac::Plane_range::iterator it = planes.begin();

  while (it != planes.end()) {

    if (Plane* plane = dynamic_cast<Plane*>(it->get()))
      {
        Vector_3 normal = plane->plane_normal();
        if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) < 1 - *(config.get_optional<double>("floor.ransac.cos_angle"))) {
          verticals.push_back(Plane_3(normal.x(), normal.y(), normal.z(), plane->d()));
        }      
        if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) < *(config.get_optional<double>("floor.ransac.cos_angle"))) {it++; continue;}
        
        // assume pass through the first supporting point
        if(normal.z() < 0)
           floors.push_back(Plane_3((points.begin() + *(plane->indices_of_assigned_points().begin()))->first, Vector_3(0,0,-1)));
        else
           floors.push_back(Plane_3((points.begin() + *(plane->indices_of_assigned_points().begin()))->first, Vector_3(0,0,1)));
        num++;
      }  
    // Proceeds with next detected shape.
    it++;
  }
  //LOG(INFO) << num << " floors are detected.";

  if(num == 1) {// only one floor
     end = clock();
     LOG(INFO) << "Floor Ransac Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s";
     return true;
  }
  else { // multi floor -> detect tiny vertical plane
    // Prints number of detected shapes and unassigned points.
    Efficient_ransac ransac2;
    Pwn_vector points2;
    for(auto &p:ransac.indices_of_unassigned_points())
       points2.push_back(points[p]);
    LOG(INFO) << "detect vertical planes from " << points2.size() << " unassigned points." ;
    ransac2.set_input(points2);
    ransac2.add_shape_factory<Plane>();// Register shapes for detection 
    Efficient_ransac::Parameters parameters2;// Sets parameters for shape detection.
    parameters2.probability = *(config.get_optional<double>("floor.ransac2.probability")); 
    parameters2.min_points = *(config.get_optional<double>("floor.ransac2.min_points"));
    parameters2.epsilon = *(config.get_optional<double>("floor.ransac2.epsilon"));  
    parameters2.normal_threshold = *(config.get_optional<double>("floor.ransac2.normal_threshold"));   
    if (*(config.get_optional<double>("floor.ransac2.cluster_epsilon")))
      parameters2.cluster_epsilon = *(config.get_optional<double>("floor.ransac2.cluster_epsilon")); 

    // Detects shapes
    ransac2.detect(parameters2);

  Efficient_ransac::Plane_range planes2 = ransac2.planes();

   // Regularize detected planes.
  CGAL::regularize_planes(points,
                          Point_map(),
                          planes,
                          CGAL::Shape_detection::Plane_map<Traits>(),
                          CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                          true,  // regularize parallelism
                          true,  // regularize orthogonality
                          false, // do not regularize coplanarity
                          true,  // regularize Z-symmetry (default)
                          5,
                           0.15);   // 10 degrees of tolerance for parallelism / orthogonality
                      

      // unsure: median results
      // record vertical points
      Pwn_vector verticals_floor;
      it = planes2.begin();
      while (it != planes2.end()) {
        if (Plane* plane = dynamic_cast<Plane*>(it->get()))
          {
            Vector_3 normal = plane->plane_normal();
            if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) < 1 - *(config.get_optional<double>("floor.ransac2.cos_angle"))) {
              // centroid points
              double x,y,z;
              x=y=z=0;
              int num = plane->indices_of_assigned_points().end() - plane->indices_of_assigned_points().begin(); 
              for(auto it2 = plane->indices_of_assigned_points().begin(); it2 != plane->indices_of_assigned_points().end(); it2++){
	                verticals_floor.push_back(points2[*it2]);
                  x+=points2[*it2].first.x(); y+=points2[*it2].first.y(); z+=points2[*it2].first.z();
	            }
              verticals.push_back(Plane_3(Point_3(x/num, y/num, z/num), Vector_3(normal.x(), normal.y(), 0)));
              // fixed: 强制转化为垂直面 normal.z() = 0
            }
          }  
        // Proceeds with next detected shape.
        it++;
      }
	  PLYPointSave(wdir + "verticals_floor.ply", verticals_floor);
    LOG(INFO) << verticals.size() << " vertical(s) are deteced.";

    if(verticals.size() == 0){
      LOG(INFO) << "maybe missing tiny vertical planes.";
      return true;
    }
    // vertical and horizantal floors intersect
    for(auto h: floors)
      for(auto v: verticals){
        if(h == *floors.begin()) continue;
        auto result = CGAL::intersection(h, v);
          if(Line_3* l = boost::get<Line_3>(&*result)){
            Line_2 line(Point_2(l->projection(Point_3(0,0,0)).x(),l->projection(Point_3(0,0,0)).y()),Vector_2(l->to_vector().x(),l->to_vector().y()));
            auto it = find(intersects.begin(), intersects.end(), line);
            auto it2 = find(intersects.begin(), intersects.end(), line.opposite());
            if(it == intersects.end() && it2 == intersects.end())
              intersects.push_back(line);            
          }
      }
   
  end = clock();
  LOG(INFO) << "Floor Ransac Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s";

  return 1; 
  }   
}


// ransac maybe multi ceilings
bool ransac_detection_ceiling(Pwn_vector& points,  vector<Plane_3>& ceilings, vector<Line_2>& intersects, const cm::Config& config, string wdir){
  // timer
  clock_t start,mid,end;
  start = clock(); 
  
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>();// Register shapes for detection 
  Efficient_ransac::Parameters parameters;// Sets parameters for shape detection.
  parameters.probability = *(config.get_optional<double>("ceiling.ransac.probability")); 
  parameters.min_points = *(config.get_optional<double>("ceiling.ransac.min_points"));
  parameters.epsilon = *(config.get_optional<double>("ceiling.ransac.epsilon"));  
  parameters.normal_threshold = *(config.get_optional<double>("ceiling.ransac.normal_threshold"));   
  if (*(config.get_optional<double>("ceiling.ransac.cluster_epsilon")))
    parameters.cluster_epsilon = *(config.get_optional<double>("ceiling.ransac.cluster_epsilon")); 

  // Detects shapes
  ransac.detect(parameters);

  
  if (ransac.shapes().end() - ransac.shapes().begin() ==0)
	return 0;

  Efficient_ransac::Plane_range planes = ransac.planes();
 
   // Regularize detected planes.
  CGAL::regularize_planes(points,
                          Point_map(),
                          planes,
                          CGAL::Shape_detection::Plane_map<Traits>(),
                          CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                          true,  // regularize parallelism
                          true,  // regularize orthogonality
                          false, // do not regularize coplanarity
                          true,  // regularize Z-symmetry (default)
                          15,
                           0.15);   // 10 degrees of tolerance for parallelism / orthogonality
                          

  Efficient_ransac::Plane_range::iterator it = planes.begin();
  // record detected celings
  Pwn_vector ceiling_points;

  while (it != planes.end()) {

    if (Plane* plane = dynamic_cast<Plane*>(it->get()))
      {
        Vector_3 normal = plane->plane_normal();
        // push points
        double x=0,y=0,z=0;
        for(auto p:plane->indices_of_assigned_points()){
          x+=points[p].first.x(); y+=points[p].first.y(); z+=points[p].first.z();
          ceiling_points.push_back(points[p]);
        }
        int num = plane->indices_of_assigned_points().end() - plane->indices_of_assigned_points().begin();
        // horizantal plane
        if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) > *(config.get_optional<double>("ceiling.ransac.cos_angle")))
        {
          // assume pass through the first supporting point
          if(normal.z() < 0)
            ceilings.push_back(Plane_3(Point_3(x/num, y/num, z/num), Vector_3(0,0,-1)));
          else
            ceilings.push_back(Plane_3(Point_3(x/num, y/num, z/num), Vector_3(0,0,1)));
        }
        // vertical plane
        else if(abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) < 1 - *(config.get_optional<double>("ceiling.ransac.cos_angle"))){
           ceilings.push_back(Plane_3(Point_3(x/num, y/num, z/num), Vector_3(normal.x(), normal.y(),0)));
        }
        else // slap plane
          ceilings.push_back(Plane_3(normal.x(), normal.y(), normal.z(), plane->d()));
         
        } 
    // Proceeds with next detected shape.
    it++;
  }
  PLYPointSave(wdir + "ceilings_ceiling.ply", ceiling_points);
  end = clock();
  LOG(INFO) << "Ceilings RANSAC Time: " << (float)(end-start)/CLOCKS_PER_SEC << "s";

  return 1;   
}
