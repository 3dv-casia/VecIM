#include "base.h"
#include <glog/logging.h>
#include <limits>

#include <CGAL/property_map.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/Regularization/regularize_planes.h>

// Type declarations
typedef CGAL::First_of_pair_property_map<Point_with_normal>     Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal>    Normal_map;
typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<K, Pwn_vector, Point_map, Normal_map>    Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>     Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<Traits>    Plane;
typedef CGAL::Shape_detection_3::Cylinder<Traits>     Cylinder;

/*!
 * \description: Project facade plane to line/segment.
 * \param plane detected facade plane by RANSAC
 * \param points facade pointcloud
 * \param num plane/line index
 * \param lines projection line set
 * \return
 */
void PtoS(Plane* plane, Pwn_vector& points, int num, vector<Line>& lines)
{
    Line line;
    line.ID = num;
    double sum_distances = 0;
    int N = 0;
    double xmin = 1e9; double xmax = -1e9; 
    double x1, y1, x2 ,y2;
    // compute line parameter
    double k = -plane->plane_normal().x()/plane->plane_normal().y();
    double b = 0;
    Vector_2 normal(plane->plane_normal().x(), plane->plane_normal().y());
    Line_2 dir(Point_2(0,0), normal.perpendicular(CGAL::POSITIVE));
    

    vector<size_t>::const_iterator index_it = plane->indices_of_assigned_points().begin();
    double vx = 0, vy = 0;
    double disx = 0; 
    while (index_it != plane->indices_of_assigned_points().end()) 
    {
        // retrieves point
        const Point_with_normal &p = *(points.begin() + (*index_it));  
        sum_distances += CGAL::sqrt(plane->squared_distance(p.first));
        Point_2 p2(Point_2(p.first.x(),p.first.y()));
        line.pointset.push_back(p.first);
        if(line.height0 < p.first.z()) // max
            line.height0 = p.first.z();
        if(line.height1 > p.first.z()) // min
            line.height1 = p.first.z();
        b += -k*p.first.x() + p.first.y();
        vx += p.second.x();
        vy += p.second.y();
        if(plane->plane_normal().y() == 0) disx += p.first.x();
        // proceeds with next point.
        index_it++;
        N++;
        // record segment endpoint    
        if (xmin > (dir.projection(p2)).x() || xmin == (dir.projection(p2)).x() && p2.y() < y1) { xmin = dir.projection(p2).x(); x1 = p2.x(); y1 = p2.y(); }
        if (xmax < (dir.projection(p2)).x() || xmax == (dir.projection(p2)).x() && p2.y() > y2) { xmax = dir.projection(p2).x(); x2 = p2.x(); y2 = p2.y(); }
    }
      if(plane->plane_normal().y() == 0) {
          disx /= N;
          line.line2 = Line_2(Point_2(disx, 1), Point_2(disx, -1));
      }
      else
          line.line2 = Line_2(k, -1, b/N);
      line.distance = sum_distances/N;
      line.s = line.line2.projection(Point_2(x1, y1));
      line.t = line.line2.projection(Point_2(x2, y2));
      line.normal = Vector_2(vx/N, vy/N);
      if(!point_2_range(line.s) || !point_2_range(line.t)){
          LOG(INFO) << "invalid line: "<< line.s << " " << line.t;
          return;
      }
      lines.push_back(line);
}

/*!
 * \description: Clean circle intersection.
 * \param circles detected cylinder projection
 * \param r radius threshold
 * \return
 */
void Clean(vector<Circle>& circles, float r)
{
     sort(circles.begin(), circles.end());
     for(int i = 0; i < circles.size(); i++){
          if(circles[i].radius > r){
            circles[i].flag = 0;
            continue;
       }
       if(circles[i].flag == 0) continue;
       for(int j = i+1; j < circles.size(); j++){
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

/*!
 * \description: whether a point is in a polygon or not (for .ply file)
 * \param points polygon endpoints
 * \param poly polygon
 * \param point query point
 * \param h point height
 * \return 
 */
bool isin(vector<Point_3> &points, Point_3& poly, Point_2 point, double& h)
{
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

/*!
 * \description: whether a point is in a polygon or not (for .off file)
 * \param points polygon endpoints
 * \param poly polygon
 * \param point query point
 * \param h point height
 * \return 
 */
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


/*!
 * \description: Octagonal fitting cylinder.
 * \param circles detected cylinder projection
 * \param wdir file save/load directory
 * \return
 */
bool Toplane(vector<Circle>& circles, string wdir)
{
    int num = 0;
    vector<Point_2> vertexes;
    // generate cylinder file
    ofstream ofss(wdir + "cylinder_final");
    for(int i = 0; i < circles.size(); i ++)
    {
        if(circles[i].flag == 0) continue;
        num++;
        // 8 points
        Vector_2 add[8] = {Vector_2(0, circles[i].radius), Vector_2(circles[i].radius/1.414,circles[i].radius/1.414), 
                          Vector_2(circles[i].radius, 0), Vector_2(circles[i].radius/1.414, -circles[i].radius/1.414), 
                          Vector_2(0, -circles[i].radius), Vector_2(-circles[i].radius/1.414, -circles[i].radius/1.414), 
                          Vector_2(-circles[i].radius, 0), Vector_2(-circles[i].radius/1.414, circles[i].radius/1.414)};
        for(int j = 0; j < 8; j++)
            vertexes.push_back(circles[i].center + add[j]);
        for(int j = 8*i + 0; j < 8*i + 8 -1; j++)
            ofss << setiosflags(ios::fixed) << setprecision(8) << vertexes[j].x() << " " << vertexes[(j+1)%(8*(i+1))].x() << " " << vertexes[j].y() << " " << vertexes[(j+1)%(8*(i+1))].y() << "\n";
        ofss << setiosflags(ios::fixed) << setprecision(8) << vertexes[8*i + 7].x() << " " << vertexes[8*i + 0].x() << " " << vertexes[8*i + 7].y() << " " << vertexes[8*i + 0].y() << "\n";
    }
    LOG(INFO) << num << " cylinders are remained."; 
    ofss.close();

    // generate height according to floor.off ceiling.off
    vector<Point_2> height(vertexes.size()); // {max, min}
    ifstream floor(wdir + "floor.off");
    ifstream ceil(wdir + "ceiling.off");
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
    for(int i = 0; i < vertexes.size(); i++){
       double hc = 1e9, hf = 1e9;
       for(int j = 0; j < polysf.size(); j++){
         if(isin(pointsf, polysf[j], vertexes[i], hf))
            break;
       }
       for(int j = 0; j < polysc.size(); j++){
         if(isin(pointsc, polysc[j], vertexes[i], hc)){
            break;
         }
       }
       if(hc == 1e9 || hf == 1e9){
            LOG(INFO) << vertexes[i] << " Cylinder is failed.";
            return false;
       }
       height[i] = Point_2(hc, hf);
    }
/*
    // generate height according to floor_vec.ply ceiling_vec.ply
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
        double hc = 1e9, hf = 1e9;
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
*/
  // obtain height for each vertex, 8 vertexes -> 1 cylinder
  vector<Point_3> points;
  vector<vector<int>> polygons;
  for(int i = 0; i < vertexes.size(); i=i+8){
    for(int j = 0; j < 8; j++){
      Point_3 p1 = Point_3(vertexes[i+j].x(), vertexes[i+j].y(), height[i+j].x());
      Point_3 p2 = Point_3(vertexes[i+j].x(), vertexes[i+j].y(), height[i+j].y());
      Point_3 p3 = Point_3(vertexes[i+(j+1)%8].x(), vertexes[i+(j+1)%8].y(), height[i+(j+1)%8].y());
      Point_3 p4 = Point_3(vertexes[i+(j+1)%8].x(), vertexes[i+(j+1)%8].y(), height[i+(j+1)%8].x());
      int id1, id2, id3, id4;
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
/*
  // generate cylinder_vec.ply
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
*/
  return true;
}

bool ransac_detection_p(Pwn_vector& points, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, vector<Line>& lines, string wdir)
{
    // timer
    clock_t start, mid, end;
    start = clock(); 
    
    Efficient_ransac ransac;
    ransac.set_input(points);
    ransac.add_shape_factory<Plane>(); // rgister shapes for detection 
    Efficient_ransac::Parameters parameters; // sets parameters for shape detection.
    parameters.probability = probability; 
    parameters.min_points = min_points;
    parameters.epsilon = epsilon;  
    parameters.normal_threshold = normal_threshold;   
    if(cluster_epsilon)
      parameters.cluster_epsilon = cluster_epsilon; 

    // detects shapes
    for(int i = 0; i < 20; i++){
        ransac.detect(parameters);
        if(ransac.shapes().end() - ransac.shapes().begin() > 0) 
            continue;
    }
    LOG(INFO) << ransac.shapes().end() - ransac.shapes().begin() << " detected planes, "
      << ransac.number_of_unassigned_points()
      << " unassigned points." ;

    if (ransac.shapes().end() - ransac.shapes().begin() == 0)
        return 0;

    Efficient_ransac::Plane_range planes = ransac.planes();
    double re_angle = cm::get_config().get<double>("facade.ransac.angle");
    LOG(INFO) << "ransac regularization angle: " << re_angle;
                    
    // regularize detected planes.
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
                            0.02);   // tolerance for parallelism / orthogonality
    int num = 0;
    Efficient_ransac::Plane_range::iterator it = planes.begin();
    while (it != planes.end()) {
        if (Plane* plane = dynamic_cast<Plane*>(it->get()))
        {
            Vector_3 normal = plane->plane_normal(); 
            if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) > cos_angle) 
            {
                it++; continue;
            } // remove non-facade (default threshold = 80)
            // store 2D projection lines and their supporting points
            PtoS(plane, points, num, lines);
            num++;
        }  
        // proceeds with next detected shape.
        it++;
    }
/*
    // the second detection
    Efficient_ransac ransac2;
    Pwn_vector points2;
    for(auto &p:ransac.indices_of_unassigned_points())
        points2.push_back(points[p]);
    ransac2.set_input(points2);
    ransac2.add_shape_factory<Plane>();
    parameters.min_points = 50;
    ransac2.detect(parameters);
    Efficient_ransac::Plane_range planes2 = ransac2.planes();
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
                            0.02);   // tolerance for parallelism / orthogonality          
    it = planes2.begin();
    while (it != planes2.end()) 
    {
        if (Plane* plane = dynamic_cast<Plane*>(it->get()))
        {
            Vector_3 normal = plane->plane_normal();
            if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) > cos_angle) {it++; continue;} // remove non-facade (theshold=80)
            PtoS(plane, points2, num, lines);
            num++;
        }  
        it++;
    }
    LOG(INFO) << ransac2.shapes().end() - ransac2.shapes().begin() << " detected planes";
*/
    LOG(INFO) << num << " lines are generated. ";
    end = clock();
    LOG(INFO) << "Facade detection and projection time: " << (float)(end-start)/CLOCKS_PER_SEC << "s";
    if(num != lines.size())
        LOG(INFO) << "The number of lines is not match. (" << num << ", " << lines.size() << ")";
    // record 2d segments and facade points
    Pwn_vector facade_points;
    ofstream file(wdir + "2dsegments", ios::out);
    for(int i = 0; i < lines.size(); i++)
    {
        file << setiosflags(ios::fixed) << setprecision(8) << lines[i].s.x() << " " << lines[i].t.x() << " "<< lines[i].s.y() << " " << lines[i].t.y() << endl;
        for(auto& p: lines[i].pointset)
            facade_points.push_back(Point_with_normal({p,Vector_3(0,1,0)}));
        
    }
    file.close();
    PLYPointSave(wdir + "facade_points_ransac.ply", facade_points);
    return 1;
}

bool ransac_detection_c(Pwn_vector& points, float radius, float probability, int min_points, float epsilon, float cluster_epsilon, float normal_threshold, float cos_angle, vector<Circle>& circles, string wdir)
{
    // timer
    clock_t start, end;
    start = clock(); 
    
    Efficient_ransac ransac;
    ransac.set_input(points);
    ransac.add_shape_factory<Cylinder>();
    Efficient_ransac::Parameters parameters;
    parameters.probability = probability; 
    parameters.min_points = min_points;
    parameters.epsilon = epsilon;  
    parameters.normal_threshold = normal_threshold;   
    if (cluster_epsilon)
        parameters.cluster_epsilon = cluster_epsilon; 

    ransac.detect(parameters);
    LOG(INFO) << ransac.shapes().end() - ransac.shapes().begin() << " detected cylinders, "
        << ransac.number_of_unassigned_points()
        << " unassigned points." ;

    if (ransac.shapes().end() - ransac.shapes().begin() == 0)
    {
        LOG(INFO) << "No cylinder is detected.";
        return 0;
    }

    Efficient_ransac::Shape_range cylinders = ransac.shapes();
    int num = 0;
    Efficient_ransac::Shape_range::iterator it = cylinders.begin();
    while (it != cylinders.end()) {
        if (Cylinder* cylinder = dynamic_cast<Cylinder*>(it->get()))
        {
            Line_3 axis = cylinder->axis();   
            if (abs(axis.to_vector().z()/sqrt(axis.to_vector().x()*axis.to_vector().x()+axis.to_vector().y()*axis.to_vector().y()+axis.to_vector().z()*axis.to_vector().z())) < cos_angle) {it++; continue;}
            double mid = 0;
            int num_c = 0;// find the middle height!
            for(auto p:cylinder->indices_of_assigned_points()){
                mid += points[p].first.z();
                num_c++;
            }    
            Plane_3 mid_p(Point_3(0,0,mid/num_c), Vector_3(0,0,1)); 
            auto result = CGAL::intersection(mid_p, axis);    
            if (result)
                if(const Point_3* pt = boost::get<Point_3>(&*result))	
                    circles.push_back(Circle(Point_2(pt->x(), pt->y()), cylinder->radius(),cylinder->indices_of_assigned_points().end()-cylinder->indices_of_assigned_points().begin(), 0, 0));
                else
                    LOG(INFO) << "Cylinder detection is failed!";       
        }  
        it++;
        num++;
    }
    LOG(INFO) << num << " cylinders are generated. ";
    
    // exclude intersection circles
    Clean(circles, radius);
    // valid circles produce eight planes
    if(!Toplane(circles, wdir))
        return 0;
    // save cylinder file
    string fname = wdir + "cylinder";
    ofstream ofs(fname);
    for(int i = 0; i < circles.size(); i++)
    {
        if (circles[i].flag == 0) continue;
        ofs << circles[i].center.x() << " " << circles[i].center.y() <<" "<< circles[i].radius << "\n";
    }
    ofs.close();

    end = clock();
    LOG(INFO) << "Cylinder reconstruction time: " << (float)(end-start)/CLOCKS_PER_SEC << "s";
    return 1;
}

bool ransac_detection_floor(Pwn_vector& points, vector<Plane_3>& floors, vector<Line_2>& intersects, const cm::Config& config, string wdir)
{
    // timer
    clock_t start, end;
    start = clock(); 
    
    vector<Plane_3> verticals; // record vertical floors
    Efficient_ransac ransac;
    ransac.set_input(points);
    ransac.add_shape_factory<Plane>();
    Efficient_ransac::Parameters parameters;
    parameters.probability = *(config.get_optional<double>("floor.ransac.probability")); 
    parameters.min_points = *(config.get_optional<double>("floor.ransac.min_points"));
    parameters.epsilon = *(config.get_optional<double>("floor.ransac.epsilon"));  
    parameters.normal_threshold = *(config.get_optional<double>("floor.ransac.normal_threshold"));   
    if (*(config.get_optional<double>("floor.ransac.cluster_epsilon")))
        parameters.cluster_epsilon = *(config.get_optional<double>("floor.ransac.cluster_epsilon")); 

    ransac.detect(parameters);
    if (ransac.shapes().end() - ransac.shapes().begin() == 0)
        return 0;
    Efficient_ransac::Plane_range planes = ransac.planes();
    CGAL::regularize_planes(points,
                            Point_map(),
                            planes,
                            CGAL::Shape_detection::Plane_map<Traits>(),
                            CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                            true, true, false, true, 10, 0.10); 

    int num=0;
    Efficient_ransac::Plane_range::iterator it = planes.begin();
    while (it != planes.end()) 
    {
        if (Plane* plane = dynamic_cast<Plane*>(it->get()))
        {
            Vector_3 normal = plane->plane_normal();
            if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) < 1 - *(config.get_optional<double>("floor.ransac.cos_angle"))) 
            {
                verticals.push_back(Plane_3(normal.x(), normal.y(), normal.z(), plane->d()));
            }      
            if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) < *(config.get_optional<double>("floor.ransac.cos_angle"))) {it++; continue;}     
            // assume that the plane passes through the first supporting point
            if(normal.z() < 0)
                floors.push_back(Plane_3((points.begin() + *(plane->indices_of_assigned_points().begin()))->first, Vector_3(0,0,-1)));
            else
                floors.push_back(Plane_3((points.begin() + *(plane->indices_of_assigned_points().begin()))->first, Vector_3(0,0,1)));
            num++;
        }  
        it++;
    }
    LOG(INFO) << num << " horizontal floor plane(s) are deteced.";
    if(num == 1) { // only one horizontal floor
        end = clock();
        LOG(INFO) << "Floor detection time: " << (float)(end-start)/CLOCKS_PER_SEC << "s";
        return true;
    }
    else // multi floor, may exist steps
    { 
        Efficient_ransac ransac2;
        Pwn_vector points2;
        for(auto &p:ransac.indices_of_unassigned_points())
            points2.push_back(points[p]);
        LOG(INFO) << "detect vertical planes from " << points2.size() << " unassigned points." ;
        ransac2.set_input(points2);
        ransac2.add_shape_factory<Plane>();
        Efficient_ransac::Parameters parameters2;
        parameters2.probability = *(config.get_optional<double>("floor.ransac2.probability")); 
        parameters2.min_points = *(config.get_optional<double>("floor.ransac2.min_points"));
        parameters2.epsilon = *(config.get_optional<double>("floor.ransac2.epsilon"));  
        parameters2.normal_threshold = *(config.get_optional<double>("floor.ransac2.normal_threshold"));   
        if (*(config.get_optional<double>("floor.ransac2.cluster_epsilon")))
            parameters2.cluster_epsilon = *(config.get_optional<double>("floor.ransac2.cluster_epsilon")); 

        ransac2.detect(parameters2);
        Efficient_ransac::Plane_range planes2 = ransac2.planes();
        CGAL::regularize_planes(points,
                                Point_map(),
                                planes,
                                CGAL::Shape_detection::Plane_map<Traits>(),
                                CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                                true, true, false, true, 5, 0.10); 
        it = planes2.begin();
        while (it != planes2.end()) 
        {
            if (Plane* plane = dynamic_cast<Plane*>(it->get()))
              {
                Vector_3 normal = plane->plane_normal();
                if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) < 1 - *(config.get_optional<double>("floor.ransac2.cos_angle"))) 
                {
                    // centroid points
                    double x,y,z;
                    x=y=z=0;
                    int num = plane->indices_of_assigned_points().end() - plane->indices_of_assigned_points().begin(); 
                    for(auto it2 = plane->indices_of_assigned_points().begin(); it2 != plane->indices_of_assigned_points().end(); it2++)
                    { 
                        x += points2[*it2].first.x(); 
                        y += points2[*it2].first.y(); 
                        z += points2[*it2].first.z();
                    }
                    verticals.push_back(Plane_3(Point_3(x/num, y/num, z/num), Vector_3(normal.x(), normal.y(), 0)));
                }
              }  
            it++;
        }
        LOG(INFO) << verticals.size() << " vertical plane(s) are deteced.";

        if(verticals.size() == 0)
        {
            return true;
        }
        // vertical and horizantal floors intersect
        for(auto& h: floors)
        {
            if(h == *floors.begin()) continue; // the begin is the blank plane
            for(auto& v: verticals)
            {
                auto result = CGAL::intersection(h, v);
                if(Line_3* l = boost::get<Line_3>(&*result))
                {
                    Line_2 line(Point_2(l->projection(Point_3(0,0,0)).x(), l->projection(Point_3(0,0,0)).y()), Vector_2(l->to_vector().x(), l->to_vector().y()));
                    auto it = find(intersects.begin(), intersects.end(), line);
                    auto it2 = find(intersects.begin(), intersects.end(), line.opposite());
                    if(it == intersects.end() && it2 == intersects.end())
                        intersects.push_back(line);            
                }
            }
        }
      end = clock();
      LOG(INFO) << "Floor detection time: " << (float)(end-start)/CLOCKS_PER_SEC << "s.";
      return 1; 
    }   
}

bool ransac_detection_ceiling(Pwn_vector& points,  vector<Plane_3>& ceilings, vector<Line_2>& intersects, const cm::Config& config, string wdir)
{
    // timer
    clock_t start, end;
    start = clock(); 
    
    Efficient_ransac ransac;
    ransac.set_input(points);
    ransac.add_shape_factory<Plane>();
    Efficient_ransac::Parameters parameters;
    parameters.probability = *(config.get_optional<double>("ceiling.ransac.probability")); 
    parameters.min_points = *(config.get_optional<double>("ceiling.ransac.min_points"));
    parameters.epsilon = *(config.get_optional<double>("ceiling.ransac.epsilon"));  
    parameters.normal_threshold = *(config.get_optional<double>("ceiling.ransac.normal_threshold"));   
    if (*(config.get_optional<double>("ceiling.ransac.cluster_epsilon")))
      parameters.cluster_epsilon = *(config.get_optional<double>("ceiling.ransac.cluster_epsilon")); 


    ransac.detect(parameters);
    if (ransac.shapes().end() - ransac.shapes().begin() ==0)
        return 0;
    Efficient_ransac::Plane_range planes = ransac.planes();
    CGAL::regularize_planes(points,
                            Point_map(),
                            planes,
                            CGAL::Shape_detection::Plane_map<Traits>(),
                            CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                            true, true, false, true, 10, 0.10); 
    Efficient_ransac::Plane_range::iterator it = planes.begin();
    // record detected celings
    Pwn_vector ceiling_points;
    while (it != planes.end()) 
    {
        if (Plane* plane = dynamic_cast<Plane*>(it->get()))
        {
            Vector_3 normal = plane->plane_normal();
            double x=0, y=0, z=0;
            for(auto p:plane->indices_of_assigned_points())
            {
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
            else if (abs(normal.z()/sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z())) < 1 - *(config.get_optional<double>("ceiling.ransac.cos_angle")))
            {
                ceilings.push_back(Plane_3(Point_3(x/num, y/num, z/num), Vector_3(normal.x(), normal.y(),0)));
            }
            else // slap plane
                ceilings.push_back(Plane_3(normal.x(), normal.y(), normal.z(), plane->d()));
            
          } 
        it++;
    }
    PLYPointSave(wdir + "ceiling_points_ransac.ply", ceiling_points);
         
    // compute intersection, blank plane (i = 0)
    for(int i = 1; i < ceilings.size(); i++)
        for(int j = i+1; j < ceilings.size(); j++)
        {
            auto result = CGAL::intersection(ceilings[i], ceilings[j]);
            if(!result) continue;
            if(const Line_3* l = boost::get<Line_3>(&*result))
            {
                Line_2 line(Point_2(l->projection(Point_3(0,0,0)).x(), l->projection(Point_3(0,0,0)).y()), Vector_2(l->to_vector().x(),l->to_vector().y()));
                auto it = find(intersects.begin(), intersects.end(), line);
                auto it2 = find(intersects.begin(), intersects.end(), line.opposite());
                if(it == intersects.end() && it2 == intersects.end())
                    intersects.push_back(line);  
            }
        }

    end = clock();
    LOG(INFO) << "Ceiling detection time: " << (float)(end-start)/CLOCKS_PER_SEC << "s.";
    return 1;   
}
