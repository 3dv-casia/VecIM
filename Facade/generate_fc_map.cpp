#include "arr.h"

#include <glog/logging.h>
#include <malloc.h>

#include <CGAL/Iso_rectangle_2.h>


bool ready_for_fc(vector<Plane_3>& planes, cv::Mat& hmap, cv::Mat& nmap, vector<Line_2>& intersect, const cm::Config& config,  Point_3 corner, double step, string wdir, int tag)
{   
    const auto file = tag == 1 ? config.get_optional<string>("file.floorname") : config.get_optional<string>("file.ceilingname");
    Pwn_vector points;
    if (!PLYPointLoad3(wdir + *file, points)){
        LOG(INFO) << "Load failed.";
        return false;
    }
    if(tag == 1)
        LOG(INFO)  << "Load finished, "<< points.size() << " floor points. ";
    else
        LOG(INFO)  << "Load finished, "<< points.size() << " ceiling points. ";

    if(tag == 1 && *config.get_optional<bool>("floor.enabled") == true || tag != 1 && *config.get_optional<bool>("ceiling.enabled") == true)
    {
        // push blank plane 
        planes.push_back(Plane_3(Point_3(0,0,corner.z()-10), Vector_3(0,0,1)));
        if(tag == 1)
        {
            // ransac detetcion, only accept horizantal/vertical floors
            // the first detection detects big shape, and the second one detects small vertical floors
            ransac_detection_floor(points, planes, intersect, config, wdir);
            LOG(INFO)  << planes.size()-1 << " floors are detected, " << intersect.size() << " intersections are generated.";
        
        }
        else
        {
            // accept horizantal and slap planes
            ransac_detection_ceiling(points, planes, intersect, config, wdir);
            LOG(INFO)  << planes.size()-1 << " ceilings are detected, "<< intersect.size() << " intersections are generated.";    
        }
        // save planes and intersection
        /* file format:
        plane_num intersection_num
        a b c d
        a b c*/ 
        string fn = tag == 1 ? "floor_intersectf" : "ceiling_intersectc";
        ofstream file(wdir + fn, ios::out);
        file << planes.size() << " " << intersect.size() << "\n"; 
        for(auto& f: planes)
            file << setiosflags(ios::fixed) << setprecision(8) << f.a() << " " << f.b() << " " << f.c() << " " << f.d() << "\n";
        for(auto& s:intersect)
            file << setiosflags(ios::fixed) << setprecision(8) << s.a() << " " << s.b() << " " << s.c() << "\n";
        file.close();
    }
    else
    { // load data from existed file
        string fn = tag == 1 ? "floor_intersectf" : "ceiling_intersectc";
        ifstream file(wdir + fn);
        if(!file.is_open())
        {
            LOG(INFO) << "Failed to load floor_intersectf/ceiling_intersectc file.";
            return false;     
        }
        int num_f, num_inf;
        double x, y, z, d;
        file >> num_f >> num_inf;
        while(num_f-->0){
          file >> x >> y >> z >> d;
          planes.push_back(Plane_3(x,y,z,d));
        }
        while(num_inf-->0){
          file >> x >> y >> z;
          intersect.push_back(Line_2(x,y,z));
        }
     }
    
    // hmap
    cv::Mat num_map = cv::Mat::zeros(hmap.rows, hmap.cols, CV_64FC2); // channels: number, z height
    for(int  i = 0; i < points.size(); i++)
    {
        int r = int((corner.y() - points[i].first.y()) / step);
        int c = int((points[i].first.x() - corner.x()) / step);
        if( r < 0 || r > hmap.rows || c < 0 || c > hmap.cols) continue;
        num_map.at<cv::Vec2d>(r, c)[0]++;
        num_map.at<cv::Vec2d>(r, c)[1] += points[i].first.z();    
    }
    for(int r = 0; r < hmap.rows; r++)
        for(int c = 0; c < hmap.cols; c++)
        {
            hmap.at<cv::Vec3d>(r, c)[0] = corner.x() + (c + 0.5)*step;
            hmap.at<cv::Vec3d>(r, c)[1] = corner.y() - (r + 0.5)*step;
            if(num_map.at<cv::Vec2d>(r, c)[0] == 0) // no points, assign blank height: zmin - 10
                hmap.at<cv::Vec3d>(r, c)[2] = corner.z() - 10;
            else
                hmap.at<cv::Vec3d>(r, c)[2] = num_map.at<cv::Vec2d>(r, c)[1] / num_map.at<cv::Vec2d>(r, c)[0];
            nmap.at<int>(r, c) = num_map.at<cv::Vec2d>(r, c)[0];                     
        }
    // save hmap to height map, scale to [0, 255]
    double zmin = 0.0, zmax = 0.0;
    vector<cv::Mat> channels;
    cv::split(hmap, channels); // split z height value
    cv::minMaxLoc(channels[2], &zmin, &zmax);
    LOG(INFO) << "height map: " << zmin << " zmin, " << zmax << " zmax.";
    const double zrange = zmax - zmin;
    cv::Mat img;
    channels[2].convertTo(img, CV_8U, 255.0 / zrange, -zmin * 255.0 / zrange);
    string fn = tag == 1 ? "floor_height.png" : "ceilings_height.png";
    cv::imwrite(wdir + fn, img);
    return true;
}

bool generate_fcf(vector<Segment_2>& segments, const Point_3 cornerz, const Point_2 length, const cm::Config& config, string wdir)
{
    I2E to_exact;
    // generate grid arr
    double step = *(config.get_optional<double>("arr.step"));
    int rows = std::ceil(length.y() / step);
    int cols = std::ceil(length.x() / step);
    const Point_2 corner(cornerz.x(),cornerz.y());
    LOG(INFO) << rows << " rows, " << cols << " cols.";

    // for building class
    cv::Mat hfmap(rows, cols, CV_64FC3); // channels: real x position, real y position, average height
    cv::Mat hcmap(rows, cols, CV_64FC3);
    cv::Mat nfmap(rows, cols, CV_32S); // point number
    cv::Mat ncmap(rows, cols, CV_32S);
    vector<Plane_3> floors;
    vector<Line_2> intersectf;
    vector<Plane_3> ceilings;
    vector<Line_2> intersectc;

    if(!ready_for_fc(floors, hfmap, nfmap, intersectf, config, cornerz, step, wdir, 1) // tag = 1 for floors
    || !ready_for_fc(ceilings, hcmap, ncmap, intersectc, config, cornerz, step, wdir, 0) // tag = 0 for ceilings
    || floors.size() == 0 || ceilings.size() == 0){ 
        LOG(INFO) << "floor/ceiling detection is failed.";
        return false;
    }

    // bbox
    const Point_2 corner2(corner.x()+length.x(), corner.y()-length.y());
    CGAL::Iso_rectangle_2<K_epec> bbox(to_exact(corner), to_exact(corner2)); 
    // push all segment together
    vector<K_epec::Segment_2> total_segs;

    // 1. push intersection (segments)
    for(auto &l:intersectf)
        if(auto result = CGAL::intersection(bbox, to_exact(l)))
            if(K_epec::Segment_2* rl = boost::get<K_epec::Segment_2>(&*result))
                total_segs.push_back(*rl);
        
    for(auto &l:intersectc)
        if(auto result = CGAL::intersection(bbox, to_exact(l)))
            if(K_epec::Segment_2* rl = boost::get<K_epec::Segment_2>(&*result))
                total_segs.push_back(*rl);
          
    // 2. push extended facade segments
    vector<K_epec::Segment_2> selected;
    double extend = *(config.get_optional<double>("arr.extend"));
    double dialog = CGAL::sqrt(length.x()*length.x() + length.y()*length.y());
    for(auto& s: segments){
        selected.push_back(to_exact(s));
        const auto mid = CGAL::midpoint(selected.back().source(), selected.back().target());
        const double len = CGAL::sqrt(s.squared_length());
        const auto dir = selected.back().to_vector()/len;
        total_segs.push_back({mid + (len/2 + extend*dialog)*dir, mid - (len/2 + extend*dialog)*dir});
    }
    
    //3. add bounding segments to close the arrangement
    total_segs.push_back(K_epec::Segment_2(bbox.vertex(0),bbox.vertex(1))); 
    total_segs.push_back(K_epec::Segment_2(bbox.vertex(1),bbox.vertex(2)));
    total_segs.push_back(K_epec::Segment_2(bbox.vertex(2),bbox.vertex(3)));
    total_segs.push_back(K_epec::Segment_2(bbox.vertex(3),bbox.vertex(0)));

    building fc(wdir, bbox, total_segs, selected, corner, step, hfmap, hcmap, nfmap, ncmap, floors, ceilings, cornerz.z() - 10); // min height = zmin - 10
    if(!fc.segment_arrangement_modeling())
        return false;
    return true;
}

bool generate_fc(vector<Segment_2> segments, double top, double bottom, string wdir)
{
   map<Point_2,int> points;
   vector<Point_2> single_points;
   vector<int> poly; // assume one object
   
   double xmin = 1e9; double ymin = 1e9;
   double xmax = -1e9; double ymax = -1e9;

   for(auto s:segments)
   {
       if (points.count(s.source()))
           points[s.source()]++;
       else
           points.insert(make_pair(s.source(),1));
       if (points.count(s.target()))
           points[s.target()]++;
       else
           points.insert(make_pair(s.target(),1));
       
   }
    int num_points = points.size();
    map<Point_2,int>::iterator iter;
    for (iter=points.begin(); iter!=points.end(); iter++)
    {
       if(iter->second == 1)
         single_points.push_back(iter->first);
       if (iter->first.x() < xmin) xmin = iter->first.x();
       if (iter->first.x() > xmax) xmax = iter->first.x();
       if (iter->first.y() < ymin) ymin = iter->first.y();
       if (iter->first.y() > ymax) ymax = iter->first.y();  
	}
    // clear
    map<Point_2,int>().swap(points);
    malloc_trim(0);

    vector<Point_2> vertexs;
    vector<vector<int>> indexs(num_points);
    int *tag = new int[num_points]();

    if(single_points.size()%2==1){
        LOG(INFO) << "Facade error.";
        return false;
    }
    vector<Point_2> xmins;	vector<Point_2> xmaxs;
	vector<Point_2> ymins; 	vector<Point_2> ymaxs;
	for(auto p: single_points){
		if (p.x() == xmin) {xmins.push_back(p);continue;}
		if (p.x() == xmax) {xmaxs.push_back(p);continue;}
		if (p.y() == ymin) {ymins.push_back(p);continue;}
		if (p.y() == ymax) {ymaxs.push_back(p);continue;}
	}

	if(xmins.size()%2 != 0 || ymins.size()%2 != 0 || xmaxs.size()%2 != 0 || ymaxs.size()%2 != 0 ){
		LOG(INFO) << "Point error.";
        return false;
	}
	for(int i = 0; i < xmins.size(); i=i+2)
       segments.push_back(Segment_2(xmins[i],xmins[i+1]));
	for(int i = 0; i < ymins.size(); i=i+2)
       segments.push_back(Segment_2(ymins[i],ymins[i+1]));
	for(int i = 0; i < xmaxs.size(); i=i+2)
       segments.push_back(Segment_2(xmaxs[i],xmaxs[i+1]));
	for(int i = 0; i < ymaxs.size(); i=i+2)
       segments.push_back(Segment_2(ymaxs[i],ymaxs[i+1]));

    for(auto s:segments)
    {
       int s1,s2;
       vector<Point_2>::iterator it = find(vertexs.begin(), vertexs.end(),s.source());
       if(it == vertexs.end()){
            vertexs.push_back(s.source());
            s1 = vertexs.size()-1;
       }
       else
           s1 = it - vertexs.begin();     

       it = find(vertexs.begin(), vertexs.end(),s.target());
       if(it == vertexs.end()){
            vertexs.push_back(s.target());
            s2 = vertexs.size()-1;
       }
       else
           s2 = it - vertexs.begin();
        if(s1 >= num_points || s1 >= num_points)
           return false;
        indexs[s1].push_back(s2);
        indexs[s2].push_back(s1);
    }
    
    int index = 0;
    poly.push_back(index);
    tag[index] = 1;
    while(tag[indexs[index][0]] == 0 || tag[indexs[index][1]] == 0){
       if(tag[indexs[index][0]] == 0){
           tag[indexs[index][0]] = 1;
           index = indexs[index][0];
           poly.push_back(index);
           continue;
       }
        if(tag[indexs[index][1]] == 0){
           tag[indexs[index][1]] = 1;
           index = indexs[index][1];
           poly.push_back(index);
           continue;
       }

    }
    delete [] tag;

    // write off files
    ofstream ofs1(wdir + "floor.off");
    ofstream ofs2(wdir + "ceiling.off");
    ofs1 << "COFF\n" << vertexs.size() << ' ' << 1 << ' ' << "0\n";
    ofs2 << "COFF\n" << vertexs.size() << ' ' << 1 << ' ' << "0\n";
    for (const auto &p : vertexs){
        ofs1 << p.x() << ' ' << p.y() << ' ' << bottom << '\n';
        ofs2 << p.x() << ' ' << p.y() << ' ' << top << '\n';
    }
    ofs1 << poly.size();
    ofs2 << poly.size();
    for (int i = 0; i < poly.size(); i++){
        ofs1 << ' ' << poly[i];
        ofs2 << ' ' << poly[poly.size()-1 -i];
    }
    ofs1 << '\n';
    ofs2 << '\n';
    return true;
}

bool generate_map(string wdir, double size)
{
    ifstream ifs(wdir + "floor.off"); 
    if(!ifs){
      LOG(INFO) << "floor.off does not exit.";
      return false;
    }
    string f;
    int np,nf,c;
    vector<Point_2> points;
    double x, y, z;
    vector<Segment_2> poly; 
    double xmin = 1e9; double ymin = 1e9;
    double xmax = -1e9; double ymax = -1e9;

    ifs >> f >> np >> nf >> c;
    for(int i = 0; i < np; i++){
       ifs >> x >> y >> z;
       points.push_back(Point_2(x,y));
       if (x < xmin) xmin = x;
       if (x > xmax) xmax = x;
       if (y < ymin) ymin = y;
       if (y > ymax) ymax = y;       
    }
    for(int i = 0; i < nf; i++){
       int nfp; 
       int s, t, ss;
       ifs >> nfp;
       ifs >> s;
       ss = s;
       for(int j = 1; j < nfp; j++){
           ifs >> t;
           poly.push_back(Segment_2(points[s], points[t]));
           s = t;
       }
       poly.push_back(Segment_2(points[s], points[ss]));
    }

    int rows = ceil((ymax - ymin)/size);
    int cols = ceil((xmax - xmin)/size);
    cv::Mat maps = cv::Mat::zeros(rows, cols, CV_8UC1);
    for(int i = 0; i < rows; i++)
      for(int j = 0; j < cols; j++){
          Point_2 p(xmin + (j+0.5)*size, ymax - (i+0.5)*size);
          int t = 0;
          for(auto s:poly) // up射线
             if((s.source().x() > p.x()) != (s.target().x() > p.x())){
              double x1 = s.source().x(); double y1 = s.source().y();
              double x2 = s.target().x(); double y2 = s.target().y();
              double l = (x1-x2)*p.y();
              double r = (y1-y2)*p.x() + x1*y2-x2*y1;
              if (x1 > x2 && l < r || x1 < x2 && l > r)
                t = !t;
             }
          if (t % 2 == 1) 
            maps.at<uchar>(i, j) = 1;
      }

    // exclude cylinder
    ifstream ifs2(wdir + "cylinder");
    if(!ifs2){
        LOG(INFO) << "cylinder does not exit.";
        return false;
    }    
    vector<Point_3> centers;
    double radious;
    while(ifs2 >> x >> y >> radious)
       centers.push_back(Point_3(x, y, radious));
    for(int i = 0; i < centers.size(); i++){
        Point_2 lc = Point_2(centers[i].x() - centers[i].z(), centers[i].y() - centers[i].z()); // left bttom corner
        Point_2 rc = Point_2(centers[i].x() + centers[i].z(), centers[i].y() + centers[i].z()); // right up corner
        int cols1 = int((lc.x() - xmin)/size);
        int rows2 = int((ymax - lc.y())/size);
        int cols2 = int((rc.x() - xmin)/size);
        int rows1 = int((ymax - rc.y())/size);
        for(int i = rows1; i <= rows2; i++)
          for(int j = cols1; j <= cols2; j++)
             maps.at<uchar>(i, j) = 0;
    }
    maps *= 255;
    const string fname = wdir + "freespace.png";
    cv::imwrite(fname, maps);
    return true;
}


bool generate_off_modeling(string wdir)
{
    vector<Point_3> points;
    vector<vector<int>> polygons;
    string str;
    int np, nf, n;
    double x, y, z;
    int num_p = 0;
    vector<Point_3>::iterator it;

    string fname[4] = {"ceiling.off", "floor.off", "facade.off", "cylinder.off"};
    for(int k = 0; k < sizeof(fname)/sizeof(fname[0]); k++){
        ifstream ifs(wdir + fname[k]);
        ifs >> str;
        ifs >> np >> nf >> n;
        int *order = new int[np];
        for(int i = 0; i < np; i++){
            ifs >> x >> y >> z;
            it = find(points.begin(), points.end(), Point_3(x, y, z));
            if(it == points.end()){
                points.push_back(Point_3(x, y, z));
                order[i] = num_p++;
            }
            else
                order[i] = it - points.begin();      
        }

        for(int i = 0; i < nf; i++){
            ifs >> n;
            int index;
            vector<int> poly;
            for(int i = 0; i < n; i++){
                ifs >> index;
                poly.push_back(order[index]);
            }   
            polygons.push_back(poly);
        }  
        ifs.close();
        delete [] order;
    }
    if(num_p != points.size())
       return false;

    // merge
    ofstream ofs(wdir + "modeling.off");
    ofs << "COFF\n" << points.size() << ' ' << polygons.size() << ' ' << "0\n";
     for (const auto &p : points)
       ofs << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
     for(auto poly: polygons){
        ofs << poly.size();
        for (const auto &p : poly){
            ofs << ' ' << p;
        }
        ofs << '\n';
     }
    return true;
}

bool generate_ply_modeling(string wdir)
{
    vector<Point_3> points;
    vector<vector<int>> polygons;
    string str;
    int np, nf, n;
    double x, y, z;
    int num_p = 0;
    vector<Point_3>::iterator it;

    string fname[4] = {"ceiling_vec.ply", "floor_vec.ply", "facade_vec.ply", "cylinder_vec.ply"};
    Mesh<Point_3> mesh;
    for(int k = 0; k < sizeof(fname)/sizeof(fname[0]); k++){
        ifstream ifs(wdir + fname[k]);
        if(!ifs) continue;
        ifs.close();
        if(k == 0) PLYMeshLoad(wdir+fname[k], mesh);
        else{
            Mesh<Point_3> sub_mesh;
            PLYMeshLoad(wdir+fname[k], sub_mesh);       
            for(auto v: sub_mesh.vertices){
               auto it = find(mesh.vertices.begin(),mesh.vertices.end(),v);
               if(it == mesh.vertices.end()) mesh.vertices.push_back(v);
            }
            for(auto f: sub_mesh.faces){
               Point_3 v0 = sub_mesh.vertices[f[0]];
               Point_3 v1 = sub_mesh.vertices[f[1]];
               Point_3 v2 = sub_mesh.vertices[f[2]];
               auto it0 = find(mesh.vertices.begin(),mesh.vertices.end(),v0);
               auto it1 = find(mesh.vertices.begin(),mesh.vertices.end(),v1);
               auto it2 = find(mesh.vertices.begin(),mesh.vertices.end(),v2);
                if(it0 == mesh.vertices.end() || it1 == mesh.vertices.end() || it2 == mesh.vertices.end())
                   LOG(INFO) << "Merge failed.";
                else{
                   mesh.faces.push_back(Point_3(it0-mesh.vertices.begin(), it1-mesh.vertices.begin(), it2-mesh.vertices.begin()));
                }
            }
        }
    }
    PLYMeshSave(wdir + "modeling.ply", mesh, 1);
    return true;
}