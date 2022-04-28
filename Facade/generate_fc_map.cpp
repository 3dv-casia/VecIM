#include "arr.h"

#include <glog/logging.h>
#include <malloc.h>

#include <CGAL/Iso_rectangle_2.h>
/*
K_epec::Segment_2 merge_seg(vector<K_epec::Segment_2>& segments, vector<int>& idxs){
    if(idxs.size() == 1)
      return segments[idxs[0]];
    double len = 0;
    double px = 0;
    double py = 0;
    for(auto& i: idxs){
        len += sqrt(CGAL::to_double(segments[i].squared_length()));
        px += CGAL::to_double((CGAL::midpoint(segments[i].min(), segments[i].max())).x()) * sqrt(CGAL::to_double(segments[i].squared_length()));
        py += CGAL::to_double((CGAL::midpoint(segments[i].min(), segments[i].max())).y()) * sqrt(CGAL::to_double(segments[i].squared_length()));
    }
    K_epec::Line_2 l(K_epec::Point_2(px/len, py/len), segments[idxs[0]].to_vector());
    K_epec::Point_2 p1 = l.projection(segments[idxs[0]].min());
    K_epec::Point_2 p2 = l.projection(segments[idxs[0]].max());
    for(int i = 1; i < idxs.size(); i++){
        p1 = p1 < l.projection(segments[idxs[i]].min())? p1 : l.projection(segments[idxs[i]].min());
        p2 = p2 > l.projection(segments[idxs[i]].max())? p2 : l.projection(segments[idxs[i]].max());
    }
    return K_epec::Segment_2(p1, p2);
}

void merge(vector<K_epec::Segment_2>& segments, const Point_2& corner){
    constexpr double EPSILON_ANGLE = 0.99984769515; // cos1 degree
    constexpr double EPSILON_DISTANCE = 1.0e-8;
    typedef CGAL::Cartesian_converter<K_epec, K> To_geom;
    To_geom to_geom;
    typedef CGAL::Cartesian_converter<K, K_epec> To_geom2;
    To_geom2 to_geom2;
    sort(segments.begin(), segments.end(), [](K_epec::Segment_2& s1, K_epec::Segment_2& s2){return CGAL::to_double(s1.squared_length()) > CGAL::to_double(s2.squared_length());});
    vector<int> visit(segments.size(), 0);
    struct parallel_dis{
        double dis = 0;
        int idx = -1;
        parallel_dis(double dis_, int idx_):dis(dis_), idx(idx_){}
        bool operator < (const parallel_dis& s) const
        {
            return dis > s.dis;
        }
    };
    // regular then cluster parallel segments into different sets
    vector<K_epec::Segment_2> re_segments;
    for(int i = 0; i < segments.size(); i++){
        if(visit[i] == 1)
           continue;
        visit[i] = 1;
        set<parallel_dis> sub_p;
        K_epec::Segment_2 s1 = segments[i];
        Point_2 p1 = to_geom(s1.min());
        Point_2 p2 = to_geom(s1.max());
        Vector_2 v1(p1, p2);
        double len1 = CGAL::sqrt(v1.squared_length());
        sub_p.insert(parallel_dis(CGAL::squared_distance(corner, Line_2(p1, p2)), i));
        for(int j = i+1; j < segments.size(); j++){
            if(visit[j] == 1)
               continue;
            K_epec::Segment_2& s2 = segments[j];
            Point_2 p3 = to_geom(s2.min());
            Point_2 p4 = to_geom(s2.max());
            Vector_2 v2(p3, p4);
            double len2 = CGAL::sqrt(v2.squared_length());
            double cos = v1 * v2 / (len1 * len2);
            double d = robust_squared_distance(p1, p2, p3, p4);
            if(cos >= EPSILON_ANGLE && d <= EPSILON_DISTANCE){ 
                Point_2 mid = CGAL::midpoint(p3, p4);
                Line_2 l(mid, v1);
                s2 = to_geom2(Segment_2(l.projection(p3), l.projection(p4))); // adjust s2 to parallel to s1
                sub_p.insert(parallel_dis(CGAL::squared_distance(corner, l), j));
                visit[j] = 1;
            }
        }
        if(sub_p.size() == 1)
          re_segments.push_back(segments[sub_p.begin()->idx]);
        else{
            for(auto it = sub_p.begin(); it != sub_p.end();){
                auto prev = it++;
                vector<int> merge_idx;
                merge_idx.push_back(prev->idx);
                while(it != sub_p.end()){
                    double d = robust_squared_distance(to_geom(segments[prev->idx].min()), to_geom(segments[prev->idx].max()), to_geom(segments[it->idx].min()), to_geom(segments[it->idx].max()));
                    if(d == 0)
                       LOG(INFO) << "TEST: " << segments[prev->idx] << " " << segments[it->idx];
                    if(d <= EPSILON_DISTANCE){
                        merge_idx.push_back(it->idx);
                        it++;
                    }
                    else
                       break;
                }
                re_segments.push_back(merge_seg(segments, merge_idx));
            }
        }
    }
    segments.clear();
    segments = re_segments;
}

void merge(vector<K_epec::Segment_2>& segments, const Point_2& corner){
    constexpr double EPSILON_ANGLE = 0.99984769515; // cos1 degree
    constexpr double EPSILON_DISTANCE = 1.0e-8;
    constexpr double EPSILON_GENERAL = 1.192092896e-07;
    bool merge = false;
    typedef CGAL::Cartesian_converter<K_epec, K> To_geom;
    To_geom to_geom;
    typedef CGAL::Cartesian_converter<K, K_epec> To_geom2;
    To_geom2 to_geom2;
    //test
    sort(segments.begin(), segments.end(), [](K_epec::Segment_2& s1, K_epec::Segment_2& s2){return CGAL::to_double(s1.squared_length()) > CGAL::to_double(s2.squared_length());});
        LOG(INFO) << segments[0].squared_length() << " " << segments[1].squared_length();
    do
    {
        merge = false;
        sort(segments.begin(), segments.end(), [](K_epec::Segment_2& s1, K_epec::Segment_2& s2){return CGAL::to_double(s1.squared_length()) > CGAL::to_double(s2.squared_length());});
        for(int i = 0; i < segments.size(); i++){
            K_epec::Segment_2 s1 = segments[i];
            Point_2 p1 = to_geom(s1.min());
            Point_2 p2 = to_geom(s1.max());
            Point_2 mid1 = CGAL::midpoint(p1, p2);
            Vector_2 v1(p1, p2);
            double len1 = CGAL::sqrt(v1.squared_length());
            for(int j = i+1; j < segments.size(); j++){
                K_epec::Segment_2 s2 = segments[j];
                Point_2 p3 = to_geom(s2.min());
                Point_2 p4 = to_geom(s2.max());
                Point_2 mid2 = CGAL::midpoint(p3, p4);
                Vector_2 v2(p3, p4);
                double len2 = CGAL::sqrt(v2.squared_length());
                double cos = v1 * v2 / (len1 * len2);
                double d = robust_squared_distance(p1, p2, p3, p4);
                double d2 = ((p4.y() - p3.y()) * (p1.x() - p2.x())) - ((p1.y() - p2.y()) * (p4.x() - p3.x()));
                if(d <= EPSILON_DISTANCE && std::abs(d2) <= EPSILON_GENERAL){ //merge
                    merge = true;
                    Point_2 mid((len1*mid1.x() + len2*mid2.x())/(len1+len2), (len1*mid1.y() + len2*mid2.y())/(len1+len2));
                    Line_2 line(mid, v1 + v2);
                    Point_2 pp1 = line.projection(p1);
                    Point_2 pp2 = line.projection(p2);
                    Point_2 pp3 = line.projection(p3);
                    Point_2 pp4 = line.projection(p4);
                    Point_2 left = pp1 < pp3 ? pp1 : pp3;
                    Point_2 right = pp2 > pp4 ? pp2: pp4; 
                //test
                if(abs(p1.x() - 14.8663) < 1e-3 && abs(p1.y() + 66.7255) < 1e-3 || abs(p3.x() - 14.8663) < 1e-3 && abs(p3.y() + 66.7255) < 1e-3){
                    LOG(INFO) << "TEST: " << s2 << " " << s1 << " " << cos << " " << d;
                    Point_2 mid = CGAL::midpoint(p1, p2);
                    Line_2 l(mid, v2);
                    Segment_2 s(l.projection(p1), l.projection(p2));
                    p1 = s.min(); p2 = s.max();
                    d = ((p4.y() - p3.y()) * (p1.x() - p2.x())) - ((p1.y() - p2.y()) * (p4.x() - p3.x()));
                    if(std::abs(d) < 1.192092896e-07)
                    LOG(INFO) << "PARALLEL " << robust_squared_distance(p1, p2, p3, p4);
                    else
                    LOG(INFO) << d;
                }   
                    segments.erase(segments.begin() + j);
                    segments.erase(segments.begin() + i);              
                    segments.push_back(K_epec::Segment_2(to_geom2(left), to_geom2(right)));
                    break;
                }
            }
            if(merge)
                break;
        }
    }while(merge);
}*/

// fixed: generate floor
bool generate_floor(vector<Plane_3>& floors, cv::Mat& hfmap, cv::Mat& nfmap, vector<Line_2>& intersectf, const cm::Config& config,  Point_3 corner, double step, string wdir){
    
    const auto floorfile = config.get_optional<string>("file.floorname");
    Pwn_vector flpoints;
    if (!PLYPointLoad3(wdir + *floorfile, flpoints)){
        LOG(INFO)  << "Load failed.";
        return false;
    }
    LOG(INFO)  << "Load finished, "<< flpoints.size() << " floor points. ";
    
    // ransac detetcion, only accept horizantal floors
    // the first detection detect big shape, the second detection detect small facade floor
    if(*config.get_optional<bool>("floor.enabled") == true){ // detection
        // push black plane 
        // assume min height
        floors.push_back(Plane_3(Point_3(0,0,corner.z()-10),Vector_3(0,0,1)));
        ransac_detection_floor(flpoints, floors, intersectf, config, wdir);
        LOG(INFO)  << floors.size()-1 << " floors are detected. " << intersectf.size() << " intersects are generated.";
        
        // save floors and intersectf
        /* floor_num intersectf_num
        a b c d
        a b c*/ 
        ofstream file(wdir + "floor_intersectf", ios::out);
        file << floors.size() << " " << intersectf.size() << "\n"; 
        for(auto f:floors)
        file << setiosflags(ios::fixed) << setprecision(8) << f.a() << " " << f.b() << " " << f.c() << " " << f.d() << "\n";
        for(auto s:intersectf)
        file << setiosflags(ios::fixed) << setprecision(8) << s.a() << " " << s.b() << " " << s.c() << "\n";
        file.close();
    }
    else{ // load from file
        ifstream file(wdir + "floor_intersectf");
        int num_f,num_inf;
        double x, y, z, d;
        file >> num_f >> num_inf;
        while(num_f-->0){
          file >> x >> y >> z >> d;
          floors.push_back(Plane_3(x,y,z,d));
        }
        while(num_inf-->0){
          file >> x >> y >> z;
          intersectf.push_back(Line_2(x,y,z));
        }
     }
    
    // hfmap
    cv::Mat num_map = cv::Mat::zeros(hfmap.rows, hfmap.cols, CV_64FC2); // {num, z}
    for(int  i = 0; i < flpoints.size(); i++){
        int r = int((corner.y() - flpoints[i].first.y()) / step);
        int c = int((flpoints[i].first.x() - corner.x()) / step);
        if( r < 0 || r > hfmap.rows || c < 0 || c > hfmap.cols) continue;
        num_map.at<cv::Vec2d>(r, c)[0]++;
        num_map.at<cv::Vec2d>(r, c)[1] += flpoints[i].first.z();    
    }
    for(int r = 0; r < hfmap.rows; r++)
      for(int c = 0; c < hfmap.cols; c++){
        hfmap.at<cv::Vec3d>(r, c)[0] = corner.x() + (c+0.5)*step;
        hfmap.at<cv::Vec3d>(r, c)[1] = corner.y() - (r+0.5)*step;
        if(num_map.at<cv::Vec2d>(r, c)[0] == 0){
          // no points, assign black height: zmin - 10
          hfmap.at<cv::Vec3d>(r, c)[2] = corner.z() - 10;
          nfmap.at<int>(r, c) = 0;
        }
        else
          hfmap.at<cv::Vec3d>(r, c)[2] = num_map.at<cv::Vec2d>(r, c)[1]/num_map.at<cv::Vec2d>(r, c)[0];
          nfmap.at<int>(r, c) = num_map.at<cv::Vec2d>(r, c)[0];                     
      }

    // save hcmap to height map, scale to [0, 255]
    double zmin = 0.0, zmax = 0.0;
    vector<cv::Mat> channels;
    cv::split(hfmap,channels); // split z height value
    cv::minMaxLoc(channels[2], &zmin, &zmax);
    LOG(INFO) << "floor height map: " << zmin << " zmin, " << zmax << " zmax.";
    const double zrange = zmax - zmin;
    cv::Mat img;
    channels[2].convertTo(img, CV_8U, 255.0 / zrange, -zmin * 255.0 / zrange);
    cv::imwrite(wdir + "floor_height.png", img);

}

// fixed: generate ceiling
bool generate_ceiling(vector<Plane_3>& ceilings, cv::Mat& hcmap, cv::Mat& ncmap, vector<Line_2>& intersectc, const cm::Config& config, Point_3 corner, double step, string wdir){

    const auto ceilingfile = config.get_optional<string>("file.ceilingname");
    Pwn_vector cpoints;
    bool tag = config.get<bool>("facade.ceiling_detection");
    if (!PLYPointLoad3(wdir + *ceilingfile, cpoints)){
        LOG(INFO)  << "Load failed.";
        return false;
    }
    LOG(INFO)  << "Load finished, "<< cpoints.size() << " ceiling points. ";
    
    if(*config.get_optional<bool>("ceiling.enabled") == true){
        // push black plane 
         // assume min height
        ceilings.push_back(Plane_3(Point_3(0,0,corner.z()-10),Vector_3(0,0,1)));
        // accept horizantal and slap big plane
        ransac_detection_ceiling(cpoints, ceilings, intersectc, config, wdir);

        LOG(INFO)  << ceilings.size()-1 << " ceilings are detected firstly. ";

        if(tag){
            // add ceilings from facade pc
            // recorded in ceiling_facade
            vector<Plane_3> ceilings_;
            map<Point_3, set<double>> c_f; // {normal, d}
            ifstream file(wdir + "ceilings_facade", ios::in);
            double x, y, z, d;
            while(file >> x >> y >> z >> d){
                if(c_f.count(Point_3(x, y, z)))
                c_f[Point_3(x, y, z)].insert(d);
                else if(c_f.count(Point_3(-x, -y, -z)))
                c_f[Point_3(-x, -y, -z)].insert(-d);
                else
                c_f.insert(make_pair<Point_3, set<double>>(Point_3(x, y, z),{d}));
            }
            file.close();
            // merge close parallel plane
            // unsure: could be solved by ransac regularization
            double dis_thred = *(config.get_optional<double>("ceiling.regularize.dis_thred"));
            for(auto it = c_f.begin(); it != c_f.end(); it++){
                for(auto it2 = it->second.begin(); it2 != it->second.end();){
                    double dn = *it2;
                    for(auto it3 = it2; it3 !=it->second.end();){
                        if(*it3 <= dn + dis_thred)
                        it3++;
                        else{
                        dn = (*it2 + *(--it3))/2;
                        ceilings_.push_back(Plane_3(it->first.x(), it->first.y(), it->first.z(), dn));
                        it2 = ++it3;
                        break;
                        }
                    }
                    dn = (*it2 + *(--it->second.end()))/2;
                    ceilings_.push_back(Plane_3(it->first.x(), it->first.y(), it->first.z(), dn));
                    it2++;
                }
            }
            c_f.clear();
            for(auto p : ceilings)
            if(p == *ceilings.begin()) //black plane
              continue;
            else
              c_f.insert(make_pair<Point_3, set<double>>(Point_3(p.a(), p.b(), p.c()),{p.d()}));
            for(auto p: ceilings_){
                // close enough plane
                auto it = c_f.find(Point_3(p.a(), p.b(), p.c()));
                if(it!=c_f.end() && p.d() < *it->second.begin() + dis_thred && p.d() > *it->second.begin() - dis_thred) continue;
                it = c_f.find(Point_3(-p.a(), -p.b(), -p.c()));
                if(it!=c_f.end() && p.d() < -(*it->second.begin()) + dis_thred && p.d() > -(*it->second.begin()) - dis_thred) continue; 
                ceilings.push_back(p);
            }

        }
        
        // intersect to form intersectc
        for(int i = 1; i < ceilings.size();i++)
        for(int j = i+1; j < ceilings.size();j++){
          auto result = CGAL::intersection(ceilings[i], ceilings[j]);
          if(!result) continue;
          if(const Line_3* l = boost::get<Line_3>(&*result)){
            Line_2 line(Point_2(l->projection(Point_3(0,0,0)).x(),l->projection(Point_3(0,0,0)).y()),Vector_2(l->to_vector().x(),l->to_vector().y()));
            auto it = find(intersectc.begin(), intersectc.end(), line);
            if(it == intersectc.end())
              intersectc.push_back(line);  
        }
        }
        LOG(INFO)  << ceilings.size()-1 << " ceilings are generated finally. "<< intersectc.size() << " intersects are generated.";
        

        // save ceilings and intersectc
        /* ceiling_num intersectc_num
        a b c d
        a b c*/ 
        ofstream file2(wdir + "ceiling_intersectc", ios::out);
        file2 << ceilings.size() << " " << intersectc.size() << "\n"; 
        for(auto f:ceilings)
        file2 << f.a() << " " << f.b() << " " << f.c() << " " << f.d() << "\n";
        for(auto s:intersectc)
        file2 << s.a() << " " << s.b() << " " << s.c() << "\n";
        file2.close();
    }

    else{ // load from file
        ifstream file(wdir + "ceiling_intersectc");
        int num_f,num_inf;
        double x, y, z, d;
        file >> num_f >> num_inf;
        while(num_f-->0){
            file >> x >> y >> z >> d;
            ceilings.push_back(Plane_3(x,y,z,d));
        }
        while(num_inf-->0){
            file >> x >> y >> z;
            intersectc.push_back(Line_2(x,y,z));
        }
    }
    
    // hcmap
    cv::Mat num_map = cv::Mat::zeros(hcmap.rows, hcmap.cols, CV_64FC2); // {num, z}
    for(int  i = 0; i < cpoints.size(); i++){
        int r = int((corner.y() - cpoints[i].first.y()) / step);
        int c = int((cpoints[i].first.x() - corner.x()) / step);
        if( r < 0 || r > hcmap.rows || c < 0 || c > hcmap.cols) continue;
        num_map.at<cv::Vec2d>(r, c)[0]++;
        num_map.at<cv::Vec2d>(r, c)[1] += cpoints[i].first.z();    
    }
    if(tag){
        // load ceiling points from ceilings_facade.ply
        Pwn_vector cfpoints;
        PLYPointLoad3(wdir + "ceilings_facade.ply", cfpoints);
        for(int  i = 0; i < cfpoints.size(); i++){
        int r = int((corner.y() - cfpoints[i].first.y()) / step);
        int c = int((cfpoints[i].first.x() - corner.x()) / step);
        if( r < 0 || r > hcmap.rows || c < 0 || c > hcmap.cols) continue;
        num_map.at<cv::Vec2d>(r, c)[0]++;
        num_map.at<cv::Vec2d>(r, c)[1] += cfpoints[i].first.z();    
        }
    }

    // compute ave z height
    for(int r = 0; r < hcmap.rows; r++)
      for(int c = 0; c < hcmap.cols; c++){
        hcmap.at<cv::Vec3d>(r, c)[0] = corner.x() + (c+0.5)*step;
        hcmap.at<cv::Vec3d>(r, c)[1] = corner.y() - (r+0.5)*step;
        if(num_map.at<cv::Vec2d>(r, c)[0] == 0){
          // no points, assign black height: zmin - 10
          hcmap.at<cv::Vec3d>(r, c)[2] = corner.z() - 10;
          ncmap.at<int>(r, c) = 0;     
        }
        else
          hcmap.at<cv::Vec3d>(r, c)[2] = num_map.at<cv::Vec2d>(r, c)[1]/num_map.at<cv::Vec2d>(r, c)[0];
          ncmap.at<int>(r, c) = num_map.at<cv::Vec2d>(r, c)[0];                        
      }

    // save hcmap to height map, scale to [0, 255]
    double zmin = 0.0, zmax = 0.0;
    vector<cv::Mat> channels;
    cv::split(hcmap,channels); // split z height value
    cv::minMaxLoc(channels[2], &zmin, &zmax);
    const double zrange = zmax - zmin;
    LOG(INFO) << "ceiling height map: " << zmin << " zmin, " << zmax << " zmax.";
    cv::Mat img;
    channels[2].convertTo(img, CV_8U, 255.0 / zrange, -zmin * 255.0 / zrange);
    cv::imwrite(wdir + "ceilings_height.png", img);

}


// non-manifold and normal condition
bool generate_fc(vector<Segment_2>& segments, const Point_3 cornerz, const Point_2 length, const cm::Config& config, string wdir){
    I2E to_exact;

    // generate grid arr
    double step = *(config.get_optional<double>("arr.step"));
    int rows = std::ceil(length.y() / step);
    int cols = std::ceil(length.x() / step);
    const Point_2 corner(cornerz.x(),cornerz.y());
    LOG(INFO) << rows << " rows, " << cols << " cols.";

    // for building class
    cv::Mat hfmap(rows, cols, CV_64FC3);
    cv::Mat hcmap(rows, cols, CV_64FC3);
    cv::Mat nfmap(rows, cols, CV_32S);
    cv::Mat ncmap(rows, cols, CV_32S);
    vector<Plane_3> floors;
    vector<Line_2> intersectf;
    vector<Plane_3> ceilings;
    vector<Line_2> intersectc;

    generate_floor(floors, hfmap, nfmap, intersectf, config, cornerz, step, wdir);
    generate_ceiling(ceilings, hcmap, ncmap, intersectc, config, cornerz, step, wdir);
    
    if(floors.size() == 0 || ceilings.size() == 0){
        LOG(INFO) << "floor, ceiling ransac detection failed.";
        return false;
    }
    // bbox
    const Point_2 corner2(corner.x()+length.x(), corner.y()-length.y());
    CGAL::Iso_rectangle_2<K_epec> bbox(to_exact(corner), to_exact(corner2)); 
    // push all segment together
    vector<K_epec::Segment_2> total_segs;

    //1. push intersection segments
    for(auto &l:intersectf)
     if(auto result = CGAL::intersection(bbox, to_exact(l)))
        if(K_epec::Segment_2* rl = boost::get<K_epec::Segment_2>(&*result))
          total_segs.push_back(*rl);
        
    for(auto &l:intersectc)
     if(auto result = CGAL::intersection(bbox, to_exact(l)))
        if(K_epec::Segment_2* rl = boost::get<K_epec::Segment_2>(&*result))
          total_segs.push_back(*rl);
          
    //2. push extended facade segments
    vector<K_epec::Segment_2> selected;
    double extend = *(config.get_optional<double>("arr.extend"));
    double dialog = CGAL::sqrt(length.x()*length.x() + length.y()*length.y());
    for(auto s: segments){
        selected.push_back(to_exact(s));
        const auto mid = CGAL::midpoint(selected.back().source(), selected.back().target());
        const double len = CGAL::sqrt(s.squared_length());
        const auto dir = selected.back().to_vector()/len;
        total_segs.push_back({mid + (len/2 + extend*dialog)*dir, mid - (len/2 + extend*dialog)*dir});
    }

   
    // TODO:merge close segments in total_segs
    //merge(total_segs, corner);

    
    //3. add bounding segments to close the arrangement
    total_segs.push_back(K_epec::Segment_2(bbox.vertex(0),bbox.vertex(1))); 
    total_segs.push_back(K_epec::Segment_2(bbox.vertex(1),bbox.vertex(2)));
    total_segs.push_back(K_epec::Segment_2(bbox.vertex(2),bbox.vertex(3)));
    total_segs.push_back(K_epec::Segment_2(bbox.vertex(3),bbox.vertex(0)));

    //LOG(INFO) << total_segs.size() - segments.size() << " segments are added.";


    building fc(wdir, bbox, total_segs, selected, corner, step, hfmap, hcmap, nfmap, ncmap, floors, ceilings, cornerz.z() - 10); // min height = zmin - 10
    if(!fc.segment_arrangement_modeling())
        return false;
    return true;

}

// two-manifold and open boundary
bool generate_fc(vector<Segment_2> segments, double top, double bottom, string wdir){
   map<Point_2,int> points;
   vector<Point_2> single_points;
   //vector< vector<int> > polgons;
   vector<int> poly; // assume one object
   
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
    int num_points = points.size();
    map<Point_2,int>::iterator iter;
    for (iter=points.begin(); iter!=points.end(); iter++){
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

    for(auto s:segments){
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

bool generate_map(string wdir, double size){

    
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
       int s,t,ss;
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
    maps*=255;
    const string fname = wdir + "freespace.png";
    cv::imwrite(fname, maps);
    return true;

}

/*
bool generate_modeling(string wdir){
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
        if(!ifs) continue;
        ifs >> str;
        ifs >> np >> nf >> n;
        int *order = new int[np];
        for(int i = 0; i < np; i++){
            ifs >> x >> y >> z;
            it = find(points.begin(), points.end(), Point_3(x, y, z));
            if(it == points.end()){
            points.push_back(Point_3(x, y, z));
            order[i] = num_p;
            num_p++;
            }
            else
                order[i] = it - points.begin();      
        }

        for(int i = 0; i < nf; i++){
            ifs >> n;
            int index;
            vector<int> poly;
            for(int i = 0; i < n ; i++){
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
*/

bool generate_modeling(string wdir){
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
    /*
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
     */
    return true;
}