#include "base.h"
#include "../basic/Config.h"
#include <glog/logging.h>

int main(int argc, char* argv[])
{ 
/*
  //test
  Pwn_vector points;
  string name = "/home/hjl/modeling/Algo3.0/Algo3.0/build_release/Facade/s3dis_area5/area5-no-facade-points.ply";
  PLYPointLoad3(name, points);
  int len = points.size();
  Pwn_vector points2;
  for(std::size_t i = 0; i < 300000; i++){
    int idx = rand()%len;
    points2.push_back(points[i]);
  }
  PLYPointSave("/home/hjl/modeling/Algo3.0/Algo3.0/build_release/Facade/s3dis_area5/area5-no-facade-points-sample.ply", points2);
  return 0;*/

  // working directory
  const std::string wdir(argv[1]);
  
  if (argc < 2)
    return 1;
    
  google::InitGoogleLogging("");
  google::SetLogDestination(google::GLOG_INFO, (wdir + "logtestInfo").c_str()); 
	google::SetStderrLogging(google::GLOG_INFO);


  // load configuration
  cm::read_config(wdir + "/config_modeling.xml");
  LOG(INFO) << "Config file loaded.";

  const cm::Config &config = cm::get_config();
 
   // timer
   clock_t startp,endp;
   startp = clock();


   // read ply file, load 3D pointcloud with normal 
  const auto facadefile = config.get_optional<string>("file.facadename");
  const auto cylinderfile = config.get_optional<string>("file.cylindername");

  Pwn_vector fapoints, cypoints;
  if (!PLYPointLoad3(wdir + *facadefile, fapoints)){
    LOG(INFO)  << "Load failed.";
    return -1;
  }
  LOG(INFO)  << "Load finished, "<< fapoints.size() << "  facade points. ";
/*
  // lidar pc segmantation
  Pwn_vector fpoints, cfpoints, cpoints, flpoints;
  for(auto p:fapoints){
    if (abs(p.second.z()/sqrt(p.second.x()*p.second.x()+p.second.y()*p.second.y()+p.second.z()*p.second.z())) < 0.5)
        fpoints.push_back(p);
    if (abs(p.second.z()/sqrt(p.second.x()*p.second.x()+p.second.y()*p.second.y()+p.second.z()*p.second.z())) > 0.5)
        cfpoints.push_back(p);
    //else
       //fpoints.push_back(p);
  }
  double h = 0;
  for(auto p:fpoints)
      h += p.first.z();
  h = h/fpoints.size();
  for(auto p:cfpoints)
     if(p.first.z() > h) cpoints.push_back(p);
     else flpoints.push_back(p);
  PLYPointSave(wdir + "wall-test.ply", fpoints);
  PLYPointSave(wdir + "ceiling-test.ply", cpoints);
  PLYPointSave(wdir + "floor-test.ply", flpoints);
  endp = clock();
  LOG(INFO) << "Segmentation time: " << (float)(endp-startp)/CLOCKS_PER_SEC << "s" ;
  return 0;
*/

   // result
   vector<Line> lines; 
   vector<IntersectPoint> interpoints;
   vector<Segment> segments;
   vector<Segment_2> selected; 
   vector<Point_2> heights;
   vector<Vector_2> oritations;
   unordered_map<Point_2, int> points_idx;


   // obtain xy bbox
   float xmin = 1e9; float xmax = -1e9; float ymin = 1e9; float ymax = -1e9; float zmin = 1e9;
   for(int i = 0; i < fapoints.size(); i++)
   {
	if (fapoints[i].first.x() < xmin) xmin = fapoints[i].first.x();
	if (fapoints[i].first.x() > xmax) xmax = fapoints[i].first.x();
	if (fapoints[i].first.y() < ymin) ymin = fapoints[i].first.y();
	if (fapoints[i].first.y() > ymax) ymax = fapoints[i].first.y();
  if (fapoints[i].first.z() < zmin) zmin = fapoints[i].first.z();
  }

   //Bbox_2 bbox = Bbox_2(xmin, ymin, xmax, ymax);
   // box egdes: a little bigger than bbox
   xmin -= 2.0; xmax += 2.0; ymin -= 2.0; ymax += 2.0;
   vector<Segment> bsegments;  
   bsegments.push_back(Segment(Segment_2(Point_2(xmin,ymin),Point_2(xmin,ymax))));
   bsegments.push_back(Segment(Segment_2(Point_2(xmin,ymin),Point_2(xmax,ymin))));
   bsegments.push_back(Segment(Segment_2(Point_2(xmax,ymin),Point_2(xmax,ymax))));
   bsegments.push_back(Segment(Segment_2(Point_2(xmax,ymax),Point_2(xmin,ymax)))); 
  
  // 1. 2D facade detection
  if(*config.get_optional<bool>("facade.enabled") == true)
  {

   // detect plane and convert to 2d line2
   LOG(INFO)  << " /***Plane detection and lines generation.***/ ";
   // load parameters
   const auto probability = config.get_optional<double>("facade.ransac.probability");
   const auto min_points = config.get_optional<int>("facade.ransac.min_points");
   const auto epsilon = config.get_optional<double>("facade.ransac.epsilon");
   const auto cluster_epsilon = config.get_optional<double>("facade.ransac.cluster_epsilon");
   const auto normal_threshold = config.get_optional<double>("facade.ransac.normal_threshold");
   const auto cos_angle = config.get_optional<double>("facade.ransac.cos_angle");

  if (!ransac_detection_p(fapoints, *probability, *min_points, *epsilon, *cluster_epsilon, *normal_threshold, *cos_angle, lines, wdir)) 
    {
	     LOG(INFO)  <<"No plane is detected.";
	     return -1;
    }
    Pwn_vector().swap(fapoints); // clear vector

   // construct hypothesis
   LOG(INFO)  << " /***Hypothesis generation.***/ ";
   // load parameters
   const auto angle_thred = config.get_optional<double>("facade.regularize.angle_thred");
   const auto num_thred = config.get_optional<double>("facade.regularize.num_thred");
   const auto radio = config.get_optional<double>("facade.regularize.extend_ratio");
   LOG(INFO) << "angle_thred: " << *angle_thred << " num_thred: " << *num_thred << " radio: " << *radio;

   if(!generation(lines, bsegments, segments, points_idx, interpoints, *angle_thred, *num_thred, *radio, wdir)){
	LOG(INFO)  << "Hypothesis generation failed!";
	return -1;
   }

   // segment selection
   LOG(INFO)  << " /***Segments Selection.***/ " ;
   // load parameters
   const auto lambda_data_fitting = config.get_optional<double>("facade.segment.lamda_data_fitting");
   const auto lambda_model_coverage = config.get_optional<double>("facade.segment.lamda_model_coverage");
   const auto lambda_model_complexity = config.get_optional<double>("facade.segment.lamda_model_complxity");
   const auto alpha = config.get_optional<double>("facade.segment.alpha");

   optimize(selected, points_idx, heights, oritations, lines, segments, interpoints, bsegments, *lambda_data_fitting, *lambda_model_coverage, *lambda_model_complexity,  (*epsilon)*30, *alpha, wdir);  
   if(selected.size()==0){
	    LOG(INFO) << "No segment is selected, segment selection failed!";
    	return -1;
   }
   else{ // save selected segments
    ofstream file(wdir + "segments_final", ios::out);
    for(int i = 0; i<selected.size(); i++){
      file << double(selected[i].source().x())<<" "<< double(selected[i].target().x())<<" "<< double(selected[i].source().y())<<" "<< double(selected[i].target().y())<< endl;
    }
	  file.close();
  }
  
  }
  else{
    ifstream f(wdir + "segments_final");
    double x1,x2,y1,y2;
    while(f >> x1 >> x2 >> y1 >> y2){
        selected.push_back({{x1, y1},{x2, y2}});
    }
  }

  // 2. floor ceiling facade reconstruction
  if(*config.get_optional<bool>("ceiling_floor") == true)
  {
    const Point_3 corner(xmin, ymax, zmin);
    const Point_2 length(xmax - xmin, ymax - ymin);
    if (!generate_fc(selected, corner, length, config, wdir)){
      LOG(INFO) << "Floor, ceiling and facade failed.";
      return -1;
    } 
    else
    {
      LOG(INFO) << "Floor, ceiling and facade done.";
    }
    

  }
  
  // 2. no ceiling/floor point clouds recontruction
  if(*config.get_optional<bool>("fheight.enabled") == true){
     Mesh planes;
     const auto top = config.get_optional<double>("fheight.top");
     const auto bottom = config.get_optional<double>("fheight.bottom");
     LOG(INFO) << "top: " << top << ", bottom: " << bottom;
 	   producemesh(selected, planes, heights, oritations, *top, *bottom);
     string filename = wdir + "facade.ply";
     PLYTriMeshSave(planes, config, filename);
     if(!generate_fc(selected, *top, *bottom, wdir)){
         LOG(INFO) << "Floor and ceiling generation failed.";
         return -1;
      }
  }

  // fixed: cylinder
  // 3. cylinder reconstruction
  if (*config.get_optional<bool>("cylinder.enabled") == true && PLYPointLoad3(wdir + *cylinderfile, cypoints)){

    vector<Circle> circles;

    // detect cylinder
   LOG(INFO)  << " /***Cylinder detection and circle generation.***/ ";
   // load parameters
   const auto probability = config.get_optional<double>("cylinder.ransac.probability");
   const auto min_points = config.get_optional<int>("cylinder.ransac.min_points");
   const auto epsilon = config.get_optional<double>("cylinder.ransac.epsilon");
   const auto cluster_epsilon = config.get_optional<double>("cylinder.ransac.cluster_epsilon");
   const auto normal_threshold = config.get_optional<double>("cylinder.ransac.normal_threshold");
   const auto cos_angle = config.get_optional<double>("cylinder.ransac.cos_angle");
   const auto radious = config.get_optional<double>("cylinder.regularize.cylinder_r");
   

   if (!ransac_detection_c(cypoints, *radious, *probability, *min_points, *epsilon, *cluster_epsilon, *normal_threshold, *cos_angle, circles, wdir)) 
	     return -1;
    else
      LOG(INFO) << "Cylinder generation down.";

 }


  // 4. merge
  if (config.get_optional<bool>("merge"))
  {
     clock_t m1 = clock();

      if(!generate_modeling(wdir)){
      LOG(INFO) << "Merge failed. ";
      return -1;
    }
    LOG(INFO) << "Merge done. ";
    clock_t m2 = clock();

    LOG(INFO) << "Merge Time: " << (float)(m2-m1)/CLOCKS_PER_SEC << "s";
  }

  endp = clock();
  LOG(INFO) << "Total time: " << (float)(endp-startp)/CLOCKS_PER_SEC << "s" ;
  
  return 0;
}
