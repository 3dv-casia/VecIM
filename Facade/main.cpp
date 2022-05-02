#include "base.h"
#include "../basic/Config.h"
#include <glog/logging.h>

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        LOG(INFO) << "Incomplete input.";
        return -1;
    }

    // working directory
    const std::string wdir(argv[1]);   
    google::InitGoogleLogging("");
    google::SetLogDestination(google::GLOG_INFO, (wdir + "logtestInfo").c_str()); 
    google::SetStderrLogging(google::GLOG_INFO);

    // load configuration
    cm::read_config(wdir + "/config_modeling.xml");
    LOG(INFO) << "Config file loaded.";
    const cm::Config &config = cm::get_config();
  
    // timer
    clock_t startp, endp;
    startp = clock();

    // read ply file, load 3D pointcloud with normal 
    const auto facadefile = config.get_optional<string>("file.facadename");
    Pwn_vector fapoints;
    if (!PLYPointLoad3(wdir + *facadefile, fapoints))
    {
        LOG(INFO)  << "Failed to load facade points file.";
        return -1;
    }
    LOG(INFO)  << "Load "<< fapoints.size() << " facade points. ";

    // record results
    vector<Line> lines;     // facade projection lines
    vector<IntersectPoint> interpoints;     // segment intersection
    vector<Segment> segments;     // facade projection segments
    vector<Segment_2> selected;     // floorplan segments
    vector<Point_2> heights;    // segment z-height
    vector<Vector_2> oritations;    // segment normals
    vector<Segment> bsegments;    // boundary segments 
    unordered_map<Point_2, int> points_idx; // point index in pointcloud

    // obtain x and y boundaries
    float xmin = 1e9; float xmax = -1e9; float ymin = 1e9; float ymax = -1e9; float zmin = 1e9;
    for(int i = 0; i < fapoints.size(); i++)
    {
        if (fapoints[i].first.x() < xmin) xmin = fapoints[i].first.x();
        if (fapoints[i].first.x() > xmax) xmax = fapoints[i].first.x();
        if (fapoints[i].first.y() < ymin) ymin = fapoints[i].first.y();
        if (fapoints[i].first.y() > ymax) ymax = fapoints[i].first.y();
        if (fapoints[i].first.z() < zmin) zmin = fapoints[i].first.z();
  }

    // bbox egdes: a little bigger than boundaries
    xmin -= 2.0; xmax += 2.0; ymin -= 2.0; ymax += 2.0;
    bsegments.push_back(Segment(Segment_2(Point_2(xmin,ymin),Point_2(xmin,ymax))));
    bsegments.push_back(Segment(Segment_2(Point_2(xmin,ymin),Point_2(xmax,ymin))));
    bsegments.push_back(Segment(Segment_2(Point_2(xmax,ymin),Point_2(xmax,ymax))));
    bsegments.push_back(Segment(Segment_2(Point_2(xmax,ymax),Point_2(xmin,ymax)))); 
    
    // Part 1. 2D facade detection
    if(*config.get_optional<bool>("facade.enabled") == true)
    {
        // detect 3D facade planes and then project them to generate 2D segments
        LOG(INFO)  << " /***3D facade plane detection and projection generation.***/ ";
        // load parameters
        const auto probability = config.get_optional<double>("facade.ransac.probability");
        const auto min_points = config.get_optional<int>("facade.ransac.min_points");
        const auto epsilon = config.get_optional<double>("facade.ransac.epsilon");
        const auto cluster_epsilon = config.get_optional<double>("facade.ransac.cluster_epsilon");
        const auto normal_threshold = config.get_optional<double>("facade.ransac.normal_threshold");
        const auto cos_angle = config.get_optional<double>("facade.ransac.cos_angle");

        if (!ransac_detection_p(fapoints, *probability, *min_points, *epsilon, *cluster_epsilon, *normal_threshold, *cos_angle, lines, wdir)) 
        {
            LOG(INFO) << "No facade plane is detected.";
            return -1;
        }
        Pwn_vector().swap(fapoints); // clear facade points vector

        // construct segment hypothesis
        LOG(INFO)  << " /***Hypothesis generation.***/ ";
        // load parameters
        const auto angle_thred = config.get_optional<double>("facade.regularize.angle_thred");
        const auto num_thred = config.get_optional<double>("facade.regularize.num_thred");
        const auto radio = config.get_optional<double>("facade.regularize.extend_ratio");
        
        if(!generation(lines, bsegments, segments, points_idx, interpoints, *angle_thred, *num_thred, *radio, wdir))
        {
            LOG(INFO)  << "Hypothesis generation is failed.";
            return -1;
        }

        // select optimal segments to generate floorplan
        LOG(INFO)  << " /***Segments Selection.***/ " ;
        // load parameters
        const auto lambda_data_fitting = config.get_optional<double>("facade.segment.lamda_data_fitting");
        const auto lambda_model_coverage = config.get_optional<double>("facade.segment.lamda_model_coverage");
        const auto lambda_model_complexity = config.get_optional<double>("facade.segment.lamda_model_complxity");
        const auto alpha = config.get_optional<double>("facade.segment.alpha");
        
        optimize(selected, points_idx, heights, oritations, lines, segments,  interpoints, bsegments, *lambda_data_fitting, *lambda_model_coverage, *lambda_model_complexity, *alpha, wdir); 
        if(selected.size() == 0)
        {
            LOG(INFO) << "No segment is selected, segment selection is failed.";
            return -1;
        }
        else
        { 
            // save selected segments
            ofstream file(wdir + "floorplan", ios::out);
            for(int i = 0; i < selected.size(); i++){
              file << setiosflags(ios::fixed) << setprecision(8) << selected[i].source().x() <<" "<< selected[i].target().x() <<" "<< selected[i].source().y() <<" "<< selected[i].target().y() << endl;
            }
            file.close();
        }
      
    }
    else
    {
        ifstream f(wdir + "floorplan");
        if(!f.is_open())
        {
            LOG(INFO) << "Failed to load floorplan file.";
            return -1;     
        }
        double x1, x2, y1, y2;
        while(f >> x1 >> x2 >> y1 >> y2){
            selected.push_back({{x1, y1}, {x2, y2}});
        }
    }

    // Part 2-1. floor, ceiling, and facade reconstruction with pointcloud
    if(*config.get_optional<bool>("ceiling_floor") == true)
    {
        const Point_3 corner(xmin, ymax, zmin);
        const Point_2 length(xmax - xmin, ymax - ymin);
        if (!generate_fcf(selected, corner, length, config, wdir))
        {
            LOG(INFO) << "Floor, ceiling and facade reconstruction is failed.";
            return -1;
        } 
        else
        {
            LOG(INFO) << "Floor, ceiling and facade reconstruction is done.";
        }
    }
    
    // Part 2-2. floor, ceiling, and facade reconstruction without pointcloud
    if(*config.get_optional<bool>("fheight.enabled") == true)
    {
        Mesh<K_epec::Point_3> facades;
        const auto top = config.get_optional<double>("fheight.top"); // horizontal ceiling height
        const auto bottom = config.get_optional<double>("fheight.bottom"); // horizontal floor height
        producemesh(selected, facades, heights, oritations, *top, *bottom);
        string filename = wdir + "facade.off";
        OFFMeshSave(filename, facades);
        if(!generate_fc(selected, *top, *bottom, wdir)){
            LOG(INFO) << "Floor and ceiling generation is failed.";
            return -1;
          }
    }

    // Part 3. cylinder reconstruction
    const auto cylinderfile = config.get_optional<string>("file.cylindername");
    Pwn_vector cypoints;
    vector<Circle> circles;
    if(*config.get_optional<bool>("cylinder.enabled") == true && PLYPointLoad3(wdir + *cylinderfile, cypoints))
    {
        // detect cylinder
        LOG(INFO)  << " /***Cylinder detection and circle generation.***/ ";
        // load parameters
        const auto probability = config.get_optional<double>("cylinder.ransac.probability");
        const auto min_points = config.get_optional<int>("cylinder.ransac.min_points");
        const auto epsilon = config.get_optional<double>("cylinder.ransac.epsilon");
        const auto cluster_epsilon = config.get_optional<double>("cylinder.ransac.cluster_epsilon");
        const auto normal_threshold = config.get_optional<double>("cylinder.ransac.normal_threshold");
        const auto cos_angle = config.get_optional<double>("cylinder.ransac.cos_angle");
        const auto radius = config.get_optional<double>("cylinder.regularize.cylinder_r");

        if (!ransac_detection_c(cypoints, *radius, *probability, *min_points, *epsilon, *cluster_epsilon, *normal_threshold, *cos_angle, circles, wdir)) 
            return -1;
        else
            LOG(INFO) << "Cylinder generation is done.";
    }

    // Part 4. merge different modeling results
    if (config.get_optional<bool>("merge"))
    {
        clock_t m1 = clock();
        if(!generate_off_modeling(wdir))
        {
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
