#include <fstream>
#include <limits.h>
#include <cmath>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/Regularization/regularize_planes.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/enum.h>
#include <list>


#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/ml.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "gco-v3.0/GCoptimization.h"

#include "arr.h"
#include "../basic/Config.h"
#include <glog/logging.h>

typedef CGAL::Triangulation_2<K>    Triangulation;
typedef Triangulation::Face_handle    Face_handle;
typedef Triangulation::Point      Point;

typedef CGAL::Partition_traits_2<K>     Traits;
typedef Traits::Polygon_2     Polygon2;
typedef std::list<Polygon2>     Polygon_list;


building::building(
    const std::string &wdir,
    const CGAL::Iso_rectangle_2<K_epec> &bbox,
    const std::vector<K_epec::Segment_2> &total_segs,
    const std::vector<K_epec::Segment_2> &facade_segs,
    const Point_2 corner,
    double step,
    const cv::Mat &hfmap,
    const cv::Mat &hcmap,
    const cv::Mat &nfmap,
    const cv::Mat &ncmap,
    const vector<Plane_3> floors,
    const vector<Plane_3> ceilings,
    const double &min_height):
    m_dir(wdir),
    m_bbox(bbox),
    m_segments(total_segs),
    m_facade_segs(facade_segs),
    m_corner(corner),
    m_step(step),
    m_hfmap(hfmap),
    m_hcmap(hcmap),
    m_nfmap(nfmap),
    m_ncmap(ncmap),
    m_fplanes(floors),
    m_cplanes(ceilings),
    m_min_height(min_height),
    m_rows(m_hfmap.rows),
    m_cols(m_hcmap.cols){}

bool building::segment_arrangement_modeling() 
{
    // timer
    clock_t start, mid, end;
    start = clock(); 
    
    const auto &config = cm::get_config();
    bool r1, r2, r3;
    if(m_segments.size() == 0 || m_facade_segs.size() == 0){
        LOG(INFO) << "Segments error, #m_segments: " << m_segments.size() << " , #m_facade_segs: " << m_facade_segs.size();
        return false;
    }

    LOG(INFO) << m_fplanes.size() << " floor planes, " << m_cplanes.size() << " ceiling planes.";

    write_all_segments();
    LOG(INFO) << "#m_segments: " << m_segments.size();

    // detected line segments to arrangement
    // obtain grid_arr m_arr overlay_arr m_arr_fdata. m_arr_cdata
    construct_segment_arrangement();
    write_arrangement();
    LOG(INFO) << "#floor data " << m_arr_fdata.size();
    LOG(INFO) << "#ceiling data " << m_arr_cdata.size();

    // recognize facade segments in m_arr edges
    // obtain m_facade facades_tag
    facade_segment_recognization();
    write_facade_segments();
    LOG(INFO) << "#facade data " << m_facade.size();
    
    // the whole MRF
    segment_arrangement_MRF_labeling_cf(m_arr_cdata, m_arr_fdata, m_fplanes, m_cplanes, m_carr, m_farr);

    // extrude and write arrangement
    {
        // extrude arrangement contour
        r1 = segment_arrangement_extrusion(m_fplanes, 0, m_farr);
        r2 = segment_arrangement_extrusion(m_cplanes, 1, m_carr);

        write_arrangement_labeling(m_fplanes, 0, m_farr);
        LOG(INFO) << "floor model extrusion done";

        write_arrangement_labeling(m_cplanes, 1, m_carr);
        LOG(INFO) << "ceiling model extrusion done";
    }

    mid = clock();
    LOG(INFO) << "Ceiling/Floor reconstruction Time: " << (float)(mid-start)/CLOCKS_PER_SEC << "s";

    // extrude facade
    write_all_facades();
    r3 = facade_extrusion();
    
    end = clock();
    LOG(INFO) << "Wall reconstruction Time: " << (float)(end-mid)/CLOCKS_PER_SEC << "s";

    return r1 && r2 && r3;
}

void building::construct_segment_arrangement() 
{
    // grid segments, in world coordinate
    vector<K_epec::Segment_2> grid_segs;
    for (int i = 0; i < m_rows + 1; ++i)
        grid_segs.push_back({
        K_epec::Point_2(m_corner.x(), m_corner.y() - i*m_step),
        K_epec::Point_2(m_corner.x() + m_cols*m_step, m_corner.y() - i*m_step)});
    for (int j = 0; j < m_cols + 1; ++j)
        grid_segs.push_back({
        K_epec::Point_2(m_corner.x() + j*m_step, m_corner.y()),
        K_epec::Point_2(m_corner.x() + j*m_step, m_corner.y() - m_rows*m_step)});

    // 1.construct sampling grid arrangement
    LOG(INFO) << "Construct sampling grid arrangement";
    LOG(INFO) << "#grid_segs: " << grid_segs.size();
    Arr_grid grid_arr;
    CGAL::insert(grid_arr, grid_segs.begin(), grid_segs.end());
    LOG(INFO) << "#f: " << grid_arr.number_of_faces() << ", #v: " << grid_arr.number_of_vertices();
    for (auto fitr = grid_arr.faces_begin(); fitr != grid_arr.faces_end(); ++fitr) 
    {
        if (fitr->is_unbounded()) {
            fitr->set_data(cv::Vec2i(-1, -1));
            continue;
        }
        auto curr = fitr->outer_ccb();
        K_epec::Point_2 ll_pt = curr->target()->point();
        do {
            if (curr->target()->point() < ll_pt)
              ll_pt = curr->target()->point();
            ++curr;
        } while (curr != fitr->outer_ccb());
        // world x, y -> image i, j
        const int ridx = int(std::round((m_corner.y() - CGAL::to_double(ll_pt.y()))/m_step)) - 1;
        const int cidx = int(std::round(CGAL::to_double((ll_pt.x()) - m_corner.x())/m_step));
        fitr->set_data(cv::Vec2i(ridx, cidx));
    }

    // 2.line segment arrangement
    LOG(INFO) << "Construct line segment arrangement";
    LOG(INFO) << "#segments: " << m_segments.size();
    m_arr.clear();
    CGAL::insert(m_arr, m_segments.begin(), m_segments.end()); 
    LOG(INFO) << "#f: " << m_arr.number_of_faces() << ", #v: " << m_arr.number_of_vertices();
    int idx = 0;
    for (auto fitr = m_arr.faces_begin(); fitr != m_arr.faces_end(); ++fitr) 
    {
        if (fitr->is_unbounded())
          fitr->set_data(-1);
        else
          fitr->set_data(idx++);
    }

    // trim arrangement
    LOG(INFO) << "Trim arrangement";
    for (auto eitr = m_arr.edges_begin(); eitr != m_arr.edges_end(); ++eitr) 
    {
        if (eitr->face()->data() == eitr->twin()->face()->data()) 
        {
            m_arr.remove_edge(eitr);
            eitr = m_arr.edges_begin();
        }
    }
    LOG(INFO) << "#f: " << m_arr.number_of_faces() << ", #v: " << m_arr.number_of_vertices();
    // merge collinear edges
    LOG(INFO) << "Merge collinear edges";
    for (auto fitr = m_arr.faces_begin(); fitr != m_arr.faces_end(); ++fitr) 
    {
        if (fitr->is_unbounded())
          continue;
        auto curr = fitr->outer_ccb(); 
        bool merged = false;
        do {
            merged = false;
            const Arr_with_index::Halfedge_handle h0 = curr;
            const Arr_with_index::Halfedge_handle h1 = curr->next();
            if (h0->target()->degree() == 2 &&
              h0->curve().line().has_on(h1->target()->point())) {
                curr = m_arr.merge_edge(h0, h1,
                K_epec::Segment_2(h0->source()->point(), h1->target()->point()));
                merged = true;
            }
            else
              ++curr;
        } while (merged || curr != fitr->outer_ccb());
    }
    LOG(INFO) << "#f: " << m_arr.number_of_faces() << ", #v: " << m_arr.number_of_vertices();

    // 3.overlay
    LOG(INFO) << "Overlaying arrangement";
    Arr_overlay overlay_arr;
    Overlay_traits overlay_traits;
    CGAL::overlay(m_arr, grid_arr, overlay_arr, overlay_traits);
    LOG(INFO) << "#f: " << overlay_arr.number_of_faces() << ", #v: " << overlay_arr.number_of_vertices();

    // config data term 
    compute_distance(m_farr, m_fplanes, m_hfmap, m_nfmap, m_arr_fdata, overlay_arr, 0); // floor
    compute_distance(m_carr, m_cplanes, m_hcmap, m_ncmap, m_arr_cdata, overlay_arr, 1); // ceiling
}

void building::compute_distance(Arr_with_index& m_fcarr, vector<Plane_3>& m_fcplanes, const cv::Mat& m_hfcmap, const cv::Mat& m_nfcmap, vector<Arr_face_data>& m_arr_fcdata, Arr_overlay& overlay_arr, int tag)
{
    m_fcarr = m_arr;
    // compute pixel distances to planes
    vector<cv::Mat> pixel_distance_to_planes;
    for (const auto &s : m_fcplanes) 
    {
        cv::Mat dis_map(m_rows, m_cols, CV_64F, cv::Scalar(0.0));
        for (int i = 0; i < m_rows; ++i) 
        {
            for (int j = 0; j < m_cols; ++j) 
            {
                const cv::Vec3d pos = m_hfcmap.at<cv::Vec3d>(i, j);
                auto result = CGAL::intersection(s,
                  Line_3(Point_3(pos[0], pos[1], pos[2]), Vector_3(0.0, 0.0, 1.0)));
                const Point_3 *intersect_pt = nullptr;
                if (result && (intersect_pt = boost::get<Point_3>(&(*result)))) 
                {
                    Point_3 p = *intersect_pt;
                    dis_map.at<double>(i, j) = std::abs(p.z() - pos[2]);             
                }
                else
                    LOG(INFO) << "No intersection.";
            }
        }
        double zmin = 0.0, zmax = 0.0;
        cv::minMaxLoc(dis_map, &zmin, &zmax);
        pixel_distance_to_planes.push_back(dis_map);
    }

    // set segment arrangement face data
    m_arr_fcdata = std::vector<Arr_face_data>(
      m_arr.number_of_faces() - m_arr.number_of_unbounded_faces(),
      Arr_face_data(m_fcplanes.size(), 0));
    // record each face_area blank pixel number
    map<int, pair<int, int>> record_pixel; // <blank, non-blank>
    for (auto fitr = overlay_arr.faces_begin(); fitr != overlay_arr.faces_end(); ++fitr) 
    {
        const int fidx = fitr->data().fidx;
        const cv::Vec2i pixel = fitr->data().pixel;
        if (fitr->is_unbounded() || fidx < 0)
          continue;
        if(record_pixel.count(fidx) == 0)
          record_pixel[fidx] = {0, 0};
        if (pixel[0] < 0 || pixel[1] < 0)
          record_pixel[fidx].first++;
        else
          m_nfcmap.at<int>(pixel) == 0 ? record_pixel[fidx].first++ : record_pixel[fidx].second++;     
    }

    for (auto fitr = overlay_arr.faces_begin(); fitr != overlay_arr.faces_end(); ++fitr) 
    {
        const int fidx = fitr->data().fidx;
        const cv::Vec2i pixel = fitr->data().pixel;
        if (fitr->is_unbounded() || fidx < 0)
          continue;
        Arr_face_data &arrd = m_arr_fcdata[fidx];
        CGAL::Polygon_2<K_epec> plg;
        auto curr = fitr->outer_ccb();
        do {
          plg.push_back(curr->target()->point());
          ++curr;
        } while (curr != fitr->outer_ccb());
        const double area = CGAL::to_double(plg.area());
        assert(area > 0.0);
        arrd.area += area;
        if (pixel[0] < 0 || pixel[1] < 0) 
        {
            // faces out side sampling grid is assigned with ground-like data
            for (auto &d : arrd.distances)
              d = 1e20;
            arrd.distances.front() = 0.0;
            arrd.pnumber = 0;           
        }
        else 
        {
            double blank_ratio = 1.0*record_pixel[fidx].first/(record_pixel[fidx].first+record_pixel[fidx].second);
            double nblank_ratio = 1.0*record_pixel[fidx].second/(record_pixel[fidx].first+record_pixel[fidx].second);
            for (std::size_t i = 0; i < pixel_distance_to_planes.size(); ++i){
                if(tag == 0) // floor
                    arrd.distances[i] += area * (m_nfcmap.at<int>(pixel) == 0 ? blank_ratio : 
                      nblank_ratio) * pixel_distance_to_planes[i].at<double>(pixel);
                else // ceiling
                    arrd.distances[i] += area * (m_nfcmap.at<int>(pixel) == 0 ? blank_ratio *
                      (i == 0 ? pixel_distance_to_planes[i].at<double>(pixel) : 
                      pixel_distance_to_planes.back().at<double>(pixel)) : 
                      nblank_ratio * pixel_distance_to_planes[i].at<double>(pixel));
            }
            arrd.pnumber += m_nfcmap.at<int>(pixel);
            arrd.bnum = record_pixel[fidx].first;
            arrd.nbnum = record_pixel[fidx].second;
        }
    }
    for (auto &arrd : m_arr_fcdata) {
        for (auto &d : arrd.distances)
          d /= arrd.area;
    }
}

void building::facade_segment_recognization()
{
    for (auto eitr = m_arr.edges_begin(); eitr != m_arr.edges_end(); ++eitr) 
    {
        Arr_with_index::Vertex_handle vs = eitr->source();
        Arr_with_index::Vertex_handle vt = eitr->target();
        int flag = 0;
        for(int i = 0; i < m_facade_segs.size(); i++)
        {
            if(m_facade_segs[i].has_on(vs->point()) && m_facade_segs[i].has_on(vt->point()))
            {
                // the edge is on the facade, exact computation
                //if(vs->point() < vt->point())
                m_facade.push_back(Facade_data(i, K_epec::Segment_2(vs->point(), vt->point())));
                //else
                //  m_facade.push_back(Facade_data(i, K_epec::Segment_2(vt->point(), vs->point())));
                flag = 1;
                facade_tags.push_back(1);
                break;
            }
        }
        if(flag == 0){
          facade_tags.push_back(0);
        }
    }
    if(facade_tags.size() != m_arr.number_of_edges())
      LOG(INFO) << "Arrangement error.";
}

// the whole MRF
void building::segment_arrangement_MRF_labeling_cf(std::vector<Arr_face_data>& m_arr_cdata, std::vector<Arr_face_data>& m_arr_fdata, std::vector<Plane_3>& m_fplanes, std::vector<Plane_3>& m_cplanes, Arr_with_index& m_carr, Arr_with_index& m_farr) 
{
  if(m_arr_cdata.size() != m_arr_fdata.size()){
    LOG(INFO) << "Arrangment error.";
    return;
  }
  
  typedef pair<int, int> Adj_facets;
  typedef map<Adj_facets, double> Adj_map;

  const auto &config = cm::get_config();
  LOG(INFO) << "MRF params: " <<
    config.get<bool>("arr.mrf.use_swap") << ' ' <<
    config.get<double>("arr.mrf.balance") << ' ' <<
    config.get<int>("arr.mrf.iterations");
  if (m_fplanes.size() < 2 || m_cplanes.size() < 2) {
    LOG(INFO) << "Insufficient shape.";
    return;
  }
  // set up arrangement graph
  // ceil: 0,..., m_arr_data.size()-1
  // label: blank, cifj
  const int nb_sites = int(m_arr_cdata.size());
  const int nb_labels = int((m_fplanes.size()-1) * (m_cplanes.size()-1) + 1);
  LOG(INFO) << nb_sites << " graph vertices. " << nb_labels << " labels.";

  // set up neighboring
  Adj_map adj_map;
  for (auto eitr = m_carr.edges_begin(); eitr != m_carr.edges_end(); ++eitr) 
  {
      if (eitr->face()->is_unbounded() || eitr->twin()->face()->is_unbounded())
        continue;
      int si = eitr->face()->data();
      int sj = eitr->twin()->face()->data();
      if (si == sj)
        continue;
      else if (si > sj)
        std::swap(si, sj);
      const double elen = std::sqrt(CGAL::to_double(
        CGAL::squared_distance(eitr->target()->point(), eitr->source()->point())));
      auto find_itr = adj_map.find({ si, sj });
      if (find_itr == adj_map.end())
        adj_map.insert({{si, sj}, elen});
      else
        find_itr->second += elen;
  }
  LOG(INFO) << "Adjacent map down.";

  // set up data term
  vector<double> data(nb_sites * nb_labels, 0.0);
  for (int i = 0; i < nb_sites; ++i)
  {
      int cnum = m_arr_cdata[i].nbnum;
      int fnum = m_arr_fdata[i].nbnum; 
      double cratio = cnum==0 && fnum==0 ? 1.0: 1.0*cnum/(cnum+fnum);    
      double fratio = cnum==0 && fnum==0 ? 1.0: 1.0*fnum/(cnum+fnum);  
      for (int j = 0; j < nb_labels; ++j) 
      { 
          if(j == 0)
          { // blank
              data[i * nb_labels + j] = 
                cratio * m_arr_cdata[i].distances[0] * m_arr_cdata[i].area + fratio * m_arr_fdata[i].distances[0] * m_arr_fdata[i].area;
                continue;
          }
          int f = ceil(1.0*j/(m_cplanes.size()-1));
          int c = j%(m_cplanes.size()-1) == 0 ? (m_cplanes.size()-1) : j%(m_cplanes.size()-1);
          data[i * nb_labels + j] = 
            cratio * m_arr_cdata[i].distances[c] * m_arr_cdata[i].area + fratio * m_arr_fdata[i].distances[f] * m_arr_fdata[i].area;
      }
  }   
  LOG(INFO) << "Data term down.";

  // set up smooth
  class SmoothFn : public GCoptimization::SmoothCostFunctor {
  public:
    SmoothFn(const size_t ceil_size_, const size_t cplane_size_) : ceil_size(ceil_size_) ,  cplane_size(cplane_size_) {}

    double compute(int s1, int s2, int l1, int l2) {
      if(l1 == l2)
        return 0.0;
      else {
        if(l1 == 0 || l2 == 0)
          return 1.0;
        int f1 = ceil(l1/(cplane_size-1));
        int c1 = l1%(cplane_size-1) == 0 ? (cplane_size-1) : l1%(cplane_size-1);
        int f2 = ceil(l2/(cplane_size-1));
        int c2 = l2%(cplane_size-1) == 0 ? (cplane_size-1) : l2%(cplane_size-1);
        if(c1==c2 || f1==f2)
          return 0.5;
        else
          return 1.0;
      }
    }
  private:
    const size_t ceil_size;
    const size_t cplane_size;
  };
 
  // MRF labeling
  try {
      // set up smooth (alternative)
      vector<double> smooth(nb_labels * nb_labels, 1.0);
      for (int i = 0; i < nb_labels; ++i)
        smooth[i * nb_labels + i] = 0.0;

      GCoptimizationGeneralGraph gc(nb_sites, nb_labels);
      gc.setDataCost(data.data());
      SmoothFn smoothfn(m_arr_cdata.size(), m_cplanes.size());
      gc.setSmoothCostFunctor(&smoothfn);
      //gc.setSmoothCost(smooth.data());
      LOG(INFO) << "Smooth term down.";

      // set neighboring (weight)
      for (const auto &adj_pair : adj_map)
        gc.setNeighbors(
          adj_pair.first.first,
          adj_pair.first.second,
          adj_pair.second * config.get<double>("arr.mrf.balance"));

      LOG(INFO) << "Before optimization energy is: " << gc.compute_energy();
      if (config.get<bool>("arr.mrf.use_swap")) {
        LOG(INFO) << "Alpha-beta swap algorithm.";
        gc.swap(config.get<int>("arr.mrf.iterations"));
      }
      else {
        LOG(INFO) << "Alpha expansion algorithm.";
        gc.expansion(config.get<int>("arr.mrf.iterations"));
      }
      LOG(INFO) << "After optimization energy is: " << gc.compute_energy();

      // label
      for (auto fitr = m_carr.faces_begin(); fitr != m_carr.faces_end(); ++fitr) 
      {
        if (fitr->data() >= 0)
        {
          if(fitr->data() > nb_sites)
            LOG(INFO) << "Site error: " << fitr->data() << " " <<  nb_sites;
          if(gc.whatLabel(fitr->data()) == 0)
            fitr->set_data(0);
          else 
            fitr->set_data(gc.whatLabel(fitr->data())%(m_cplanes.size()-1) == 0 ? (m_cplanes.size()-1) : gc.whatLabel(fitr->data())%(m_cplanes.size()-1));
        }
      }
      for (auto fitr = m_farr.faces_begin(); fitr != m_farr.faces_end(); ++fitr) 
      {
        if (fitr->data() >= 0)
        {
          if(fitr->data() > nb_sites)
            LOG(INFO) << "Site error: " << fitr->data() << " " <<  nb_sites;
          if(gc.whatLabel(fitr->data()) == 0)
            fitr->set_data(0);
          else 
            fitr->set_data(ceil(gc.whatLabel(fitr->data())/(m_cplanes.size()-1.0)));
        }
      }
  }
  catch (GCException e) {
    e.Report();
  }
} 


// record height differences 
bool building::segment_arrangement_extrusion(std::vector<Plane_3>& m_planes, int tag, Arr_with_index& m_arrs) 
{
    E2I to_inexact;
    I2E to_exact;
    vector<K_epec::Point_3> points;
    vector<vector<std::size_t>> polygons;
/* save in adaptive triangle mesh
    Mesh mesh;
    ofstream file(m_dir + "non_sample_polygon");
    for (auto fitr = m_arrs.faces_begin(); fitr != m_arrs.faces_end(); ++fitr) {
        const int label = fitr->data();
        if (label <= 0 || label > m_planes.size()) 
          continue;  
        // extrude ceiling/floor planes
        // points
        const std::size_t offset = points.size();
        std::size_t nb_vertices = 0;
        auto curr = fitr->outer_ccb();
        // polygon partitioning
        vector<Point_2> pp;
        do{
          pp.push_back(to_geom2(curr->target()->point()));
          ++curr;
        }while(curr != fitr->outer_ccb());
        auto it = pp.begin();
        auto prev(it++);
        // delete the same neighbor points
        for(; it != pp.end();){
          if(*it == *prev){
            prev = pp.erase(prev);
            it = prev+1;
          }
          else{
            ++it,++prev;
          }
        }
        if(*(pp.begin()) == *prev)
          prev = pp.erase(prev);
        if(pp.size() < 3)
          continue;
        // guarantee1
        if(!CGAL::is_simple_2(pp.begin(), pp.end())){
          LOG(INFO) << "Non_simple polygon.";
          for(int i = 1; i < pp.size(); i++){
            file << pp[i-1].x() << " " << pp[i].x() << " " << pp[i-1].y() << " " << pp[i].y();
          }
          file.close();
          return false;
          continue;
        }
        // guarantee2
        if(CGAL::orientation_2(pp.begin(), pp.end()) != CGAL::COUNTERCLOCKWISE){
          std::reverse(pp.begin(), pp.end());
        }
        Polygon_list polys;
        CGAL::approx_convex_partition_2(pp.begin(),
                                      pp.end(),
                                      std::back_inserter(polys));
        for (Polygon_list::iterator it = polys.begin(); it != polys.end(); ++it)
        {
            Triangulation t;
            t.insert(it->vertices_begin(), it->vertices_end());
              for(auto fit : t.finite_face_handles()){
                int fidx[3];
                vector<std::size_t> faces;
                for(int i = 0; i < 3; i++){
                  K_epec::Point_3 p(fit->vertex(i)->point().x(),
                    fit->vertex(i)->point().y(),
                      0.0);
                  auto result = CGAL::intersection(
                    to_geom(m_planes[label]),
                    K_epec::Line_3(p, K_epec::Vector_3(0.0, 0.0, 1.0)));
                  const K_epec::Point_3 *intersect_pt = nullptr;
                  if (result && (intersect_pt = boost::get<K_epec::Point_3>(&(*result))))
                    p = *intersect_pt;
                  else
                    LOG(INFO) << "no intersection";
                  auto it = find(points.begin(),points.end(),p);
                    faces.push_back(it-points.begin());
                  if (it == points.end()){
                    points.push_back(p); 
                    mesh.vertices.push_back(to_geom2(p));
                  }
                }
                Triangle_3 tri(to_geom2(points[faces[0]]), to_geom2(points[faces[1]]), to_geom2(points[faces[2]]));
                if (tri.is_degenerate())
                    continue;
                if(tag == 0)
                  mesh.faces.push_back({faces[0],faces[1],faces[2]});
                else
                  mesh.faces.push_back({faces[1],faces[0],faces[2]});
              }
        }
    }
    file.close();
    
    // .ply
    string fname;
    if(tag == 0)
     fname = m_dir + "floor_vec.ply";
    else
     fname = m_dir + "ceiling_vec.ply";
    const auto &config = cm::get_config();
    if(!PLYTriMeshSave(mesh, config, fname))
      LOG(INFO) << "Ceiling/Floor failed to save.";
*/
     
    if(tag == 0)
      LOG(INFO) << "Generate floor...";
    else
      LOG(INFO) << "Generate ceiling...";
    // save in .off
    for (auto fitr = m_arrs.faces_begin(); fitr != m_arrs.faces_end(); ++fitr) 
    {
        const int label = fitr->data();
        if (label <= 0 || label > m_planes.size())
          continue;
        // extrude ceiling/floor planes
        // points
        const std::size_t offset = points.size();
        std::size_t nb_vertices = 0;
        // polygons
        std::vector<std::size_t> roof_plg;
        auto curr = fitr->outer_ccb();
        do {
          K_epec::Point_3 p(curr->target()->point().x(),
            curr->target()->point().y(),
              0.0);
          auto result = CGAL::intersection(
            to_exact(m_planes[label]),
            K_epec::Line_3(p, K_epec::Vector_3(0.0, 0.0, 1.0)));
          const K_epec::Point_3 *intersect_pt = nullptr;
          if (result && (intersect_pt = boost::get<K_epec::Point_3>(&(*result))))
            p = *intersect_pt;
          else
            LOG(INFO) << "No intersection.";
          auto it = find(points.begin(),points.end(),p);
          roof_plg.push_back(it-points.begin());
          if (it == points.end())
              points.push_back(p);
          ++curr;
          ++nb_vertices;
        } while (curr != fitr->outer_ccb());
        polygons.push_back(roof_plg);
    }
    std::ofstream ofs;
    if(tag == 0)
      ofs.open(m_dir + "floor.off");
    else
      ofs.open(m_dir + "ceiling.off");
    ofs << "COFF\n" << points.size() << ' ' << polygons.size() << ' ' << "0\n";
    for (const auto &p : points)
      ofs << setiosflags(ios::fixed) << setprecision(8) << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
    for (const auto &plg : polygons) 
    {
      ofs << plg.size();
      for (int i = 0; i < plg.size(); i++)
        if(tag == 0)
          ofs << ' ' << plg[i];
        else
          ofs << ' ' << plg[plg.size()-1-i];
        ofs << '\n';
    }
    ofs.close();
    if(tag == 0)
      LOG(INFO) << "Generate floor down.";
    else
      LOG(INFO) << "Generate ceiling down";

    // deal with inconsistant height (c-c, f-f) and store them in <facades> which is the final wall model
    //some non-facade edges cause height inconsistance
    int idx = -1; // index in m_facade
    for (auto eitr = m_arrs.edges_begin(); eitr != m_arrs.edges_end(); ++eitr) 
    {
        int l1 = eitr->face()->data();
        int l2 = eitr->twin()->face()->data(); // l1 != l2
        K_epec::Point_3 ps(eitr->source()->point().x(),
          eitr->source()->point().y(),
          0.0);
        K_epec::Point_3 pt(eitr->target()->point().x(),
          eitr->target()->point().y(),
          0.0);
        // non-facade height inconsistance due to different ceiling/floor labels
        if(facade_tags[std::distance(m_arrs.edges_begin(), eitr)] == 0)
        { // non-facade
            if(l1 <= 0 || l1 > m_planes.size() || l2 <= 0 || l2 > m_planes.size())
                continue;
            if(l1 != l2)
            {
                K_epec::Point_3 *p1s = nullptr, *p1t = nullptr, *p2s = nullptr, *p2t = nullptr;                             
                Facade f;
                f.seg = K_epec::Segment_2(eitr->source()->point(), eitr->target()->point());
                auto result1 = CGAL::intersection(
                  to_exact(m_planes[l1]),
                  K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result1 && (p1s = boost::get<K_epec::Point_3>(&(*result1)))){
                    f.heights.push_back(CGAL::to_double(p1s->z()));
                }
                else{
                    LOG(INFO) << "Extrusion error2.";
                    return false;    
                }
                auto result2 = CGAL::intersection(
                  to_exact(m_planes[l2]),
                  K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result2 && (p2s = boost::get<K_epec::Point_3>(&(*result2)))){
                    f.heights.push_back(CGAL::to_double(p2s->z()));                
                }
                else{
                    LOG(INFO) << "Extrusion error2.";
                    return false;    
                }          
                auto result3 = CGAL::intersection(
                  to_exact(m_planes[l1]),
                  K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result3 && (p1t = boost::get<K_epec::Point_3>(&(*result3)))){
                    f.heights.push_back(CGAL::to_double(p1t->z()));  
                }
                else{
                    LOG(INFO) << "Extrusion error2.";
                    return false;    
                }      
                auto result4 = CGAL::intersection(
                  to_exact(m_planes[l2]),
                  K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result4 && (p2t = boost::get<K_epec::Point_3>(&(*result4)))){
                    f.heights.push_back(CGAL::to_double(p2t->z()));                 
                }
                else{
                    LOG(INFO) << "Extrusion error2.";
                    return false;    
                }
                if(f.heights[0] < f.heights[1]){ // height: sh sl th tl
                    std::swap(f.heights[0], f.heights[1]);
                    std::swap(f.heights[2], f.heights[3]);
                }
                Plane_3 n;
                // decide the normal
                if(tag == 0) // floor
                  n = Plane_3(to_inexact(*p2t), to_inexact(*p1t), to_inexact(*p1s));
                else // ceiling
                  n = Plane_3(to_inexact(*p2t), to_inexact(*p2s), to_inexact(*p1s));  
                f.normal = Vector_2(n.orthogonal_vector().x(), n.orthogonal_vector().y());
                facades.push_back(f);
             }
        }
        else
        { // facade 
            idx++;
            // push facade height value
            if((l1 <= 0 || l1 > m_planes.size()) && (l2 <= 0 || l2 > m_planes.size()))
              continue;
            K_epec::Point_3 *p1s = nullptr, *p1t = nullptr, *p2s = nullptr, *p2t = nullptr;
              
            if(!(l1 <= 0 || l1 > m_planes.size()))
            {
                auto result = CGAL::intersection(
                  to_exact(m_planes[l1]),
                  K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result && (p1s = boost::get<K_epec::Point_3>(&(*result))))
                {
                    auto it = find(points.begin(), points.end(),*p1s);
                    if(it==points.end()){
                        points.push_back(*p1s);
                    }                       
                }
                auto result2 = CGAL::intersection(
                  to_exact(m_planes[l1]),
                  K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));                
                if (result2 && (p1t = boost::get<K_epec::Point_3>(&(*result2))))
                {
                    auto it = find(points.begin(), points.end(),*p1t);
                    if(it==points.end()){
                        points.push_back(*p1t);
                    }                       
                }
                Point_2 z(CGAL::to_double((*p1s).z()), CGAL::to_double((*p1t).z())); //{source, target}
                if(tag == 0) {// floor
                    m_facade[idx].bottoms.push_back(z);
                }
                else // ceiling
                    m_facade[idx].ups.push_back(z);
            }
            
            if(!(l2 <= 0 || l2 > m_planes.size()))
            {
                auto result = CGAL::intersection(
                  to_exact(m_planes[l2]),
                  K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result && (p2s = boost::get<K_epec::Point_3>(&(*result))))
                {
                    auto it = find(points.begin(), points.end(),*p2s);
                    if(it==points.end()){                     
                        points.push_back(*p2s);
                    }                       
                }
                auto result2 = CGAL::intersection(
                  to_exact(m_planes[l2]),
                  K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result2 && (p2t = boost::get<K_epec::Point_3>(&(*result2))))
                {
                    auto it = find(points.begin(), points.end(),*p2t);
                    if(it==points.end()){                       
                        points.push_back(*p2t);
                    }                       
                }
                Point_2 z(CGAL::to_double((*p2s).z()), CGAL::to_double((*p2t).z()));
                if(tag == 0) {// floor
                  auto it = find(m_facade[idx].bottoms.begin(), m_facade[idx].bottoms.end(), z);
                  if(it == m_facade[idx].bottoms.end())
                    m_facade[idx].bottoms.push_back(z);
                }
                else{
                  auto it = find(m_facade[idx].ups.begin(), m_facade[idx].ups.end(), z);
                  if(it == m_facade[idx].ups.end())
                    m_facade[idx].ups.push_back(z);
                } // ceiling  
                if(l1 > 0 && l1 < m_planes.size() && l1 == l2)
                  m_facade[idx].tag = 1;
            }

        }
    }
    return true;
}


// extrude facade -> m_facade
bool building::facade_extrusion()
{

    E2I to_inexact;
    I2E to_exact;
    const auto &config = cm::get_config();
    ifstream file(m_dir + "3dpointsindex", ios::in);
    int segs_num; 
    int points_num;
    file >> segs_num;
    double x,y,z;
    if(segs_num != m_facade_segs.size()) {
      LOG(INFO) << "facade record error.";
      //return false;
    } 
    vector<vector<Point_3>> points(segs_num);
    vector<Vector_2> normals;
    for(auto &p_set: points){
        file >> points_num;
        file >> x >> y;
        normals.push_back(Vector_2(x,y));     
        for(int i = 0; i < points_num; i++){
            file >> x >> y >> z;
            p_set.push_back(Point_3(x,y,z));
        }
    } 
   
    // for linear program
    vector<Facade> potential_facades;
    set<string> non_overlap;
    vector<int> suppoints_num;
    vector<double> area_ratio;
    double area_total = 0;

    for(int i = 0; i < m_facade.size(); i++)
    {
        auto facade = m_facade[i];
        if(facade.bottoms.size() == 0 || facade.ups.size() == 0)
          continue;

        facade.normal = normals[facade.index];
        vector<Point_2> height;
        for(auto p:facade.ups)
          height.push_back(p); // height: {source, target}
        for(auto p:facade.bottoms)
          height.push_back(p);

        std::sort(height.begin(), height.end(), 
          [](Point_2 p1, Point_2 p2){ return (p1.x()+p1.y())/2 > (p2.x()+p2.y())/2;}); // decreasingly sort by average z       
        bool c2 = height.size() == 2 && facade.tag == 0; // boundary facade
        bool c3 = facade.ups.size() == 1 && facade.bottoms.size() == 2;
        bool c4 = facade.ups.size() == 2 && facade.bottoms.size() == 1;
        bool c5 = facade.ups.size() == 2 && facade.bottoms.size() == 2;
        // guarantee manifold models
        for(int i = 0; i < height.size() - 1; i++)
        {  
          if(c2 || (c3 && i == 1) || (c4 && i == 0) || (c5 && (i == 0 || i == 2))){ // push these height differences into the wall solution           
               facades.push_back(Facade(K_epec::Segment_2(facade.seg.source(), facade.seg.target()), {height[i].x(), height[i+1].x(), height[i].y(), height[i+1].y()}));
               facades[facades.size()-1].normal = facade.normal;            
            }
          else
          { // unsure height differences
              Point_2 margin(std::max(height[i].x(), height[i].y()), std::min(height[i+1].x(), height[i+1].y())); // {big, small}
              int margin_num = 0;
              unordered_set<Point_2> counts;
              K_epec::Point_2 lc = facade.seg.min();
              K_epec::Point_2 rc = facade.seg.max();
              double len = sqrt(CGAL::to_double((rc-lc).squared_length()));
              if(len < 1e-3){ // too small segments, just for connection
                potential_facades.push_back(Facade(K_epec::Segment_2(facade.seg.source(), facade.seg.target()), {height[i].x(), height[i+1].x(), height[i].y(), height[i+1].y()}));
                potential_facades[potential_facades.size()-1].normal = facade.normal;
                suppoints_num.push_back(0);
                area_ratio.push_back(0);    
                continue;       
              }
              int cols = std::ceil(len/m_step);
              int rows = std::ceil((margin.x() - margin.y())/m_step);
              Line_2 l1 = to_inexact(facade.seg.supporting_line().perpendicular(facade.seg.min()));
              Line_2 l2 = to_inexact(facade.seg.supporting_line().perpendicular(facade.seg.max()));

              for(auto& p : points[facade.index])
              {
                  K_epec::Point_2 pp = facade.seg.supporting_line().projection(K_epec::Point_2(p.x(), p.y()));
                  double pl_len = sqrt(CGAL::to_double((pp - lc).squared_length()));
                  if(pp < lc && pl_len > 1e-4 || pp > rc && pl_len > 1e-4) // accepted error
                    continue;
                  if(pp == facade.seg.source() && CGAL::to_double(p.z()) >= height[i+1].x() && CGAL::to_double(p.z()) <= height[i].x()
                    || pp == facade.seg.target() && CGAL::to_double(p.z()) >= height[i+1].y() && CGAL::to_double(p.z()) <= height[i].y())
                  {
                      margin_num++;
                      int c = std::ceil(pl_len/m_step);
                      int r = std::ceil(CGAL::to_double(p.z()) - margin.y()/m_step);
                      if(margin.y() >= CGAL::to_double(p.z()) && margin.y() - CGAL::to_double(p.z()) < 1e-2) // accepted error
                        r = 1;
                      if(CGAL::to_double(p.z()) >= margin.x() && CGAL::to_double(p.z()) - margin.x() < 1e-2)
                        r = rows;
                      if(c > 0 && c <= cols && r > 0 && r <= rows){
                        counts.insert({c, r});
                      }
                  }
                  else if(pp == facade.seg.source() || pp == facade.seg.target())
                    continue;
                  else
                  {               
                      Segment_2 seg1(Point_2(0, height[i+1].x()),Point_2(len, height[i+1].y())); // low
                      Segment_2 seg2(Point_2(0, height[i].x()), Point_2(len, height[i].y())); // high
                      Ray_2 ray(Point_2(pl_len, p.z()), Vector_2(0, 1));                      
                      int tag = 0;
                      auto result = CGAL::intersection(seg1, ray);
                      const Point_2 *intersect_pt = nullptr;
                      if(result && (intersect_pt = boost::get<Point_2>(&(*result))))
                      {
                          if(intersect_pt->y() - p.z() < 1e-3){
                            tag = 1;
                          }
                          else
                            continue;
                      }
                      else
                      {
                          result = CGAL::intersection(seg2, ray);
                          if(result){
                            tag = 1;
                          }
                          else 
                            continue;
                      }
                      if(tag == 1)
                      { // within valid area
                          margin_num++;
                          int c = std::ceil(pl_len/m_step);
                          int r = std::ceil((CGAL::to_double(p.z()) - margin.y())/m_step);
                          if(margin.y() >= CGAL::to_double(p.z()) && margin.y() - CGAL::to_double(p.z()) < 1e-2) // accepted error
                            r = 1;
                          if(CGAL::to_double(p.z()) >= margin.x() && CGAL::to_double(p.z()) - margin.x() < 1e-2)
                            r = rows;
                          if(c > 0 && c <= cols && r > 0 && r <= rows){
                            counts.insert({c, r});
                          }
                      }
                  }
                }

                potential_facades.push_back(Facade(K_epec::Segment_2(facade.seg.source(), facade.seg.target()), {height[i].x(), height[i+1].x(), height[i].y(), height[i+1].y()}));
                potential_facades[potential_facades.size()-1].normal = facade.normal;
                suppoints_num.push_back(margin_num);
                int r1 = std::ceil(std::abs(height[i].x() - height[i].y()) / m_step);
                int r2 = std::ceil(std::abs(height[i+1].x() - height[i+1].y()) / m_step);
                area_ratio.push_back(1-(counts.size()*1.0)/(rows*cols - r1*cols/2 - r2*cols/2));
                area_total += area_ratio[area_ratio.size()-1];
            }        
        }
      }
      vector<vector<Point_3>>().swap(points);
    
    for(auto &a:area_ratio)
       a /= area_total;

  {
    // non consistant non-facade segments are also added into <facades>
    auto e2 = m_carr.edges_begin();
    for(auto e1 = m_farr.edges_begin(); e1 != m_farr.edges_end(); e1++, e2++)
    {
        const int l11 = e1->face()->data();
        const int l12 = e1->twin()->face()->data();
        const int l21 = e2->face()->data();
        const int l22 = e2->twin()->face()->data();
        if(facade_tags[std::distance(m_farr.edges_begin(),e1)] == 1)  continue; // ignore facade
        if(l11 <= 0  && l12 > 0 && l12 < m_fplanes.size() || l12 <= 0 && l11 > 0 && l11 < m_fplanes.size())
        {
            Facade f;
            // judge normal
            if(l11 <= 0 && l12 > 0 && l12 < m_fplanes.size()){
              K_epec::Vector_2 v = (K_epec::Line_2(e1->twin()->source()->point(), e1->twin()->target()->point())).perpendicular(K_epec::Point_2(0,0)).to_vector();
              f.normal = Vector_2(CGAL::to_double(v.x()),CGAL::to_double(v.y()));
            }
            else{
              K_epec::Vector_2 v = (K_epec::Line_2(e1->source()->point(), e1->target()->point())).perpendicular(K_epec::Point_2(0,0)).to_vector();
              f.normal = Vector_2(CGAL::to_double(v.x()),CGAL::to_double(v.y()));
            }
            f.seg = K_epec::Segment_2(e1->source()->point(), e1->target()->point());
            K_epec::Point_3 ps(e1->source()->point().x(),
              e1->source()->point().y(),
                0.0);
            K_epec::Point_3 pt(e1->target()->point().x(),
              e1->target()->point().y(),
                0.0);
            K_epec::Plane_3 p1, p2;
            if(l11 > 0)
              p1 = to_exact(m_fplanes[l11]);
            else
              p1 = to_exact(m_fplanes[l12]);
            if(l21 > 0)
              p2 = to_exact(m_cplanes[l21]);
            else
              p2 = to_exact(m_cplanes[l22]);          
            K_epec::Point_3 *p1s = nullptr, *p1t = nullptr, *p2s = nullptr, *p2t = nullptr;
            auto result = CGAL::intersection(
              p1,
              K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
            if (!(result && (p1s = boost::get<K_epec::Point_3>(&(*result))))){
              LOG(INFO) << "Extrusion error3_1.";
              return false;    
            }
            else
            {
              f.heights.push_back(CGAL::to_double(p1s->z()));
            }

            result = CGAL::intersection(
              p2,
              K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
            if (!(result && (p2s = boost::get<K_epec::Point_3>(&(*result))))){
              LOG(INFO) << "Extrusion error3_2.";
              return false;    
            } 
            else{
              f.heights.push_back(CGAL::to_double(p2s->z()));
            }
                              
            result = CGAL::intersection(
              p1,
              K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));               
            if (!(result && (p1t = boost::get<K_epec::Point_3>(&(*result))))){
              LOG(INFO) << "Extrusion error3_3.";
              return false;    
            }  
            else{
              f.heights.push_back(CGAL::to_double(p1t->z()));
            }  
                  
            result = CGAL::intersection(
              p2,
              K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));               
            if (!(result && (p2t = boost::get<K_epec::Point_3>(&(*result))))){
              LOG(INFO) << "Extrusion error3_4.";
              return false;    
            }     
            else{
              f.heights.push_back(CGAL::to_double(p2t->z()));
            }       
            if(f.heights[0] < f.heights[1]){ // height: sh sl th tl
              std::swap(f.heights[0], f.heights[1]);
              std::swap(f.heights[2], f.heights[3]);
            }
            facades.push_back(f); // arr.cpp     
          }
    }
  }

    /*
    facades: must be in the solution, guarantee the manifold model
    potential_facades: maybe in the solution, solved by lp
    */
    vector<IntersectPoint> inter_points;  
    vector<pair<K_epec::Point_2, K_epec::Point_2>> sub_points;
    vector<Point_2> points_height; // {high, low}
    for(int i = 0; i < potential_facades.size(); i++)
    {
        auto f = potential_facades[i];
        auto it = find(sub_points.begin(), sub_points.end(), 
          make_pair(f.seg.source(), K_epec::Point_2(f.heights[0], f.heights[1])));
        if(it == sub_points.end()){
          sub_points.push_back(make_pair(f.seg.source(), K_epec::Point_2(f.heights[0], f.heights[1]))); // s s-high s-low
          inter_points.push_back(
            IntersectPoint(to_inexact(f.seg.source()), Point_2(f.heights[0], f.heights[1]), {i}));
        }
        else{
          inter_points[std::distance(sub_points.begin(), it)].IDs.push_back(i);
        }
        it = find(sub_points.begin(), sub_points.end(),
          make_pair(f.seg.target(), K_epec::Point_2(f.heights[2], f.heights[3])));
        if(it == sub_points.end()){
          sub_points.push_back(make_pair(f.seg.target(), K_epec::Point_2(f.heights[2], f.heights[3]))); // t t-high t-low
          inter_points.push_back(
            IntersectPoint(to_inexact(f.seg.target()), Point_2(f.heights[2], f.heights[3]), {i}));
        }
        else{
          inter_points[std::distance(sub_points.begin(), it)].IDs.push_back(i);
        }
    }
    vector<pair<K_epec::Point_2, K_epec::Point_2>>().swap(sub_points);

    // judge inter_point.flag
    for(auto &ip:inter_points)
      for(auto &s:facades){
        if(point_segment_squared_distance(ip.p.x(), ip.p.y(), 
          CGAL::to_double(s.seg.source().x()), CGAL::to_double(s.seg.source().y()), 
          CGAL::to_double(s.seg.target().x()), CGAL::to_double(s.seg.target().y())) < 1e-8)
          {         
            K_epec::Point_2 p = s.seg.supporting_line().projection(to_exact(ip.p));
            int idx_max = s.heights[0] > s.heights[2] ? 0 : 1;
            int idx_min = s.heights[0] > s.heights[2] ? 1 : 0;
            double len_max = sqrt(CGAL::to_double(s.seg.squared_length()));
            double len_min = sqrt(CGAL::to_double((s.seg[idx_min] - p).squared_length())); 
            double hmax = s.heights[2*idx_min] + (s.heights[2*idx_max] - s.heights[2*idx_min]) * len_min / len_max;
            idx_max = s.heights[1] > s.heights[3] ? 0 : 1;
            idx_min = s.heights[1] > s.heights[3] ? 1 : 0;
            double hmin = s.heights[2*idx_min+1] + (s.heights[2*idx_min+1] - s.heights[2*idx_min+1]) * len_min / len_max;
            if(abs(ip.heights.x() - hmax) < 1e-3 && abs(ip.heights.y() - hmin) < 1e-3){
              ip.flag = 1; // boundary point
              break;
            }
        }    
      }

    // record inter_points
    ofstream fb(m_dir + "boundary_points");  
    ofstream fnb(m_dir + "non-boundary_points");  
    int boundarynums=0;
    for(auto ip:inter_points){
      if(ip.flag == 1){
        boundarynums++;
        fb << ip.p << "\n";
      }
      else
        fnb << ip.p << "\n";
    }
    fb.close();
    fnb.close();
    LOG(INFO) << boundarynums << " boundary points. " << inter_points.size()-boundarynums << " non-boundary points.";

    // lp
    double lambda_data_fitting = cm::get_config().get<double>("facade.height.lamda_data_fitting");
    double lambda_model_coverage = cm::get_config().get<double>("facade.height.lamda_model_coverage");
    double lambda_model_complexity = cm::get_config().get<double>("facade.height.lamda_model_complexity");
    optimize2(potential_facades, inter_points, suppoints_num, area_ratio, lambda_data_fitting, lambda_model_coverage, lambda_model_complexity, m_dir);
    
    // record results
    ofstream f1(m_dir + "potential_facades");
    ofstream f2(m_dir + "selected_facades");  
    for(auto f:potential_facades){
      f1 << setiosflags(ios::fixed) << setprecision(8) << f.seg.source().x() << " " << f.seg.target().x() << " " << f.seg.source().y() << " " << f.seg.target().y() << "\n";
      if(f.flag == 1)
      f2 << setiosflags(ios::fixed) << setprecision(8) << f.seg.source().x() << " " << f.seg.target().x() << " " << f.seg.source().y() << " " << f.seg.target().y() << "\n";
    }
    f1.close();
    f2.close(); 

    ofstream f3(m_dir + "pre_facades");
    for(auto f:facades)
      f3 << setiosflags(ios::fixed) << setprecision(8) << f.seg.source().x() << " " <<  f.seg.target().x() << " " << f.seg.source().y() << " " << f.seg.target().y() << "\n";
    f3.close();  

    // extrude facade
    std::vector<K_epec::Point_3> fpoints;
    std::vector<std::vector<std::size_t>> polygons;
    std::size_t id[4];
    K_epec::Point_3 p[4];
    for(auto f:facades){
      p[0] = K_epec::Point_3(f.seg.source().x(), f.seg.source().y(), f.heights[0]);  //sh
      p[1] = K_epec::Point_3(f.seg.source().x(), f.seg.source().y(), f.heights[1]);  // sl
      p[2] = K_epec::Point_3(f.seg.target().x(), f.seg.target().y(), f.heights[2]);  // th 
      p[3] = K_epec::Point_3(f.seg.target().x(), f.seg.target().y(), f.heights[3]);  // tl
      for(int i = 0; i < 4; i++){
        fpoints.push_back(p[i]);
        id[i] = fpoints.size()-1;
      }
      // accordng to normal
      K_epec::Plane_3 n(p[0],p[1],p[3]);
      K_epec::Vector_2 n2(n.orthogonal_vector().x(), n.orthogonal_vector().y());
      if(n2 * to_exact(f.normal) >= 0 )
        polygons.push_back({id[0], id[1], id[3], id[2]});
      else
        polygons.push_back({id[0], id[2], id[3], id[1]});

    }

    for(auto f:potential_facades){
      if(f.flag != 1) continue;
      p[0] = K_epec::Point_3(f.seg.source().x(), f.seg.source().y(), f.heights[0]);  //sh
      p[1] = K_epec::Point_3(f.seg.source().x(), f.seg.source().y(), f.heights[1]);  // sl
      p[2] = K_epec::Point_3(f.seg.target().x(), f.seg.target().y(), f.heights[2]);  // th 
      p[3] = K_epec::Point_3(f.seg.target().x(), f.seg.target().y(), f.heights[3]);  // tl
      for(int i = 0; i < 4; i++){
        fpoints.push_back(p[i]);
        id[i] = fpoints.size()-1;
      }
      // accordng to normal
      K_epec::Plane_3 n(p[0],p[1],p[3]);
      K_epec::Vector_2 n2(n.orthogonal_vector().x(), n.orthogonal_vector().y());
      if(n2 * to_exact(f.normal) >= 0 )
         polygons.push_back({id[0], id[1], id[3], id[2]});
      else
         polygons.push_back({id[0], id[2], id[3], id[1]});
    } 

    string ofs = m_dir + "facade.off";
    Mesh<K_epec::Point_3> mesh;
    for(auto& p: fpoints)
      mesh.vertices.push_back(p);
    for(auto& poly: polygons){
      for(int i = 2; i < poly.size(); i++){
          mesh.faces.push_back(K_epec::Point_3(poly[0], poly[i-1], poly[i]));
      }
    }
    OFFMeshSave(ofs, mesh);
/*  save in adaptive triangle mesh
    string ofs = m_dir + "facade_vec.ply";
    Mesh mesh;
    for(auto& p: fpoints)
      mesh.vertices.push_back(p);
    for(auto& poly: polygons){
      for(int i = 2; i < poly.size(); i++){
          mesh.faces.push_back({poly[0], poly[i-1], poly[i]});
      }
    }
    if(!PLYTriMeshSave(mesh, config, ofs))
       LOG(INFO) << "Facade saves failed.";
*/
    return true;
}

void building::write_all_facades()
{
    E2I to_inexact;
    cv::Mat img_out(m_rows, m_cols, CV_8UC3, cv::Scalar(255, 255, 255));
    for (auto it = m_farr.edges_begin(); it != m_farr.edges_end(); it++) 
    {
      if(facade_tags[std::distance(m_farr.edges_begin(), it)]==0) continue;
      Segment_2 s(to_inexact(it->source()->point()), to_inexact(it->target()->point()));
      const cv::Point2i ps((int)((CGAL::to_double(s.source().x()) - m_corner.x())/m_step), (int)((m_corner.y() - CGAL::to_double(s.source().y()))/m_step));
      const cv::Point2i pt((int)((CGAL::to_double(s.target().x()) - m_corner.x())/m_step), (int)((m_corner.y() - CGAL::to_double(s.target().y()))/m_step));
      cv::line(img_out, ps, pt, cv::Scalar(0, 0, 0));
    }
    const string fname = m_dir + "facades_all.png";
    if (!cv::imwrite(fname, img_out))
      throw std::runtime_error("Failed to write file " + fname);
}

void building::write_all_segments() 
{
  cv::Mat img_out(m_rows, m_cols, CV_8UC3, cv::Scalar(255, 255, 255));
  for (const auto &s : m_segments) {
    const cv::Point2i ps((int)((CGAL::to_double(s.source().x()) - m_corner.x())/m_step), (int)((m_corner.y() - CGAL::to_double(s.source().y()))/m_step));
    const cv::Point2i pt((int)((CGAL::to_double(s.target().x()) - m_corner.x())/m_step), (int)((m_corner.y() - CGAL::to_double(s.target().y()))/m_step));
    cv::line(img_out, ps, pt, cv::Scalar(0, 0, 0));
  }
  const std::string fname = m_dir + "segments_all.png";
  if (!cv::imwrite(fname, img_out))
    throw std::runtime_error("Failed to write file " + fname);
}


void building::write_arrangement() 
{
  cv::Mat img_out(m_rows, m_cols, CV_8UC3, cv::Scalar(255.0, 255.0, 255.0));
  for (auto eitr = m_arr.edges_begin(); eitr != m_arr.edges_end(); ++eitr) 
  {
    const cv::Point2i ps(
      int((CGAL::to_double(eitr->source()->point().x())-m_corner.x())/m_step),
      int((m_corner.y()-CGAL::to_double(eitr->source()->point().y()))/m_step));
    const cv::Point2i pt(
      int((CGAL::to_double(eitr->target()->point().x())-m_corner.x())/m_step),
      int((m_corner.y()-CGAL::to_double(eitr->target()->point().y()))/m_step));
    cv::line(img_out, ps, pt, cv::Scalar(0.0, 0.0, 0.0));
  }
  const std::string fname = m_dir + "arrangment.png";
  if (!cv::imwrite(fname, img_out))
    throw std::runtime_error("Failed to write file " + fname);
}

void building::write_facade_segments()
{
    ofstream file(m_dir + "facade_segments");
    for(auto facade:m_facade)
       file << setiosflags(ios::fixed) << setprecision(8) << facade.seg.source().x() << " " << facade.seg.target().x() << " " << facade.seg.source().y() << " " << facade.seg.target().y() << "\n";
    file.close();file.close();
}

void building::write_arrangement_labeling(vector<Plane_3>& m_planes, int tag, Arr_with_index& m_arrs) 
{
  cv::Mat img_out = cv::Mat::zeros(m_rows, m_cols, CV_8U);
  for (auto fitr = m_arrs.faces_begin(); fitr != m_arrs.faces_end(); ++fitr) 
  {
    if (fitr->data() <= 0)
      continue;
    vector<cv::Point2i> plg;
    auto curr = fitr->outer_ccb();
    do {
      plg.push_back({
        int((CGAL::to_double(curr->target()->point().x())-m_corner.x())/m_step),
        int((m_corner.y() - CGAL::to_double(curr->target()->point().y()))/m_step)});
      ++curr;
    } while (curr != fitr->outer_ccb());
    const cv::Point2i *pts[] = {plg.data()};
    const int npts[] = {int(plg.size())};
    cv::fillPoly(img_out, pts, npts, 1, cv::Scalar(255.0 * fitr->data() / m_planes.size()));
  }
  cv::applyColorMap(img_out, img_out, cv::COLORMAP_JET);
  string fname;
  if(tag == 0)
     fname = m_dir + "floor_mrf.png";
  else
     fname = m_dir + "ceiling_mrf.png";
  if (!cv::imwrite(fname, img_out))
    throw std::runtime_error("Failed to write file " + fname);
}