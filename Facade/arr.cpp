#include <fstream>
#include <limits.h>
#include <cmath>

//#include <CGAL/Point_set_3.h>
//#include <CGAL/Point_set_3/IO.h>
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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K> Triangulation;
typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Point  Point;

typedef CGAL::Partition_traits_2<K>                         Traits;
typedef Traits::Polygon_2                                   Polygon2;
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

bool building::segment_arrangement_modeling() {
  // timer
  clock_t start, mid, end;
  start = clock(); 
  
  const auto &config = cm::get_config();
  bool r1, r2, r3;
  if(m_segments.size() == 0 || m_facade_segs.size() == 0){
    LOG(INFO) << "Segments error.";
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
  
  // floor MRF tag = 0, ceiling MRF tag = 1, adjust to uniform labels
  {
    // MRF labeling
    //segment_arrangement_MRF_labeling(m_arr_fdata, m_fplanes, m_farr);
    //segment_arrangement_MRF_labeling(m_arr_cdata, m_cplanes, m_carr);
    // label adjust
    //arrangement_adjust();
    
    // the whole MRF
    segment_arrangement_MRF_labeling_cf(m_arr_cdata, m_arr_fdata, m_fplanes, m_cplanes, m_carr, m_farr);
  }

  

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

void building::construct_segment_arrangement() {
  
  // grid_arr and segments_arr ->  overlay_arr data 
  // fidx: segment_arr polygon. pixel: grid index
  struct Overlay_data {
    Overlay_data() : fidx(0), pixel(0, 0) {}
    Overlay_data(const int &f, const cv::Vec2i &p) : fidx(f), pixel(p) {}
    Overlay_data(const Overlay_data &obj) : fidx(obj.fidx), pixel(obj.pixel) {}
    int fidx;
    cv::Vec2i pixel;
  };
  struct Overlay_label {
    Overlay_data operator() (const int &f, const cv::Vec2i &p) const {
      return {f, p};
    }
  };
  typedef CGAL::Arr_face_extended_dcel<Arr_traits, cv::Vec2i> Dcel_grid;
  typedef CGAL::Arrangement_2<Arr_traits, Dcel_grid> Arr_grid;
  typedef CGAL::Arr_face_extended_dcel<Arr_traits, Overlay_data> Dcel_overlay;
  typedef CGAL::Arrangement_2<Arr_traits, Dcel_overlay> Arr_overlay;
  typedef CGAL::Arr_face_overlay_traits<
    Arr_with_index, Arr_grid, Arr_overlay, Overlay_label> Overlay_traits;

  // grid segments, in world coordinate
  std::vector<K_epec::Segment_2> grid_segs;
  for (int i = 0; i < m_rows + 1; ++i)
    grid_segs.push_back({
      K_epec::Point_2(m_corner.x(), m_corner.y() - i*m_step),
      K_epec::Point_2(m_corner.x() + m_cols*m_step, m_corner.y() - i*m_step)});
  for (int j = 0; j < m_cols + 1; ++j)
    grid_segs.push_back({
      K_epec::Point_2(m_corner.x() + j*m_step, m_corner.y()),
      K_epec::Point_2(m_corner.x() + j*m_step, m_corner.y() - m_rows*m_step)});

  // construct sampling grid arrangement
  LOG(INFO) << "Construct sampling grid arrangement";
  LOG(INFO) << "#grid_segs: " << grid_segs.size();
  Arr_grid grid_arr;
  CGAL::insert(grid_arr, grid_segs.begin(), grid_segs.end());
  LOG(INFO) << "#f: " << grid_arr.number_of_faces()
    << ", #v: " << grid_arr.number_of_vertices();
  for (auto fitr = grid_arr.faces_begin(); fitr != grid_arr.faces_end(); ++fitr) {
    if (fitr->is_unbounded()) {
      fitr->set_data(cv::Vec2i(-1, -1));
      continue;
    }
    // find lower left corner point in world coordinate
    // it is the upper left corner in image coordinate
    auto curr = fitr->outer_ccb();
    K_epec::Point_2 ll_pt = curr->target()->point();
    do {
      if (curr->target()->point() < ll_pt)
        ll_pt = curr->target()->point();
      ++curr;
    } while (curr != fitr->outer_ccb());
    // set face row / column in image: world x, y -> image x, y -> image i, j
    const int ridx = int(std::round((m_corner.y() - CGAL::to_double(ll_pt.y()))/m_step)) - 1;
    const int cidx = int(std::round(CGAL::to_double((ll_pt.x()) - m_corner.x())/m_step));
    fitr->set_data(cv::Vec2i(ridx, cidx));
  }

  // line segment arrangement
  LOG(INFO) << "Construct line segment arrangement";
  LOG(INFO) << "#segments: " << m_segments.size();
  m_arr.clear();
  CGAL::insert(m_arr, m_segments.begin(), m_segments.end()); 
  LOG(INFO) << "#f: " << m_arr.number_of_faces()
    << ", #v: " << m_arr.number_of_vertices();

  int idx = 0;
  for (auto fitr = m_arr.faces_begin(); fitr != m_arr.faces_end(); ++fitr) {
    if (fitr->is_unbounded())
      fitr->set_data(-1);
    else
      fitr->set_data(idx++);
  }

  // trim arrangement
  LOG(INFO) << "Trim arrangement";
  for (auto eitr = m_arr.edges_begin(); eitr != m_arr.edges_end(); ++eitr) {
    if (eitr->face()->data() == eitr->twin()->face()->data()) {
      m_arr.remove_edge(eitr);
      eitr = m_arr.edges_begin();
    }
  }
  LOG(INFO) << "#f: " << m_arr.number_of_faces()
    << ", #v: " << m_arr.number_of_vertices();
  // merge collinear edges
  LOG(INFO) << "Merge collinear edges";
  for (auto fitr = m_arr.faces_begin(); fitr != m_arr.faces_end(); ++fitr) {
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
  LOG(INFO) << "#f: " << m_arr.number_of_faces()
    << ", #v: " << m_arr.number_of_vertices();

  // overlay
  LOG(INFO) << "Overlaying arrangement";
  Arr_overlay overlay_arr;
  Overlay_traits overlay_traits;
  CGAL::overlay(m_arr, grid_arr, overlay_arr, overlay_traits);
  LOG(INFO) << "#f: " << overlay_arr.number_of_faces()
    << ", #v: " << overlay_arr.number_of_vertices();

  // handle with floor ->  m_arr_fdata
  {
    // fixed: 是否正确赋值
   m_farr = m_arr;
  // compute pixel distances to planes
  std::vector<cv::Mat> pixel_distance_to_planes;
  for (const auto &s : m_fplanes) {
    cv::Mat dis_map(m_rows, m_cols, CV_64F, cv::Scalar(0.0));
    for (int i = 0; i < m_rows; ++i) {
      for (int j = 0; j < m_cols; ++j) {
        const cv::Vec3d pos = m_hfmap.at<cv::Vec3d>(i, j);
        auto result = CGAL::intersection(
          s,
          Line_3(Point_3(pos[0], pos[1], pos[2]), Vector_3(0.0, 0.0, 1.0)));
        const Point_3 *intersect_pt = nullptr;
        if (result && (intersect_pt = boost::get<Point_3>(&(*result)))) {
          Point_3 p = *intersect_pt;
          dis_map.at<double>(i, j) = std::abs( p.z() - pos[2] );
            //std::sqrt(CGAL::squared_distance(p, Point_3(pos[0], pos[1], pos[2])));
        }

        else
          LOG(INFO) << "no intersection";
      }
    }
    double zmin=0.0,zmax=0.0;
    cv::minMaxLoc(dis_map,&zmin,&zmax);
    LOG(INFO) << "distance to plane " << pixel_distance_to_planes.size()+1 << " zmin: " << zmin << " zmax: " << zmax;
    pixel_distance_to_planes.push_back(dis_map);
  }

  // set segment arrangement face data
  m_arr_fdata = std::vector<Arr_face_data>(
    m_arr.number_of_faces() - m_arr.number_of_unbounded_faces(),
    Arr_face_data(m_fplanes.size(),0));
  // record each face_area black pixel number
  map< int, pair<int, int> > record_pixel; // <black, non-black>
  for (auto fitr = overlay_arr.faces_begin(); fitr != overlay_arr.faces_end(); ++fitr) {
    const int fidx = fitr->data().fidx;
    const cv::Vec2i pixel = fitr->data().pixel;
    if (fitr->is_unbounded() || fidx < 0)
      continue;
    if(record_pixel.count(fidx) == 0)
       record_pixel[fidx] = {0, 0};
    if (pixel[0] < 0 || pixel[1] < 0)
        record_pixel[fidx].first++;
    else
        m_nfmap.at<int>(pixel) == 0 ? record_pixel[fidx].first++ : record_pixel[fidx].second++;
     
  }

  for (auto fitr = overlay_arr.faces_begin(); fitr != overlay_arr.faces_end(); ++fitr) {
    const int fidx = fitr->data().fidx;
    const cv::Vec2i pixel = fitr->data().pixel;
    if (fitr->is_unbounded() || fidx < 0)
      continue;
    Arr_face_data &arrd = m_arr_fdata[fidx];

    CGAL::Polygon_2<K_epec> plg;
    auto curr = fitr->outer_ccb();
    do {
      plg.push_back(curr->target()->point());
      ++curr;
    } while (curr != fitr->outer_ccb());
    const double area = CGAL::to_double(plg.area());
    assert(area > 0.0);
    arrd.area += area;

    if (pixel[0] < 0 || pixel[1] < 0) {
      // faces out side sampling grid is assigned with black-like data
      for (auto &d : arrd.distances)
        d = 1e20;
      arrd.distances.front() = 0.0;
      arrd.pnumber = 0;
      //arrd.normal = cv::Vec3d(0.0, 0.0, 1.0);
    }
    else {
      double black_ratio = 1.0*record_pixel[fidx].first/(record_pixel[fidx].first+record_pixel[fidx].second);
      double nblack_ratio = 1.0*record_pixel[fidx].second/(record_pixel[fidx].first+record_pixel[fidx].second);
      for (std::size_t i = 0; i < pixel_distance_to_planes.size(); ++i)
        arrd.distances[i] += area * (m_nfmap.at<int>(pixel) == 0 ? black_ratio : nblack_ratio) * pixel_distance_to_planes[i].at<double>(pixel);
      arrd.pnumber += m_nfmap.at<int>(pixel);
      arrd.bnum = record_pixel[fidx].first;
      arrd.nbnum =record_pixel[fidx].second;
    }
  }
    for (auto &arrd : m_arr_fdata) {
        for (auto &d : arrd.distances)
        d /= arrd.area;
    }

  // record distance to planes
  ofstream f(m_dir + "arr2fplane_info");
  for(auto arr:m_arr_fdata){
     f << arr.area << " ";
  for(auto dis:arr.distances)
     f << dis << " ";
     f << "\n";
  }
  }



// handle with ceiling ->  m_arr_cdata
  {
    m_carr = m_arr; // fixed: 浅拷贝
  // compute pixel distances to planes
  std::vector<cv::Mat> pixel_distance_to_planes;
  for (const auto &s : m_cplanes) {
    cv::Mat dis_map(m_rows, m_cols, CV_64F, cv::Scalar(0.0));
    for (int i = 0; i < m_rows; ++i) {
      for (int j = 0; j < m_cols; ++j) {
        const cv::Vec3d pos = m_hcmap.at<cv::Vec3d>(i, j);
        auto result = CGAL::intersection(
          s,
          Line_3(Point_3(pos[0], pos[1], pos[2]), Vector_3(0.0, 0.0, 1.0)));
        const Point_3 *intersect_pt = nullptr;
        if (result && (intersect_pt = boost::get<Point_3>(&(*result)))) {
          Point_3 p = *intersect_pt;
          dis_map.at<double>(i, j) = std::abs( p.z() - pos[2] );
            //std::sqrt(CGAL::squared_distance(p, Point_3(pos[0], pos[1], pos[2])));
        }

        else
          LOG(INFO) << "no intersection";
      }
    }
    pixel_distance_to_planes.push_back(dis_map);
  }

  // set segment arrangement face data
  m_arr_cdata = std::vector<Arr_face_data>(
    m_arr.number_of_faces() - m_arr.number_of_unbounded_faces(),
    Arr_face_data(m_cplanes.size(), 0));
  // record each face_area black pixel number
  map< int, pair<int, int> > record_pixel; // <black, non-black>
  for (auto fitr = overlay_arr.faces_begin(); fitr != overlay_arr.faces_end(); ++fitr) {
    const int fidx = fitr->data().fidx;
    const cv::Vec2i pixel = fitr->data().pixel;
    if (fitr->is_unbounded() || fidx < 0)
      continue;
    if(record_pixel.count(fidx) == 0)
       record_pixel[fidx] = {0, 0};
    if (pixel[0] < 0 || pixel[1] < 0)
        record_pixel[fidx].first++;
    else
        m_ncmap.at<int>(pixel) == 0 ? record_pixel[fidx].first++ : record_pixel[fidx].second++;
   
  }
  for (auto fitr = overlay_arr.faces_begin(); fitr != overlay_arr.faces_end(); ++fitr) {
    const int fidx = fitr->data().fidx;
    const cv::Vec2i pixel = fitr->data().pixel;
    if (fitr->is_unbounded() || fidx < 0)
      continue;
    Arr_face_data &arrd = m_arr_cdata[fidx];

    CGAL::Polygon_2<K_epec> plg;
    auto curr = fitr->outer_ccb();

    do {
      plg.push_back(curr->target()->point());
      ++curr;
    } while (curr != fitr->outer_ccb());
    const double area = CGAL::to_double(plg.area());
    assert(area > 0.0);
    arrd.area += area;

    if (pixel[0] < 0 || pixel[1] < 0) {
      // faces out side sampling grid is assigned with ground-like data
      for (auto &d : arrd.distances)
        d = 1e20;
      arrd.distances.front() = 0.0;
      arrd.pnumber = 0;
      //arrd.normal = cv::Vec3d(0.0, 0.0, 1.0);
    }
    else {
      double black_ratio = 1.0*record_pixel[fidx].first/(record_pixel[fidx].first+record_pixel[fidx].second);
      double nblack_ratio = 1.0*record_pixel[fidx].second/(record_pixel[fidx].first+record_pixel[fidx].second);
      for (std::size_t i = 0; i < pixel_distance_to_planes.size(); ++i)
        arrd.distances[i] += area * (m_ncmap.at<int>(pixel) == 0 ? black_ratio *  (i == 0 ? pixel_distance_to_planes[i].at<double>(pixel) : pixel_distance_to_planes[pixel_distance_to_planes.size()-1].at<double>(pixel)) : nblack_ratio *  pixel_distance_to_planes[i].at<double>(pixel));
      arrd.pnumber += m_ncmap.at<int>(pixel); // unsure: 将空pixel到非black平面的距离设置为一样的值，对缺失更鲁棒?
      arrd.bnum = record_pixel[fidx].first;
      arrd.nbnum =record_pixel[fidx].second;
    }


  }
  for (auto &arrd : m_arr_cdata) {
    for (auto &d : arrd.distances)
      d /= arrd.area;
  }

  }
}

void building::facade_segment_recognization(){
    constexpr double EPSILON_GENERAL = 1.0e-4; 
    typedef CGAL::Cartesian_converter<K_epec, K> To_geom;
    To_geom to_geom;
    for (auto eitr = m_arr.edges_begin(); eitr != m_arr.edges_end(); ++eitr) {
        Arr_with_index::Vertex_handle vs = eitr->source();
        Arr_with_index::Vertex_handle vt = eitr->target();
        Point_2 s = to_geom(vs->point());
        Point_2 t = to_geom(vt->point());

        double s1 = s.x() < t.x() ? s.x():t.x();
        double s2 = s.x() > t.x() ? s.x():t.x();

        int flag = 0;
        double dis1 = 1, dis2 = 1;
        for(int i = 0; i < m_facade_segs.size(); i++){
           //fixed: avoid overlap            
            auto seg = to_geom(m_facade_segs[i]);
             dis1 = min(dis1, point_segment_squared_distance(s.x(), s.y(), seg.source().x(), seg.source().y(), seg.target().x(), seg.target().y()));
             dis2 = min(dis2, point_segment_squared_distance(t.x(), t.y(), seg.source().x(), seg.source().y(), seg.target().x(), seg.target().y()));
           if(point_segment_squared_distance(s.x(), s.y(), seg.source().x(), seg.source().y(), seg.target().x(), seg.target().y()) < EPSILON_GENERAL
            && point_segment_squared_distance(t.x(), t.y(), seg.source().x(), seg.source().y(), seg.target().x(), seg.target().y()) < EPSILON_GENERAL){
                // edge is on the facade
                if(s < t){
                m_facade.push_back(Facade_data(i, K_epec::Segment_2(vs->point(),vt->point())));
               
                }
                else{
                m_facade.push_back(Facade_data(i, K_epec::Segment_2(vt->point(),vs->point())));
     
                }
                flag = 1;
                facade_tags.push_back(1);
                break;
            }
        }
        if(flag == 0)
          facade_tags.push_back(0);

    }

    if(facade_tags.size() != m_arr.number_of_edges())
      LOG(INFO) << "Arrangement error.";

   

}

void building::segment_arrangement_MRF_labeling(std::vector<Arr_face_data>& m_arr_data, std::vector<Plane_3>& m_planes, Arr_with_index& m_arrs) {

  // two faces may be separated by different segments, we need mapping
  typedef std::pair<int, int> Adj_facets;
  typedef std::map<Adj_facets, double> Adj_map;

  const auto &config = cm::get_config();
  LOG(INFO) << "MRF params: " <<
    config.get<bool>("arr.mrf.use_swap") << ' ' <<
    config.get<double>("arr.mrf.balance") << ' ' <<
    config.get<int>("arr.mrf.iterations");
  if (m_planes.size() < 2) {
    LOG(INFO) << "Insufficient shape.";
    return;
  }

  // set up neighboring
  Adj_map adj_map;
  for (auto eitr = m_arrs.edges_begin(); eitr != m_arrs.edges_end(); ++eitr) {
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

  // set up arrangement graph
  const int nb_sites = int(m_arr_data.size());
  const int nb_labels = int(m_planes.size());

  // setup data term
  std::vector<double> data(nb_sites * nb_labels, 0.0);
  for (int i = 0; i < nb_sites; ++i)
    for (int j = 0; j < nb_labels; ++j) {
      data[i * nb_labels + j] =
        m_arr_data[i].distances[j] * m_arr_data[i].area;
        LOG(INFO) << "TAST: " << data[i * nb_labels + j];
    }


  // setup smooth
  std::vector<double> smooth(nb_labels * nb_labels, 1.0);
  for (int i = 0; i < nb_labels; ++i)
    smooth[i * nb_labels + i] = 0.0;

  // MRF labeling
  try {
    GCoptimizationGeneralGraph gc(nb_sites, nb_labels);
    gc.setDataCost(data.data());
    gc.setSmoothCost(smooth.data());

    // set neighboring
    // arrangement face area and pixel length are both in image coordinate
    // const double pixel_length = 0.2;
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

    // unsure: 浅拷贝 or 引用
    for (auto fitr = m_arrs.faces_begin(); fitr != m_arrs.faces_end(); ++fitr) {
      if (fitr->data() >= 0)
        fitr->set_data(gc.whatLabel(fitr->data()));
    }

  }
  catch (GCException e) {
    e.Report();
  }
} 

// fixed: label adjust
// unsure: merge facades of label0, and then adjust
void building::arrangement_adjust(){
    const auto &config = cm::get_config();
    constexpr double EPSILON_GENERAL = 1.0e-10;
    double lamda = config.get<double>("arr.mrf.balance");

    // bbox edges
    // fixed:bbox boundary face label->0
    vector<K_epec::Segment_2> boundary;
    boundary.push_back(K_epec::Segment_2(m_bbox.vertex(0),m_bbox.vertex(1))); 
    boundary.push_back(K_epec::Segment_2(m_bbox.vertex(1),m_bbox.vertex(2)));
    boundary.push_back(K_epec::Segment_2(m_bbox.vertex(2),m_bbox.vertex(3)));
    boundary.push_back(K_epec::Segment_2(m_bbox.vertex(3),m_bbox.vertex(0)));

    auto eitr2 = m_carr.edges_begin();
    for(auto eitr = m_farr.edges_begin(); eitr != m_farr.edges_end(); eitr++){
      if(facade_tags[std::distance(m_farr.edges_begin(), eitr)] == 0){
        int tag = 0;
        for(int i = 0; i < 4; i++){
          if(point_segment_squared_distance(CGAL::to_double(eitr->source()->point().x()), CGAL::to_double(eitr->source()->point().y()), CGAL::to_double(boundary[i].source().x()), CGAL::to_double(boundary[i].source().y()), CGAL::to_double(boundary[i].target().x()), CGAL::to_double(boundary[i].target().y())) < EPSILON_GENERAL
            && point_segment_squared_distance(CGAL::to_double(eitr->target()->point().x()), CGAL::to_double(eitr->target()->point().y()), CGAL::to_double(boundary[i].source().x()), CGAL::to_double(boundary[i].source().y()), CGAL::to_double(boundary[i].target().x()), CGAL::to_double(boundary[i].target().y())) < EPSILON_GENERAL){
              tag = 1;
              break;
            }
        }
        if(tag == 1) {
           eitr->face()->data() = 0;
           eitr->twin()->face()->data() = 0;
           eitr2->face()->data() = 0;
           eitr2->twin()->face()->data() = 0;
        }// boundary non-facade face         
      } // non-facade
      eitr2++;
    }
    int time = 0;
    while(time < 10){ // todo: merge same label
    time++;
    int tag = 1;
    auto arr2 = m_carr.faces_begin();
    for(auto arr=m_farr.faces_begin(); arr!= m_farr.faces_end(); arr++){
          if(arr->is_unbounded()){
            arr2++;
            continue;
          }
          if(arr->data() == 0 && arr2->data() != 0 || arr2->data() == 0 && arr->data() != 0){
            if(arr->data() == 0){
              vector<double> len(m_fplanes.size()-1,0);
              auto curr = arr->outer_ccb();
              double len_all = 0;
              do {
                if(curr->twin()->face()->data()){
                  const double elen = std::sqrt(CGAL::to_double(
                  CGAL::squared_distance(curr->target()->point(),curr->source()->point())));
                  len_all += elen;                 
                  len[curr->twin()->face()->data()-1] += elen;
                }
                ++curr;
              } while (curr != arr->outer_ccb());
               if(len_all > 1e-4){ // do not need merge
                  // only consider nerghbor labels
                  for(int i = 0; i < len.size() ; i++){
                    len[i] = lamda*(len_all - len[i]);// + m_arr_fdata[std::distance(m_farr.faces_begin(),arr)].distances[i+1];
                  }
                  int mini = 0; double minl = len[0];
                  for(int i = 1; i < len.size(); i++){  
                  if(minl > len[i]) {minl = len[i]; mini = i;}
                  }
                  arr->set_data(mini+1);
               }
               else
               {
                 tag = 0;
               }
               
             
            }
            if(arr2->data() == 0){
              vector<int> len(m_cplanes.size()-1,0);
              auto curr = arr2->outer_ccb();
              double len_all = 0;
              do {
                if(curr->twin()->face()->data()){
                    const double elen = std::sqrt(CGAL::to_double(
                    CGAL::squared_distance(curr->target()->point(),curr->source()->point())));
                    len_all += elen;
                    len[curr->twin()->face()->data()-1] += elen;
                }

                ++curr;
              } while (curr != arr2->outer_ccb());
              if(len_all > 1e-4){ // DO NOT need merge
                 for(int i = 0; i < len.size(); i++){
                    len[i] = lamda*(len_all - len[i]) + m_arr_cdata[std::distance(m_carr.faces_begin(),arr2)].distances[i+1];
                  }
                  int mini = 0; double minl = len[0];
                  for(int i = 1; i < len.size(); i++)
                  if(minl > len[i]) {minl = len[i]; mini = i;}
                  arr2->set_data(mini+1);
              }
              else
              {
                tag = 0;
              }
            }
          }
          arr2++;
    }
    if(tag == 1)
     break;
    }

    return;

}

// todo: the whole MRF
void building::segment_arrangement_MRF_labeling_cf(std::vector<Arr_face_data>& m_arr_cdata, std::vector<Arr_face_data>& m_arr_fdata, std::vector<Plane_3>& m_fplanes, std::vector<Plane_3>& m_cplanes, Arr_with_index& m_carr, Arr_with_index& m_farr) {
  // two faces may be separated by different segments, we need mapping
  if(m_arr_cdata.size() != m_arr_fdata.size()){
    LOG(INFO) << "Arrangment error.";
    return;
  }
  
  typedef std::pair<int, int> Adj_facets;
  typedef std::map<Adj_facets, double> Adj_map;

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
  // label: black, cifj
  const int nb_sites = int(m_arr_cdata.size());
  const int nb_labels = int((m_fplanes.size()-1) * (m_cplanes.size()-1) + 1);
  LOG(INFO) << nb_sites << " graph vertices. " << nb_labels << " labels.";

  // set up neighboring
  Adj_map adj_map;
  for (auto eitr = m_carr.edges_begin(); eitr != m_carr.edges_end(); ++eitr) {
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


  // setup data term
  std::vector<double> data(nb_sites * nb_labels, 0.0);
  for (int i = 0; i < nb_sites; ++i){
      int cnum = m_arr_cdata[i].nbnum;
      int fnum = m_arr_fdata[i].nbnum; 
      double cratio = cnum==0 && fnum==0 ? 1.0: 1.0*cnum/(cnum+fnum);    
      double fratio = cnum==0 && fnum==0 ? 1.0: 1.0*fnum/(cnum+fnum);  
    for (int j = 0; j < nb_labels; ++j) { // maybe different data weight
      if(j == 0){ // black
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
    /*
  for (int i = 0; i < nb_sites; ++i){
    for (int j = 0; j < nb_labels; ++j) {
        LOG(INFO) << data[i * nb_labels + j] << " ";
    }
      LOG(INFO) << "\n";
  }
   */
    
  LOG(INFO) << "Data term down.";
  

  // setup smooth
  class SmoothFn : public GCoptimization::SmoothCostFunctor {
  public:
    SmoothFn(const size_t ceil_size_, const size_t cplane_size_) : ceil_size(ceil_size_) ,  cplane_size(cplane_size_) {}

    double compute(int s1, int s2, int l1, int l2) {
      if(l1==l2)
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
  // setup smooth
  std::vector<double> smooth(nb_labels * nb_labels, 1.0);
  for (int i = 0; i < nb_labels; ++i)
    smooth[i * nb_labels + i] = 0.0;

    GCoptimizationGeneralGraph gc(nb_sites, nb_labels);
    gc.setDataCost(data.data());
    SmoothFn smoothfn(m_arr_cdata.size(), m_cplanes.size());
    gc.setSmoothCostFunctor(&smoothfn);
    //gc.setSmoothCost(smooth.data());
    LOG(INFO) << "Smooth term down.";

    // set neighboring
    // arrangement face area and pixel length are both in image coordinate
    // const double pixel_length = 0.2;
    for (const auto &adj_pair : adj_map)
      gc.setNeighbors(
        adj_pair.first.first,
        adj_pair.first.second,
        adj_pair.second * config.get<double>("arr.mrf.balance"));

    LOG(INFO) << "Before optimization energy is: " << gc.compute_energy();
    LOG(INFO) << gc.giveDataEnergy() << " " << gc.giveSmoothEnergy() << " " << gc.giveLabelEnergy();
    if (config.get<bool>("arr.mrf.use_swap")) {
      LOG(INFO) << "Alpha-beta swap algorithm.";
      gc.swap(config.get<int>("arr.mrf.iterations"));
    }
    else {
      LOG(INFO) << "Alpha expansion algorithm.";
      gc.expansion(config.get<int>("arr.mrf.iterations"));
    }
    LOG(INFO) << "After optimization energy is: " << gc.compute_energy();


    // 引用
    for (auto fitr = m_carr.faces_begin(); fitr != m_carr.faces_end(); ++fitr) {
      if (fitr->data() >= 0){
        if(fitr->data() > nb_sites)
          LOG(INFO) << "Site error: " << fitr->data() << " " <<  nb_sites;
        if(gc.whatLabel(fitr->data()) == 0)
           fitr->set_data(0);
        else 
          fitr->set_data(gc.whatLabel(fitr->data())%(m_cplanes.size()-1) == 0 ? (m_cplanes.size()-1) : gc.whatLabel(fitr->data())%(m_cplanes.size()-1));
        }

    }
    for (auto fitr = m_farr.faces_begin(); fitr != m_farr.faces_end(); ++fitr) {
      if (fitr->data() >= 0){
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


// fixed: 立面差的高度补全 记录facade高度 增加数据结构 Facade_data
bool building::segment_arrangement_extrusion(std::vector<Plane_3>& m_planes, int tag, Arr_with_index& m_arrs) {

  // lod2 floor/ceiling extrusion, with color

    typedef CGAL::Cartesian_converter<K, K_epec> To_geom;
    To_geom to_geom;
    typedef CGAL::Cartesian_converter<K_epec, K> To_geom2;
    To_geom2 to_geom2;
    std::vector<K_epec::Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;
    Mesh mesh;

    for (auto fitr = m_arrs.faces_begin(); fitr != m_arrs.faces_end(); ++fitr) {
        const int label = fitr->data();
        if (label < 0 || label > m_planes.size()) // include faces label = 0
          continue;
        // fixed: 退化为线段的faces处理
        /*
        CGAL::Polygon_2<K_epec> plg;
        auto currs = fitr->outer_ccb();
        do { 
          plg.push_back(currs->target()->point());
          ++currs;
        } while (currs != fitr->outer_ccb());
        double area = CGAL::to_double(plg.area());
        if(area < 5e-3){  // unsure: filter out small faces, pick the longest edges
          currs = fitr->outer_ccb();

          auto currs_longest = currs;
          double len_max = 0;
          do { 
          double len = sqrt(CGAL::to_double((currs->source()->point()-currs->target()->point()).squared_length()));
          if(len > 5e-3){
            if(len_max < len){
                len_max = len;
                currs_longest = currs;
            }
          }
          ++currs;
        } while (currs != fitr->outer_ccb());

          CGAL::Polygon_2<K_epec> plg3; 
          if(currs_longest->twin()->face()->number_of_outer_ccbs() == 1){  // case: torus has two outer_ccbs, precondition: number_of_outer_ccbs() == 1
          auto currs3 = currs_longest->twin()->face()->outer_ccb();
          do { 
            plg3.push_back(currs3->target()->point());
            ++currs3;
            } while (currs3 != currs->twin()->face()->outer_ccb());
            double area3 = CGAL::to_double(plg3.area());
          if(area3 > 5e-3){
            fitr->set_data(currs_longest->twin()->face()->data());// degenerated faces
            break;}
          }
          continue;
        }
        */
        CGAL::Polygon_2<K_epec> plg;
        auto currs = fitr->outer_ccb();
        do { 
          plg.push_back(currs->target()->point());
          ++currs;
        } while (currs != fitr->outer_ccb());
        double area = CGAL::to_double(plg.area());
        if(area < 5e-3){  // unsure: filter out small faces, pick the longest edges
          currs = fitr->outer_ccb();
          double max_len = 5e-3;
          do { 
          double len = sqrt(CGAL::to_double((currs->source()->point()-currs->target()->point()).squared_length()));
          if(len > max_len && currs->twin()->face()->number_of_outer_ccbs() == 1){
              fitr->set_data(currs->twin()->face()->data());// degenerated faces
              max_len = len;
          }
          ++currs;
        } while (currs != fitr->outer_ccb());
        }
        

        if(label == 0)
          continue;
      
        // extrude ceiling/floor planes
        // points
        const std::size_t offset = points.size();
        std::size_t nb_vertices = 0;
        auto curr = fitr->outer_ccb();
        // 多边形划分
        vector<Point_2> pp;
        do{
          pp.push_back(to_geom2(curr->target()->point()));
          ++curr;
        }while (curr != fitr->outer_ccb());
        auto it = pp.begin();
        auto prev(it++);
        // delete the same neighbor points
        for(; it != pp.end();){
          Point_2 p1 = Point_2(static_cast<int>(it->x()*1000)/1000.0, (static_cast<int>((*it).y()*1000))/1000.0);
          Point_2 p2 = Point_2((static_cast<int>((*prev).x()*1000))/1000.0, (static_cast<int>((*prev).y()*1000))/1000.0);
          if(p1 == p2){
            prev = pp.erase(prev);
            it = prev+1;
          }
          else{
            ++it,++prev;
          }
        }
        Point_2 p1 = Point_2((static_cast<int>(pp.begin()->x()*1000))/1000.0, (static_cast<int>(pp.begin()->y()*1000))/1000.0);
        Point_2 p2 = Point_2((static_cast<int>(prev->x()*1000))/1000.0, (static_cast<int>(prev->y()*1000))/1000.0);
        if(p1 == p2)
          prev = pp.erase(prev);
        if(pp.size() < 3)
          continue;
        // guarantee1
        if(!CGAL::is_simple_2(pp.begin(), pp.end())){
          LOG(INFO) << "Non_simple polygon.";
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
        for (Polygon_list::iterator it = polys.begin(); it != polys.end();++it)
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

    // .ply
    string fname;
    if(tag == 0)
     fname = m_dir + "floor_vec.ply";
    else
     fname = m_dir + "ceiling_vec.ply";

    //OFFMeshSave(fname, mesh);
    const auto &config = cm::get_config();
    if(!PLYTriMeshSave(mesh, config, fname))
      LOG(INFO) << "Ceiling/Floor saves failed.";
    //PLYMeshSave(fname, mesh, 1);

    // fixed:solve inconsistant height, store in <facades> sturcture
    //some non-facade edges cause height inconsistance
    int idx = -1; // index in m_facade
    unordered_map<string, int> re_edges; // todo: record close edges and index
     for (auto eitr = m_arrs.edges_begin(); eitr != m_arrs.edges_end(); ++eitr) {
        int l1 = eitr->face()->data();
        int l2 = eitr->twin()->face()->data();
        K_epec::Point_3 ps(eitr->source()->point().x(),
        eitr->source()->point().y(),
        0.0);
        K_epec::Point_3 pt(eitr->target()->point().x(),
        eitr->target()->point().y(),
        0.0);


         // fixed: non-facade height inconsistance
         if(facade_tags[std::distance(m_arrs.edges_begin(), eitr)] == 0){ // non-facade
            // unsure:test
            if(l1 <=0 || l1 > m_planes.size() || l2 <=0 || l2 > m_planes.size())
            continue;
             if(l1 != l2){
                K_epec::Point_3 *p1s = nullptr, *p1t = nullptr, *p2s = nullptr, *p2t = nullptr;                             
                Facade f;
                f.seg = K_epec::Segment_2(eitr->source()->point(), eitr->target()->point());
                auto result1 = CGAL::intersection(
                to_geom(m_planes[l1]),
                K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result1 && (p1s = boost::get<K_epec::Point_3>(&(*result1)))){
                    f.heights.push_back(CGAL::to_double(p1s->z()));
                }
                else{
                    LOG(INFO) << "Extrusion error2.";
                    return false;    
                }

                auto result2 = CGAL::intersection(
                to_geom(m_planes[l2]),
                K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result2 && (p2s = boost::get<K_epec::Point_3>(&(*result2)))){
                    f.heights.push_back(CGAL::to_double(p2s->z()));                
                }
                else{
                    LOG(INFO) << "Extrusion error2.";
                    return false;    
                }
                
                auto result3 = CGAL::intersection(
                to_geom(m_planes[l1]),
                K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result3 && (p1t = boost::get<K_epec::Point_3>(&(*result3)))){
                    f.heights.push_back(CGAL::to_double(p1t->z()));  
                }
                else{
                    LOG(INFO) << "Extrusion error2.";
                    return false;    
                }
                
                auto result4 = CGAL::intersection(
                to_geom(m_planes[l2]),
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
            
               // fixed: decide the normal
                if(tag == 0) // floor
                  n = Plane_3(to_geom2(*p2t), to_geom2(*p1t), to_geom2(*p1s));
                else // ceiling
                  n = Plane_3(to_geom2(*p2t), to_geom2(*p2s), to_geom2(*p1s));  

                f.normal = Vector_2(n.orthogonal_vector().x(), n.orthogonal_vector().y());
                facades.push_back(f); // arr.cpp
             }
         }
         else{
             idx++;
             // facade 
             // fixed: push facade height value
            if((l1 <=0 || l1 > m_planes.size()) && (l2 <=0 || l2 > m_planes.size()))
            continue;
            K_epec::Point_3 *p1s = nullptr, *p1t = nullptr, *p2s = nullptr, *p2t = nullptr;
              
            if(!(l1 <=0 || l1 > m_planes.size())){
                auto result = CGAL::intersection(
                    to_geom(m_planes[l1]),
                    K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result && (p1s = boost::get<K_epec::Point_3>(&(*result)))){
                    auto it = find(points.begin(), points.end(),*p1s);
                    if(it==points.end()){
                        //LOG(INFO) << "Facade height assign error.";
                        //return false;
                        points.push_back(*p1s);
                    }                       
                    }
                auto result2 = CGAL::intersection(
                    to_geom(m_planes[l1]),
                    K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));                
                if (result2 && (p1t = boost::get<K_epec::Point_3>(&(*result2)))){
                    auto it = find(points.begin(), points.end(),*p1t);
                    if(it==points.end()){
                        //LOG(INFO) << "Facade height assign error.";
                        //return false;
                        points.push_back(*p1t);
                    }                       
                    }
                Point_2 z = *p1s < *p1t ? Point_2(CGAL::to_double((*p1s).z()),CGAL::to_double((*p1t).z())) : Point_2(CGAL::to_double((*p1t).z()),CGAL::to_double((*p1s).z()));
                if(tag == 0) {// floor
                   m_facade[idx].bottoms.push_back(z);
                }
                else // ceiling
                m_facade[idx].ups.push_back(z);

            }
            
            if(!(l2 <=0 || l2 > m_planes.size())){
                 auto result = CGAL::intersection(
                    to_geom(m_planes[l2]),
                    K_epec::Line_3(ps, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result && (p2s = boost::get<K_epec::Point_3>(&(*result)))){
                    auto it = find(points.begin(), points.end(),*p2s);
                    if(it==points.end()){
                        //LOG(INFO) << "Facade height assign error.";
                        //return false;
                        points.push_back(*p2s);
                    }                       
                    }
                auto result2 = CGAL::intersection(
                    to_geom(m_planes[l2]),
                    K_epec::Line_3(pt, K_epec::Vector_3(0.0, 0.0, 1.0)));               
                if (result2 && (p2t = boost::get<K_epec::Point_3>(&(*result2)))){
                    auto it = find(points.begin(), points.end(),*p2t);
                    if(it==points.end()){
                        //LOG(INFO) << "Facade height assign error.";
                        //return false;
                        points.push_back(*p2t);
                    }                       
                    }
                Point_2 z = *p2s < *p2t ? Point_2(CGAL::to_double((*p2s).z()),CGAL::to_double((*p2t).z())) : Point_2(CGAL::to_double((*p2t).z()),CGAL::to_double((*p2s).z()));
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
                  m_facade[idx].tag =1;
            }

         }
     }

    return true;

}


// fixed: version2: extrude facade -> m_facade
bool building::facade_extrusion(){

    typedef CGAL::Cartesian_converter<K_epec, K> To_geom;
    To_geom to_geom;
    typedef CGAL::Cartesian_converter<K, K_epec> To_geom2;
    To_geom2 to_geom2;

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

    //fixed: m_facade sort, then delete repeated edges
   sort(m_facade.begin(), m_facade.end()); // todo:facade sort
   vector<int> nonrepeated(m_facade.size(), 0);
   set<string> re_facades;
   for(int i = 0; i < m_facade.size(); ){
      if(CGAL::to_double((m_facade[i].seg.source()-m_facade[i].seg.target()).squared_length()) < 1e-6){
         i++;
         continue; 
      }// filter out too small segments
      nonrepeated[i] = 1;
      K_epec::Segment_2 seg = m_facade[i].seg;
      Point_2 s1 = to_geom(seg.min());
      Point_2 t1 = to_geom(seg.max());
      std::string str1 = std::to_string((int)std::round(s1.x()*1000)) + std::to_string((int)std::round(s1.y()*1000)) + std::to_string((int)std::round(t1.x()*1000)) + std::to_string((int)std::round(t1.y()*1000));
      re_facades.insert(str1);
      int j = i+1;
      for(; j < m_facade.size(); j++){        
        Point_2 s2 = to_geom(m_facade[j].seg.min());
        Point_2 t2 = to_geom(m_facade[j].seg.max());// compute error
        std::string str2 = std::to_string((int)std::round(s2.x()*1000)) + std::to_string((int)std::round(s2.y()*1000)) + std::to_string((int)std::round(t2.x()*1000)) + std::to_string((int)std::round(t2.y()*1000));
        if(re_facades.count(str2)){ // repeated facade exist
        //test
        LOG(INFO) << s1 << " " << t1 << ", " << s2 << " " << t2;
        if(m_facade[j].tag == 0)
          m_facade[i].tag = 0; //  tag update
        for(auto t:m_facade[j].ups)
          if(m_facade[i].ups.size()==0 || find(m_facade[i].ups.begin(),m_facade[i].ups.end(),t)==m_facade[i].ups.end())
            m_facade[i].ups.push_back(t);
        for(auto t:m_facade[j].bottoms)
          if(m_facade[i].bottoms.size()==0 || find(m_facade[i].bottoms.begin(),m_facade[i].bottoms.end(),t)==m_facade[i].bottoms.end())
            m_facade[i].bottoms.push_back(t);
         }
         else
         {
           i=j;break;
         }
         
       }
       i=j;
   }

    for(int i = 0; i < m_facade.size(); i++){
        auto facade = m_facade[i];
        Point_2 s(CGAL::to_double(facade.seg.source().x()), CGAL::to_double(facade.seg.source().y()));
        Point_2 t(CGAL::to_double(facade.seg.target().x()), CGAL::to_double(facade.seg.target().y()));
        double s1 = s.x() < t.x()? s.x() : t.x();
        double s2 = s.x() > t.x()? s.x() : t.x();
       
        if(facade.bottoms.size() == 0 || facade.ups.size() == 0 || nonrepeated[i]==0) // unsure: || or &&
          continue;
        facade.normal = normals[facade.index];
        vector<Point_2> height;
        for(auto p:facade.ups)
           height.push_back(p);
        for(auto p:facade.bottoms)
           height.push_back(p);

        vector<Point_2> margin(height.size()-1);
        vector<int> margin_num(height.size()-1,0);
        vector<set<Point_2>> counts(height.size()-1);
        vector<int> rows_grid(height.size()-1);
        // fixed: bbox grid for points 

        std::sort(height.begin(), height.end(), [](Point_2 p1, Point_2 p2){ return (p1.x()+p1.y())/2 > (p2.x()+p2.y())/2;});// sort by avrage z
        for(int i = 0; i < height.size()-1; i++)
           margin[i] = Point_2((height[i].x()+height[i].y())/2, (height[i+1].x()+height[i+1].y())/2); // {big,small}
                         
        K_epec::Point_2 lc = facade.seg.source();
        K_epec::Point_2 rc = facade.seg.target();
        int cols = std::ceil(sqrt(CGAL::to_double((lc-rc).squared_length()))/m_step);
        for(int i = 0; i < margin_num.size(); i++){
               int rows = std::ceil((margin[i].x() - margin[i].y())/m_step);
               rows_grid[i] = rows;
        }
          
        Line_2 l1 = to_geom(facade.seg.supporting_line().perpendicular(facade.seg.source()));
        Line_2 l2 =  to_geom(facade.seg.supporting_line().perpendicular(facade.seg.target()));
        for(auto p : points[facade.index]){
            if(l1.has_on_positive_side(Point_2(p.x(), p.y())) != l2.has_on_positive_side(Point_2(p.x(), p.y()))){
                for(int i = 0; i < margin.size(); i++){
                    if(p.z() <= margin[i].x() && p.z() >= margin[i].y()){
                        margin_num[i]++;
                        int c = sqrt(CGAL::to_double(((facade.seg.supporting_line().projection(K_epec::Point_2(p.x(), p.y()))) - lc).squared_length()))/m_step;
                        int r = (margin[i].x() - p.z())/m_step;
                        if(c > 0 && c < cols && r > 0 && r < rows_grid[i]){
                          counts[i].insert({c, r});
                        }
                        break;
                    }
                }
            }
        }


        Point_3 p1, p2, p3, p4;
        bool c2 = margin.size() == 1 && facade.tag == 0; // unsure:处理连接面只有一个的边的情况,立面不能确定，非立面强制确定
        bool c3 = facade.ups.size() == 1 && facade.bottoms.size() == 2;
        bool c4 = facade.ups.size() == 2 && facade.bottoms.size() == 1;
        bool c5 = facade.ups.size() == 2 && facade.bottoms.size() == 2;
        // unsure: guarantee manifold models
        bool c6 = facade.ups.size() == 1;
        for(int i = 0; i < margin_num.size(); i++){  
           if( c2 || (c3 && i == 1) || (c4 && i == 0) || (c5 && (i == 0 || i == 2))){ // push these height differences into the wall solution           
               facades.push_back(Facade(K_epec::Segment_2(facade.seg.source(), facade.seg.target()), {height[i].x(), height[i+1].x(),height[i].y(), height[i+1].y()}));
               facades[facades.size()-1].normal = facade.normal;            
            }
            else{ // unsure height differences
              //if(facade.seg.squared_length() < 1e-8) // unsure:filter out single points (really small segment)
                 //continue;
              double sx = to_geom(facade.seg.source().x()); double sy = to_geom(facade.seg.source().y());
              double tx = to_geom(facade.seg.target().x()); double ty = to_geom(facade.seg.target().y());
              string temp1 = std::to_string(floor(sx*100000)/100000) + std::to_string(floor(sy*100000)/100000) + std::to_string(floor(tx*100000)/100000) + std::to_string(floor(ty*100000)/100000);
              string temp2 = std::to_string(floor(tx*100000)/100000) + std::to_string(floor(ty*100000)/100000) + std::to_string(floor(sx*100000)/100000) + std::to_string(floor(sy*100000)/100000);
             
              if(non_overlap.count(temp1)==0 && non_overlap.count(temp2)==0){
                 potential_facades.push_back(Facade(K_epec::Segment_2(facade.seg.source(), facade.seg.target()), {height[i].x(), height[i+1].x(),height[i].y(), height[i+1].y()}));
                 potential_facades[potential_facades.size()-1].normal = facade.normal;
                 non_overlap.insert(temp1);
                 non_overlap.insert(temp2);
              
               suppoints_num.push_back(margin_num[i]);
               area_ratio.push_back(1-(counts[i].size()*1.0)/(rows_grid[i] * cols));
               area_total += area_ratio[area_ratio.size()-1];
               //LOG(INFO) << area_ratio[area_ratio.size()-1] << " " << counts[i].size() << " " << rows_grid[i] * cols;
              }
            }
              
        }
        vector<Point_2>().swap(margin); 
        vector<int>().swap(margin_num);
        vector<set<Point_2>>().swap(counts);
        vector<int>().swap(rows_grid);
      }
      vector<vector<Point_3>>().swap(points);
    
    for(auto &a:area_ratio)
       a /= area_total;

    // fixed: 当立面的支撑点空时，若该立面为连接顶地面的唯一桥梁则强制恢复该立面->假设。

    // fixed: non consistant non-facade segments are also added into <facades> structure
    auto e2 = m_carr.edges_begin();
    for(auto e1 = m_farr.edges_begin(); e1 != m_farr.edges_end(); e1++,e2++){
      const int l11 = e1->face()->data();
      const int l12 = e1->twin()->face()->data();
      const int l21 = e2->face()->data();
      const int l22 = e2->twin()->face()->data();
      if(facade_tags[std::distance(m_farr.edges_begin(),e1)] == 1)  continue;// for non facade
      // fixed: boundary edges -> face index = -1
      if(l11 <=0  && l12 >0 && l12 < m_fplanes.size() || l12 <= 0 && l11 >0 && l11 < m_fplanes.size())
      {
          Facade f;
          // judge normal
          if(l11 <= 0 && l12 >0 && l12 < m_fplanes.size()){
            K_epec::Vector_2 v = (K_epec::Line_2(e1->twin()->source()->point(),e1->twin()->target()->point())).perpendicular(K_epec::Point_2(0,0)).to_vector();
            f.normal = Vector_2(CGAL::to_double(v.x()),CGAL::to_double(v.y()));
          }
          else{
            K_epec::Vector_2 v = (K_epec::Line_2(e1->source()->point(),e1->target()->point())).perpendicular(K_epec::Point_2(0,0)).to_vector();
            f.normal = Vector_2(CGAL::to_double(v.x()),CGAL::to_double(v.y()));
          }
          f.seg = K_epec::Segment_2(e1->source()->point(), e1->target()->point());
           K_epec::Point_3 ps(e1->source()->point().x(),
          e1->source()->point().y(),
            0.0);
          K_epec::Point_3 pt(e1->target()->point().x(),
          e1->target()->point().y(),
          0.0);
          K_epec::Plane_3 p1,p2;
          if(l11 > 0 )
             p1 = to_geom2(m_fplanes[l11]);
          else
             p1 = to_geom2(m_fplanes[l12]);
          if(l21 > 0 )
             p2 = to_geom2(m_cplanes[l21]);
          else
             p2 = to_geom2(m_cplanes[l22]);          
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

    /*
    facades: must be in the solution, guarantee the manifold models
    potential_facades: maybe in the solution, solved by lp
    */
    vector<IntersectPoint> inter_points;  
    vector<pair<K_epec::Point_2, K_epec::Point_2>> sub_points;
    vector<Point_2> points_height; // {high, low}
    for(int i = 0; i < potential_facades.size(); i++){
      auto f = potential_facades[i];
      auto it = find(sub_points.begin(), sub_points.end(),make_pair(f.seg.source(), K_epec::Point_2(f.heights[0], f.heights[1])));
      if(it == sub_points.end()){
        sub_points.push_back(make_pair(f.seg.source(), K_epec::Point_2(f.heights[0], f.heights[1]))); // s s-high s-low
        inter_points.push_back(IntersectPoint(to_geom(f.seg.source()), Point_2(f.heights[0], f.heights[1]), {i}));
      }
      else{
        inter_points[std::distance(sub_points.begin(), it)].IDs.push_back(i);
      }
      it = find(sub_points.begin(), sub_points.end(),make_pair(f.seg.target(), K_epec::Point_2(f.heights[2], f.heights[3])));
      if(it == sub_points.end()){
        sub_points.push_back(make_pair(f.seg.target(), K_epec::Point_2(f.heights[2], f.heights[3]))); // t t-high t-low
        inter_points.push_back(IntersectPoint(to_geom(f.seg.target()), Point_2(f.heights[2], f.heights[3]), {i}));
      }
      else{
        inter_points[std::distance(sub_points.begin(), it)].IDs.push_back(i);
      }
    }
    vector<pair<K_epec::Point_2, K_epec::Point_2>>().swap(sub_points);

    // judge inter_point.flag
    for(auto &ip:inter_points)
      for(auto &s:facades){
        if(point_segment_squared_distance(ip.p.x(), ip.p.y(), CGAL::to_double(s.seg.source().x()), CGAL::to_double(s.seg.source().y()), CGAL::to_double(s.seg.target().x()), CGAL::to_double(s.seg.target().y())) < 1e-8){
          
           K_epec::Point_2 p = s.seg.supporting_line().projection(to_geom2(ip.p));
           int idx_max = s.heights[0] > s.heights[2] ? 0 : 1;
           int idx_min = s.heights[0] > s.heights[2] ? 1 : 0;
           double len_max = sqrt(CGAL::to_double(s.seg.squared_length()));
           double len_min = sqrt(CGAL::to_double((s.seg[idx_min] - p).squared_length())); 
           double hmax = s.heights[2*idx_min] + (s.heights[2*idx_max] - s.heights[2*idx_min]) * len_min / len_max;
           idx_max = s.heights[1] > s.heights[3] ? 0 : 1;
           idx_min = s.heights[1] > s.heights[3] ? 1 : 0;
           double hmin = s.heights[2*idx_min+1] + (s.heights[2*idx_min+1] - s.heights[2*idx_min+1]) * len_min / len_max;
           if(abs(ip.heights.x() - hmax) < 1e-3 && abs(ip.heights.y() - hmin) < 1e-3){
             ip.flag = 1; // fixed: boundary point
             break;
           }

        }
      
      }
 
    // deal with close segments and points
    vector<IntersectPoint> re_inter_points;
    vector<int> re_points(inter_points.size(), 0);
    // todo: flag isolate small segments which is in just one segment
    int k = -1;
    for(int i = 0; i < inter_points.size(); i++){
         if(re_points[i] == 1)
             continue;
         re_inter_points.push_back(inter_points[i]);
         k++;
         int j = i+1;
         for(; j < inter_points.size(); j++){
            if(re_points[j] == 1)
             continue;
            auto p1 = inter_points[i].p;
            auto p2 = inter_points[j].p;
            auto h1 = inter_points[i].heights;
            auto h2 = inter_points[j].heights;
            //TODO: dela with small segments
            if((p1 - p2).squared_length() < 25e-6){ // && (h1 - h2).squared_length() < 25e-6){ //5e-3 //unsure: ignore heights
                re_points[j] = 1;
                if(inter_points[j].flag == 1) // boundary point
                  re_inter_points[k].flag = 1;
                re_inter_points[k].IDs.insert(re_inter_points[k].IDs.end(), inter_points[j].IDs.begin(), inter_points[j].IDs.end()); 
                sort(re_inter_points[k].IDs.begin(), re_inter_points[k].IDs.end());
                re_inter_points[k].IDs.erase(unique(re_inter_points[k].IDs.begin(), re_inter_points[k].IDs.end()), re_inter_points[k].IDs.end());
            }
         }
    }

    // unsure: potential facades delete small segments
    vector<int> isolate_seg(potential_facades.size(), 0);
    for(auto r: re_inter_points){
         for(auto i: r.IDs){
           isolate_seg[i]++;
         }
    }
    // unsure: potential facades merge close segments, consider heights
    vector<int> re_segs(potential_facades.size(), 0);    
    for(int i = 0; i < potential_facades.size();i++){
        if(re_segs[i] == 1 || isolate_seg[i] == 1) continue; // consider
        auto p = potential_facades[i];
        auto s = p.seg.source() > p.seg.target() ? p.seg.target() : p.seg.source();
        auto t = p.seg.source() < p.seg.target() ? p.seg.target() : p.seg.source();
        double hs = potential_facades[i].heights[0] - potential_facades[i].heights[1];
        double ht = potential_facades[i].heights[2] - potential_facades[i].heights[3];
        for(int j = i+1; j < potential_facades.size(); j++){
           if(re_segs[j] == 1) continue;
           auto p2 = potential_facades[j];
           auto s2 = p2.seg.source() > p2.seg.target() ? p2.seg.target() : p2.seg.source();
           auto t2 = p2.seg.source() < p2.seg.target() ? p2.seg.target() : p2.seg.source();
           double hs2 = potential_facades[j].heights[0] - potential_facades[j].heights[1];
           double ht2 = potential_facades[j].heights[2] - potential_facades[j].heights[3];
           if(CGAL::to_double((s-s2).squared_length()) < 25e-6 && CGAL::to_double((t-t2).squared_length()) < 25e-6 && abs(hs - hs2) < 5e-3 && abs(ht - ht2) < 5e-3){
              isolate_seg[j] = -1;
              re_segs[j] = 1;
           }
        }
    }
    
    // record inter_points
    ofstream fb(m_dir + "boundary_points");  
    ofstream fnb(m_dir + "non-boundary_points");  
    int boundarynums=0;
    for(auto ip:re_inter_points){
      if(ip.flag == 1){
        boundarynums++;
        fnb << ip.p << "\n";
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
    optimize2(potential_facades, re_inter_points, isolate_seg, suppoints_num, area_ratio, lambda_data_fitting, lambda_model_coverage, lambda_model_complexity, m_dir);
    
    // record results
    ofstream f1(m_dir + "potential_facades");
    ofstream f2(m_dir + "selected_facades");  
    for(auto f:potential_facades){
      f1 << CGAL::to_double(f.seg.source().x()) << " " << CGAL::to_double(f.seg.target().x()) << " " << CGAL::to_double(f.seg.source().y()) << " " << CGAL::to_double(f.seg.target().y()) << "\n";
      if(f.flag == 1)
      f2 << CGAL::to_double(f.seg.source().x()) << " " << CGAL::to_double(f.seg.target().x()) << " " << CGAL::to_double(f.seg.source().y()) << " " << CGAL::to_double(f.seg.target().y()) << "\n";
    }
    f1.close();
    f2.close(); 

    ofstream f3(m_dir + "pre_facades");
    for(auto f:facades)
      f3 << CGAL::to_double(f.seg.source().x()) << " " <<  CGAL::to_double(f.seg.target().x()) << " " << CGAL::to_double(f.seg.source().y()) << " " << CGAL::to_double(f.seg.target().y()) << "\n";
    f3.close();  

    // extrude facade off
    std::vector<Point_3> fpoints;
    std::vector<std::vector<std::size_t>> polygons;
    std::size_t id[4];
    Point_3 p[4];
    for(auto f:facades){
      p[0] = Point_3(CGAL::to_double(f.seg.source().x()), CGAL::to_double(f.seg.source().y()), f.heights[0]);  //sh
      p[1] = Point_3(CGAL::to_double(f.seg.source().x()), CGAL::to_double(f.seg.source().y()), f.heights[1]);  // sl
      p[2] = Point_3(CGAL::to_double(f.seg.target().x()), CGAL::to_double(f.seg.target().y()), f.heights[2]);  // th 
      p[3] = Point_3(CGAL::to_double(f.seg.target().x()), CGAL::to_double(f.seg.target().y()), f.heights[3]);  // tl
      //TODO: filter too small facades
      if(f.seg.squared_length() < 1e-6)
        continue;
      for(int i = 0; i < 4; i++){
        //auto it = find(fpoints.begin(), fpoints.end(), p[i]);
        //id[i] = it - fpoints.begin();
        //if(it == fpoints.end())
        fpoints.push_back(p[i]);
        id[i] = fpoints.size()-1;
      }
      // accordng to normal
      Plane_3 n(p[0],p[1],p[3]);
      Vector_2 n2(n.orthogonal_vector().x(), n.orthogonal_vector().y());
      //LOG(INFO) << f.seg.source() << " " << f.seg.target() << " " << n2 << " " << f.normal;
      if(n2 * f.normal >= 0 )
         polygons.push_back({id[0], id[1], id[3], id[2]});
      else
         polygons.push_back({id[0], id[2], id[3], id[1]});

    }

    for(auto f:potential_facades){
      if(f.flag != 1) continue;
      p[0] = Point_3(CGAL::to_double(f.seg.source().x()), CGAL::to_double(f.seg.source().y()), f.heights[0]);  //sh
      p[1] = Point_3(CGAL::to_double(f.seg.source().x()), CGAL::to_double(f.seg.source().y()), f.heights[1]);  // sl
      p[2] = Point_3(CGAL::to_double(f.seg.target().x()), CGAL::to_double(f.seg.target().y()), f.heights[2]);  // th 
      p[3] = Point_3(CGAL::to_double(f.seg.target().x()), CGAL::to_double(f.seg.target().y()), f.heights[3]);  // tl
      //TODO: filter too small facades
      if(f.seg.squared_length() < 1e-6)
        continue;
      for(int i = 0; i < 4; i++){
        //auto it = find(fpoints.begin(), fpoints.end(), p[i]);
        //id[i] = it - fpoints.begin();
        //if(it == fpoints.end())
        fpoints.push_back(p[i]);
        id[i] = fpoints.size()-1;
      }
       // accordng to normal
      Plane_3 n(p[0],p[1],p[3]);
      Vector_2 n2(n.orthogonal_vector().x(), n.orthogonal_vector().y());
      if(n2 * f.normal >= 0 )
         polygons.push_back({id[0], id[1], id[3], id[2]});
      else
         polygons.push_back({id[0], id[2], id[3], id[1]});
    }
 
    string ofs = m_dir + "facade_vec.ply";
    Mesh mesh;
    for(auto& p: fpoints)
      mesh.vertices.push_back(p);
    for(auto& poly: polygons){
      for(int i = 2; i < poly.size(); i++){
          mesh.faces.push_back({poly[0], poly[i-1], poly[i]});
      }
    }
    //OFFMeshSave(ofs, mesh);
    if(!PLYTriMeshSave(mesh, config, ofs))
       LOG(INFO) << "Facade saves failed.";
    return true;
}


void building::write_all_facades(){
  typedef CGAL::Cartesian_converter<K_epec, K> To_geom;
  To_geom to_geom;
  cv::Mat img_out(m_rows, m_cols, CV_8UC3, cv::Scalar(255, 255, 255));
  for (auto it = m_farr.edges_begin(); it != m_farr.edges_end(); it++) {
    if(facade_tags[std::distance(m_farr.edges_begin(), it)]==0) continue;
    Segment_2 s(to_geom(it->source()->point()), to_geom(it->target()->point()));
    const cv::Point2i ps((int)((CGAL::to_double(s.source().x()) - m_corner.x())/m_step), (int)((m_corner.y() - CGAL::to_double(s.source().y()))/m_step));
    const cv::Point2i pt((int)((CGAL::to_double(s.target().x()) - m_corner.x())/m_step), (int)((m_corner.y() - CGAL::to_double(s.target().y()))/m_step));
    cv::line(img_out, ps, pt, cv::Scalar(0, 0, 0));
  }

  const std::string fname = m_dir + "facades_all.png";
  if (!cv::imwrite(fname, img_out))
    throw std::runtime_error("Failed to write file " + fname);
}




// fixed: write segments/arrangments to images
void building::write_all_segments() {
  // write on white canvas
  //test
  ofstream file(m_dir + "segments_all", ios::out);
  cv::Mat img_out(m_rows, m_cols, CV_8UC3, cv::Scalar(255, 255, 255));
  for (const auto &s : m_segments) {
    const cv::Point2i ps((int)((CGAL::to_double(s.source().x()) - m_corner.x())/m_step), (int)((m_corner.y() - CGAL::to_double(s.source().y()))/m_step));
    const cv::Point2i pt((int)((CGAL::to_double(s.target().x()) - m_corner.x())/m_step), (int)((m_corner.y() - CGAL::to_double(s.target().y()))/m_step));
    cv::line(img_out, ps, pt, cv::Scalar(0, 0, 0));
  file << CGAL::to_double(s.source().x())<< " " <<CGAL::to_double(s.target().x())<<" " << CGAL::to_double(s.source().y()) <<" " << CGAL::to_double(s.target().y())<<endl;

  }
  file.close();
  const std::string fname = m_dir + "segments_all.png";
  if (!cv::imwrite(fname, img_out))
    throw std::runtime_error("Failed to write file " + fname);
}

void building::write_arrangement() {
  cv::Mat img_out(m_rows, m_cols, CV_8UC3, cv::Scalar(255.0, 255.0, 255.0));
   ofstream file(m_dir + "arrangement", ios::out);
  for (auto eitr = m_arr.edges_begin(); eitr != m_arr.edges_end(); ++eitr) {
    const cv::Point2i ps(
      int((CGAL::to_double(eitr->source()->point().x())-m_corner.x())/m_step),
      int((m_corner.y()-CGAL::to_double(eitr->source()->point().y()))/m_step));
    const cv::Point2i pt(
      int((CGAL::to_double(eitr->target()->point().x())-m_corner.x())/m_step),
      int((m_corner.y()-CGAL::to_double(eitr->target()->point().y()))/m_step));
    cv::line(img_out, ps, pt, cv::Scalar(0.0, 0.0, 0.0));
    file << CGAL::to_double(eitr->source()->point().x())<<" "<<CGAL::to_double(eitr->target()->point().x())<<" "<<CGAL::to_double(eitr->source()->point().y())<<" "<<CGAL::to_double(eitr->target()->point().y())<<endl;
  }
  file.close();
  const std::string fname = m_dir + "arrangment.png";
  if (!cv::imwrite(fname, img_out))
    throw std::runtime_error("Failed to write file " + fname);
}

void building::write_facade_segments(){
    ofstream file(m_dir + "facade_segments");
    for(auto facade:m_facade)
       file << facade.seg.source().x() << " " << facade.seg.target().x() << " " << facade.seg.source().y() << " " << facade.seg.target().y();
    file.close();file.close();
}

void building::write_arrangement_labeling(std::vector<Plane_3>& m_planes, int tag, Arr_with_index& m_arrs) {
  cv::Mat img_out = cv::Mat::zeros(m_rows, m_cols, CV_8U);
  for (auto fitr = m_arrs.faces_begin(); fitr != m_arrs.faces_end(); ++fitr) {
    // do not draw black polygon, it can have huge holes
    // treat it as huge background
    if (fitr->data() <= 0)
      continue;
    std::vector<cv::Point2i> plg;
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
  std::string fname;
  if(tag == 0 )
     fname = m_dir + "floor_mrf.png";
  else
     fname = m_dir + "ceiling_mrf.png";
  if (!cv::imwrite(fname, img_out))
    throw std::runtime_error("Failed to write file " + fname);
}

void building::write_pixel_label(){
    
}