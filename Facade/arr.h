#include "base.h"

#include <opencv2/core.hpp>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>



/*!
 * \brief Single building arrangement.
 */
class building {
  typedef CGAL::Arr_segment_traits_2<K_epec> Arr_traits;
  typedef CGAL::Arr_face_extended_dcel<Arr_traits, int> Dcel_with_index;
  typedef CGAL::Arrangement_2<Arr_traits, Dcel_with_index> Arr_with_index;

  // segment arrangement attached data
  struct Arr_face_data {
    Arr_face_data(const std::size_t &n, const std::size_t &num) :
      distances(n, 0.0),
      area(0.0),
      normal(0.0, 0.0, 0.0),
      pnumber(num),  
      bnum(0),
      nbnum(0){}
    // averaged face distance to detected planes
    std::vector<double> distances;
    // area
    double area;
    // averaged normal
    cv::Vec3d normal;
    // contained 3D points number
    int pnumber;
    // black pixel number
    int bnum;
    // non-black pixel number
    int nbnum;
    // possible future data
    bool operator < (const Arr_face_data& f) const{
        return area < f.area || area == f.area && pnumber < f.pnumber;
    }
  };
  
  // facade
  struct Facade_data{
    Facade_data(const std::size_t i, const K_epec::Segment_2 s):
    index(i),
    seg(s){}
    Facade_data(vector<vector<Point_3>> poly):polygons(poly){}
    int index; // segment(line) index
    K_epec::Segment_2 seg;
    vector<Point_2> bottoms; // {small, big}
    vector<Point_2> ups;
    vector<vector<Point_3>> polygons;
    Vector_2 normal;
    int tag = 0; //if floor and ceiling both one side black and the other non-black, then set it to 0.
   
    bool operator < (const Facade_data& f) const{
      /*
         Point_2 s1 = {(int)std::round((CGAL::to_double(f.seg.source().x())*1000)),(int)std::round((CGAL::to_double(f.seg.source().y())*1000))};
         Point_2 s2 = {(int)std::round((CGAL::to_double(seg.source().x())*1000)),(int)std::round((CGAL::to_double(seg.source().y())*1000))};         
         Point_2 t1 = {(int)std::round((CGAL::to_double(f.seg.target().x())*1000)),(int)std::round((CGAL::to_double(f.seg.target().y())*1000))};
         Point_2 t2 = {(int)std::round((CGAL::to_double(seg.target().x())*1000)),(int)std::round((CGAL::to_double(seg.target().y())*1000))};
      */
        return f.seg.min() > seg.min() || f.seg.min() == seg.min() && f.seg.max() > seg.max();
    }
  };

public:
  building(
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
    const double &min_height);

  /*!
   * \brief Line segment arrangement modeling.
   */
  bool segment_arrangement_modeling();

private:

  /*!
   * \brief Construct regularized segment arrangement.
   * Overlay grid arrangement onto it to compute data.
   */
  void construct_segment_arrangement();

  /*!
   * \brief recognize smaller facade segments in m_facade_segs by m_arr.
   */
  void facade_segment_recognization(); 

  /*!
   * \brief Segment arrangement labeling.
   * Label set: detected planes.
   * Graph: segment arrangement faces.
   */
  void segment_arrangement_MRF_labeling(std::vector<Arr_face_data>& m_arr_data, std::vector<Plane_3>& m_planes, Arr_with_index& m_arrs);

  /*!
   * \brief Arrangement adjust.
   * floor and ceiling both have labels or not.
   */
  void arrangement_adjust();

  /*!
   * \brief Segment arrangement labeling.
   * Label set: detected ceiling and floor planes.
   * Graph: segment arrangement faces * 2.
   */
  void segment_arrangement_MRF_labeling_cf(std::vector<Arr_face_data>& m_arr_cdata, std::vector<Arr_face_data>& m_arr_fdata, std::vector<Plane_3>& m_fplanes, std::vector<Plane_3>& m_cplanes, Arr_with_index& m_carr, Arr_with_index& m_farr);


  /*!
   * \brief Extrude labeled arrangement faces to their respective planes.
   */
  bool segment_arrangement_extrusion(std::vector<Plane_3>& m_planes, int tag, Arr_with_index& m_arrs);
  
  /*!
   * \brief Extrude facade up and bottom according to m_facade and supporting points.
   */
 bool facade_extrusion();

  /*!
   * \brief Write all segments.
   */
  void write_all_segments();

  /*!
   * \brief Write all facades in m_arr.
   */
  void write_all_facades();

  /*!
   * \brief Write arrangement.
   */
  void write_arrangement();

/*!
   * \brief Write facade segments.
   */
  void write_facade_segments();

  /*!
   * \brief Write arrangement MRF labeling.
   */
  void write_arrangement_labeling(std::vector<Plane_3>& m_planes, int tag, Arr_with_index& m_arrs);

  /*!
   * \brief Write pixel label.
   */
  void write_pixel_label();


private:
  // building data naming convention prefix
  const std::string m_dir;


  // building pixel height floor and ceiling map, point number map
  const cv::Mat m_hfmap;
  const cv::Mat m_hcmap;
  const cv::Mat m_nfmap;
  const cv::Mat m_ncmap;

  // min height
  const double m_min_height;

  // building bbox
  const CGAL::Iso_rectangle_2<K_epec> m_bbox;

  const int m_rows;
  const int m_cols;
  const double m_step;
  const Point_2 m_corner;

  // segments
  const std::vector<K_epec::Segment_2> m_segments;
  const std::vector<K_epec::Segment_2> m_facade_segs; // segments_final
  std::vector<Facade_data> m_facade; // smaller segments
  std::vector<int> facade_tags; //0-nonfacade 1-facade, record m_arr_edges

  // detected shapes
  std::vector<Plane_3> m_fplanes;
  std::vector<Plane_3> m_cplanes;

  // segment arrangement
  Arr_with_index m_arr;
  Arr_with_index m_farr;
  Arr_with_index m_carr;

  // segment arrangement data
  std::vector<Arr_face_data> m_arr_fdata;
  std::vector<Arr_face_data> m_arr_cdata;

  //facade solution
  std::vector<Facade> facades;

};

