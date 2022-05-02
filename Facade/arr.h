#ifndef _ARR_H_ 
#define _ARR_H_

#include "base.h"
#include <opencv2/core.hpp>
#include <iomanip>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/squared_distance_2.h>

/*!
 * \description: Single building arrangement.
 */
class building 
{
    typedef CGAL::Arr_segment_traits_2<K_epec> Arr_traits;
    typedef CGAL::Arr_face_extended_dcel<Arr_traits, int> Dcel_with_index;
    typedef CGAL::Arrangement_2<Arr_traits, Dcel_with_index> Arr_with_index;

    // segment arrangement attached data
    struct Arr_face_data 
    {
        Arr_face_data(const std::size_t &n, const std::size_t &num) :
            distances(n, 0.0),
            area(0.0),
            normal(0.0, 0.0, 0.0),
            pnumber(num),  
            bnum(0),
            nbnum(0){}
        std::vector<double> distances; // averaged face distance to detected planes    
        double area; // area
        cv::Vec3d normal; // averaged normal
        int pnumber; // contained 3D points number
        int bnum; // blank pixel number
        int nbnum; // non-blank pixel number
        bool operator < (const Arr_face_data& f) const{
            return area < f.area || area == f.area && pnumber < f.pnumber;
        }
    };
  
    // smaller facade segments
    struct Facade_data 
    {
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
        int tag = 0; //if floor and ceiling both one side blank and the other non-blank, then set it to 0.
        bool operator < (const Facade_data& f) const {
            return f.seg.min() > seg.min() || f.seg.min() == seg.min() && f.seg.max() > seg.max();
        }
    };

    // grid_arr and segments_arr -> overlay_arr data 
    struct Overlay_data 
    {
      Overlay_data() : fidx(0), pixel(0, 0) {}
      Overlay_data(const int &f, const cv::Vec2i &p) : fidx(f), pixel(p) {}
      Overlay_data(const Overlay_data &obj) : fidx(obj.fidx), pixel(obj.pixel) {}
      int fidx; // segments_arr index
      cv::Vec2i pixel; // grid_arr index
    };
    struct Overlay_label 
    {
      Overlay_data operator() (const int &f, const cv::Vec2i &p) const {
        return {f, p};
      }
    };
    typedef CGAL::Arr_face_extended_dcel<Arr_traits, cv::Vec2i> Dcel_grid;
    typedef CGAL::Arrangement_2<Arr_traits, Dcel_grid> Arr_grid;
    typedef CGAL::Arr_face_extended_dcel<Arr_traits, Overlay_data> Dcel_overlay;
    typedef CGAL::Arrangement_2<Arr_traits, Dcel_overlay> Arr_overlay;
    typedef CGAL::Arr_face_overlay_traits<Arr_with_index, Arr_grid, Arr_overlay, Overlay_label> Overlay_traits;

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
     * \description: Line segment arrangement modeling.
     */
    bool segment_arrangement_modeling();

private:

    /*!
     * \description: Construct regularized segment arrangement.
     */
    void construct_segment_arrangement();

    /*!
     * \description: Compute distances between cells to detected planes.
     * \param m_fcarr arrangement
     * \param m_fcplanes detected planes
     * \param m_hfcmap height map
     * \param m_arr_fcdata arrangement attached data
     * \param overlay_arr basic arrangement unit
     * \param tag deal with ceiling (1) or floor (0)
     * \return
     */    
    void compute_distance(Arr_with_index& m_fcarr, vector<Plane_3>& m_fcplanes, const cv::Mat& m_hfcmap, const cv::Mat& m_nfcmap, vector<Arr_face_data>& m_arr_fcdata, Arr_overlay& overlay_arr, int tag);

    /*!
     * \description: Recognize smaller facade segments.
     */
    void facade_segment_recognization(); 

    /*!
     * \description: Segment arrangement labeling.
     * \param m_arr_cdata ceiling arrangement attached data
     * \param m_arr_fdata floor arrangement attached data
     * \param m_fplanes detected floor planes (label set)
     * \param m_cplanes detected ceiling planes (label set)
     * \param m_carr ceiling arrangement
     * \param m_farr floor arrangement
     */
    void segment_arrangement_MRF_labeling_cf(std::vector<Arr_face_data>& m_arr_cdata, std::vector<Arr_face_data>& m_arr_fdata, std::vector<Plane_3>& m_fplanes, std::vector<Plane_3>& m_cplanes, Arr_with_index& m_carr, Arr_with_index& m_farr);

    /*!
     * \description: Extrude labeled arrangement faces to their respective planes.
     * \param m_planes detected planes
     * \param tag deal with ceiling (1) or floor (0)
     * \param m_arrs arrangement
     * \return
     */     
    bool segment_arrangement_extrusion(std::vector<Plane_3>& m_planes, int tag, Arr_with_index& m_arrs);
    
    /*!
     * \description: Extrude facade height according to m_facade and supporting points.
     */
    bool facade_extrusion();

    /*!
     * \description: Write all segments.
     */
    void write_all_segments();

    /*!
     * \description: Write all facades to image.
     */
    void write_all_facades();

    /*!
     * \description: Write arrangement.
     */
    void write_arrangement();

    /*!
     * \description: Write facade segments to file.
     */
    void write_facade_segments();

    /*!
     * \description: Write arrangement MRF labeling.
     * \param m_planes detected planes
     * \param tag deal with ceiling or floor
     * \param m_arrs arrangement
     * \return
     */     
    void write_arrangement_labeling(std::vector<Plane_3>& m_planes, int tag, Arr_with_index& m_arrs);

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
    const vector<K_epec::Segment_2> m_segments; // total segments
    const vector<K_epec::Segment_2> m_facade_segs; // floorplan segments
    vector<Facade_data> m_facade; // smaller facade segments in m_arr_edges
    vector<int> facade_tags; // m_arr_edges: 0-non_facade 1-facade

    // detected shapes
    vector<Plane_3> m_fplanes;
    vector<Plane_3> m_cplanes;

    // segment arrangement
    Arr_with_index m_arr;
    Arr_with_index m_farr;
    Arr_with_index m_carr;

    // segment arrangement attached data
    vector<Arr_face_data> m_arr_fdata;
    vector<Arr_face_data> m_arr_cdata;

    //facade solution
    vector<Facade> facades;
};

#endif // _ARR_H_
