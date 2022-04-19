#ifndef APPROXIMATE_SQUARED_DISTANCE_H
#define APPROXIMATE_SQUARED_DISTANCE_H

// minimum segment squared length
#define EPSILON_MINIMUM_SQUARED_LENGTH 1.0e-8

/*!
 * \brief Approximate point to point squared distance.
 * \note Possible underflow and overflow of float number type.
 * \param p1x first point x
 * \param p1y first point y
 * \param p2x second point x
 * \param p2y second point y
 * \return squared distance
 */
inline double point_point_squared_distance(
  const double &p1x, const double &p1y,
  const double &p2x, const double &p2y) {
  // squared distance of two points
  const double dx = p2x - p1x;
  const double dy = p2y - p1y;
  return (dx * dx) + (dy * dy);
}

/*!
 * \brief Approximate point to segment squared distance.
 * \note Possible underflow and overflow of float number type.
 * \param px point x
 * \param py point y
 * \param p1x segment source point x
 * \param p1y segment source point y
 * \param p2x segment target point x
 * \param p2y segment target point y
 * \return squared distance
 */
inline double point_segment_squared_distance(
  const double &px, const double &py,
  const double &p1x, const double &p1y,
  const double &p2x, const double &p2y) {
  // tiny epsilon
  constexpr double EPSILON_TINY = 1.0e-14;

  // segment is a point
  const double seg_len_squared = point_point_squared_distance(p1x, p1y, p2x, p2y);
  if (seg_len_squared < EPSILON_MINIMUM_SQUARED_LENGTH)
    return point_point_squared_distance(px, py, p1x, p1y);

  // orthogonal projection intersection
  const double dx = p2x - p1x;
  const double dy = p2y - p1y;
  const double dp1x = px - p1x;
  const double dp1y = py - p1y;
  const double t = ((dp1x * dx) + (dp1y * dy)) / seg_len_squared;

  // the point is 'outside' of the segment
  if (t <= EPSILON_TINY)
    return point_point_squared_distance(px, py, p1x, p1y);
  else if (t >= (1.0 - EPSILON_TINY))
    return point_point_squared_distance(px, py, p2x, p2y);

  // the point is 'within' the bounds of the segment
  const double qx = ((1.0 - t) * p1x) + (t * p2x);
  const double qy = ((1.0 - t) * p1y) + (t * p2y);
  return point_point_squared_distance(px, py, qx, qy);
}

/*!
 * \brief Approximate segment to segment squared distance.
 * \note Possible underflow and overflow of float number type.
 * \param p1x first segment source point x
 * \param p1y first segment source point y
 * \param p2x first segment target point x
 * \param p2y first segment target point y
 * \param p3x second segment source point x
 * \param p3y second segment source point y
 * \param p4x second segment target point x
 * \param p4y second segment target point y
 * \return squared distance
 */
inline double segment_segment_squared_distance(
  const double &p1x, const double &p1y,
  const double &p2x, const double &p2y,
  const double &p3x, const double &p3y,
  const double &p4x, const double &p4y) {
  // general epsilon for places where you need more "slop" than EPSILON_TINY
  constexpr double EPSILON_GENERAL = 1.192092896e-07;

  // treat short segments as a single point
  const double seg12_squared_length = point_point_squared_distance(p1x, p1y, p2x, p2y);
  const double seg34_squared_length = point_point_squared_distance(p3x, p3y, p4x, p4y);
  if (seg12_squared_length <= EPSILON_MINIMUM_SQUARED_LENGTH) {
    if (seg34_squared_length <= EPSILON_MINIMUM_SQUARED_LENGTH)
      return point_point_squared_distance(p1x, p1y, p3x, p3y);
    else
      return point_segment_squared_distance(p1x, p1y, p3x, p3y, p4x, p4y);
  }
  else if (seg34_squared_length <= EPSILON_MINIMUM_SQUARED_LENGTH)
    return point_segment_squared_distance(p3x, p3y, p1x, p1y, p2x, p2y);

  // check if parallel by the cross product of vector34 (p3 -> p4) and vector21 (p2 -> p1)
  const double d = ((p4y - p3y) * (p1x - p2x)) - ((p1y - p2y) * (p4x - p3x));
  const bool is_parallel = std::abs(d) < EPSILON_GENERAL;

  if (!is_parallel) {
    // check for intersection
    // vector24 (p2 -> p4) cross product with vector34 (p3 -> p4)
    // t = 1.0 - (vector24.cross_product(vector34) / d);
    const double t = 1.0 - ((((p4y - p3y) * (p4x - p2x)) - ((p4y - p2y) * (p4x - p3x))) / d);
    const double seg12_eps = std::sqrt(EPSILON_MINIMUM_SQUARED_LENGTH / seg12_squared_length);
    if (t >= -seg12_eps && t <= (1.0 + seg12_eps)) {
      // inside [p1,p2], segments may intersect.
      // double s = 1.0 - (vector21.cross_product(vector24) / d);
      const double s = 1.0 - ((((p4y - p2y) * (p1x - p2x)) - ((p1y - p2y) * (p4x - p2x))) / d);
      const double seg34_eps = std::sqrt(EPSILON_MINIMUM_SQUARED_LENGTH / seg34_squared_length);
      // segments intersect
      if (s >= -seg34_eps && s <= (1.0 + seg34_eps))
        return 0.0;
    }
  }

  // segments parallel, collinear or do not intersect
  // the minimum distance of those 4 tests is the closest point
  const double dist1 = point_segment_squared_distance(p3x, p3y, p1x, p1y, p2x, p2y);
  const double dist2 = point_segment_squared_distance(p4x, p4y, p1x, p1y, p2x, p2y);
  const double dist3 = point_segment_squared_distance(p1x, p1y, p3x, p3y, p4x, p4y);
  const double dist4 = point_segment_squared_distance(p2x, p2y, p3x, p3y, p4x, p4y);

  return std::min(std::min(dist1, dist2), std::min(dist3, dist4));
}

#undef EPSILON_MINIMUM_SQUARED_LENGTH

#endif // APPROXIMATE_SQUARED_DISTANCE_H
