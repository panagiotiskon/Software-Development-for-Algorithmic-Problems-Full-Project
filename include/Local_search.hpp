#ifndef LOCAL_HPP
#define LOCAL_HPP

#include "Area_maximization_minimization.hpp"
#include "Incremental.hpp"

int construct_polygon(int i, int j, Point_2 v, Segment_2 u, segments chain, int polygon_area, int l, Points tp);
segments final_polygon(int i, int j, Point_2 v, Segment_2 u, segments chain, int l);
segments change_direction(segments chain_seg);
int find_blue_edge(Segment_2 k, Points convex_hull, segments chain, int mid);
bool point_of_segment(Segment_2 s, Point_2 p);


#endif