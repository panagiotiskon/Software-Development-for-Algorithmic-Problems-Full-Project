#include "Area_maximization_minimization.hpp"

Segment_2 edge_exists(Point_2 , Point_2 , segments );
segments create_segments(Points );
int find_red_segments(Segment_2 , Points , segments , int );
segments incremental_min(Points , Points , segments , segments ,Segment_2);
segments incremental_max(Points , Points , segments , segments,Segment_2);
segments incremental(Points , Points , segments, segments,Segment_2 ); 