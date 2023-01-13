#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <sys/types.h>

#include "include/Area_maximization_minimization.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std::chrono;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;               // Point_2 object
typedef K::Segment_2 Segment_2;           // Segment_2 object
typedef CGAL::Polygon_2<K> Polygon_2;     // Polygon_2 object
typedef std::vector<Point_2> Points;      // vector with Point_2 objects
typedef std::vector<Segment_2> segments;  // vector with Segment_2 objects
typedef std::vector<Polygon_2> Polygon_v; // vector with Polygon_2 objects

typedef CGAL::Search_traits_2<K> T;
typedef CGAL::Fuzzy_iso_box<T> box;
typedef CGAL::Kd_tree<T> tree;

typedef std::vector<double> dist; // vector with distances from a point to an edge
typedef std::vector<int> areas;   // vector with polugon areas
typedef std::vector<int> findd;   // vector with position of visible edges from an interior point(position in polygon chain)
typedef std::vector<int> Areas;
typedef std::vector<Point_2>::iterator pveciterator;  // iterator gia vector apo points
typedef std::vector<Segment_2>::iterator segiterator; // iterator gia segments
typedef std::vector<double> distance;

int glob_flag=0;
int flag_algo = -1;    // algorithm we choose
int flag_min_max = -1; // minimization or maximixation
int option = -1;       // global local or subdivision
int L = -1;            // path length
int flagalgo = -1;
int flaginit = -1;
int flagedge = -1;
double threshold = -1; // threshold
Polygon_2 p;           // polygon
segments chain;        // chain
Points points;         // points

int main(int argc, char *argv[])
{
    polygonization_application(argv);
}

