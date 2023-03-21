#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <sys/types.h>

#include "include/Pol_app.hpp"

#include <cstdlib>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>


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

