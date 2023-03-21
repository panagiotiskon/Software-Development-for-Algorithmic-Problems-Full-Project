#include "Area_maximization_minimization.hpp"
#include "Incremental.hpp"
#include "Local_search.hpp"

using namespace std::chrono;

static Segment_2 error = Segment_2(Point_2(-1, -1), Point_2(-1, -1));

void handle_input_p1(char **argv)
{
  if (strcmp("incremental", argv[18]) == 0)
  {
    flagalgo = 1; // if -algorithm=incremental
    if (strcmp("1a", argv[22]) == 0)
      flaginit = 1;
    else if (strcmp("1b", argv[22]) == 0)
      flaginit = 2;
    else if (strcmp("2a", argv[22]) == 0)
      flaginit = 3;
    else if (strcmp("2b", argv[22]) == 0)
      flaginit = 4;
  }
  else if (strcmp("convex_hull", argv[18]) == 0)
  {
    flagalgo = 2; // if -algorithm=convex_hull
  }
  if (strcmp("1", argv[20]) == 0)
  {
    flagedge = 1; // if -edgeselection=random
  }
  else if (strcmp("2", argv[20]) == 0)
  {
    flagedge = 2; // if -edgeselection=min
  }
  else
  {
    flagedge = 3; // if -edgeselection=max
  }
}
void handle_input(char **argv)
{
  std::ifstream in(argv[2]); // input file
  std::ofstream outfile;
  outfile.open(argv[4]); // output file
  if (strcmp("local_search", argv[6]) == 0)
  {
    flag_algo = 1; // if -algorithm=local_search
    threshold = atof(argv[11]);
  }
  else if (strcmp("simulated_annealing", argv[6]) == 0)
  {
    flag_algo = 2; // if -algorithm=simulated_annealing
    if (strcmp("local", argv[11]) == 0)
    {
      option = 1; // option=local
    }
    else if (strcmp("global", argv[11]) == 0)
    {
      option = 2; // option=global
    }
    else if (strcmp("subdivision", argv[11]) == 0)
    {
      option = 3;
    }
  }
  L = atoi(argv[8]);
  if (strcmp("-min", argv[9]) == 0)
  {
    flag_min_max = 1; // if -min
  }
  else if (strcmp("-max", argv[9]) == 0)
  {
    flag_min_max = 2; // if -max
  }
}
int create_polygon(char *arg)
{
  int how_many_points;

  std::string input_name;
  input_name.assign(arg, arg + 24); // finding how many points we have
  while (input_name.at(0) < '1' || input_name.at(0) > '9')
  {
    input_name.erase(input_name.begin());
  }

  int i = 0;
  while (input_name.at(i) >= '0' && input_name.at(i) <= '9')
  {
    i++;
  }
  input_name.erase(input_name.begin() + i, input_name.end());

  how_many_points = stoi(input_name);

  std::string line;
  std::ifstream poly("output_polygon.txt");
  std::getline(poly, line); // skip first line
  i = 0;
  char *temp;
  int pos;
  int x;
  int y;
  std::string sub1;
  std::string sub2;
  while (i < how_many_points)
  { // extracting the points of the polygon from my output and making the polygon we will use
    std::getline(poly, line);
    pos = line.find(" ");
    sub1 = line.substr(0, pos);
    sub2 = line.substr(pos + 1);
    x = stoi(sub1);
    y = stoi(sub2);
    p.push_back(Point_2(x, y));
    i++;
  }
  return how_many_points;
}
Polygon_2 create_polygon_2(Points subset_polygon_points) // EDW NA GINEI TO DEBUG
{
  std::ofstream ofs("output_polygon_subdiv.txt"); // create a new text file for its subset polygon
  int j;
  for (j = 0; j < subset_polygon_points.size(); j++)
  {
    std::cout << subset_polygon_points[j] << std::endl;
  }
  if (flagalgo == 2)
  {
    convex_hull_fun(subset_polygon_points, ofs); // create a polygon using the selected method with subset points
  }
  else if (flagalgo == 1)
  {
    incremental_fun(subset_polygon_points, ofs); // create a polygon using the selected method with subset points
  }
  std::ifstream ifs;
  std::string s;
  ifs.open("output_polygon_subdiv.txt"); // open the file we wrote the polygon so we can read the points
  getline(ifs, s);                       // skip first line
  int i;
  int pos;
  int x;
  int y;
  std::string sub1;
  std::string sub2;
  Polygon_2 temp_polygon;
  while (i < subset_polygon_points.size())
  { // extracting the points of the polygon from my output and making the polygon we will use
    std::getline(ifs, s);
    pos = s.find(" ");
    sub1 = s.substr(0, pos);
    sub2 = s.substr(pos + 1);
    x = stoi(sub1);
    y = stoi(sub2);
    temp_polygon.push_back(Point_2(x, y));
    i++;
  } // create the polygon from the outfile we created
  return temp_polygon;
}

void create_chain(int how_many_points)
{
  int i = 0;
  while (i < how_many_points)
  {
    chain.push_back(Segment_2(p[i], p[i + 1]));
    if (i == how_many_points - 2)
    {
      break;
    }
    i++;
  }
  chain.push_back(Segment_2(p[i + 1], p[0])); // made the segments of the polygon
}
void get_points(int how_many_points)
{
  int i;
  for (i = 0; i < how_many_points; i++)
  {
    points.push_back(p[i]);
  }
}
void local_search(void)
{
  // pana
}
double simulated_annealing(int how_many_points, std::ofstream &outp)
{
  double previous_energy;
  double convex_hull_area;
  double polygon_area;
  double polygon_area2;
  srand((unsigned)time(NULL));
  Points convex_hull_points;
  CGAL::convex_hull_2(p.begin(), p.end(), std::back_inserter(convex_hull_points));           // find the convex hull
  CGAL::area_2(convex_hull_points.begin(), convex_hull_points.end(), convex_hull_area, K()); // compute convex hull area
  CGAL::area_2(p.begin(), p.end(), polygon_area, K());                                       // compute polygons area

  if (flag_min_max == 1)
  {                                                                                  // min
    previous_energy = how_many_points * (abs(polygon_area) / abs(convex_hull_area)); // finding initial energy for minimization
    if (option == 1)
    {
      polygon_area2 = sa_local(previous_energy, how_many_points);
    }
    else if (option == 2)
    { // global
      polygon_area2 = sa_global(previous_energy, how_many_points);
    }
    else if (option == 3)
    { // subdivision
      sa_subdiv(how_many_points);
    }
  }
  else if (flag_min_max == 2)
  {                                                                                        // max
    previous_energy = how_many_points * (1 - (abs(polygon_area) / abs(convex_hull_area))); // finding energy for maximization
    if (option == 1)
    { // local
      polygon_area2 = sa_local(previous_energy, how_many_points);
    }
    else if (option == 2)
    {
      polygon_area2 = sa_global(previous_energy, how_many_points);
    }
    else if (option == 3)
    { // subdivision
      sa_subdiv(how_many_points);
    }
  }

  chain.clear();
  create_chain(how_many_points);
  double ratio=print_result(how_many_points, convex_hull_area, polygon_area, polygon_area2, outp);
  return ratio;
}
int find_intersection(Segment_2 candicate_edge, Segment_2 pr, Segment_2 qs)
{
  int flag;
  flag = 0;
  auto result = CGAL::intersection(candicate_edge, pr); // checking if the candicate edge intersects with pr edge
  if (result)
  {
    if (const Segment_2 *s = boost::get<Segment_2>(&*result)) // if they intersect to a segment then we have an intersection
    {
      flag = 1;
    }
    else
    {
      Point_2 *p = boost::get<Point_2>(&*result); // if they intersect to p point we dont have an intersection
      if (*p == pr[0])
      {
        flag = 0;
      }
      else // if its any other point we have an intersect
      {
        flag = 1;
      }
    }
  }
  if (flag == 1)
  {
    return 1;
  }
  // same process for the qs edge
  result = CGAL::intersection(candicate_edge, qs);
  if (result)
  {
    if (Segment_2 *s = boost::get<Segment_2>(&*result))
    {
      flag = 1;
    }
    else
    {
      Point_2 *p = boost::get<Point_2>(&*result);
      if (*p == qs[1]) // if they intersect to s point we dont have an intersection
      {
        flag = 0;
      }
      else
      {
        flag = 1;
      }
    }
  }

  if (flag == 1)
  {
    return 1; // we have an intersection
  }
  else
  {
    return 0;
  }
}
int find_intersection_1(Segment_2 segm1, Segment_2 segm2)
{
  int flag;
  flag = 0;
  auto result = CGAL::intersection(segm1, segm2); // checking if the candicate edge intersects with pr edge
  if (result)
  {
    if (const Segment_2 *s = boost::get<Segment_2>(&*result)) // if they intersect to a segment then we have an intersection
    {
      flag = 1;
    }
    else
    {
      if (const Point_2 *p = boost::get<Point_2>(&*result))
      {
        flag = 1;
      }
    }
  }
  return flag;
}
void create_new_polygon()
{
  p.clear();
  int i;
  p.push_back(chain[0][0]);
  for (i = 0; i < chain.size() - 1; i++)
  {
    p.push_back(chain[i][1]);
  }
}
int find_intersection_2(Segment_2 segm1)
{
  int i;
  int flag;
  for (i = 0; i < chain.size(); i++)
  {
    flag = 0;
    if (chain[i] != segm1)
    {
      auto result = CGAL::intersection(segm1, chain[i]); // checking if the candicate edge intersects with pr edge
      if (result)
      {
        if (const Segment_2 *s = boost::get<Segment_2>(&*result)) // if they intersect to a segment then we have an intersection
        {
          flag = 1;
        }
        else
        {
          Point_2 *p = boost::get<Point_2>(&*result);
          if (*p == segm1[0] || *p == segm1[1])
          {
            flag = 0;
          }
          else
          {
            flag = 1;
          }
        }
      }
    }
    else
    {
      continue;
    }
    if (flag == 1)
    {
      return 1;
    }
  }
  return 0;
}
//----------------------------------------------------------------------------------


Points handleinput(std::ifstream &in, Points result2)
{
  std::string line;
  std::string temp;
  char ccc;
  int i;
  int n;
  int n1;
  int length;
  while (in.peek() != EOF)
  {
    // get each line
    std::getline(in, line);
    i = 0;
    ccc = line.at(i);
    while (ccc != '\t')
    {
      i++;
      ccc = line.at(i);
    }
    line.erase(0, i + 1); // finding the first characters that i dont want 0,1,2,3,4,5 etc and erase them
    temp = line;
    i = 0;
    ccc = line.at(i);
    while (ccc != '\t')
    {
      i++;
      ccc = line.at(i);
    }
    line.erase(line.begin() + i, line.end()); // finding the first point of the line
    n = stoi(line);                           // converting it to int
    line = temp;
    length = trunc(log10(n)) + 1;
    line.erase(0, length + 1);         // finding the second point of the line
    n1 = stoi(line);                   // converting it to int
    result2.push_back(Point_2(n, n1)); // pushing each point to the vector
  }
  return result2;
}
Point_2 pointdistance1(Points interior, segments ch, dist d, int count)
{
  // finding the count closest interior point from an edge
  int i;
  Point_2 temppoint;
  d.clear();
  for (i = 0; i < interior.size(); i++)
  {
    d.push_back(CGAL::squared_distance(ch.front(), interior[i])); // for every point find the distance from an edge
  }
  std::sort(d.begin(), d.end()); // sort the distances
  for (i = 0; i < interior.size(); i++)
  {
    if (CGAL::squared_distance(ch.front(), interior[i]) == d[count]) // finding the point i want
    {
      temppoint = interior[i];
    }
  }
  return temppoint;
}
Point_2 pointdistance(Points interior, segments ch, dist d)
{ // same as pointdistance 1 but finds the closest point from an edge
  int i;
  Point_2 temppoint;
  d.clear();
  for (i = 0; i < interior.size(); i++)
  {
    d.push_back(CGAL::squared_distance(ch.front(), interior[i]));
  }
  std::sort(d.begin(), d.end());
  for (i = 0; i < interior.size(); i++)
  {
    if (CGAL::squared_distance(ch.front(), interior[i]) == d.front())
    {
      temppoint = interior[i];
    }
  }
  return temppoint;
}
segments findvisible(Point_2 t, segments temp, segments chain1, segments visible)
{
  // takes one edge and two edges that conect 1 interior point with this edge and test for every edge if they intersect with any other edge
  int c = 0;
  int flag;
  while (c < chain1.size())
  {
    temp.push_back(Segment_2(t, chain1[c][0]));                   // this is an edge from interior point to a peak of a polygons edge
    temp.push_back(Segment_2(t, chain1[c][1]));                   // this is an edge from interior point to the other peak of a polygons edge
    flag = findintersection(temp[0], temp[1], chain1, chain1[c]); // check if its visible
    if (flag == 0)                                                // no intersections so its a visible edge
    {
      visible.push_back(chain1[c]);
    }
    temp.clear();
    c++;
  }
  return visible;
}
int findintersection(Segment_2 interioredge1, Segment_2 interioredge2, segments pchain, Segment_2 chainedge)
{
  // finds if two edges intersect to a point or to an edge
  int c = 0;
  int flag;
  while (c < pchain.size())
  {
    flag = 0;
    auto result = CGAL::intersection(interioredge1, pchain[c]); // giving one edge from a point to a peak of a polygons edge
    if (result)
    {
      if (const Segment_2 *s = boost::get<Segment_2>(&*result)) // if they intersect to a segment then we have an intersection
      {
        flag = 1;
      }
      else
      {
        Point_2 *p = boost::get<Point_2>(&*result); // if they intersect to a point and the point is the peak of the edge then we dont have an intersection
        if (pchain[c] == chainedge && *p == interioredge1[1])
        {
          flag = 0;
        }
        else // if its any other point we have an intersect
        {
          flag = 1;
        }
      }
    }
    // same process for the edge from interior point to the other peak of the edge
    result = CGAL::intersection(interioredge2, pchain[c]);
    if (result)
    {
      if (Segment_2 *s = boost::get<Segment_2>(&*result))
      {
        flag = 1;
      }
      else
      {
        Point_2 *p = boost::get<Point_2>(&*result);
        if (pchain[c] == chainedge && *p == interioredge2[1])
        {
          flag = 0;
        }
        else
        {
          flag = 1;
        }
      }
    }
    c++;
  }
  return flag;
}
int check_inside(Point_2 pt, Point_2 *pgn_begin, Point_2 *pgn_end, K traits)
{
  // checks if the polygon surrounds every point or not
  int flag = 0;
  switch (CGAL::bounded_side_2(pgn_begin, pgn_end, pt, traits))
  {
  case CGAL::ON_BOUNDED_SIDE:
    break;
  case CGAL::ON_BOUNDARY:
    break;
  case CGAL::ON_UNBOUNDED_SIDE:
    flag = 1;
    break;
  }
  return flag;
}

bool comp1a(Point_2 pt1, Point_2 pt2)
{
  if (pt1.x() == pt2.x())
  {
    return pt1.y() < pt2.y();
  }
  else
    return (pt1.x() < pt2.x());
}
bool comp1b(Point_2 pt1, Point_2 pt2)
{
  if (pt1.x() == pt2.x())
  {
    return pt1.y() > pt2.y();
  }
  else
    return (pt1.x() > pt2.x());
}
bool comp2a(Point_2 pt1, Point_2 pt2)
{
  if (pt1.y() == pt2.y())
  {
    return pt1.x() < pt2.x();
  }
  else
    return (pt1.y() < pt2.y());
}
bool comp2b(Point_2 pt1, Point_2 pt2)
{
  if (pt1.y() == pt2.y())
  {
    return pt1.x() > pt2.x();
  }
  else
    return (pt1.y() > pt2.y());
}

Points init_1a(Points p)
{
  std::sort(p.begin(), p.end(), comp1a);
  return p;
}

Points init_1b(Points p)
{
  std::sort(p.begin(), p.end(), comp1b);
  return p;
}
Points init_2a(Points p)
{
  std::sort(p.begin(), p.end(), comp2a);
  return p;
}
Points init_2b(Points p)
{
  std::sort(p.begin(), p.end(), comp2b);
  return p;
}

void convex_hull_fun(Points result2, std::ofstream &outfile)
{

  Points result1; // init vector with convex hull's points
  Points result;  // init vector with convex hull's points
  double conarea; // convec hull area
  int j;
  int k;
  int i;
  int numberofsegments; // number of segments
  segments chain1;      // init vector that keeps the polygon chain
  Polygon_2 poly;       // init my polygon
  int pointcount;       // if its 0 then i found a visible edge for this point if its more than zero then i went to the next closest point of an edge to searh visible edges
  dist dis;             // init vector that keeps the distances between interior point and edges
  int flag1;            // value 0 if the interior point is inside the polygon and 1 if outside the polygon
  segments visible;     // init vector that keeps the visible edges
  int random;           // random number i create for random edge selection
  segments temp;        // init vector with temporary edges needed
  int e;
  int pos;      // position of an edge in the polygon chain
  int counterr; // using this counter to check each visible edge in min max edge selection
  Points keep;  // init vector that keeps the temporary pollygon i make
  areas areas;  // init vector with area of polygons
  findd finde;  // init vector with position of visible edges from an interior point(position in polygon chain)
  int index;    // position of an edge in the polygon chain
  double re;    // area of a polygon
  double area;  // area of a polygon

  std::copy(result2.begin(), result2.end(), std::back_inserter(result1));          // creating a copy of vector result2
  CGAL::convex_hull_2(result2.begin(), result2.end(), std::back_inserter(result)); // geting the convex hull from all the points i have
  CGAL::area_2(result.begin(), result.end(), conarea, K());
  for (j = 0; j < result.size(); j++)
  {
    result1.erase(std::remove(result1.begin(), result1.end(), result[j]), result1.end()); // erase the points from the convex hull so i can keep the interior points
  }
  numberofsegments = result.size();
  for (i = 0; i < numberofsegments - 1; i++)
  {
    chain1.push_back(Segment_2(result[i], result[i + 1])); // initialize the polugon chain,polugon chain should be equal with the convex hull at first
  }
  chain1.push_back(Segment_2(result[i], result[0]));
  for (i = 0; i < result.size(); i++)
  {
    poly.push_back(result[i]); // initializing polygon
  }
  pointcount = 0;             // we are looking for the closest interior point from an edge
  while (result1.size() != 0) // while my polygon does not contain every interior point
  {
    Point_2 t; // t is the closest point from an edge
    i = 0;
    if (pointcount != 0)
    {
      t = pointdistance1(result1, chain1, dis, pointcount); // couldnt find visible edge from this point find the next closest point
    }
    else
    {
      t = pointdistance(result1, chain1, dis); // function returns closest point from an edge
    }
    if (flagedge == 1) // if edge selection random is selected
    {
      flag1 = 0;
      visible = findvisible(t, temp, chain1, visible); // function returns visible edges from the interior point t
      while (visible.size() != 0)                      // check every visible edge
      {
        random = rand() % visible.size(); // find a ranodm number
        e = 0;
        while (e < chain1.size())
        {
          if (visible[random] == chain1[e])
          {
            pos = e; // find visible edge in the polygon chain1
            break;
          }
          e++;
        }
        Segment_2 tempppp = chain1[pos];                                   // temp store the visible edge
        Point_2 temppp = Point_2(chain1[pos][1]);                          // temp store the second point of the edge
        chain1.insert(chain1.begin() + pos, Segment_2(chain1[pos][0], t)); // insert in chain1 at the position that the visible edge was found the new edge connecting with the interior point
        chain1.insert(chain1.begin() + pos + 1, Segment_2(t, temppp));
        chain1.erase(chain1.begin() + pos + 2);
        result1.erase(std::find(result1.begin(), result1.end(), t)); // delete the interior point
        poly.clear();
        keep.clear();
        e = 0;
        poly.push_back(chain1[e][0]); // initialize the new polygon
        keep.push_back(chain1[e][0]); // make a copy of this polygon
        poly.push_back(chain1[e][1]);
        keep.push_back(chain1[e][1]);
        e++;
        while (e < chain1.size() - 1)
        {
          poly.push_back(chain1[e][1]);
          keep.push_back(chain1[e][1]);
          e++;
        }
        e = 0;
        while (e < result1.size())
        {
          flag1 = check_inside(result1[e], &*keep.begin(), &*keep.end(), K()); // check if with the new edges i made the polygon surrrounds every point
          e++;
        }
        if (flag1 == 1 || poly.is_simple() == 0) // if with the new edges the polugon is not simple or the polugon does not surround every point i am backtracking
        {
          chain1.erase(chain1.begin() + pos);
          chain1.erase(chain1.begin() + pos);
          chain1.insert(chain1.begin() + pos, tempppp);                      // backtracking deleting the edges i created and bringing back the previous edge
          result1.push_back(t);                                              // pushing back again the interior point
          visible.erase(std::find(visible.begin(), visible.end(), tempppp)); // deleting the visible edges because it doesnt meet the criteria
          if (visible.size() == 0)
          {
            pointcount++; // if there are not visible edges that meet the criteria need to search visible edge for the next close point
          }
        }
        else
        {
          pointcount = 0;
          break;
        }
      }
      visible.clear();
    }
    else if (flagedge == 2 || flagedge == 3) // if edge selection is min or max
    {
      flag1 = 0;
      counterr = 0;
      visible = findvisible(t, temp, chain1, visible); // finds visible edges
      while (visible.size() > counterr)                // checks every visible edge
      {
        random = counterr;
        e = 0;
        while (e < chain1.size())
        {
          if (visible[random] == chain1[e])
          {
            pos = e;
            break;
          }
          e++;
        } // finds visible edge position in chain1
        // same prosess as for random edge selection
        Segment_2 tempppp = chain1[pos];
        Point_2 temppp = Point_2(chain1[pos][1]);
        chain1.insert(chain1.begin() + pos, Segment_2(chain1[pos][0], t));
        chain1.insert(chain1.begin() + pos + 1, Segment_2(t, temppp));
        chain1.erase(chain1.begin() + pos + 2);
        result1.erase(std::find(result1.begin(), result1.end(), t));
        poly.clear();
        keep.clear();
        e = 0;
        poly.push_back(chain1[e][0]);
        keep.push_back(chain1[e][0]);
        poly.push_back(chain1[e][1]);
        keep.push_back(chain1[e][1]);
        e++;
        while (e < chain1.size() - 1)
        {
          poly.push_back(chain1[e][1]);
          keep.push_back(chain1[e][1]);
          e++;
        }
        e = 0;
        while (e < result1.size())
        {
          flag1 = check_inside(result1[e], &*keep.begin(), &*keep.end(), K());
          e++;
        }
        // same process
        if (flag1 == 1 || poly.is_simple() == 0)
        {
          chain1.erase(chain1.begin() + pos);
          chain1.erase(chain1.begin() + pos);
          chain1.insert(chain1.begin() + pos, tempppp);
          result1.push_back(t);
          if (visible.size() == 0)
          {
            pointcount++;
          }
        }
        else
        {
          // found visible edge that meets the criteria
          pointcount = 0;
          CGAL::area_2(poly.begin(), poly.end(), re, K()); // compute polugons area with the current edges
          areas.push_back(re);
          finde.push_back(counterr); // store the position of the visible edge that meets the criteria
          chain1.erase(chain1.begin() + pos);
          chain1.erase(chain1.begin() + pos);
          chain1.insert(chain1.begin() + pos, tempppp); // backtrack for more visible edges that meet the criteria
          result1.push_back(t);
        }
        counterr++; // go to the next edge
      }
      if (pointcount == 0) // found visible edge
      {
        if (flagedge == 2) // if edge selection =min
        {
          index = std::distance(std::begin(areas), std::min_element(std::begin(areas), std::end(areas))); // find min area
        }
        else if (flagedge == 3) // if edge selection =max
        {
          index = std::distance(std::begin(areas), std::max_element(std::begin(areas), std::end(areas))); // find max area
        }
        random = finde[index]; // this is the index of the visible edge that makes a polygon with min/max area
        e = 0;
        while (e < chain1.size()) // finding the edge in polygon chain1
        {
          if (visible[random] == chain1[e])
          {
            pos = e;
            break;
          }
          e++;
        }
        // same process creating the edges that make a polygon with min/max area
        Segment_2 tempppp = chain1[pos];
        Point_2 temppp = Point_2(chain1[pos][1]);
        chain1.insert(chain1.begin() + pos, Segment_2(chain1[pos][0], t));
        chain1.insert(chain1.begin() + pos + 1, Segment_2(t, temppp)); // ta kanw ola auta gia na exw mia swsth seira me ta edges gia na ftiaxw eukola to polugwno
        chain1.erase(chain1.begin() + pos + 2);
        result1.erase(std::find(result1.begin(), result1.end(), t));
        poly.clear();
        keep.clear();
        e = 0;
        poly.push_back(chain1[e][0]);
        keep.push_back(chain1[e][0]);
        poly.push_back(chain1[e][1]);
        keep.push_back(chain1[e][1]);
        e++;
        while (e < chain1.size() - 1)
        {
          poly.push_back(chain1[e][1]);
          keep.push_back(chain1[e][1]);
          e++;
        }
      }
      finde.clear();
      areas.clear();
      visible.clear();
    }
  }
  // printing the resutlts
  k = 0;
  outfile << "Polygonization" << std::endl;
  while (k < poly.size())
  {
    outfile << poly[k] << std::endl;
    k++;
  }
  k = 0;
  while (k < chain1.size())
  {
    outfile << chain1[k][0] << "<----" << chain1[k][1] << std::endl;
    k++;
  }
  outfile << "Algorithm: convex_hull edge selection " << flagedge << std::endl;
  CGAL::area_2(poly.begin(), poly.end(), area, K());
  outfile << "Calculated area: " << abs(area) << std::endl;
  outfile << "ratio: " << conarea / abs(area) << std::endl;
  outfile << "ratio: " << abs(area) / abs(conarea) << std::endl;
}

void incremental_fun(Points result2, std::ofstream &outfile)
{
  Points result; // init vector with convex hull's points
  Points temppoint;
  pveciterator p1;
  pveciterator p2;
  Polygon_2 poly;  // init my polygon
  segments chain1; // init vector that keeps the polygon chain
  segments convex_seg;
  Points convex_hull;
  Segment_2 dont_touch;
  int k;                                                                 // variables for for and while                                              // output file
  std::copy(result2.begin(), result2.end(), std::back_inserter(result)); // kanw copy sto vector result1 ta points                                                // paw xana sthn arxh tou arxeiou
  if (flaginit == 1)
    result = init_1a(result);
  else if (flaginit == 2)
  {
    result = init_2a(result);
  }
  else if (flaginit == 3)
  {
    result = init_1b(result);
  }
  else if (flaginit == 4)
  {
    result = init_2b(result);
  }

  if (option == 3)
  {
    dont_touch = Segment_2(result[result.size() - 2], result[result.size() - 1]);

  }
  if (result.size() >= 3)
  { // build first triangle to temp
    for (int i = 0; i < 3; i++)
    {
      Point_2 t = result[i];
      temppoint.push_back(t);
    }
  }

  p1 = result.begin(); // delete three points from result
  p2 = result.begin() + 3;
  result.erase(p1, p2);

  int numberofsegments = temppoint.size();
  int i;

  for (i = 0; i < numberofsegments; i++) // initialise triangle
  {
    if (i == numberofsegments - 1)
    {
      chain1.push_back(Segment_2(temppoint[numberofsegments - 1], temppoint[0]));
      break;
    }
    chain1.push_back(Segment_2(temppoint[i], temppoint[i + 1]));
  }

  CGAL::convex_hull_2(temppoint.begin(), temppoint.end(), std::back_inserter(convex_hull));
  int d = convex_hull.size();
  for (int i = 0; i < convex_hull.size(); i++)
  {
    poly.push_back(convex_hull[i]);
  }

  convex_seg = create_segments(temppoint);
  segments chain11;
  if (flagedge == 1)
  {
    chain11 = incremental(result, temppoint, convex_seg, chain1, dont_touch);
  }
  else if (flagedge == 2)
  {
    chain11 = incremental_min(result, temppoint, convex_seg, chain1, dont_touch);
  }
  else if (flagedge == 3)
  {
    chain11 = incremental_max(result, temppoint, convex_seg, chain1, dont_touch);
  }
  poly.clear();
  int e = 0;
  poly.push_back(chain11[e][0]); // initialize the new polygon
  poly.push_back(chain11[e][1]);
  e++;
  while (e < chain11.size() - 1)
  {
    poly.push_back(chain11[e][1]);
    e++;
  }
  if (poly.is_simple() == 0)
  {
    std::cout << "ERROR" << std::endl;
  }
  outfile << "Polygonization" << std::endl;
  k = 0;
  while (k < poly.size())
  {
    outfile << poly[k] << std::endl;
    k++;
  }
  k = 0;
  while (k < chain11.size())
  {
    outfile << chain11[k][0] << "<----" << chain11[k][1] << std::endl;
    k++;
  }
  outfile << "Algorithm: incremental edge selection " << flagedge << " initialization " << flaginit << std::endl;
  double re;
  CGAL::area_2(poly.begin(), poly.end(), re, K());
  outfile << "Calculated area " << abs(re) << std::endl;
  Points convex_hull_area;
  CGAL::convex_hull_2(result2.begin(), result2.end(), std::back_inserter(convex_hull_area));
  double convex_hull_re;
  CGAL::area_2(convex_hull_area.begin(), convex_hull_area.end(), convex_hull_re, K());
  outfile << "ratio: " << abs(re) / abs(convex_hull_re) << std::endl;
}
double sa_local(double previous_energy, int how_many_points)
{ // local
  int cut_off=how_many_points*500;
  unsigned long int time_c=0;
  int i;
  int random_position; // position of a random point
  int random_position_1;
  int counter = 0;
  double sequent_energy;
  double convex_hull_area;
  double polygon_area;
  double temperature = 1; // initialize T as 1
  double metropolis_criterion;
  double DE;              // energy rate of change
  double R;               // random number between 0,1
  bool do_they_intersect; // find if they intersect
  srand((unsigned)time(NULL));
  Point_2 random_point;      // a random point
  Point_2 successor_point;   // next point from the random_point
  Point_2 d_successor_point; // next x2 point from the random_point
  Point_2 predecessor_point; // previous point from the random_point
  Points temp;               // holds the points of the fuzzy iso box
  Points tempp;
  Points convex_hull_points;
  Polygon_2 temp_polygon; // keeps a temp polygon
  tree tr;                // kd tree
  while (temperature >= 0)
  {
    high_resolution_clock::time_point start = high_resolution_clock::now();
    tempp.clear();
    temp.clear();
    convex_hull_points.clear();
    tr.clear();
    R = (float)rand() / RAND_MAX;
    temp_polygon = p;
    random_position = 1 + rand() % (how_many_points - 3); // find a random position of a point cant be the first or last point or before last point
    random_point = p[random_position];                    // find a random point
    successor_point = p[random_position + 1];             // find its successor
    d_successor_point = p[random_position + 2];           // fits its x2 successor
    predecessor_point = p[random_position - 1];           // find its predecessor
    p[random_position] = successor_point;
    p[random_position + 1] = random_point; // swap positions,random point with its successor
    random_point = p[random_position];
    successor_point = p[random_position + 1];                                                                                           // new values for succesor point and the random point
    chain.clear(); 
                                                                                                                     // clear the chain because we have a new polygon
    create_chain(how_many_points);                                                                                                      // create the new chain
    do_they_intersect = find_intersection_1(Segment_2(predecessor_point, random_point), Segment_2(successor_point, d_successor_point)); // check if pr,qs intersect // check if pr,qs intersect
    if (do_they_intersect == 1)
    { // the new edges intersect each other
      p.clear();
      p = temp_polygon;
      chain.clear();
      create_chain(how_many_points);
      continue; // they intersect so we find another point
    }
    else if (do_they_intersect == 0)
    {
      for (i = 0; i < how_many_points; i++)
      {

        tr.insert(p[i]); // construct a kd tree with every point
      }
      tr.build(); // build the tree
      tempp.push_back(predecessor_point);
      tempp.push_back(random_point);
      tempp.push_back(d_successor_point);
      tempp.push_back(successor_point);

      sort(tempp.begin(), tempp.end());
      int min_x = tempp[0].x();
      int max_x = tempp[3].x();
      int arr[] = {(int)tempp[0].y(), (int)tempp[1].y(), (int)tempp[2].y(), (int)tempp[3].y()};
      int *min_y = std::min_element(std::begin(arr), std::end(arr));
      int *max_y = std::max_element(std::begin(arr), std::end(arr));
      box search_box(Point_2(min_x, *min_y), Point_2(max_x, *max_y)); // construct a fuzzy box with the range of the points you need
      tr.search(std::back_inserter(temp), search_box);                // get the points inside the fuzzy box
      int j;
      int flag;
      int flag1;

      int result = 0;

      int result1 = 0;

      Segment_2 temp1;
      Segment_2 temp2;

      for (i = 0; i < temp.size(); i++)
      {

        result = 0;
        result1 = 0;

        flag = 0;
        flag1 = 0;
        if (temp[i] == random_point || temp[i] == successor_point)
        
        { // dont need to chekc for points q and r(successor and random) because their edges dont intersect

          continue;
        }
        else
        { // for any other point
          for (j = 0; j < chain.size(); j++)
          {

            if (chain[j][0] == temp[i] && chain[j][1] != random_point && chain[j][1] != successor_point)
            { // if i find an edge which does not have q or r as a point i need to check if it intersects
              flag = 1;
              temp1 = chain[j];
            }
            else if (chain[j][1] == temp[i] && chain[j][0] != random_point && chain[j][0] != successor_point)
            {
              flag1 = 1;
              temp2 = chain[j];
            }
          }

        }
        if (flag == 1)
        {
          result = find_intersection(temp1, Segment_2(predecessor_point, random_point), Segment_2(successor_point, d_successor_point));
        }
        if (flag1 == 1)
        {
          result1 = find_intersection(temp2, Segment_2(predecessor_point, random_point), Segment_2(successor_point, d_successor_point));
        }

        if (result == 1 || result1 == 1)
        {
          p.clear();
          p = temp_polygon;
          chain.clear();
          create_chain(how_many_points);
          continue; // they intersect so we find another point
        }
      }
    }

    // if we have a valid change
    CGAL::convex_hull_2(p.begin(), p.end(), std::back_inserter(convex_hull_points));           // find the convex hull
    CGAL::area_2(convex_hull_points.begin(), convex_hull_points.end(), convex_hull_area, K()); // compute convex hull area
    CGAL::area_2(p.begin(), p.end(), polygon_area, K());                                       // compute polygons area
    if (flag_min_max == 1)
    {
      // if min
      sequent_energy = how_many_points * (abs(polygon_area) / abs(convex_hull_area));
    }
    else if (flag_min_max == 2)
    {
      // if max
      sequent_energy = how_many_points * (1 - (abs(polygon_area) / abs(convex_hull_area))); // new energy
    }
    DE = sequent_energy - previous_energy; // DE
    metropolis_criterion = exp(-DE / temperature);
    if (DE < 0 || metropolis_criterion >= R)
    {
      temperature = temperature - (1 / L);
      previous_energy = sequent_energy;
    }
    else if (counter < L)
    {
      counter++;
      p.clear();
      p = temp_polygon;
      chain.clear();
      create_chain(how_many_points);
    }
    else if (counter == L)
    {
      break;
    }
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start).count();
    time_c+=duration;
    if(time_c>cut_off){
      if(flag_min_max==1){
        glob_flag=1;
      }
      else{
        glob_flag=2;
      }
      break;
    }
  }
  return polygon_area;
}
double sa_global(double previous_energy, int how_many_points)
{
  
  int i;
  int random_position; // position of a random point
  int random_position_1;
  int counter = 0;
  int valid_check_counter = 0; // for 10 points for example we may not find a valid move so we stop the program
  double sequent_energy;
  double convex_hull_area;
  double polygon_area;
  double temperature = 1; // initialize T as 1
  double metropolis_criterion;
  double DE; // energy rate of change
  int cut_off=how_many_points*500;
  unsigned long int time_c=0;
  double R;  // random number between 0,1
  srand((unsigned)time(NULL));
  int cant_find=0;
  Point_2 random_point;      // a random point
  Point_2 random_point_1;    // 2nd random point (global)
  Point_2 successor_point;   // next point from the random_point
  Point_2 predecessor_point; // previous point from the random_point
  Point_2 successor_point_1; // successor for the 2nd random point(global)
  Points convex_hull_points;
  Polygon_2 temp_polygon; // keeps a temp polygon
  while (temperature >= 0)
  {
    high_resolution_clock::time_point start = high_resolution_clock::now();
    // global
    if (how_many_points == 10 && valid_check_counter > 1000)
    {
      std::cout << "NO valid moves" << std::endl;
      break; // we didnt find a valid move because we have only 10 points
    }
    convex_hull_points.clear();
    R = (float)rand() / RAND_MAX;
    if (option == 3)
    {
      random_position = 2 + rand() % ((how_many_points) -2); // find a random position of a point cant be the first or second or before last or last
      int kl;
      int maxx=-1;
      Point_2 found_p;
      for(kl=0;kl<p.size();kl++){
        if(p[kl].x()>maxx){
          maxx=p[kl].x();
          found_p=p[kl];
        }
      }
      while(p[random_position]==found_p || p[random_position+1]==found_p)
        random_position = 2 + rand() % ((how_many_points) - 2);

    }                                                             // if we call from sudivision
    else if (option != 3)
    {
      random_position = 1 + rand() % ((how_many_points)-2); // find a random position of a point cant be the first or last point
    }
    temp_polygon = p;
    random_point = p[random_position];          // find a random point
    successor_point = p[random_position + 1];   // find its successor
    predecessor_point = p[random_position - 1]; // find its predecessor
    if (option == 3)
    {
      int kl;
      int maxx=-1;
      Point_2 found_p;
      for(kl=0;kl<p.size();kl++){
        if(p[kl].x()>maxx){
          maxx=p[kl].x();
          found_p=p[kl];
        }
      }
      random_position_1 = 1 + rand() % (how_many_points-1); // find random position for te other point.Cant be the first point because we use subdivision
      while(p[random_position_1+1]==found_p)
        random_position_1 = 1 + rand() % (how_many_points-1);
    }
    else if (option != 3)
    {
      random_position_1 = rand() % (how_many_points-1); // find random position for te other point.Can be the first point because we only get its successor cant be the last point
    }
    while (random_position_1 == random_position || random_position_1 == random_position + 1 || random_position_1 == random_position - 1){
      cant_find++;
      if(cant_find>how_many_points*1000){
        if(flag_min_max==1){
          glob_flag=1;
        }
        else{
          glob_flag=2;
        }
        return 0;
      }
      // new random point cant be same as the previous random point or its successor or its predecessor or its predecessor -1 because his succesor will be equal to predecessor of the previous point
      if (option == 3)
      {
        int kl;
        int maxx = -1;
        Point_2 found_p;
        for (kl = 0; kl < p.size(); kl++)
        {
          if (p[kl].x() > maxx)
          {
            maxx = p[kl].x();
            found_p = p[kl];
          }
        }
        random_position_1 = 1 + rand() % (how_many_points-1); // find random position for te other point.Cant be the first point because we use subdivision
        while (p[random_position_1 + 1] == found_p)
          random_position_1 = 1 + rand() % (how_many_points-1);
      }
      else if (option != 3)
      {
        random_position_1 = rand() % (how_many_points-1); // find random position for te other point.Can be the first point because we only get its successor cant be the last point
      }                                                 // find another point
    }
    random_point_1 = p[random_position_1];        // find random point 2
    successor_point_1 = p[random_position_1 + 1]; // find its successor
    chain.erase(std::find(chain.begin(), chain.end(), Segment_2(random_point, successor_point))); // erase segment qp
    for (i = 0; i < chain.size(); i++)
    {
      if (chain[i][0] == Point_2(predecessor_point))
      {
        chain[i] = Segment_2(Point_2(predecessor_point), Point_2(successor_point)); // make the new edge pr
      }
    }
    int pos;
    for (i = 0; i < chain.size(); i++)
    {
      if (Segment_2(random_point_1, successor_point_1) == chain[i])
      {
        chain[i] = Segment_2(Point_2(random_point_1), Point_2(random_point)); // connect s and q
        pos = i;
        break;
      }
    }
    chain.insert(chain.begin() + pos + 1, Segment_2(random_point, successor_point_1)); // connect q and t
    create_new_polygon();
    int result;
    int result3;
    //std::cout << successor_point << std::endl;
    result = find_intersection_1(Segment_2(predecessor_point, successor_point), Segment_2(random_point_1, random_point));
    result3 = find_intersection_1(Segment_2(predecessor_point, successor_point), Segment_2(random_point, successor_point_1));
    // check if pr intersects sq qt
    if (result == 1 || result3 == 1 || p.is_simple()==0)
    {
      p.clear();
      p = temp_polygon;
      chain.clear();
      create_chain(how_many_points);
      valid_check_counter++;
      continue; // they intersect so we find another point
    }
    else
    {
      int result1;
      int result2;
      result = find_intersection_2(Segment_2(predecessor_point, successor_point));
      result1 = find_intersection_2(Segment_2(random_point_1, random_point));
      result2 = find_intersection_2(Segment_2(random_point, successor_point_1));
      if (result == 1 || result1 == 1 || result2 == 1 || p.is_simple()==0)
      {
        p.clear();
        p = temp_polygon;
        chain.clear();
        create_chain(how_many_points);
        valid_check_counter++;
        continue; // they intersect so we find another point
      }
    }
    valid_check_counter = 0;
    CGAL::convex_hull_2(p.begin(), p.end(), std::back_inserter(convex_hull_points));           // find the convex hull
    CGAL::area_2(convex_hull_points.begin(), convex_hull_points.end(), convex_hull_area, K()); // compute convex hull area
    CGAL::area_2(p.begin(), p.end(), polygon_area, K());                                       // compute polygons area
    if (flag_min_max == 1)
    {
      // if min
      sequent_energy = how_many_points * (abs(polygon_area) / abs(convex_hull_area));
    }
    else if (flag_min_max == 2)
    {
      // if max
      sequent_energy = how_many_points * (1 - (abs(polygon_area) / abs(convex_hull_area))); // new energy
    }
    DE = sequent_energy - previous_energy; // DE
    metropolis_criterion = exp(-DE / temperature);
    if (DE < 0 || metropolis_criterion >= R)
    {
      temperature = temperature - (1 / L);
      previous_energy = sequent_energy;
    }

    else if (counter < L)
    {
      counter++;
      p.clear();
      p = temp_polygon;
      chain.clear();
      create_chain(how_many_points);
    }
    else if (counter == L)
    {
      break;
    }
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start).count();
    time_c+=duration;
    if(time_c>cut_off){
      if(flag_min_max==1){
        glob_flag=1;
      }
      else{
        glob_flag=2;
      }
      return 0;
    }

  }
  return polygon_area;
}
void sa_subdiv(int how_many_points)
{
  int k;
  int m = 100; // define m for subvision
  int i;
  int point_counter = 0;
  double previous_energy;
  double convex_hull_area;
  double polygon_area;
  p.clear();
  chain.clear();
  int end_it = 0;
  int must_end = 0;
  sort(points.begin(), points.end()); // sort points by x
  for (k = 0; k < points.size(); k++)
  {
    std::cout << points[k] << std::endl;
  }
  std::cout << std::endl;
  Points keep_that;
  Points temp_points;
  int j;
  int flag_merge = 0;
  k = (how_many_points) / (m); // k subsets
  int rest = how_many_points % m;
  if (rest != 0)
  {
    k = k + 1;
  }
  std::cout << k << std::endl;
  Polygon_v polygon_v;
  Polygon_v polygon_sorted_v;
  Polygon_2 temp_polygon1;
  Points shared_point;
  for (i = 0; i < k; i++)
  { // for every subset
    if (i == 0)
    {
      copy(points.begin() + point_counter, points.begin() + m + point_counter, back_inserter(temp_points));
    }
    else if (i > 0 && i < k - 1)
    {
      copy(points.begin() + point_counter - 1, points.begin() + point_counter + m, back_inserter(temp_points));
    }
    if (i != k - 1)
    {
      point_counter = point_counter + m;
    }
    if (i == 0)
    {
      int important_flag = 0;
      int important_counter = 0;
      while (important_flag != 1)
      {
        Point_2 last_point = temp_points[m - 1 + important_counter];        // get last point
        Point_2 before_last_point = temp_points[m - 2 + important_counter]; // get before last point
        Point_2 next_point = points[point_counter + important_counter];     // next point
        int last_p_y = (int)last_point.y();                                 // get last point y
        int before_p_y = (int)before_last_point.y();                        // get before last point y
        int after_p_y = (int)next_point.y();
        if (last_p_y > before_p_y && after_p_y < last_p_y)
        {                     // if we have monotone increasing
          important_flag = 1; // we have a good subset continue
        }
        else
        {
          point_counter++;
          important_counter++;
          temp_points.push_back(points[(m - 1) + important_counter]); // add one more point till we have monotone increasing edge(thats for first subset)
        }
      }
    } // swsto
    else if (i <= k - 2)
    {
      int important_flag = 0;
      int important_counter = 0;
      while (important_flag != 1)
      {
        Point_2 last_point = temp_points[m + important_counter];            // get last point
        Point_2 before_last_point = temp_points[m - 1 + important_counter]; // get before last point
        Point_2 next_point = points[point_counter + important_counter];     // next point
        int last_p_y = (int)last_point.y();                                 // get last point y
        int before_p_y = (int)before_last_point.y();                        // get before last point y
        int after_p_y = (int)next_point.y();
        if (last_p_y > before_p_y && after_p_y < last_p_y)
        {                     // if we have monotone increasing
          important_flag = 1; // we have a good subset continue
        }
        else
        {
          important_counter++;
          if (point_counter - 1 + important_counter > how_many_points - 1)
          {
            must_end = 1;
          }
          else
          {
            temp_points.push_back(points[point_counter - 1 + important_counter]); // add one more point till we have monotone increasing edge(thats for first subset)
          }
        }
      }
      point_counter = point_counter + important_counter;
      if (how_many_points - point_counter < m)
      {
        int l;
        for (l = point_counter; l < how_many_points; l++)
        {
          temp_points.push_back(points[l]);
        }
        end_it = 1;
      }
    }
    else if (i == k - 1)
    {
      if (points.size() - point_counter < m)
      {
        for (j = point_counter; j < points.size(); j++)
        {
          keep_that.push_back(points[j]); // add the rest points to the last subset

          flag_merge = 1;
        }
      }
      else
      {
        for (j = point_counter; j < points.size(); j++)
        {
          temp_points.push_back(points[j]); // create new subset for rest points
          flag_merge = 0;
        }
      }
    }
    if (must_end == 0 && flag_merge == 0)
    {
      temp_polygon1 = create_polygon_2(temp_points);
      polygon_v.push_back(temp_polygon1); // push back to the vector the polygon with subsets points
      if (end_it == 1)
      {
        break;
      }
    }
    else if (must_end == 0 && flag_merge == 1)
    {
      temp_polygon1 = create_polygon_2(keep_that);
      polygon_v.push_back(temp_polygon1);
    }
    else if (must_end == 1 && flag_merge == 0)
    {
      temp_polygon1 = create_polygon_2(temp_points);
      polygon_v.push_back(temp_polygon1); // push back to the vector the polygon with subsets points
      break;
    }
    else if (must_end == 1 && flag_merge == 1)
    {
      temp_polygon1 = create_polygon_2(keep_that);
      polygon_v.push_back(temp_polygon1);
      break;
    }
    shared_point.push_back(temp_points.back());
    temp_points.clear();
  }
  for(i=0;i<polygon_v.size();i++){
    for(j=0;j<polygon_v[i].size();j++){
      std::cout << polygon_v[i][j] << std::endl;
      
    }
    std::cout << "_______________________" << std::endl;
  }
  std::cout << "_______________________" << std::endl;
  p.clear();
  for (i = 0; i < polygon_v.size(); i++)
  {
    //for its polygon compute the previous energy and run global
    Points convex_hull_points;
    CGAL::convex_hull_2(polygon_v[i].begin(), polygon_v[i].end(), std::back_inserter(convex_hull_points));   // find the convex hull
    CGAL::area_2(convex_hull_points.begin(), convex_hull_points.end(), convex_hull_area, K()); // compute convex hull area
    CGAL::area_2(polygon_v[i].begin(), polygon_v[i].end(), polygon_area, K());                                       // compute polygons area
    if(flag_min_max==1){
      //min
      previous_energy = polygon_v[i].size() * (abs(polygon_area) / abs(convex_hull_area)); // finding initial energy for minimization
    }
    else{
      previous_energy = polygon_v[i].size() * (1 - (abs(polygon_area) / abs(convex_hull_area))); // finding energy for maximization
    }
    for(j=0;j<polygon_v[i].size();j++){
      p.push_back(polygon_v[i][j]);
    }
    create_chain(polygon_v[i].size());
    std::cout << abs(previous_energy) << std::endl;
    sa_global(abs(previous_energy),polygon_v[i].size());
    polygon_sorted_v.push_back(p);//store the minimized/maximized polygon
    p.clear();
    chain.clear();
  }
  std::cout << "sorted " << std::endl;
  p.clear();
  for (i = 0; i < polygon_sorted_v.size(); i++)
  {
    for(j=0;j<polygon_sorted_v[i].size();j++){
      std::cout << polygon_sorted_v[i][j] << std::endl;
      p.push_back(polygon_sorted_v[i][j]);
    }
    std::cout << "______" << std::endl;
    if(polygon_sorted_v[i].is_simple()==0){
      std::cout << "er" << std::endl;
    }
    p.clear();

  }
  for (i = 0; i < polygon_sorted_v.size() - 1; i++)
  {
    auto it = find(polygon_sorted_v[i].begin(), polygon_sorted_v[i].end(), shared_point[i]);
    int index1=it-polygon_sorted_v[i].begin()-1;//position of p
    int j=0;
    Point_2 temp_p=polygon_sorted_v[i][j];
    while(temp_p!=polygon_sorted_v[i][index1]){
      p.push_back(temp_p);
      j++;
      temp_p=polygon_sorted_v[i][j];
    }
    p.push_back(temp_p);
    j++;
    j++;
    int k;
    for(k=1;k<polygon_sorted_v[i+1].size();k++){
      p.push_back(polygon_sorted_v[i+1][k]);
    }
    p.push_back(polygon_sorted_v[i+1][0]);
    for(k=j;k<polygon_sorted_v[i].size();k++){
      p.push_back(polygon_sorted_v[i][k]);
    }
    p.clear();
  }
}

double print_result(int how_many_points, double convex_hull_area, double polygon_area_1, double polygon_area_2, std::ofstream &outf)
{
  int i;
  outf << "Optimal Area Polygonization" << std::endl;
  for (i = 0; i < how_many_points; i++)
  {
    outf << p[i] << std::endl;
  }
  for (i = 0; i < how_many_points; i++)
  {
    outf << chain[i][0] << " <---- " << chain[i][1] << std::endl;
  }
  if (flag_algo == 1)
  {
    outf << "local search" << std::endl;
  }
  else if (flag_algo == 2)
  {
    outf << "simulated annealing" << std::endl;
  }
  if (flag_min_max == 1)
  {
    outf << "min" << std::endl;
  }
  else
  {
    outf << "max" << std::endl;
  }
  outf << "area initial: " << abs(polygon_area_1) << std::endl;
  outf << "area: " << abs(polygon_area_2) << std::endl;
  outf << "ratio initial: " << abs(polygon_area_1) / abs(convex_hull_area) << std::endl;
  outf << "ratio: " << abs(polygon_area_2) / abs(convex_hull_area) << std::endl;
  return abs(polygon_area_2) / abs(convex_hull_area);
}
