#include "Incremental.hpp"

static Segment_2 error = Segment_2(Point_2(-1, -1), Point_2(-1, -1));


int find_red_segments(Segment_2 k, Points convex_hull, segments chain1, int mid)
{ // returns 0 if there is no intersection of k with any of the segments of the chain1
  int flag;
  int counter = 0;

  for (int i = 0; i < chain1.size(); i++)
  {
    auto result = CGAL::intersection(k, chain1[i]);
    if (result)
    {
      if (const Segment_2 *s = boost::get<Segment_2>(&*result))
      {
        flag = 1;
        counter++;
        break;
      }
      else
      {
        Point_2 *p = boost::get<Point_2>(&*result);
        if (std::find(convex_hull.begin(), convex_hull.end(), *p) != convex_hull.end() && mid == 0) // check if p is an edge of the chain1
        {
          flag = 0;
        }
        else if (mid == 1 && *p == k[1])
        { // if we check the middle point of a segment the only point of intersection between that and the new point should be only itself
          flag = 0;
        }
        else
        {
          flag = 1;
          counter++;
          break;
        }
      }
    }
  }

  return counter == 0 ? 0 : 1;
}

Segment_2 edge_exists(Point_2 a, Point_2 b, segments chain1)
{ // returns segment if seg(a,b) exists in chain1 else returns error segment
  if (std::find(chain1.begin(), chain1.end(), Segment_2(a, b)) != chain1.end())
    return Segment_2(a, b);
  else if (std::find(chain1.begin(), chain1.end(), Segment_2(b, a)) != chain1.end())
    return Segment_2(b, a);
  else
    return error;
}

segments create_segments(Points p)
{ // give ordered set of points to create vector of segments
  segments s;
  for (int i = 0; i < p.size(); i++)
  {
    if (i != p.size() - 1)
      s.push_back(Segment_2(p[i], p[i + 1]));
  }
  s.push_back(Segment_2(p[p.size() - 1], p[0]));

  return s;
}




segments incremental_min(Points result, Points curr_points, segments convex_seg, segments chain_seg, Segment_2 dont_touch)
{
  Points convex_hull;
  Points red_points;
  segments red_segments_final;
  Polygon_2 p1;
  Polygon_2 keep;
  Areas areas;
  Areas find;

  while (result.size() != 0)
  {
    red_points.clear();
    convex_hull.clear();
    convex_seg.clear();
    red_segments_final.clear();

    CGAL::convex_hull_2(curr_points.begin(), curr_points.end(), std::back_inserter(convex_hull)); // create convex hull for curr points

    int numberofsegments = convex_hull.size();
    int i;

    for (i = 0; i < numberofsegments; i++) // initialise triangle
    {
      if (i == numberofsegments - 1)
      {
        convex_seg.push_back(Segment_2(convex_hull[numberofsegments - 1], convex_hull[0]));
        break;
      }

      convex_seg.push_back(Segment_2(convex_hull[i], convex_hull[i + 1]));
    }
    Point_2 k = result[0]; // find new point k

    for (int i = 0; i < convex_hull.size(); i++)
    {
      int res = find_red_segments(Segment_2(k, convex_hull[i]), convex_hull, convex_seg, 0); // find red segments
      if (!res)
      {
        red_points.push_back(convex_hull[i]);
      }
    }

    if (!red_points.empty()) // if some edges of the convex hull are red try to find if the match any side of the polygon
    {
      for (int i = 0; i < red_points.size(); i++)
      {
        for (int j = 1; j < red_points.size(); j++)
        {
          Segment_2 temp = edge_exists(red_points[i], red_points[j], convex_seg);
          if (edge_exists(temp[0], temp[1], chain_seg) != error)
          {
            if (std::find(red_segments_final.begin(), red_segments_final.end(), edge_exists(temp[0], temp[1], chain_seg)) == red_segments_final.end())
            {                                                                                                                            // check if segment of 2 red points exists in polygon
              int res = find_red_segments(Segment_2(k, midpoint(edge_exists(temp[0], temp[1], chain_seg))), convex_hull, convex_seg, 1); // check if middle of the red segment is visible
              if (!res)
              {
                red_segments_final.push_back(edge_exists(temp[0], temp[1], chain_seg));
              }
            }
          }
        }
      }
    }

    red_points.clear();
    Points tempp = curr_points;
    // if no visible edges in convex hull then search in polygon chain
    if (red_segments_final.size() == 0)
    {
      for (int i = 0; i < curr_points.size(); i++)
      {
        int res = find_red_segments(Segment_2(k, curr_points[i]), curr_points, chain_seg, 0);
        if (!res)
        {
          red_points.push_back(curr_points[i]); // find red edges
        }
      }
      if (!red_points.empty()) // for every combination of these red edges find if they match a segment of the polygon chain
      {
        for (int i = 0; i < red_points.size(); i++)
        {
          for (int j = 1; j < red_points.size(); j++)
          {
            Segment_2 temp = edge_exists(red_points[i], red_points[j], chain_seg);
            if (temp != error)
            {
              if (std::find(red_segments_final.begin(), red_segments_final.end(), temp) == red_segments_final.end())
              {
                int res = find_red_segments(Segment_2(k, midpoint(edge_exists(temp[0], temp[1], chain_seg))), curr_points, chain_seg, 1); // check if middle of the segment is visible to the new point
                if (!res)
                {
                  red_segments_final.push_back(temp);
                }
              }
            }
          }
        }
      }
    }

    int min = 0;
    int pos;
    int min_index = 0;
    while (min < red_segments_final.size())
    {

      int i = 0;
      while (i < chain_seg.size())
      {
        if (red_segments_final[min] == chain_seg[i])
        {
          pos = i;
          break;
        }
        i++;
      }
      Segment_2 tempppp = chain_seg[pos];
      Point_2 temppp = Point_2(chain_seg[pos][1]);
      curr_points.push_back(k);
      chain_seg.insert(chain_seg.begin() + pos, Segment_2(chain_seg[pos][0], k));
      // insert in chain at the position that the visible edge was found the new edge connecting with the interior point
      chain_seg.insert(chain_seg.begin() + pos + 1, Segment_2(k, temppp));
      chain_seg.erase(chain_seg.begin() + pos + 2);
      result.erase(result.begin());

      p1.clear();
      keep.clear();
      int e = 0;
      p1.push_back(chain_seg[e][0]);   // initialize the new polygon
      keep.push_back(chain_seg[e][0]); // make a copy of this polygon
      p1.push_back(chain_seg[e][1]);
      keep.push_back(chain_seg[e][1]);
      e++;
      while (e < chain_seg.size() - 1)
      {
        p1.push_back(chain_seg[e][1]);
        keep.push_back(chain_seg[e][1]);
        e++;
      }
      if (option == 3 && k == dont_touch[0] && std::find(chain_seg.begin(), chain_seg.end(), Segment_2(dont_touch[0], dont_touch[1])) != chain_seg.end())
      {
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.insert(chain_seg.begin() + pos, tempppp); // backtracking deleting the edges i created and bringing back the previous edge
        result.insert(result.begin(), k);
        curr_points.pop_back(); // pushing back again the interior point
      }
      if (p1.is_simple() == 0) // if with the new edges the polugon is not simple or the polugon does not surround every point i am backtracking
      {
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.insert(chain_seg.begin() + pos, tempppp); // backtracking deleting the edges i created and bringing back the previous edge
        result.insert(result.begin(), k);
        curr_points.pop_back(); // pushing back again the interior point
      }
      else
      {
        double re;
        CGAL::area_2(p1.begin(), p1.end(), re, K());
        areas.push_back(re);
        find.push_back(min); // save the indexes of the red segments
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.insert(chain_seg.begin() + pos, tempppp); // backtracking deleting the edges i created and bringing back the previous edge
        result.insert(result.begin(), k);
        curr_points.pop_back();
      }
      min++;
    }

    i = 0;
    int index = std::distance(std::begin(areas), std::min_element(std::begin(areas), std::end(areas))); // find the of the segment which creates the minimum area
    min_index = find[index];
    while (i < chain_seg.size())
    {
      if (red_segments_final[min_index] == chain_seg[i])
      {
        pos = i;
        break;
      }
      i++;
    }
    Segment_2 tempppp = chain_seg[pos];

    Point_2 temppp = Point_2(chain_seg[pos][1]);

    curr_points.push_back(k);
    chain_seg.insert(chain_seg.begin() + pos, Segment_2(chain_seg[pos][0], k));
    // insert in chain at the position that the visible edge was found the new edge connecting with the interior point
    chain_seg.insert(chain_seg.begin() + pos + 1, Segment_2(k, temppp));
    chain_seg.erase(chain_seg.begin() + pos + 2);
    result.erase(result.begin());

    p1.clear();
    keep.clear();
    int e = 0;
    p1.push_back(chain_seg[e][0]);   // initialize the new polygon
    keep.push_back(chain_seg[e][0]); // make a copy of this polygon
    p1.push_back(chain_seg[e][1]);
    keep.push_back(chain_seg[e][1]);
    e++;
    while (e < chain_seg.size() - 1)
    {
      p1.push_back(chain_seg[e][1]);
      keep.push_back(chain_seg[e][1]);
      e++;
    }
    find.clear();
    areas.clear();
  }
  return chain_seg;
}

segments incremental_max(Points result, Points curr_points, segments convex_seg, segments chain_seg, Segment_2 dont_touch)
{ // same exact function as the above but creates the polygon with the max area
  Points convex_hull;
  Points red_points;
  segments red_segments_final;
  Polygon_2 p1;
  Polygon_2 keep;

  Areas areas;
  Areas find;

  while (result.size() != 0)
  {
    red_points.clear();
    convex_hull.clear();
    convex_seg.clear();
    red_segments_final.clear();

    CGAL::convex_hull_2(curr_points.begin(), curr_points.end(), std::back_inserter(convex_hull)); // create convex hull for curr points

    int numberofsegments = convex_hull.size();
    int i;

    for (i = 0; i < numberofsegments; i++) // arxikopoiw to convex_seg me ta 3 segments
    {
      if (i == numberofsegments - 1)
      {
        convex_seg.push_back(Segment_2(convex_hull[numberofsegments - 1], convex_hull[0]));
        break;
      }

      convex_seg.push_back(Segment_2(convex_hull[i], convex_hull[i + 1]));
    }

    Point_2 k = result[0]; // brhskw to neo shmeio k

    for (int i = 0; i < convex_hull.size(); i++)
    {
      int res = find_red_segments(Segment_2(k, convex_hull[i]), convex_hull, convex_seg, 0);

      if (!res)
      {
        red_points.push_back(convex_hull[i]);
      }
    }

    if (!red_points.empty())
    {
      for (int i = 0; i < red_points.size(); i++)
      {
        for (int j = 1; j < red_points.size(); j++)
        {
          Segment_2 temp = edge_exists(red_points[i], red_points[j], convex_seg);
          if (edge_exists(temp[0], temp[1], chain_seg) != error)
          {
            if (std::find(red_segments_final.begin(), red_segments_final.end(), edge_exists(temp[0], temp[1], chain_seg)) == red_segments_final.end())
            {
              int res = find_red_segments(Segment_2(k, midpoint(edge_exists(temp[0], temp[1], chain_seg))), convex_hull, convex_seg, 1);
              if (!res)
              {
                red_segments_final.push_back(edge_exists(temp[0], temp[1], chain_seg));
              }
            }
          }
        }
      }
    }

    red_points.clear();
    Points tempp = curr_points;
    // if no visible edges in convex hull then search in polygon chain
    if (red_segments_final.size() == 0)
    {
      for (int i = 0; i < curr_points.size(); i++)
      {
        int res = find_red_segments(Segment_2(k, curr_points[i]), curr_points, chain_seg, 0);
        if (!res)
        {
          red_points.push_back(curr_points[i]);
        }
      }
      if (!red_points.empty())
      {
        for (int i = 0; i < red_points.size(); i++)
        {
          for (int j = 1; j < red_points.size(); j++)
          {
            Segment_2 temp = edge_exists(red_points[i], red_points[j], chain_seg);
            if (temp != error)
            { // for every combination of red edges check if they match a segment in polygon chain
              if (std::find(red_segments_final.begin(), red_segments_final.end(), temp) == red_segments_final.end())
              {
                int res = find_red_segments(Segment_2(k, midpoint(edge_exists(temp[0], temp[1], chain_seg))), curr_points, chain_seg, 1);
                if (!res)
                {
                  red_segments_final.push_back(temp);
                }
              }
            }
          }
        }
      }
    }
    int max = 0;
    int pos;
    int max_index = 0;

    while (max < red_segments_final.size())
    {
      int i = 0;
      while (i < chain_seg.size())
      {
        if (red_segments_final[max] == chain_seg[i])
        {
          pos = i;
          break;
        }
        i++;
      }
      Segment_2 tempppp = chain_seg[pos];
      Point_2 temppp = Point_2(chain_seg[pos][1]);
      curr_points.push_back(k);
      chain_seg.insert(chain_seg.begin() + pos, Segment_2(chain_seg[pos][0], k));
      // insert in chain at the position that the visible edge was found the new edge connecting with the interior point
      chain_seg.insert(chain_seg.begin() + pos + 1, Segment_2(k, temppp));
      chain_seg.erase(chain_seg.begin() + pos + 2);
      result.erase(result.begin());

      p1.clear();
      keep.clear();
      int e = 0;
      p1.push_back(chain_seg[e][0]);   // initialize the new polygon
      keep.push_back(chain_seg[e][0]); // make a copy of this polygon
      p1.push_back(chain_seg[e][1]);
      keep.push_back(chain_seg[e][1]);
      e++;
      while (e < chain_seg.size() - 1)
      {
        p1.push_back(chain_seg[e][1]);
        keep.push_back(chain_seg[e][1]);
        e++;
      }
      if (option == 3 && k == dont_touch[0] && std::find(chain_seg.begin(), chain_seg.end(), Segment_2(dont_touch[0], dont_touch[1])) != chain_seg.end())
      {
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.insert(chain_seg.begin() + pos, tempppp); // backtracking deleting the edges i created and bringing back the previous edge
        result.insert(result.begin(), k);
        curr_points.pop_back(); // pushing back again the interior point
      }
      if (p1.is_simple() == 0) // if with the new edges the polugon is not simple or the polugon does not surround every point i am backtracking
      {
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.insert(chain_seg.begin() + pos, tempppp); // backtracking deleting the edges i created and bringing back the previous edge
        result.insert(result.begin(), k);
        curr_points.pop_back(); // pushing back again the interior point
      }
      else
      {
        double re;
        CGAL::area_2(p1.begin(), p1.end(), re, K());
        areas.push_back(re);
        find.push_back(max);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.insert(chain_seg.begin() + pos, tempppp); // backtracking deleting the edges i created and bringing back the previous edge
        result.insert(result.begin(), k);
        curr_points.pop_back();
      }
      max++;
    }
    int index = std::distance(std::begin(areas), std::max_element(std::begin(areas), std::end(areas)));
    max_index = find[index];
    i = 0;
    while (i < chain_seg.size())
    {
      if (red_segments_final[max_index] == chain_seg[i])
      {
        pos = i;
        break;
      }
      i++;
    }
    Segment_2 tempppp = chain_seg[pos];
    Point_2 temppp = Point_2(chain_seg[pos][1]);

    curr_points.push_back(k);
    chain_seg.insert(chain_seg.begin() + pos, Segment_2(chain_seg[pos][0], k));
    // insert in chain at the position that the visible edge was found the new edge connecting with the interior point
    chain_seg.insert(chain_seg.begin() + pos + 1, Segment_2(k, temppp));
    chain_seg.erase(chain_seg.begin() + pos + 2);
    result.erase(result.begin());

    p1.clear();
    keep.clear();
    int e = 0;
    p1.push_back(chain_seg[e][0]);   // initialize the new polygon
    keep.push_back(chain_seg[e][0]); // make a copy of this polygon
    p1.push_back(chain_seg[e][1]);
    keep.push_back(chain_seg[e][1]);
    e++;
    while (e < chain_seg.size() - 1)
    {
      p1.push_back(chain_seg[e][1]);
      keep.push_back(chain_seg[e][1]);
      e++;
    }
    find.clear();
    areas.clear();
  }

  return chain_seg;
}

segments incremental(Points result, Points curr_points, segments convex_seg, segments chain_seg, Segment_2 dont_touch)
{ // random icremental
  srand((unsigned)time(NULL));

  Points convex_hull;
  Points red_points;

  segments red_segments_convex;
  segments red_segments_chain;
  segments red_segments_final;
  Polygon_2 p1;
  Polygon_2 keep;
  int global_counter = 1;
  while (result.size() != 0)
  {
    red_points.clear();
    convex_hull.clear();
    convex_seg.clear();
    red_segments_final.clear();

    CGAL::convex_hull_2(curr_points.begin(), curr_points.end(), std::back_inserter(convex_hull)); // create convex hull for curr points

    int numberofsegments = convex_hull.size();
    int i;

    for (i = 0; i < numberofsegments; i++) // initialise triangle
    {
      if (i == numberofsegments - 1)
      {
        convex_seg.push_back(Segment_2(convex_hull[numberofsegments - 1], convex_hull[0]));
        break;
      }

      convex_seg.push_back(Segment_2(convex_hull[i], convex_hull[i + 1]));
    }

    Point_2 k = result[0]; // new point k
    for (int i = 0; i < convex_hull.size(); i++)
    {
      int res = find_red_segments(Segment_2(k, convex_hull[i]), convex_hull, convex_seg, 0);

      if (!res)
      {
        red_points.push_back(convex_hull[i]);
      }
    }

    if (!red_points.empty())
    {
      for (int i = 0; i < red_points.size(); i++)
      {
        for (int j = 1; j < red_points.size(); j++)
        {
          Segment_2 temp = edge_exists(red_points[i], red_points[j], convex_seg);
          if (edge_exists(temp[0], temp[1], chain_seg) != error)
          {
            if (std::find(red_segments_final.begin(), red_segments_final.end(), edge_exists(temp[0], temp[1], chain_seg)) == red_segments_final.end())
            {
              int res = find_red_segments(Segment_2(k, midpoint(edge_exists(temp[0], temp[1], chain_seg))), convex_hull, convex_seg, 1);
              if (!res)
              {
                red_segments_final.push_back(edge_exists(temp[0], temp[1], chain_seg));
              }
            }
          }
        }
      }
    }

    red_points.clear();
    Points tempp = curr_points;
    // if no visible edges in convex hull then search in polygon chain
    if (red_segments_final.size() == 0)
    {
      for (int i = 0; i < curr_points.size(); i++)
      {
        int res = find_red_segments(Segment_2(k, curr_points[i]), curr_points, chain_seg, 0);
        if (!res)
        {
          red_points.push_back(curr_points[i]);
        }
      }
      if (!red_points.empty())
      {
        for (int i = 0; i < red_points.size(); i++)
        {
          for (int j = 1; j < red_points.size(); j++)
          {
            Segment_2 temp = edge_exists(red_points[i], red_points[j], chain_seg);
            if (temp != error)
            { // for every combination of red edges check if they match a segment in CH
              if (std::find(red_segments_final.begin(), red_segments_final.end(), temp) == red_segments_final.end())
              {
                int res = find_red_segments(Segment_2(k, midpoint(edge_exists(temp[0], temp[1], chain_seg))), curr_points, chain_seg, 1);
                if (!res)
                {
                  red_segments_final.push_back(temp);
                }
              }
            }
          }
        }
      }
    }
    while (red_segments_final.size() > 0)
    {
      int random = rand() % red_segments_final.size(); // find a random segment

      int i = 0;
      int pos;
      while (i < chain_seg.size())
      {
        if (red_segments_final[random] == chain_seg[i])
        {
          pos = i;
          break;
        }
        i++;
      }
      Segment_2 tempppp = chain_seg[pos];
      Point_2 temppp = Point_2(chain_seg[pos][1]);
      curr_points.push_back(k);
      chain_seg.insert(chain_seg.begin() + pos, Segment_2(chain_seg[pos][0], k));
      // insert in chain at the position that the visible edge was found the new edge connecting with the interior point
      chain_seg.insert(chain_seg.begin() + pos + 1, Segment_2(k, temppp));
      chain_seg.erase(chain_seg.begin() + pos + 2);
      result.erase(result.begin());
      p1.clear();
      keep.clear();
      int e = 0;
      p1.push_back(chain_seg[e][0]);   // initialize the new polygon
      keep.push_back(chain_seg[e][0]); // make a copy of this polygon
      p1.push_back(chain_seg[e][1]);
      keep.push_back(chain_seg[e][1]);
      e++;
      while (e < chain_seg.size() - 1)
      {
        p1.push_back(chain_seg[e][1]);
        keep.push_back(chain_seg[e][1]);
        e++;
      }
      if (option == 3 && k == dont_touch[0] && std::count(chain_seg.begin(), chain_seg.end(), Segment_2(dont_touch[0], dont_touch[1]))!=0)
      {
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.insert(chain_seg.begin() + pos, tempppp); // backtracking deleting the edges i created and bringing back the previous edge
        result.insert(result.begin(), k);
        curr_points.pop_back();                                                                             // pushing back again the interior point
        red_segments_final.erase(std::find(red_segments_final.begin(), red_segments_final.end(), tempppp)); // deleting the visible edges because it doesnt meet the criteria
      }
      if (p1.is_simple() == 0) // if with the new edges the polugon is not simple or the polugon does not surround every point i am backtracking
      {
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.erase(chain_seg.begin() + pos);
        chain_seg.insert(chain_seg.begin() + pos, tempppp); // backtracking deleting the edges i created and bringing back the previous edge
        result.insert(result.begin(), k);
        curr_points.pop_back();                                                                             // pushing back again the interior point
        red_segments_final.erase(std::find(red_segments_final.begin(), red_segments_final.end(), tempppp)); // deleting the visible edges because it doesnt meet the criteria
      }
      else
      {
        break;
      }
    }
  }
  return chain_seg;
}
