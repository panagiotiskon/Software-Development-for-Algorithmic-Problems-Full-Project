#include "Local_search.hpp"

using namespace std::chrono;

static Segment_2 error = Segment_2(Point_2(-1, -1), Point_2(-1, -1));

double local_search(int min_or_max, int cut_off,std::ofstream &outf)
{
  unsigned long int time=0;
  double convex_hull_area;
  double polygon_area;
  double temp_area;
  double DA = 100;
  List la;
  List li;
  List lj;
  List lk;
  Points temp_points;
  Points tp;
  Points convex_hull_points;
  Polygon_2 temp_polygon; // keeps a temp polygon

  std::vector<Segment_2>::iterator seg_it;

  segments tchain;

  std::copy(chain.begin(), chain.end(), std::back_inserter(tchain)); // copy initial chain

  for (int i = 0; i < p.size(); i++)
  {
    temp_polygon.push_back(p[i]); // copy initial polygon
  }

  CGAL::convex_hull_2(p.begin(), p.end(), std::back_inserter(convex_hull_points));

  std::copy(p.begin(), p.end(), std::back_inserter(temp_points));

  std::copy(p.begin(), p.end(), std::back_inserter(tp));

  CGAL::area_2(convex_hull_points.begin(), convex_hull_points.end(), convex_hull_area, K()); // compute convex hull area

  CGAL::area_2(tp.begin(), tp.end(), temp_area, K()); // compute polygons area
  double area_before = temp_area;
  double conv_area = convex_hull_area;
  while (DA >= threshold)
  {
    high_resolution_clock::time_point start = high_resolution_clock::now();
    for (int j = 0; j < tchain.size(); j++)
    {
      int k = 1;
      while (k <= L)
      {
        for (int i = 0; i < tp.size(); i++)
        {
          if (k == 1)
          {
            if (!point_of_segment(tchain[j], tp[i]))
            { // if the selected point is not part of the edge
              if (!find_blue_edge(Segment_2(tp[i], tchain[j][0]), temp_points, tchain, 0) || edge_exists(tp[i], tchain[j][0], tchain) != error)
              { // if the edge from p[i] to chain[j][0] is visible
                if (!find_blue_edge(Segment_2(tp[i], tchain[j][1]), temp_points, tchain, 0) || edge_exists(tp[i], tchain[j][1], tchain) != error)
                { // if the edge is visible
                  if (i == 0)
                  {
                    // if the edge from the two vertices intersects with the polygon
                    if (!find_blue_edge(Segment_2(tp[i + 1], tp[tp.size() - 1]), temp_points, tchain, 0))
                    { // if new edge from the two vertices of the v is not intersecting
                      int a = construct_polygon(i, j, tp[i], tchain[j], tchain, temp_area, 1, tp);
                      if (a != -1)
                      {
                        li.push_back(i);
                        lj.push_back(j);
                        la.push_back(a);
                        lk.push_back(1);
                      }
                    }
                  }
                  if (i == tp.size() - 1)
                  { // if the edge from the two vertices intersects with the polygon
                    if (!find_blue_edge(Segment_2(tp[0], tp[i - 1]), temp_points, tchain, 0))
                    {
                      int a = construct_polygon(i, j, tp[i], tchain[j], tchain, temp_area, 1, tp);
                      if (a != -1)
                      {
                        li.push_back(i);
                        lj.push_back(j);
                        la.push_back(a);
                        lk.push_back(1);
                      }
                    }
                  }
                  else
                  {
                    if (!find_blue_edge(Segment_2(tp[i - 1], tp[i + 1]), temp_points, tchain, 0))
                    {
                      int a = construct_polygon(i, j, tp[i], tchain[j], tchain, temp_area, 1, tp);
                      if (a != -1)
                      {
                        li.push_back(i);
                        lj.push_back(j);
                        la.push_back(a);
                        lk.push_back(1);
                      }
                    }
                  }
                }
              }
            }
          }

          else
          {
            if (i + k <= tp.size())
            {
              bool belongs = false;
              for (int bl = 0; bl < k; bl++)
              { // check if any points of v belongs to chain[j]
                if (point_of_segment(tchain[j], tp[i + bl]) == true)
                {
                  belongs = true;
                  break;
                }
              }

              if (belongs == false)
              { // if none of the points of v belongs to chain[j]
                if (!find_blue_edge(Segment_2(tp[i], tchain[j][0]), temp_points, tchain, 0))
                {
                  if (!find_blue_edge(Segment_2(tp[i + k], tchain[j][1]), temp_points, tchain, 0)) // we want to check the edges from the start point and the final point to the chain edge
                  {                                                                                // if the edge is visible
                    if (i == 0)
                    {
                      if (!find_blue_edge(Segment_2(tp[tp.size() - 1], tp[i + k + 1]), temp_points, tchain, 0)) // check seg (i-1, i+k+1)
                      {
                        int a = construct_polygon(i, j, tp[i], tchain[j], tchain, temp_area, k, tp);
                        if (a != -1)
                        {
                          li.push_back(i);
                          lj.push_back(j);
                          la.push_back(a);
                          lk.push_back(k);
                        }
                      }
                    }

                    else if (i + k < tp.size())
                    {
                      if (!find_blue_edge(Segment_2(tp[i - 1], tp[i + k + 1]), temp_points, tchain, 0))
                      {
                        int a = construct_polygon(i, j, tp[i], tchain[j], tchain, temp_area, k, tp);
                        if (a != -1)
                        {
                          li.push_back(i);
                          lj.push_back(j);
                          la.push_back(a);
                          lk.push_back(k);
                        }
                      }
                    }
                    else if ((i + k) == tp.size())
                    {
                      if (!find_blue_edge(Segment_2(tp[i - 1], tp[0]), temp_points, tchain, 0))
                      {
                        int a = construct_polygon(i, j, tp[i], tchain[j], tchain, temp_area, k, tp);
                        if (a != -1)
                        {
                          li.push_back(i);
                          lj.push_back(j);
                          la.push_back(a);
                          lk.push_back(k);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }

        k++;
      }
    }

    if (la.size() == 0)
    {
      break;
    }
    else
    {

      int a = temp_area;
      int pos = -1;

      if (flag_min_max == 2)
      { // for max polygon
        for (int i = 0; i < la.size(); i++)
        {
          if (a < la[i])
          {
            a = la[i];
            pos = i;
          }
        }
      }
      else if (flag_min_max == 1)
      { // for min polygon
        for (int i = 0; i < la.size(); i++)
        {
          if (a > la[i])
          {
            a = la[i];
            pos = i;
          }
        }
      }

      if (pos == -1)
      {
        break;
      }

      segments chain2;
      chain2 = final_polygon(li[pos], lj[pos], tp[li[pos]], tchain[lj[pos]], tchain, lk[pos]);
      tchain.clear();
      tp.clear();

      std::copy(chain2.begin(), chain2.end(), std::back_inserter(tchain));

      temp_polygon.clear();
      int e = 0;
      temp_polygon.push_back(tchain[e][0]); // initialize the new polygon
      tp.push_back(tchain[e][0]);
      temp_polygon.push_back(tchain[e][1]);
      tp.push_back(tchain[e][1]);
      e++;
      while (e < tchain.size() - 1)
      {
        temp_polygon.push_back(tchain[e][1]);
        tp.push_back(tchain[e][1]);
        e++;
      }
      double re;
      CGAL::area_2(temp_polygon.begin(), temp_polygon.end(), re, K()); // compute polygons area
      DA = abs(re - temp_area);                                        // check the threshold
      temp_area = re;                                                  // assign new area

      la.clear();
      li.clear();
      lj.clear();
      lk.clear();
    }
    high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
    auto duration1 =duration_cast<std::chrono::milliseconds>(stop-start).count();
    time+=duration1;
    if(time>cut_off){
      if(min_or_max==1){
        return 1;
      }
      else{
        return 0;
      }
    }
  }
  p.clear();
  int x;
  p.push_back(tchain[0][0]);
  for (x = 0; x < tchain.size() - 1; x++)
  {
    p.push_back(tchain[x][1]);
  }
  chain.clear();
  create_chain(p.size());
  double ratio = print_result(p.size(), conv_area, area_before, temp_area, outf);
  return ratio;
}

int construct_polygon(int i, int j, Point_2 v, Segment_2 u, segments tchain, int polygon_area, int k, Points tp)
{
    segments temp;
    segments temp_k;
    segments f_temp;
    Polygon_2 t_p;
    Points temp_points;
    int pos = -1;
    std::copy(tchain.begin(), tchain.end(), std::back_inserter(temp)); // copy chain
    std::copy(tp.begin(), tp.end(), std::back_inserter(temp_points));
    if (k == 1)
    {
        temp.erase(temp.begin() + j);
        temp.insert(temp.begin() + j, Segment_2(u[0], v));
        temp.insert(temp.begin() + j + 1, Segment_2(v, u[1]));

        for (int i = 0; i < tchain.size(); i++)
        {
            if (tchain[i][0] == v)
            {
                pos = i;
                break;
            }
        }

        int ind = pos;

        temp.erase(temp.begin() + ind);

        if (i != 0 && i != (tchain.size() - 1))
        {
            temp.erase(temp.begin() + ind);
            temp.insert(temp.begin() + ind, Segment_2(tchain[ind - 1][0], tchain[0][1]));
        }
        else if (i == (tchain.size() - 1))
        {
            temp.erase(temp.begin() + ind);
            temp.insert(temp.begin() + ind, Segment_2(tchain[ind - 1][0], tchain[0][0]));
        }
        else
        {
            temp.erase(temp.begin() + tchain.size() - 1);
            temp.insert(temp.begin(), Segment_2(tchain[tchain.size() - 1][0], tchain[ind][1]));
        }
    }

    else
    {

        if (i != 0 && i != (tchain.size() - 1))
        {
            for (int b = 0; b < tchain.size(); b++) // find the position of i in the chain
            {
                if (tchain[b][0] == v)
                {
                    pos = b;
                    break;
                }
            }

            std::copy(tchain.begin() + pos, tchain.begin() + pos + k - 1, std::back_inserter(temp_k)); // add to temp_k the segments to be rotated [pos, pos+k-1]

            temp_k = change_direction(temp_k); // change the direction

            temp_k.insert(temp_k.begin(), Segment_2(u[0], temp_k[0][0])); // add the first segment  to the changed part

            temp_k.push_back(Segment_2(temp_k[temp_k.size() - 1][1], u[1])); // add the second segment to the changed part

            Segment_2 new_segment = Segment_2(tchain[pos - 1][0], tchain[pos + k - 1][1]); // the segment between the two edges that connect the part to be changed

            if (i > j)
            {

                for (int a = 0; a < j; a++)
                {
                    f_temp.push_back(temp[a]); // add [0, j-1] to f_temp
                }

                std::copy(temp_k.begin(), temp_k.end(), std::back_inserter(f_temp)); // rotated part

                for (int a = j + 1; a < pos - 1; a++) //[j+1, pos-2]
                {
                    f_temp.push_back(temp[a]);
                }

                f_temp.push_back(new_segment); // push new segment

                for (int a = pos + 2; a < temp.size(); a++) //[pos+2, end]
                {
                    f_temp.push_back(temp[a]);
                }
            }

            else if (i < j)
            {

                for (int a = 0; a < pos - 1; a++) //[0,pos-2]
                {
                    f_temp.push_back(temp[a]);
                }

                f_temp.push_back(new_segment); // push new segment

                for (int a = pos + 2; a < j; a++) //[pos+2, j-1]
                {
                    f_temp.push_back(temp[a]);
                }

                std::copy(temp_k.begin(), temp_k.end(), std::back_inserter(f_temp)); // rotated part

                for (int a = j + 1; a < temp.size(); a++) //[j+1, end]
                {
                    f_temp.push_back(temp[a]);
                }
            }
            temp.clear();
            std::copy(f_temp.begin(), f_temp.end(), std::back_inserter(temp)); // finaly add to temp
        }
    }

    int e = 0;

    t_p.push_back(temp[e][0]); // create polygon from temp chain
    t_p.push_back(temp[e][1]);
    e++;
    while (e < temp.size() - 1)
    {
        t_p.push_back(temp[e][1]);
        e++;
    }

    if (t_p.is_simple()) // if the newly created polygon is simple
    {
        double temp_area;
        CGAL::area_2(t_p.begin(), t_p.end(), temp_area, K());

        if (flag_min_max == 2)
        { // for max polygon
            if (abs(temp_area) > abs(polygon_area))
            {
                return abs(temp_area);
            }
        }
        else if (flag_min_max == 1)
        { // for min polygon
            if (abs(temp_area) < abs(polygon_area))
            {
                return abs(temp_area);
            }
        }
    }
    t_p.clear();
    temp.clear();
    temp_k.clear();
    f_temp.clear();
    return -1;
}

segments final_polygon(int i, int j, Point_2 v, Segment_2 u, segments tchain, int k)
{ // constructs new polygon and returns its segments

    segments temp;
    segments temp_k;
    segments f_temp;
    Polygon_2 t_p;

    int pos = -1;
    std::copy(tchain.begin(), tchain.end(), std::back_inserter(temp)); // copy chain
    if (k == 1)
    {
        temp.erase(temp.begin() + j);
        temp.insert(temp.begin() + j, Segment_2(u[0], v));
        temp.insert(temp.begin() + j + 1, Segment_2(v, u[1]));

        for (int i = 0; i < tchain.size(); i++)
        {
            if (tchain[i][0] == v)
            {
                pos = i;
                break;
            }
        }
        int ind = pos;
        temp.erase(temp.begin() + ind);
        if (i != 0 && i != (tchain.size() - 1))
        {
            temp.erase(temp.begin() + ind);
            temp.insert(temp.begin() + ind, Segment_2(tchain[ind - 1][0], tchain[0][1]));
        }
        else if (i == (tchain.size() - 1))
        {
            temp.erase(temp.begin() + ind);
            temp.insert(temp.begin() + ind, Segment_2(tchain[ind - 1][0], tchain[0][0]));
        }
        else
        {
            temp.erase(temp.begin() + tchain.size() - 1);
            temp.insert(temp.begin(), Segment_2(tchain[tchain.size() - 1][0], tchain[ind][1]));
        }

        return temp;
    }
    else
    {

        if (i != 0 && i != (tchain.size() - 1))
        {
            for (int b = 0; b < tchain.size(); b++)
            {
                if (tchain[b][0] == v)
                {
                    pos = b;
                    break;
                }
            }

            std::copy(tchain.begin() + pos, tchain.begin() + pos + k - 1, std::back_inserter(temp_k)); // add to temp_k the segments to be rotated

            temp_k = change_direction(temp_k); // change the direction

            temp_k.insert(temp_k.begin(), Segment_2(u[0], temp_k[0][0])); // add the first segment  to the changed part

            temp_k.push_back(Segment_2(temp_k[temp_k.size() - 1][1], u[1])); // add the second segment to the changed part

            Segment_2 new_segment = Segment_2(tchain[pos - 1][0], tchain[pos + k - 1][1]); // the segment between the two edges that connect the part to be changed

            if (i > j)
            {

                for (int a = 0; a < j; a++)
                {
                    f_temp.push_back(temp[a]); // add [0, j-1] to f_temp
                }

                std::copy(temp_k.begin(), temp_k.end(), std::back_inserter(f_temp)); // rotated part

                for (int a = j + 1; a < pos - 1; a++) //[j+1, pos-2]
                {
                    f_temp.push_back(temp[a]);
                }

                f_temp.push_back(new_segment); // push new segment

                for (int a = pos + 2; a < temp.size(); a++)
                {
                    f_temp.push_back(temp[a]);
                }
            }

            else if (i < j)
            {
                for (int a = 0; a < pos - 1; a++) //[0, pos-2]
                {
                    f_temp.push_back(temp[a]);
                }

                f_temp.push_back(new_segment); // push new segment

                for (int a = pos + 2; a < j; a++) //[pos+2, j+1]
                {
                    f_temp.push_back(temp[a]);
                }

                std::copy(temp_k.begin(), temp_k.end(), std::back_inserter(f_temp)); // rotated part

                for (int a = j + 1; a < temp.size(); a++)
                {
                    f_temp.push_back(temp[a]);
                }
            }

            temp.clear();                                                      // clear temp
            std::copy(f_temp.begin(), f_temp.end(), std::back_inserter(temp)); // finaly add to temp
            return temp;
        }
    }

    return temp;
}
segments change_direction(segments chain_seg)
{ // you give segments and it return them in the opposite direction
    segments temp2;
    for (int j; j < chain_seg.size(); j++)
    {
        temp2.insert(temp2.begin() + j, Segment_2(chain_seg[chain_seg.size() - j - 1][1], chain_seg[chain_seg.size() - j - 1][0]));
    }
    return temp2;
}
int find_blue_edge(Segment_2 k, Points convex_hull, segments tchain, int mid)
{
  int flag;
  int counter = 0;
  for (int i = 0; i < tchain.size(); i++)
  {
    auto result = CGAL::intersection(k, tchain[i]);
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
        if (std::find(convex_hull.begin(), convex_hull.end(), *p) != convex_hull.end() && mid == 0) // check if p is an edge of the chain
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

bool point_of_segment(Segment_2 s, Point_2 p)
{ // returns if a point belongs to edge
  if (s[1] == p)
  {
    return true;
  }
  else if (s[0] == p)
  {
    return true;
  }
  else
    return false;
}
