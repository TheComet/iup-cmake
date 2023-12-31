/** \file
 * \brief Image Analysis
 *
 * See Copyright Notice in im_lib.h
 */


#include <im.h>
#include <im_util.h>
#include <im_math.h>

#include "im_process_counter.h"
#include "im_process_ana.h"
#include "im_process_pnt.h"

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>


#define MAX_COUNT 65536  // maximum number of regions

/* adjust the alias table to be a remap table (final step) */
static void alias_update(imushort* alias_table, int &region_count)
{
  int i, real_count = region_count;

  for (i = 0; i < region_count; i++)
  {
    if (alias_table[i])
    {
      // search for the first alias
      imushort prev = alias_table[i];
      while (alias_table[prev])
        prev = alias_table[prev];

      alias_table[i] = prev;
      real_count--;  // decrement aliases from the region count
    }
  }

  // now all the aliases in the same group point to only one alias
  // transform the alias table into a remap table

  alias_table[0] = 0;
  alias_table[1] = 0;  // border is mapped to background

  int r = 1;
  for (i = 2; i < region_count; i++)
  {
    if (!alias_table[i])
    {
      alias_table[i] = (imushort)r; // only non alias get real values
      r++;
    }
    else
      alias_table[i] = (imushort)(alias_table[alias_table[i]]);
  }

  region_count = real_count-2; // remove the regions (background,border) from the count 
}

/* find the smallest region number to be set as alias. */
static void alias_getmin(imushort* alias_table, imushort region, imushort &min)
{
  while (alias_table[region])
  {
    if (min > alias_table[region])
      min = alias_table[region];

    region = alias_table[region];
  }
}

/* replace all the aliases of a region by its smallest value. */
static void alias_setmin(imushort* alias_table, imushort region, imushort min)
{
  while (alias_table[region])
  {
    imushort next_region = alias_table[region];
    alias_table[region] = min;
    region = next_region;
  }

  if (region != min)
    alias_table[region] = min;
}

/* set a region number to be an alias of another */
static void alias_set(imushort* alias_table, imushort region1, imushort region2)
{
  if (region1 == region2)
    return;

  imushort min = region1<region2? region1: region2;

  alias_getmin(alias_table, region1, min);
  alias_getmin(alias_table, region2, min);

  if (region1 != min && alias_table[region1] != min)
    alias_setmin(alias_table, region1, min);
  if (region2 != min && alias_table[region2] != min)
    alias_setmin(alias_table, region2, min);
}

static int DoAnalyzeFindRegions(int width, int height, imbyte* map, imushort* new_map, int connect, int *region_count, int counter)
{
  int i, j;

  // mark the pixels that touch the border
  // if a region touch the border, is the invalid region 1

  imCounterTotal(counter, height, "Analyzing...");

  imbyte* pmap = map;
  imushort* new_pmap = new_map;
  for (j = 0; j < width; j++)     // first line
  {
    if (pmap[j])
      new_pmap[j] = 1;
  }
  pmap += width;
  new_pmap += width;

  for (i = 1; i < height-1; i++)  // first column
  {
    if (pmap[0])
      new_pmap[0] = 1;

    pmap += width;
    new_pmap += width;
  }

  if (!imCounterInc(counter)) 
    return 0;

  // find and connect the regions

  imbyte* pmap1 = map;         // previous line (line 0)
  imushort* new_pmap1 = new_map; 

  pmap = map + width;          // current line (line 1)
  new_pmap = new_map + width;

  *region_count = 2;  // 0- background, 1-border
  imushort* alias_table = new imushort [MAX_COUNT];
  memset(alias_table, 0, MAX_COUNT*sizeof(imushort)); // aliases are all zero at start (not used)

  for (i = 1; i < height; i++)
  {
    for (j = 1; j < width; j++)
    {
      int has_j1 = j < width-1? 1: 0;
      if (pmap[j])
      {
        if (pmap[j-1] || pmap1[j] || 
            (connect == 8 && (pmap1[j-1] || (has_j1&&pmap1[j+1])))) // 4 or 8 connected to the previous neighbors
        {
          imushort region = 0;
          if (i == height-1 || j == width-1)
          {
            region = new_pmap[j] = 1;
          }

          if (pmap[j-1])
          {
            if (!region)
              region = new_pmap[j-1];  // horizontal neighbor  -00
            else                       //                      X1
            {
              // this is a right border pixel that connects to an horizontal neighbor

              // this pixel can connect two different regions
              alias_set(alias_table, region, new_pmap[j-1]);
            }
          }

          if (pmap1[j])    // vertical neighbor
          {
            if (!region)
              region = new_pmap1[j];  // isolated vertical neighbor  -X-
            else                      //                             01
            {
              // an horizontal neighbor connects to a vertical neighbor  -X-
              //                                                         X1

              // this pixel can connect two different regions
              alias_set(alias_table, region, new_pmap1[j]);
            }
          }
          else if (region && connect==8 && (has_j1&&pmap1[j+1]))
          {
            // an horizontal neighbor connects to a right corner neighbor   00X
            //                                                              X1

            // this pixel can connect two different regions
            alias_set(alias_table, region, new_pmap1[j+1]);
          }

          if (connect == 8 && (pmap1[j-1] || (has_j1&&pmap1[j+1])) && !region) // isolated corner
          {
            // a left corner neighbor or a right corner neighbor  X0X
            //                                                    01

            if (pmap1[j-1])  // left corner
              region = new_pmap1[j-1];

            if (pmap1[j+1])  // right corner
            {
              if (!region) // isolated right corner
                region = new_pmap1[j+1];
              else
              {
                // this pixel can connect two different regions
                alias_set(alias_table, new_pmap1[j-1], new_pmap1[j+1]);
              }
            }
          }

          new_pmap[j] = region;
        }
        else
        {
          // this pixel touches no pixels

          if (i == height-1 || j == width-1)
            new_pmap[j] = 1;
          else
          {
            // create a new region  000
            //                      01
            new_pmap[j] = (imushort)*region_count;
            (*region_count)++;

            if (*region_count > MAX_COUNT)
            {
              delete [] alias_table;
              return 0;
            }
          }
        }
      }
    }

    pmap1 = pmap;
    new_pmap1 = new_pmap;
    pmap += width;
    new_pmap += width;

    if (!imCounterInc(counter))
    {
      delete[] alias_table;
      return 0;
    }
  }

  // now all pixels are marked, 
  // but some marks are aliases to others

  // adjust the alias table to be a remap table
  // and return the real region count
  alias_update(alias_table, *region_count);

  int count = width*height;
  for (i = 0; i < count; i++)
  {
    new_map[i] = alias_table[new_map[i]];
  }

  delete [] alias_table;

  return 1;
}

static int DoAnalyzeFindRegionsBorder(int width, int height, imbyte* map, imushort* new_map, int connect, int *region_count, int counter)
{
  int i, j;

  imbyte* pmap1 = map - width;         // previous line (line -1 = invalid)
  imushort* new_pmap1 = new_map - width; 

  imbyte* pmap = map;                  // current line (line 0)
  imushort* new_pmap = new_map;

  imCounterTotal(counter, height, "Analyzing...");

  *region_count = 2;  // still consider: 0- background, 1-border
  imushort* alias_table = new imushort [MAX_COUNT];
  memset(alias_table, 0, MAX_COUNT*sizeof(imushort)); // aliases are all zero at start (not used)

  for (i = 0; i < height; i++)
  {
    for (j = 0; j < width; j++)
    {
      if (pmap[j])
      {
        int b01 = j > 0? 1: 0; // valid for pmap[j-1]
        int b10 = i > 0? 1: 0; // valid for pmap1[j]
        int b11 = i > 0 && j > 0? 1: 0; // valid for pmap1[j-1]
        int b12 = i > 0 && j < width-1? 1: 0; // valid for pmap1[j+1]

        if ((b01&&pmap[j-1]) || (b10&&pmap1[j]) || 
            (connect == 8 && ((b11&&pmap1[j-1]) || (b12&&pmap1[j+1])))) // 4 or 8 connected to the previous neighbors
        {
          imushort region = 0;

          if (b01&&pmap[j-1])
          {
            if (!region)
              region = new_pmap[j-1];  // horizontal neighbor  -00
            else                       //                      X1
            {
              // this is a right border pixel that connects to an horizontal neighbor

              // this pixel can connect two different regions
              alias_set(alias_table, region, new_pmap[j-1]);
            }
          }

          if (b10&&pmap1[j])    // vertical neighbor
          {
            if (!region)
              region = new_pmap1[j];  // isolated vertical neighbor  -X-
            else                      //                             01
            {
              // an horizontal neighbor connects to a vertical neighbor  -X-
              //                                                         X1

              // this pixel can connect two different regions
              alias_set(alias_table, region, new_pmap1[j]);
            }
          }
          else if (region && connect == 8 && (b12&&pmap1[j+1]))
          {
            // an horizontal neighbor connects to a right corner neighbor   00X
            //                                                              X1

            // this pixel can connect two different regions
            alias_set(alias_table, region, new_pmap1[j+1]);
          }

          if (connect == 8 && ((b11&&pmap1[j-1]) || (b12&&pmap1[j+1])) && !region) // isolated corner
          {
            // a left corner neighbor or a right corner neighbor  X0X
            //                                                    01

            if (b11&&pmap1[j-1])  // left corner
              region = new_pmap1[j-1];

            if (b12&&pmap1[j+1])  // right corner
            {
              if (!region) // isolated right corner
                region = new_pmap1[j+1];
              else
              {
                // this pixel can connect two different regions
                alias_set(alias_table, new_pmap1[j-1], new_pmap1[j+1]);
              }
            }
          }

          new_pmap[j] = region;
        }
        else
        {
          // this pixel touches no pixels

          // create a new region  000
          //                      01
          new_pmap[j] = (imushort)*region_count;
          (*region_count)++;

          if (*region_count > MAX_COUNT)
          {
            delete [] alias_table;
            return 0;
          }
        }
      }
    }

    pmap1 = pmap;
    new_pmap1 = new_pmap;
    pmap += width;
    new_pmap += width;

    if (!imCounterInc(counter))
    {
      delete[] alias_table;
      return 0;
    }
  }

  // now all pixels are marked, 
  // but some marks are aliases to others

  // adjust the alias table to be a remap table
  // and return the real region count
  alias_update(alias_table, *region_count);

  int count = width*height;
  for (i = 0; i < count; i++)
  {
    new_map[i] = alias_table[new_map[i]];
  }

  delete [] alias_table;

  return 1;
}

int imAnalyzeFindRegions(const imImage* src_image, imImage* dst_image, int connect, int touch_border, int *region_count)
{
  int ret = 0;
  int counter = imCounterBegin("FindRegions");

  imImageSetAttribute(dst_image, "REGION_CONNECT", IM_BYTE, 1, connect == 4 ? "4" : "8");
  if (touch_border)
    ret = DoAnalyzeFindRegionsBorder(src_image->width, src_image->height, (imbyte*)src_image->data[0], (imushort*)dst_image->data[0], connect, region_count, counter);
  else
    ret = DoAnalyzeFindRegions(src_image->width, src_image->height, (imbyte*)src_image->data[0], (imushort*)dst_image->data[0], connect, region_count, counter);

  imProcessCounterEnd(counter);
  return ret;
}

int imAnalyzeMeasureArea(const imImage* image, int* data_area, int region_count)
{
  imushort* img_data = (imushort*)image->data[0];

  int counter = imProcessCounterBegin("MeasureArea");
  imCounterTotal(counter, image->height, "Analyzing...");

  memset(data_area, 0, region_count*sizeof(int));

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(image->count))
#endif
  for (int i = 0; i < image->count; i++)
  {
    if (i % image->width == 0)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_BEGIN_PROCESSING;

    if (img_data[i])
    {
      int index = img_data[i] - 1;
#ifdef _OPENMP
#pragma omp atomic
#endif
      data_area[index]++;
    }

    if (i % image->width == 0)
    {
      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_END_PROCESSING;
  }

  imProcessCounterEnd(counter);
  return processing;
}

int imAnalyzeMeasureCentroid(const imImage* image, const int* data_area, int region_count, double* data_cx, double* data_cy)
{
  imushort* img_data = (imushort*)image->data[0];
  int* local_data_area = 0;
  int ret;

  int counter = imProcessCounterBegin("MeasureCentroid");
  imCounterTotal(counter, image->height, "Analyzing...");

  if (!data_area)
  {
    local_data_area = (int*)malloc(region_count*sizeof(int));
    ret = imAnalyzeMeasureArea(image, local_data_area, region_count);
    data_area = (const int*)local_data_area;

    if (!ret)
    {
      free(local_data_area);
      imProcessCounterEnd(counter);
      return 0;
    }
  }

  if (data_cx) memset(data_cx, 0, region_count*sizeof(double));
  if (data_cy) memset(data_cy, 0, region_count*sizeof(double));

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(image->height))
#endif
  for (int y = 0; y < image->height; y++) 
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int offset = y*image->width;

    for (int x = 0; x < image->width; x++)
    {
      int region_index = img_data[offset+x];
      if (region_index)
      {
        int ri = region_index-1;
        if (data_cx) 
        {
#ifdef _OPENMP
#pragma omp atomic
#endif
          data_cx[ri] += (double)x;
        }
        if (data_cy) 
        {
#ifdef _OPENMP
#pragma omp atomic
#endif
          data_cy[ri] += (double)y;
        }
      }
    }

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  for (int i = 0; i < region_count; i++) 
  {
    if (data_cx) data_cx[i] /= (double)data_area[i];
    if (data_cy) data_cy[i] /= (double)data_area[i];
  }

  if (local_data_area)
    free(local_data_area);

  imProcessCounterEnd(counter);
  return processing;
}

static inline double ipow(double x, int j)
{
	double r = 1.0;
	for (int i = 0; i < j; i++) 
    r *= x;
	return r;
}

static int iCalcMoment(double* cm, int px, int py, const imImage* image, const double* cx, const double* cy, int region_count, int counter)
{
  imushort* img_data = (imushort*)image->data[0];

  memset(cm, 0, region_count*sizeof(double));

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(image->height))
#endif
  for (int y = 0; y < image->height; y++) 
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int offset = y*image->width;

    for (int x = 0; x < image->width; x++)
    {
      int region_index = img_data[offset+x];
      if (region_index)
      {
        int ri = region_index-1;

        if (px == 0)
        {
          double inc = ipow(y-cy[ri],py);
#ifdef _OPENMP
#pragma omp atomic
#endif
          cm[ri] += inc;
        }
        else if (py == 0)
        {
          double inc = ipow(x-cx[ri],px);
#ifdef _OPENMP
#pragma omp atomic
#endif
          cm[ri] += inc;
        }
        else
        {
          double inc = ipow(x-cx[ri],px)*ipow(y-cy[ri],py);
#ifdef _OPENMP
#pragma omp atomic
#endif
          cm[ri] += inc;
        }
      }
    }

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  return processing;
}

template<class T>
static inline int IsPerimeterPoint(T* map, int width, int height, int x, int y)
{
  // map here points to the start of the line, even if its an invalid line.

  // if outside the image, then is not a perimeter line.
  if (x == -1 || x == width ||
      y == -1 || y == height)
    return 0;

  T v = map[x]; // here v is image(x,y)
  if (!v)
    return 0;

  // if touches the border, then is a perimeter line.
  if (x == 0 || x == width-1 ||
      y == 0 || y == height-1)
    return 1;

  // if has 4 connected neighbors, then is a perimeter line.
  if (map[width+x] != v ||
      map[x+1] != v ||
      map[x-1] != v ||
      map[-width+x] != v)
    return 1;

  return 0;
}

int imAnalyzeMeasurePrincipalAxis(const imImage* image, const int* data_area, const double* data_cx, const double* data_cy,
                                   const int region_count, double* major_slope, double* major_length, 
                                                           double* minor_slope, double* minor_length)
{
  int *local_data_area = 0;
  double *local_data_cx = 0, 
         *local_data_cy = 0;
  int ret = 0;

  int counter = imProcessCounterBegin("PrincipalAxis");

  if (!data_area)
  {
    local_data_area = (int*)malloc(region_count*sizeof(int));
    ret = imAnalyzeMeasureArea(image, local_data_area, region_count);
    data_area = (const int*)local_data_area;

    if (!ret)
    {
      free(local_data_area);
      imProcessCounterEnd(counter);
      return 0;
    }
  }

  if (!data_cx || !data_cy)
  {
    if (!data_cx)
    {
      local_data_cx = (double*)malloc(region_count*sizeof(double));
      data_cx = (const double*)local_data_cx;
    }

    if (!data_cy)
    {
      local_data_cy = (double*)malloc(region_count*sizeof(double));
      data_cy = (const double*)local_data_cy;
    }

    ret = 1;
    if (local_data_cx && local_data_cy)
      ret = imAnalyzeMeasureCentroid(image, data_area, region_count, local_data_cx, local_data_cy);
    else if (local_data_cx)
      ret = imAnalyzeMeasureCentroid(image, data_area, region_count, local_data_cx, NULL);
    else if (local_data_cy)
      ret = imAnalyzeMeasureCentroid(image, data_area, region_count, NULL, local_data_cy);

    if (!ret)
    {
      if (local_data_cx) 
        free(local_data_cx);
      if (local_data_cy) 
        free(local_data_cy);
      imProcessCounterEnd(counter);
      return 0;
    }
  }

  imCounterTotal(counter, 4 * image->height + 2 * region_count, "Analyzing...");

  // additional moments
  double* cm20 = (double*)malloc(region_count*sizeof(double));
  double* cm02 = (double*)malloc(region_count*sizeof(double));
  double* cm11 = (double*)malloc(region_count*sizeof(double));
  
  ret = iCalcMoment(cm20, 2, 0, image, data_cx, data_cy, region_count, counter);
  if (!ret)
  {
    if (local_data_area) 
      free(local_data_area);
    if (local_data_cx) 
      free(local_data_cx);
    if (local_data_cy) 
      free(local_data_cy);
    free(cm20); 
    free(cm02); 
    free(cm11);
    imProcessCounterEnd(counter);
    return 0;
  }
  ret = iCalcMoment(cm02, 0, 2, image, data_cx, data_cy, region_count, counter);
  if (!ret)
  {
    if (local_data_area) 
      free(local_data_area);
    if (local_data_cx) 
      free(local_data_cx);
    if (local_data_cy) 
      free(local_data_cy);
    free(cm20); 
    free(cm02); 
    free(cm11);
    imProcessCounterEnd(counter);
    return 0;
  }
  ret = iCalcMoment(cm11, 1, 1, image, data_cx, data_cy, region_count, counter);
  if (!ret)
  {
    if (local_data_area) 
      free(local_data_area);
    if (local_data_cx) 
      free(local_data_cx);
    if (local_data_cy) 
      free(local_data_cy);
    free(cm20); 
    free(cm02); 
    free(cm11);
    imProcessCounterEnd(counter);
    return 0;
  }

  double *local_major_slope = 0, *local_minor_slope = 0;
  if (!major_slope)
  {
    local_major_slope = (double*)malloc(region_count*sizeof(double));
    major_slope = local_major_slope;
  }
  if (!minor_slope)
  {
    local_minor_slope = (double*)malloc(region_count*sizeof(double));
    minor_slope = local_minor_slope;
  }

#define RAD2DEG  57.295779513

  // We are going to find 2 axis parameters.
  // Axis 1 are located in quadrants 1-3
  // Axis 2 are located in quadrants 2-4

  // Quadrants
  //    2 | 1
  //    -----
  //    3 | 4

  // line coefficients for lines that belongs to axis 1 and 2
  double* A1 = (double*)malloc(region_count*sizeof(double));
  double* A2 = (double*)malloc(region_count*sizeof(double));
  double* C1 = (double*)malloc(region_count*sizeof(double));
  double* C2 = (double*)malloc(region_count*sizeof(double));

  double *slope1 = major_slope; // Use major_slope as a storage place, 
  double *slope2 = minor_slope; // and create an alias to make code clear.

  for (int i = 0; i < region_count; i++) 
  {
    if (cm11[i] == 0)
    {
      slope1[i] = 0;
      slope2[i] = 90;

      // These should not be used
      A1[i] = 0; 
      A2[i] = 0;  // infinite
      C1[i] = 0;  // data_cy[i]
      C2[i] = 0;  
    }
    else
    {
      double b = (cm20[i] - cm02[i])/cm11[i];
      double delta = sqrt(b*b + 4.0);
      double r1 = (-b-delta)/2.0;
      double r2 = (-b+delta)/2.0;
      double a1 = (double)(atan(r1)*RAD2DEG + 90);  // to avoid negative results
      double a2 = (double)(atan(r2)*RAD2DEG + 90);

      if (a1 == 180) a1 = 0;
      if (a2 == 180) a2 = 0;

      if (a1 < 90)             // a1 is quadrants q1-q3
      {                        
        slope1[i] = a1;   
        slope2[i] = a2;   
        A1[i] = (double)r1;
        A2[i] = (double)r2;
      }
      else                     // a2 is quadrants q1-q3
      {
        slope1[i] = a2;
        slope2[i] = a1;
        A1[i] = (double)r2;
        A2[i] = (double)r1;
      }

      C1[i] = data_cy[i] - A1[i] * data_cx[i];
      C2[i] = data_cy[i] - A2[i] * data_cx[i];
    }

    if (!imCounterInc(counter))
    {
      if (local_major_slope) 
        free(local_major_slope);
      if (local_minor_slope) 
        free(local_minor_slope);
      free(A1);
      free(A2);
      free(C1);
      free(C2);

      if (local_data_area) 
        free(local_data_area);
      if (local_data_cx) 
        free(local_data_cx);
      if (local_data_cy) 
        free(local_data_cy);
      free(cm20); 
      free(cm02); 
      free(cm11);

      imProcessCounterEnd(counter);
      return 0;
    }
  }

  // moments are not necessary anymore
  free(cm20); 
  free(cm02); 
  free(cm11);
  cm20 = 0; cm02 = 0; cm11 = 0;

  // maximum distance from a point in the perimeter to an axis in each side of the axis
  // D1 is distance to axis 1, a and b are sides
  double* D1a = (double*)malloc(region_count*sizeof(double));
  double* D1b = (double*)malloc(region_count*sizeof(double));
  double* D2a = (double*)malloc(region_count*sizeof(double));
  double* D2b = (double*)malloc(region_count*sizeof(double));
  memset(D1a, 0, region_count*sizeof(double));
  memset(D1b, 0, region_count*sizeof(double));
  memset(D2a, 0, region_count*sizeof(double));
  memset(D2b, 0, region_count*sizeof(double));

  imushort* img_data = (imushort*)image->data[0];
  int width = image->width;
  int height = image->height;
  for (int y = 0; y < height; y++) 
  {
    int offset = y*width;

    for (int x = 0; x < width; x++)
    {
      if (IsPerimeterPoint(img_data+offset, width, height, x, y))
      {
        int index = img_data[offset+x] - 1;

        double d1, d2;
        if (slope2[index] == 90)
        {
          d2 = y - data_cy[index];   // I checked this many times, looks odd but it is correct.
          d1 = x - data_cx[index];
        }
        else
        {
          d1 = A1[index]*x - y + C1[index];
          d2 = A2[index]*x - y + C2[index];
        }

        if (d1 < 0)
        {
          d1 = (double)fabs(d1);
          if (d1 > D1a[index])         
            D1a[index] = d1;
        }
        else
        {
          if (d1 > D1b[index])
            D1b[index] = d1;
        }

        if (d2 < 0)
        {
          d2 = (double)fabs(d2);
          if (d2 > D2a[index])         
            D2a[index] = d2;
        }
        else
        {
          if (d2 > D2b[index])
            D2b[index] = d2;
        }
      }
    }

    if (!imCounterInc(counter))
    {
      ret = 0;
      break;
    }
  }

  for (int i = 0; i < region_count && ret != 0; i++) 
  {
    double AB1 = (double)sqrt(A1[i]*A1[i] + 1);
    double AB2 = (double)sqrt(A2[i]*A2[i] + 1);

    double D1 = (D1a[i] + D1b[i]) / AB1; 
    double D2 = (D2a[i] + D2b[i]) / AB2;

    if (D1 < D2) // Major Axis in 2-4 quadrants
    {
      // now remember that we did an alias before
      // slope1 -> major_slope
      // slope2 -> minor_slope

      double tmp = major_slope[i];
      major_slope[i] = minor_slope[i];
      minor_slope[i] = tmp;

      if (minor_length) minor_length[i] = D1;
      if (major_length) major_length[i] = D2;
    }
    else
    {
      if (minor_length) minor_length[i] = D2;
      if (major_length) major_length[i] = D1;
    }

    if (!imCounterInc(counter))
    {
      ret = 0;
      break;
    }
  }

  free(D1b);
  free(D2b);
  free(D1a);
  free(D2a);

  if (local_major_slope) free(local_major_slope);
  if (local_minor_slope) free(local_minor_slope);
  free(A1);  
  free(A2);  
  free(C1);  
  free(C2);

  if (local_data_area) 
    free(local_data_area);
  if (local_data_cx) 
    free(local_data_cx);
  if (local_data_cy) 
    free(local_data_cy);

  imProcessCounterEnd(counter);
  return ret;
}

int imAnalyzeMeasureHoles(const imImage* image, int connect, int region_count, int* count_data, int* area_data, double* perim_data)
{
  int counter = imProcessCounterBegin("MeasureHoles");
  int ret;

  int i;
  imImage *inv_image = imImageCreate(image->width, image->height, IM_BINARY, IM_BYTE);
  imbyte* inv_data = (imbyte*)inv_image->data[0];
  imushort* img_data = (imushort*)image->data[0];

  if (!inv_image)
  {
    imProcessCounterEnd(counter);
    return 0;
  }

  memset(count_data, 0, region_count*sizeof(int));
  memset(area_data, 0, region_count*sizeof(int));
  memset(perim_data, 0, region_count*sizeof(double));

  // finds the holes in the inverted image
  for (i = 0; i < image->count; i++)
  {
    if (img_data[i])
      inv_data[i] = 0;
    else
      inv_data[i] = 1;
  }

  imImage *holes_image = imImageClone(image);
  if (!holes_image)
  {
    imImageDestroy(inv_image);
    imProcessCounterEnd(counter);
    return 0;
  }

  int holes_count = 0;
  if (!imAnalyzeFindRegions(inv_image, holes_image, connect, 0, &holes_count))
  {
    imImageDestroy(inv_image);
    imImageDestroy(holes_image);
    imProcessCounterEnd(counter);
    return 0;
  }

  imImageDestroy(inv_image);

  if (!holes_count)
  {
    imImageDestroy(holes_image);
    imProcessCounterEnd(counter);
    return 1;
  }

  // measure each holes area
  int* holes_area = (int*)malloc(holes_count*sizeof(int));
  ret = imAnalyzeMeasureArea(holes_image, holes_area, holes_count);
  if (!ret)
  {
    free(holes_area);
    imImageDestroy(holes_image);
    imProcessCounterEnd(counter);
    return 0;
  }

  double* holes_perim = 0;
  if (perim_data) 
  {
    holes_perim = (double*)malloc(holes_count*sizeof(int));
    ret = imAnalyzeMeasurePerimeter(holes_image, holes_perim, holes_count);

    if (!ret)
    {
      free(holes_perim);
      free(holes_area);
      imImageDestroy(holes_image);
      imProcessCounterEnd(counter);
      return 0;
    }
  }

  ret = 1;
  imCounterTotal(counter, image->height - 2, "Analyzing...");

  imushort* holes_data = (imushort*)holes_image->data[0];
  img_data = (imushort*)image->data[0];

  // holes do not touch the border
  for (int y = 1; y < image->height-1; y++) 
  {
    int offset_up = (y+1)*image->width;
    int offset = y*image->width;
    int offset_dw = (y-1)*image->width;

    for (int x = 1; x < image->width-1; x++)
    {
      int hole_index = holes_data[offset+x];

      if (hole_index && holes_area[hole_index-1]) // a hole not yet used
      {
        // if the hole has not been used, 
        // it is the first time we encounter a pixel of this hole.
        // then it is a pixel from the hole border.
        // now find which region this hole is inside.
        // a 4 connected neighbor is necessarily a valid region or 0.

        int region_index = 0;
        if (img_data[offset_up + x]) region_index = img_data[offset_up + x];
        else if (img_data[offset + x+1]) region_index = img_data[offset + x+1];
        else if (img_data[offset + x-1]) region_index = img_data[offset + x-1]; 
        else if (img_data[offset_dw+x]) region_index = img_data[offset_dw+x];

        if (region_index) 
        {
          if (count_data) 
            count_data[region_index-1]++;
          if (area_data) 
            area_data[region_index-1] += holes_area[hole_index-1];
          if (perim_data) 
            perim_data[region_index-1] += holes_perim[hole_index-1];

          holes_area[hole_index-1] = 0; // mark hole as used
        }
      }
    }

    if (!imCounterInc(counter))
    {
      ret = 0;
      break;
    }
  }

  if (holes_perim) 
    free(holes_perim);
  free(holes_area);
  imImageDestroy(holes_image);

  imProcessCounterEnd(counter);
  return ret;
}

template<class T>
static int DoPerimeterLine(T* map, T* new_map, int width, int height, int counter)
{
  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for (int y = 0; y < height; y++) 
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int offset = y*width;

    for (int x = 0; x < width; x++)
    {
      if (IsPerimeterPoint(map+offset, width, height, x, y))
        new_map[offset+x] = map[offset+x];
      else
        new_map[offset+x] = 0;
    }

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  return processing;
}

int imProcessPerimeterLine(const imImage* src_image, imImage* dst_image)
{
  int ret = 0;
  int counter = imProcessCounterBegin("PerimeterLine");
  imCounterTotal(counter, src_image->height, "Processing...");

  switch (src_image->data_type)
  {
  case IM_BYTE:
    ret = DoPerimeterLine((imbyte*)src_image->data[0], (imbyte*)dst_image->data[0], src_image->width, src_image->height, counter);
    break;                                                                                
  case IM_SHORT:
    ret = DoPerimeterLine((short*)src_image->data[0], (short*)dst_image->data[0], src_image->width, src_image->height, counter);
    break;                                                                                
  case IM_USHORT:
    ret = DoPerimeterLine((imushort*)src_image->data[0], (imushort*)dst_image->data[0], src_image->width, src_image->height, counter);
    break;                                                                                
  case IM_INT:                                                                           
    ret = DoPerimeterLine((int*)src_image->data[0], (int*)dst_image->data[0], src_image->width, src_image->height, counter);
    break;                                                                                
  }

  imProcessCounterEnd(counter);
  return ret;
}

/* Perimeter Templates idea based in
   Parker, Practical Computer Vision Using C

For 1.414 (sqrt(2)/2 + sqrt(2)/2) [1]:
     1 0 0   0 0 1   1 0 0   0 0 1   0 0 0   1 0 1
     0 x 0   0 x 0   0 x 0   0 x 0   0 x 0   0 x 0
     0 0 1   1 0 0   1 0 0   0 0 1   1 0 1   0 0 0
      129      36     132      33      5      160

For 1.207 (sqrt(2)/2 + 1.0/2) [2]:
     0 0 0   0 0 1   0 1 0   0 1 0   1 0 0   0 0 1   0 0 0   1 0 0
     1 x 0   1 x 0   0 x 0   0 x 0   0 x 0   0 x 0   0 x 1   0 x 1
     0 0 1   0 0 0   1 0 0   0 0 1   0 1 0   0 1 0   1 0 0   0 0 0
       17      48      68      65     130      34       12     136

     0 0 0   1 0 0   1 1 0   0 1 1   0 0 1   0 0 0   0 0 0   0 0 0
     1 x 0   1 x 0   0 x 0   0 x 0   0 x 1   0 x 1   0 x 0   0 x 0
     1 0 0   0 0 0   0 0 0   0 0 0   0 0 0   0 0 1   0 1 1   1 1 0
       20     144     192      96      40      9       3       6

For 1.0 (1.0/2 + 1.0/2) [0]:
     0 0 0   0 1 0   0 0 0   0 0 0   0 1 0   0 1 0
     1 x 1   0 x 0   1 x 0   0 x 1   1 x 0   0 x 1
     0 0 0   0 1 0   0 1 0   0 1 0   0 0 0   0 0 0
       24      66      18      10      80      72

For 0.707 (sqrt(2)/2) [3]:
     1 0 0   0 0 1   0 0 0   0 0 0
     0 x 0   0 x 0   0 x 0   0 x 0         (For Line Length)
     0 0 0   0 0 0   0 0 1   1 0 0
      128      32      1       4

For 0.5 (1.0/2) [4]:
     0 1 0   0 0 0   0 0 0   0 0 0
     0 x 0   0 x 1   0 x 0   1 x 0         (For Line Length)
     0 0 0   0 0 0   0 1 0   0 0 0
       64      8       2      16

*/
static void iInitPerimTemplate(imbyte *templ, double *v)
{
  memset(templ, 0, 256);

  templ[129] = 1;
  templ[36]  = 1;
  templ[132] = 1;
  templ[33]  = 1;
  templ[5]   = 1;
  templ[160] = 1;

  templ[17]  = 2;
  templ[48]  = 2;
  templ[68]  = 2;
  templ[65]  = 2;
  templ[130] = 2;
  templ[34]  = 2;
  templ[12]  = 2;
  templ[136] = 2;
  templ[20]  = 2;
  templ[144] = 2;
  templ[192] = 2;
  templ[96]  = 2;
  templ[40]  = 2;
  templ[9]   = 2;
  templ[3]   = 2;
  templ[6]   = 2;

  templ[24] = 0;
  templ[66] = 0;
  templ[18] = 0;
  templ[10] = 0;
  templ[80] = 0;
  templ[72] = 0;

  templ[128] = 3;
  templ[32]  = 3;
  templ[1]   = 3;
  templ[4]   = 3;

  templ[64] = 4;
  templ[8]  = 4;
  templ[2]  = 4;
  templ[16] = 4;

const double DT_SQRT2   = 1.414213562373;
const double DT_SQRT2D2 = 0.707106781187;

  v[1] = DT_SQRT2;   
  v[2] = DT_SQRT2D2 + 0.5;   
  v[0] = 1.0;
  v[3] = DT_SQRT2D2;
  v[4] = 0.5;
}

int imAnalyzeMeasurePerimeter(const imImage* image, double* perim_data, int region_count)
{
  static imbyte templ[256];
  static double vt[5];
  static int first = 1;
  if (first)
  {
    iInitPerimTemplate(templ, vt);
    first = 0;
  }

  int counter = imProcessCounterBegin("MeasurePerimeter");
  imCounterTotal(counter, image->height, "Analyzing...");

  imushort* map = (imushort*)image->data[0];

  memset(perim_data, 0, region_count*sizeof(double));

  int width = image->width;
  int height = image->height;

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for (int y = 0; y < height; y++) 
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int offset = y*image->width;

    for (int x = 0; x < width; x++)
    {
      if (IsPerimeterPoint(map+offset, width, height, x, y))
      {
        int T = 0;

        // check the 8 neighbors if they belong to the perimeter
        if (IsPerimeterPoint(map+offset+width, width, height, x-1, y+1))
          T |= 0x01;
        if (IsPerimeterPoint(map+offset+width, width, height, x, y+1))
          T |= 0x02;
        if (IsPerimeterPoint(map+offset+width, width, height, x+1, y+1))
          T |= 0x04;

        if (IsPerimeterPoint(map+offset, width, height, x-1, y))
          T |= 0x08;
        if (IsPerimeterPoint(map+offset, width, height, x+1, y))
          T |= 0x10;

        if (IsPerimeterPoint(map+offset-width, width, height, x-1, y-1))
          T |= 0x20;
        if (IsPerimeterPoint(map+offset-width, width, height, x, y-1))
          T |= 0x40;
        if (IsPerimeterPoint(map+offset-width, width, height, x+1, y-1))
          T |= 0x80;

        if (T)
        {
          int index = map[offset+x] - 1;
          double inc = vt[templ[T]];
#ifdef _OPENMP
#pragma omp atomic
#endif
          perim_data[index] += inc;
        }
      }
    }

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  imProcessCounterEnd(counter);
  return processing;
}

/* Perimeter Area Templates

For "1.0" (0):

     1 1 1
     1 x 1
     1 1 1
      255

For "0.75" (1):

     1 1 1   1 1 1   0 1 1   1 1 0   1 1 1   1 1 1   1 1 1   1 0 1
     1 x 1   1 x 1   1 x 1   1 x 1   0 x 1   1 x 0   1 x 1   1 x 1
     0 1 1   1 1 0   1 1 1   1 1 1   1 1 1   1 1 1   1 0 1   1 1 1
      251     254     127     223     239     247     253     191

For "0.625" (2):

     1 1 1   0 0 1   0 1 1   1 1 0   1 1 1   1 1 1   1 1 1   1 0 0
     1 x 1   1 x 1   0 x 1   1 x 0   0 x 1   1 x 0   1 x 1   1 x 1
     0 0 1   1 1 1   1 1 1   1 1 1   0 1 1   1 1 0   1 0 0   1 1 1
      249     63      111     215     235     246     252     159

For "0.5" (3):

     0 0 0   0 1 1   1 1 1   1 1 0   1 1 1   0 0 1   1 0 0   1 1 1  
     1 x 1   0 x 1   1 x 1   1 x 0   0 x 1   0 x 1   1 x 0   1 x 0  
     1 1 1   0 1 1   0 0 0   1 1 0   0 0 1   1 1 1   1 1 1   1 0 0  
      31      107     248     214     233     47      151     244

For "0.375" (4):

     0 0 0   1 1 1   1 1 0   0 1 1   1 0 0   0 0 1   0 0 0   1 1 1
     1 x 0   1 x 0   1 x 0   0 x 1   1 x 0   0 x 1   0 x 1   0 x 1
     1 1 1   0 0 0   1 0 0   0 0 1   1 1 0   0 1 1   1 1 1   0 0 0
      23      240     212     105     150     43      15      232

For "0.25" (5):

     0 0 0   0 0 0   1 1 0   0 1 1   1 0 0   0 0 1   0 0 0   1 1 1
     1 x 0   0 x 1   1 x 0   0 x 1   1 x 0   0 x 1   0 x 0   0 x 0
     1 1 0   0 1 1   0 0 0   0 0 0   1 0 0   0 0 1   1 1 1   0 0 0
      22      11      208     104     148     41       7      224

For "0.125" (6):

     0 0 0   0 0 0   1 1 0   0 0 1   1 0 0   0 0 0   0 0 0   0 1 1
     1 x 0   0 x 0   0 x 0   0 x 1   1 x 0   0 x 1   0 x 0   0 x 0
     1 0 0   0 1 1   0 0 0   0 0 0   0 0 0   0 0 1   1 1 0   0 0 0
      20       3      192      40     144      9       6       96

*/
static void iInitPerimAreaTemplate(imbyte *templ, double *v)
{
  memset(templ, 0, 256);

  templ[255] = 0;

  templ[251] = 1;
  templ[254] = 1;
  templ[127] = 1;
  templ[223] = 1;
  templ[239] = 1;
  templ[247] = 1;
  templ[253] = 1;
  templ[191] = 1;
        
  templ[249] = 2;
  templ[63] = 2;
  templ[111] = 2;
  templ[215] = 2;
  templ[235] = 2;
  templ[246] = 2;
  templ[252] = 2;
  templ[159] = 2;
        
  templ[31] = 3;
  templ[107] = 3;
  templ[248] = 3;
  templ[214] = 3;
  templ[233] = 3;
  templ[47] = 3;
  templ[151] = 3;
  templ[244] = 3;
        
  templ[23] = 4;
  templ[240] = 4;
  templ[212] = 4;
  templ[105] = 4;
  templ[150] = 4;
  templ[43] = 4;
  templ[15] = 4;
  templ[232] = 4;
        
  templ[22] = 5;
  templ[11] = 5;
  templ[208] = 5;
  templ[104] = 5;
  templ[148] = 5;
  templ[41] = 5;
  templ[7] = 5;
  templ[224] = 5;
        
  templ[20] = 6;
  templ[3] = 6;
  templ[192] = 6;
  templ[40] = 6;
  templ[144] = 6;
  templ[9] = 6;
  templ[6] = 6;
  templ[96] = 6;

  v[0] = 1.0;
  v[1] = 0.75;  
  v[2] = 0.625;  
  v[3] = 0.5;
  v[4] = 0.375;
  v[5] = 0.25;
  v[6] = 0.125;
}

int imAnalyzeMeasurePerimArea(const imImage* image, double* perimarea_data, int region_count)
{
  static imbyte templ[256];
  static double vt[7];
  static int first = 1;
  if (first)
  {
    iInitPerimAreaTemplate(templ, vt);
    first = 0;
  }

  imushort* map = (imushort*)image->data[0];

  memset(perimarea_data, 0, region_count*sizeof(double));

  int width = image->width;
  int height = image->height;

  int counter = imProcessCounterBegin("PerimArea");
  imCounterTotal(counter, image->height, "Analyzing...");

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for (int y = 0; y < height; y++) 
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int offset_up = (y + 1)*width;
    int offset = y*width;
    int offset_dw = (y-1)*width;

    for (int x = 0; x < width; x++)
    {
      imushort v = map[offset+x];
      if (v)
      {
        int T = 0;
        if (x>0 && y<height-1 &&       map[offset_up + x-1] == v) T |= 0x01;
        if (y<height-1 &&              map[offset_up + x  ] == v) T |= 0x02;
        if (x<width-1 && y<height-1 && map[offset_up + x+1] == v) T |= 0x04;
        if (x>0 &&                     map[offset    + x-1] == v) T |= 0x08;
        if (x<width-1 &&               map[offset    + x+1] == v) T |= 0x10; 
        if (x>0 && y>0 &&              map[offset_dw + x-1] == v) T |= 0x20;
        if (y>0 &&                     map[offset_dw + x  ] == v) T |= 0x40;
        if (x<width-1 && y>0 &&        map[offset_dw + x+1] == v) T |= 0x80;

        if (T)
        {
          int index = v-1;
          double inc = vt[templ[T]];
#ifdef _OPENMP
#pragma omp atomic
#endif
          perimarea_data[index] += inc;
        }
      }
    }

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  imProcessCounterEnd(counter);
  return processing;
}

int imProcessRemoveByArea(const imImage* src_image, imImage* dst_image, int connect, int start_size, int end_size, int inside)
{
  int counter = imProcessCounterBegin("RemoveByArea");

  imImage *region_image = imImageCreate(src_image->width, src_image->height, IM_GRAY, IM_USHORT);
  if (!region_image)
  {
    imProcessCounterEnd(counter);
    return 0;
  }

  int region_count = 0;

  int ret = imAnalyzeFindRegions(src_image, region_image, connect, 1, &region_count);
  if (!region_count || !ret)
  {
    if (ret)
      imImageClear(dst_image);

    imImageDestroy(region_image);
    imProcessCounterEnd(counter);
    return ret;
  }

  if (end_size == 0)
    end_size = src_image->width*src_image->height;

  int outside;
  if (inside)
  {
    /* remove from inside */
    inside = 0;
    outside = 1;
  }
  else
  {
    /* remove from outside */
    inside = 1;
    outside = 0;
  }

  int* area_data = (int*)malloc(region_count*sizeof(int));
  ret = imAnalyzeMeasureArea(region_image, area_data, region_count);
  if (!ret)
  {
    free(area_data);
    imImageDestroy(region_image);
    imProcessCounterEnd(counter);
    return ret;
  }

  imushort* region_data = (imushort*)region_image->data[0];
  imbyte* img_data = (imbyte*)dst_image->data[0];

  imCounterTotal(counter, src_image->height, "Processing...");

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(src_image->count))
#endif
  for (int i = 0; i < src_image->count; i++)
  {
    if (i % src_image->width == 0)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_BEGIN_PROCESSING;

    if (region_data[i])
    {
      int area = area_data[region_data[i] - 1];
      if (area < start_size || area > end_size)
        img_data[i] = (imbyte)outside;
      else
        img_data[i] = (imbyte)inside;
    }
    else
      img_data[i] = 0;

    if (i % src_image->width == 0)
    {
      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_END_PROCESSING;
  }

  free(area_data);
  imImageDestroy(region_image);

  imProcessCounterEnd(counter);
  return processing;
}

int imProcessFillHoles(const imImage* src_image, imImage* dst_image, int connect)
{
  int counter = imProcessCounterBegin("FillHoles");

  // finding regions in the inverted src_image will isolate only the holes.
  imProcessNegative(src_image, dst_image);

  imImage *region_image = imImageCreate(src_image->width, src_image->height, IM_GRAY, IM_USHORT);
  if (!region_image)
  {
    imProcessCounterEnd(counter);
    return 0;
  }

  int holes_count = 0;
  
  int ret = imAnalyzeFindRegions(dst_image, region_image, connect, 0, &holes_count);
  if (!holes_count || !ret)
  {
    imImageCopy(src_image, dst_image);
    imImageDestroy(region_image);
    imProcessCounterEnd(counter);
    return ret;
  }

  imushort* region_data = (imushort*)region_image->data[0];
  imbyte* dst_data = (imbyte*)dst_image->data[0];

  imCounterTotal(counter, src_image->height, "Processing...");

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(src_image->count))
#endif
  for (int i = 0; i < src_image->count; i++)
  {
    if (i % src_image->width == 0)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_BEGIN_PROCESSING;

    if (region_data[i])
      dst_data[i] = 1;
    else
      dst_data[i] = !(dst_data[i]);  // Fix negative data.

    if (i % src_image->width == 0)
    {
      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_END_PROCESSING;
  }

  imImageDestroy(region_image);

  imProcessCounterEnd(counter);
  return processing;
}
