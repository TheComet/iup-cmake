/** \file
 * \brief Distance Transform
 *
 * See Copyright Notice in im_lib.h
 */

#include <im.h>
#include <im_util.h>

#include "im_process_counter.h"
#include "im_process_glo.h"

#include <stdlib.h>
#include <memory.h>
#include <math.h>

const double DT_ONE    = 1.0;             // 1x0
const double DT_SQRT2 = 1.414213562373;  // 1x1
const double DT_SQRT5 = 2.2360679775;    // 2x1
const double DT_SQRT10 = 3.1622776601684; // 3x1
const double DT_SQRT13 = 3.605551275464;  // 3x2
const double DT_SQRT17 = 4.12310562562;   // 4x1
const double DT_SQRT25 = 5.0;             // 4x3

template<class T>
static inline void setValue(int r, int r1, int r2, int r3, int r4, T *image_data, int f)
{
  T v;
  T minv = image_data[r];        // (x,y)

  if (f)
    v = image_data[r - 1] + (T)DT_ONE;      // (x-1,y)
  else
    v = image_data[r + 1] + (T)DT_ONE;      // (x+1,y)
  if (v < minv) minv = v;

  v = image_data[r1] + (T)DT_ONE;          // (x,y-1)           (x,y+1)
  if (v < minv) minv = v; 

  if (minv < (T)DT_SQRT2)
    goto min_attrib;

  v = image_data[r1 - 1] + (T)DT_SQRT2;    // (x-1,y-1)         (x-1,y+1)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r1 + 1] + (T)DT_SQRT2;    // (x+1,y-1)         (x+1,y+1)
  if (v < minv) minv = v;                                      

  if (minv < (T)DT_SQRT5)
    goto min_attrib;

  v = image_data[r1 + 2] + (T)DT_SQRT5;    // (x+2,y-1)         (x+2,y+1)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r1 - 2] + (T)DT_SQRT5;    // (x-2,y-1)         (x-2,y+1)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r2 - 1] + (T)DT_SQRT5;    // (x-1,y-2)         (x-1,y+2)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r2 + 1] + (T)DT_SQRT5;    // (x+1,y-2)         (x+1,y+2)
  if (v < minv) minv = v;                                      

  if (minv < (T)DT_SQRT10)
    goto min_attrib;
                                                             
  v = image_data[r1 + 3] + (T)DT_SQRT10;   // (x+3,y-1)         (x+3,y+1)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r1 - 3] + (T)DT_SQRT10;   // (x-3,y-1)         (x-3,y+1)
  if (v < minv) minv = v;                                      

  v = image_data[r3 - 1] + (T)DT_SQRT10;   // (x-1,y-3)         (x-1,y+3)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r3 + 1] + (T)DT_SQRT10;   // (x+1,y-3)         (x+1,y+3)
  if (v < minv) minv = v;                                      

  if (minv < (T)DT_SQRT13)
    goto min_attrib;

  v = image_data[r2 - 3] + (T)DT_SQRT13;   // (x-3,y-2)         (x-3,y+2)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r2 + 3] + (T)DT_SQRT13;   // (x+3,y-2)         (x+3,y+2)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r3 + 2] + (T)DT_SQRT13;   // (x+2,y-3)         (x+2,y+3)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r3 - 2] + (T)DT_SQRT13;   // (x-2,y-3)         (x-2,y+3)
  if (v < minv) minv = v;

  if (minv < (T)DT_SQRT17)
    goto min_attrib;
                                                             
  v = image_data[r1 + 4] + (T)DT_SQRT17;   // (x+4,y-1)         (x+4,y+1)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r1 - 4] + (T)DT_SQRT17;   // (x-4,y-1)         (x-4,y+1)
  if (v < minv) minv = v;                                      

  v = image_data[r4 - 1] + (T)DT_SQRT17;   // (x-1,y-4)         (x-1,y+4)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r4 + 1] + (T)DT_SQRT17;   // (x+1,y-4)         (x+1,y+4)
  if (v < minv) minv = v;                                      

  if (minv < (T)DT_SQRT25)
    goto min_attrib;

  v = image_data[r3 - 4] + (T)DT_SQRT25;   // (x-4,y-3)         (x-4,y+3)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r3 + 4] + (T)DT_SQRT25;   // (x+4,y-3)         (x+4,y+3)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r4 + 3] + (T)DT_SQRT25;   // (x+3,y-4)         (x+3,y+4)
  if (v < minv) minv = v;                                      
                                                             
  v = image_data[r4 - 3] + (T)DT_SQRT25;   // (x-3,y-4)         (x-3,y+4)
  if (v < minv) minv = v;

min_attrib:
  image_data[r] = minv;
}

template<class T>
static inline void setValueForwardEdge(int r, int r1, int r2, int width, int x, int y, T *image_data)
{
  T v;
  T minv = image_data[r];        // (x,y)

  if (y > 0)
  {
    v = image_data[r1] + (T)DT_ONE;         // (x,y-1)
    if (v < minv) minv = v;
  }

  if (x > 0)
  {
    v = image_data[r - 1] + (T)DT_ONE;      // (x-1,y)
    if (v < minv) minv = v;
  }

  if (x > 0 && y > 0)
  {
    v = image_data[r1 - 1] + (T)DT_SQRT2;   // (x-1,y-1)
    if (v < minv) minv = v;
  }

  if (x < width-2 && y > 0)
  {
    v = image_data[r1 + 1] + (T)DT_SQRT2;   // (x+1,y-1)
    if (v < minv) minv = v;
  }

  if (x > 0 && y > 1)
  {
    v = image_data[r2 - 1] + (T)DT_SQRT5;   // (x-1,y-2)
    if (v < minv) minv = v;
  }

  if (x < width-2 && y > 1)
  {
    v = image_data[r2 + 1] + (T)DT_SQRT5;   // (x+1,y-2)
    if (v < minv) minv = v;
  }

  if (x < width-3 && y > 0)
  {
    v = image_data[r1 + 2] + (T)DT_SQRT5;   // (x+2,y-1)
    if (v < minv) minv = v;
  }

  if (x > 1 && y > 0)
  {
    v = image_data[r1 - 2] + (T)DT_SQRT5;   // (x-2,y-1)
    if (v < minv) minv = v;
  }

  image_data[r] = minv;
}

template<class T>
static inline void setValueBackwardEdge(int r, int r1, int r2, int width, int height, int x, int y, T *image_data)
{
  T v;
  T minv = image_data[r];        // (x,y)

  if (x < width-2)
  {
    v = image_data[r + 1] + (T)DT_ONE;      // (x+1,y)
    if (v < minv) minv = v;
  }

  if (y < height-2)
  {
    v = image_data[r1] + (T)DT_ONE;         // (x,y+1)
    if (v < minv) minv = v;
  }

  if (y < height-2 && x > 0)
  {
    v = image_data[r1 - 1] + (T)DT_SQRT2;   // (x-1,y+1)
    if (v < minv) minv = v;
  }

  if (y < height-2 && x < width-2)
  {
    v = image_data[r1 + 1] + (T)DT_SQRT2;   // (x+1,y+1)
    if (v < minv) minv = v;
  }

  if (y < height-2 && x < width-3)
  {
    v = image_data[r1 + 2] + (T)DT_SQRT5;   // (x+2,y+1)
    if (v < minv) minv = v;
  }

  if (y < height-3 && x < width-2)
  {
    v = image_data[r2 + 1] + (T)DT_SQRT5;   // (x+1,y+2)
    if (v < minv) minv = v;
  }

  if (y < height-3 && x > 0)
  {
    v = image_data[r2 - 1] + (T)DT_SQRT5;   // (x-1,y+2)
    if (v < minv) minv = v;
  }

  if (y < height-2 && x > 1)
  {
    v = image_data[r1 - 2] + (T)DT_SQRT5;   // (x-2,y+1)
    if (v < minv) minv = v;
  }

  image_data[r] = minv;
}

template<class T>
static void iDoDistanceTransform(int width, int height, imbyte* src_data, T* dst_data)
{
  int count  = width*height;
  double max_dist = sqrt(double(width*width + height*height));

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    // if pixel is background, then distance is zero.
    if (src_data[i])
      dst_data[i] = (T)max_dist;
  }

  /* down->top, left->right */
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for (int y = 0; y < height; y++)
  {
    int offset = y * width;
    int offset1 = offset - width;
    int offset2 = offset - 2 * width;
    int offset3 = offset - 3 * width;
    int offset4 = offset - 4 * width;

    for (int x = 0; x < width; x++)
    {
      if (src_data[offset])
      {
        if (x < 4 || x > width - 5 || y < 4 || y > height - 5)
          setValueForwardEdge(offset, offset1, offset2, width, x, y, dst_data);
        else
          setValue(offset, offset1, offset2, offset3, offset4, dst_data, 1);
      }

      offset++;
      offset1++;
      offset2++;
      offset3++;
      offset4++;
    }
  }

  /* top->down, right->left */
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for (int y = height - 1; y >= 0; y--)
  {
    int offset = y * width + width - 1;
    int offset1 = offset + width;
    int offset2 = offset + 2 * width;
    int offset3 = offset + 3 * width;
    int offset4 = offset + 4 * width;

    for (int x = width - 1; x >= 0; x--)
    {
      if (src_data[offset])
      {
        if (x < 4 || x > width - 5 || y < 4 || y > height - 5)
          setValueBackwardEdge(offset, offset1, offset2, width, height, x, y, dst_data);
        else
          setValue(offset, offset1, offset2, offset3, offset4, dst_data, 0);
      }

      offset--;
      offset1--;
      offset2--;
      offset3--;
      offset4--;
    }
  }
}

void imProcessDistanceTransform(const imImage* src_image, imImage* dst_image)
{
  if (dst_image->data_type == IM_FLOAT)
    iDoDistanceTransform(src_image->width, src_image->height, (imbyte*)src_image->data[0], (float*)dst_image->data[0]);
  else
    iDoDistanceTransform(src_image->width, src_image->height, (imbyte*)src_image->data[0], (double*)dst_image->data[0]);
}

static void iFillValue(imbyte* img_data, int x, int y, int width, int value)
{
  int r = y * width + x;
  int r1a = r - width;
  int r1b = r + width;
  int v;

  int old_value = img_data[r];
  img_data[r] = (imbyte)value;

  v = img_data[r1a];        // (x,y-1)
  if (v == old_value) 
    iFillValue(img_data, x, y-1, width, value);

  v = img_data[r - 1];      // (x-1,y)
  if (v == old_value) 
    iFillValue(img_data, x-1, y, width, value);

  v = img_data[r1a - 1];    // (x-1,y-1)
  if (v == old_value) 
    iFillValue(img_data, x-1, y-1, width, value);

  v = img_data[r1a + 1];    // (x+1,y-1)
  if (v == old_value) 
    iFillValue(img_data, x+1, y-1, width, value);

  v = img_data[r + 1];      // (x+1,y)
  if (v == old_value) 
    iFillValue(img_data, x+1, y, width, value);

  v = img_data[r1b];        // (x,y+1)
  if (v == old_value) 
    iFillValue(img_data, x, y+1, width, value);

  v = img_data[r1b - 1];    // (x-1,y+1)
  if (v == old_value) 
    iFillValue(img_data, x-1, y+1, width, value);

  v = img_data[r1b + 1];    // (x+1,y+1)
  if (v == old_value) 
    iFillValue(img_data, x+1, y+1, width, value);
}

template<class T>
static inline int iCheckFalseMaximum(int r, int r2a, int r2b, int width, T *src_data)
{
  /* we are ignoring valeys of 1 pixel. */
  /* this is not 100% fail proof */
  T v;
  T maxv = src_data[r];  // (x,y)
  int r1a = r - width;
  int r1b = r + width;

  v = src_data[r2a - 1];    // (x-1,y-2)
  if (v > maxv) return 1;

  v = src_data[r2a];        // (x,y-2)
  if (v > maxv) return 1;

  v = src_data[r2a + 1];    // (x+1,y-2)
  if (v > maxv) return 1;

  v = src_data[r2b - 1];    // (x-1,y+2)
  if (v > maxv) return 1;

  v = src_data[r2b];        // (x,y+2)
  if (v > maxv) return 1;

  v = src_data[r2b + 1];    // (x+1,y+2)
  if (v > maxv) return 1;


  v = src_data[r2b - 2];    // (x-2,y+2)
  if (v > maxv) return 1;

  v = src_data[r1b - 2];    // (x-2,y+1)
  if (v > maxv) return 1;

  v = src_data[r - 2];      // (x-2,y)
  if (v > maxv) return 1;

  v = src_data[r1a - 2];    // (x-2,y-1)
  if (v > maxv) return 1;

  v = src_data[r2a - 2];    // (x-2,y-2)
  if (v > maxv) return 1;


  v = src_data[r2a + 2];    // (x+2,y-2)
  if (v > maxv) return 1;

  v = src_data[r1a + 2];    // (x+2,y-1)
  if (v > maxv) return 1;

  v = src_data[r + 2];      // (x+2,y)
  if (v > maxv) return 1;

  v = src_data[r1b + 2];    // (x+2,y+1)
  if (v > maxv) return 1;

  v = src_data[r2b + 2];    // (x+2,y+2)
  if (v > maxv) return 1;

  return 0;
}

template<class T>
static inline void iCheckMaximum(int r, int r1a, int r1b, T *src_data, imbyte* dst_data)
{
  int unique = 1;
  T v;
  T maxv = src_data[r];  // (x,y)

  v = src_data[r1a];        // (x,y-1)
  if (v >= maxv) { maxv = v; unique = 0; }

  v = src_data[r - 1];      // (x-1,y)
  if (v >= maxv) { maxv = v; unique = 0; }

  v = src_data[r1a - 1];    // (x-1,y-1)
  if (v >= maxv) { maxv = v; unique = 0; }

  v = src_data[r1a + 1];    // (x+1,y-1)
  if (v >= maxv) { maxv = v; unique = 0; }

  v = src_data[r + 1];      // (x+1,y)
  if (v >= maxv) { maxv = v; unique = 0; }

  v = src_data[r1b];        // (x,y+1)
  if (v >= maxv) { maxv = v; unique = 0; }

  v = src_data[r1b - 1];    // (x-1,y+1)
  if (v >= maxv) { maxv = v; unique = 0; }

  v = src_data[r1b + 1];    // (x+1,y+1)
  if (v >= maxv) { maxv = v; unique = 0; }

  if (src_data[r] < maxv)   // not a maximum
    dst_data[r] = 0;
  else
  {
    if (unique)            // unique maximum
      dst_data[r] = 1;
    else                   // can be maximum
      dst_data[r] = 2;
  }
}

template<class T>
static void iDoRegionalMaximum(int width, int height, T* src_data, imbyte* dst_data)
{
  int count = width*height;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for (int y = 1; y < height - 1; y++)
  {
    int offset = y * width + 1;
    int offsetA = offset - width;
    int offsetB = offset + width;

    for (int x = 1; x < width - 1; x++)
    {
      if (src_data[offset])
        iCheckMaximum(offset, offsetA, offsetB, src_data, dst_data);

      offset++;
      offsetA++;
      offsetB++;
    }
  }

  // remove false maximum
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for (int y = 2; y < height - 2; y++)
  {
    int offset = y * width + 2;
    int offsetA = offset - 2 * width;
    int offsetB = offset + 2 * width;

    for (int x = 2; x < width - 2; x++)
    {
      if (dst_data[offset] == 2)
      {
        if (iCheckFalseMaximum(offset, offsetA, offsetB, width, src_data))
          iFillValue(dst_data, x, y, width, 0);
      }

      offset++;
      offsetA++;
      offsetB++;
    }
  }

  // update target with remaining maximum
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    if (dst_data[i] == 2)
      dst_data[i] = 1;
  }
}

void imProcessRegionalMaximum(const imImage* src_image, imImage* dst_image)
{
  if (src_image->data_type == IM_FLOAT)
    iDoRegionalMaximum(src_image->width, src_image->height, (float*)src_image->data[0], (imbyte*)dst_image->data[0]);
  else
    iDoRegionalMaximum(src_image->width, src_image->height, (double*)src_image->data[0], (imbyte*)dst_image->data[0]);
}
