/** \file
 * \brief Geometric Operations
 *
 * See Copyright Notice in im_lib.h
 */


#include <im.h>
#include <im_util.h>

#include "im_process_counter.h"
#include "im_process_loc.h"
#include "im_math_op.h"

#include <stdlib.h>
#include <memory.h>


static inline void imRect2Polar(double xr, double yr, double *radius, double *theta)
{
  *radius = sqrt(xr*xr + yr*yr);
  *theta = atan2(yr, xr);
}

static inline void imPolar2Rect(double radius, double theta, double *xr, double *yr)
{
  *xr = radius * cos(theta);
  *yr = radius * sin(theta);
}

static inline void swirl_invtransf(int x, int y, double *xl, double *yl, double k, double xc, double yc)
{
  double radius, theta;
  double xr = x + 0.5 - xc;
  double yr = y + 0.5 - yc;

  imRect2Polar(xr, yr, &radius, &theta);

  theta += k * radius;

  imPolar2Rect(radius, theta, xl, yl);

  *xl += xc;
  *yl += yc;
}

template <class DT, class DTU> 
static int Swirl(int width, int height, DT *src_map, DT *dst_map, 
                         double k, int counter, DTU Dummy, int order)
{
  double xc = double(width/2.0);
  double yc = double(height/2.0);

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

    int line_offset = y*width;

    for (int x = 0; x < width; x++)
    {
      double xl, yl;
      swirl_invtransf(x, y, &xl, &yl, k, xc, yc);
                   
      // if inside the original image broad area
      if (xl > 0.0 && yl > 0.0 && xl < width && yl < height)
      {
        if (order == 1)
          dst_map[line_offset+x] = imBilinearInterpolation(width, height, src_map, xl, yl);
        else if (order == 3)
          dst_map[line_offset+x] = imBicubicInterpolation(width, height, src_map, xl, yl, Dummy);
        else
          dst_map[line_offset+x] = imZeroOrderInterpolation(width, height, src_map, xl, yl);
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

static inline void lensdistort_invtransf(int x, int y, double *xl, double *yl, double a, double b, double c, double d, double off_x, double off_y, double xc, double yc)
{
  double aux;
  double xr = (x + 0.5) - xc;
  double yr = (y + 0.5) - yc;

  double r = sqrt(xr*xr + yr*yr);
  aux = ((a * r + b)*r + c)*r + d;

  *xl = xr*aux;
  *yl = yr*aux;

  *xl += xc;
  *yl += yc;

  *xl -= off_x;
  *yl -= off_y;
}

template <class DT, class DTU>
static int LensDistort(int src_width, int src_height, DT *src_map, DT *dst_map,
                       double a, double b, double c, int dst_width, int dst_height, int counter, DTU Dummy, int order)
{
  double off_x = double(dst_width - src_width) / 2.0;
  double off_y = double(dst_height - src_height) / 2.0;
  double xc = double(dst_width / 2.0);
  double yc = double(dst_height / 2.0);
  // normalized by half of the smaller image dimension
  double norm = min_op(dst_width, dst_height) / 2.0;

  // Normalize the coefficient instead of each coordinate
  a /= norm*norm*norm;
  b /= norm*norm;
  c /= norm;
  double d = 1.0 - a - b - c;

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(dst_height))
#endif
  for (int y = 0; y < dst_height; y++)
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int line_offset = y*dst_width;

    for (int x = 0; x < dst_width; x++)
    {
      double xl, yl;
      lensdistort_invtransf(x, y, &xl, &yl, a, b, c, d, off_x, off_y, xc, yc);

      // if inside the original image broad area
      if (xl > 0.0 && yl > 0.0 && xl < src_width && yl < src_height)
      {
        if (order == 1)
          dst_map[line_offset + x] = imBilinearInterpolation(src_width, src_height, src_map, xl, yl);
        else if (order == 3)
          dst_map[line_offset + x] = imBicubicInterpolation(src_width, src_height, src_map, xl, yl, Dummy);
        else
          dst_map[line_offset + x] = imZeroOrderInterpolation(src_width, src_height, src_map, xl, yl);
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

static inline void radial_invtransf(int x, int y, double *xl, double *yl, double k1, double xc, double yc)
{
  double aux;
  double xr = x + 0.5 - xc;
  double yr = y + 0.5 - yc;

  aux = 1.0 + k1*(xr*xr + yr*yr);

  *xl = xr*aux;
  *yl = yr*aux;

  *xl += xc;
  *yl += yc;
}

template <class DT, class DTU> 
static int Radial(int width, int height, DT *src_map, DT *dst_map, 
                         double k1, int counter, DTU Dummy, int order)
{
  double xc = double(width/2.0);
  double yc = double(height/2.0);
  int diag = (int)sqrt(double(width*width + height*height));

  // Normalize the coefficient instead of each coordinate
  k1 /= (diag * diag);
         
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

    int line_offset = y*width;

    for (int x = 0; x < width; x++)
    {
      double xl, yl;
      radial_invtransf(x, y, &xl, &yl, k1, xc, yc);
                   
      // if inside the original image broad area
      if (xl > 0.0 && yl > 0.0 && xl < width && yl < height)
      {
        if (order == 1)
          dst_map[line_offset+x] = imBilinearInterpolation(width, height, src_map, xl, yl);
        else if (order == 3)
          dst_map[line_offset+x] = imBicubicInterpolation(width, height, src_map, xl, yl, Dummy);
        else
          dst_map[line_offset+x] = imZeroOrderInterpolation(width, height, src_map, xl, yl);
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

//*******************************************************************************************
//rotate_invtransf
//   shift the center to the origin of the target image
//   rotates centrered in the origin
//   shift the origin back to the center of the original image
//*******************************************************************************************

inline void rotate_invtransf(int x, int y, double *xl, double *yl, double cos0, double sin0, double dcx, double dcy, double scx, double scy)
{
  double xr = x + 0.5 - dcx;
  double yr = y + 0.5 - dcy;
  *xl = double(xr * cos0 - yr * sin0 + scx);
  *yl = double(xr * sin0 + yr * cos0 + scy);
}

template <class DT, class DTU> 
static int RotateCenter(int src_width, int src_height, DT *src_map, 
                        int dst_width, int dst_height, DT *dst_map, 
                        double cos0, double sin0, int counter, DTU Dummy, int order)
{
  double dcx = double(dst_width/2.0);
  double dcy = double(dst_height/2.0);
  double scx = double(src_width/2.0);
  double scy = double(src_height/2.0);

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(dst_height))
#endif
  for (int y = 0; y < dst_height; y++)
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int line_offset = y*dst_width;

    for (int x = 0; x < dst_width; x++)
    {
      double xl, yl;
      rotate_invtransf(x, y, &xl, &yl, cos0, sin0, dcx, dcy, scx, scy);
                   
      // if inside the original image broad area
      if (xl > 0.0 && yl > 0.0 && xl < src_width && yl < src_height)
      {
        if (order == 1)
          dst_map[line_offset+x] = imBilinearInterpolation(src_width, src_height, src_map, xl, yl);
        else if (order == 3)
          dst_map[line_offset+x] = imBicubicInterpolation(src_width, src_height, src_map, xl, yl, Dummy);
        else
          dst_map[line_offset+x] = imZeroOrderInterpolation(src_width, src_height, src_map, xl, yl);
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

template <class DT, class DTU> 
static int Rotate(int src_width, int src_height, DT *src_map, 
                  int dst_width, int dst_height, DT *dst_map, 
                  double cos0, double sin0, int ref_x, int ref_y, int to_origin, 
                  int counter, DTU Dummy, int order)
{
  double sx = double(ref_x);
  double sy = double(ref_y);
  double dx = sx;
  double dy = sy;
  if (to_origin)
  {
    dx = 0;
    dy = 0;
  }

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(dst_height))
#endif
  for (int y = 0; y < dst_height; y++)
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int line_offset = y*dst_width;

    for (int x = 0; x < dst_width; x++)
    {
      double xl, yl;
      rotate_invtransf(x, y, &xl, &yl, cos0, sin0, dx, dy, sx, sy);
                   
      // if inside the original image broad area
      if (xl > 0.0 && yl > 0.0 && xl < src_width && yl < src_height)
      {
        if (order == 1)
          dst_map[line_offset+x] = imBilinearInterpolation(src_width, src_height, src_map, xl, yl);
        else if (order == 3)
          dst_map[line_offset+x] = imBicubicInterpolation(src_width, src_height, src_map, xl, yl, Dummy);
        else
          dst_map[line_offset+x] = imZeroOrderInterpolation(src_width, src_height, src_map, xl, yl);
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


/********************************************************************************/


template <class DT> 
static int Rotate90(int src_width, 
                   int src_height, 
                   DT *src_map, 
                   DT *dst_map, 
                   int dir, 
                   int counter)
{
  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(src_height))
#endif
  for(int y = 0; y < src_height; y++)
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int yd, xd;

    // dir = clockwise (1) or counter clockwise (-1).
    // dst_width  = src_height
    // dst_height = src_width

    if (dir == 1)  
      xd = y;
    else
      xd = src_height-1 - y;

    int line_offset = y*src_width;

    for(int x = 0; x < src_width; x++)
    {
      if (dir == 1)
        yd = src_width-1 - x;
      else
        yd = x;

      dst_map[yd*src_height + xd] = src_map[line_offset + x];
    }        

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  return processing;
}

template <class DT> 
static int Rotate180(int width, 
                   int height, 
                   DT *src_map, 
                   DT *dst_map, 
                   int counter)
{
  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for(int y = 0; y < height; y++)
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int yd = height - 1 - y;

    int src_line_offset = y*width;
    int dst_line_offset = yd*width;

    for(int x = 0; x < width; x++)
    {
      int xd = width-1 - x;
      dst_map[dst_line_offset + xd] = src_map[src_line_offset + x];
    }        

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  return processing;
}

template <class DT> 
static int Mirror(int width, 
                   int height, 
                   DT *src_map, 
                   DT *dst_map, 
                   int counter)
{
  IM_INT_PROCESSING;

  if (src_map == dst_map) // check of in-place operation
  {
    int half_width = width/2;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
    for(int y = 0 ; y < height; y++)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
      IM_BEGIN_PROCESSING;

      int line_offset = y*width;

      for(int x = 0 ; x < half_width; x++)
      {
        int xd = width-1 - x;
        DT temp_value = src_map[line_offset + xd];
        src_map[line_offset + xd] = src_map[line_offset + x];
        src_map[line_offset + x] = temp_value;
        xd--;
      }        

      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
      IM_END_PROCESSING;
    }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
    for(int y = 0 ; y < height; y++)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
      IM_BEGIN_PROCESSING;

      int line_offset = y*width;

      for(int x = 0 ; x < width; x++)
      {
        int xd = width-1 - x;
        dst_map[line_offset + xd] = src_map[line_offset + x];
      }        

      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
      IM_END_PROCESSING;
    }
  }

  return processing;
}

template <class DT> 
static int Flip(int width, 
                   int height, 
                   DT *src_map, 
                   DT *dst_map, 
                   int counter)
{
  IM_INT_PROCESSING;

  if (src_map == dst_map) // check of in-place operation
  {
    DT* temp_line = (DT*)malloc(width*sizeof(DT));
    int half_height = height/2;

    // Can NOT run in parallel
    for(int y = 0 ; y < half_height; y++)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
      IM_BEGIN_PROCESSING;

      int yd = height - 1 - y;
      memcpy(temp_line, dst_map + yd*width, width*sizeof(DT));
      memcpy(dst_map + yd*width, src_map + y*width, width*sizeof(DT));
      memcpy(src_map + y*width, temp_line, width*sizeof(DT));

      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
      IM_END_PROCESSING;
    }

    free(temp_line);
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
    for(int y = 0 ; y < height; y++)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
      IM_BEGIN_PROCESSING;

      int yd = height - 1 - y;
      memcpy(dst_map + yd*width, src_map + y*width,width * sizeof(DT));

      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
      IM_END_PROCESSING;
    }
  }

  return processing;
}

template <class DT> 
static int InterlaceSplit(int width, 
                   int height, 
                   DT *src_map, 
                   DT *dst_map1,
                   DT *dst_map2, 
                   int counter)
{
  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for(int y = 0; y < height; y++)
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int yd = y / 2;
    if (y%2)
      memcpy(dst_map2 + yd*width, src_map + y*width, width*sizeof(DT));
    else
      memcpy(dst_map1 + yd*width, src_map + y*width, width*sizeof(DT));

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  return processing;
}


/********************************************************************************/


int imProcessRotate90(const imImage* src_image, imImage* dst_image, int dir)
{
  int ret = 0;

  int src_depth = src_image->has_alpha && dst_image->has_alpha ? src_image->depth + 1 : src_image->depth;
  int counter = imProcessCounterBegin("Rotate90");
  imCounterTotal(counter, src_depth*src_image->height, "Processing...");  /* size of the source image */

  for (int i = 0; i < src_depth; i++)
  {
    switch(src_image->data_type)
    {
    case IM_BYTE:
      ret = Rotate90(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)dst_image->data[i], dir, counter);
      break;
    case IM_SHORT:
      ret = Rotate90(src_image->width, src_image->height, (short*)src_image->data[i], (short*)dst_image->data[i], dir, counter);
      break;
    case IM_USHORT:
      ret = Rotate90(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)dst_image->data[i], dir, counter);
      break;
    case IM_INT:
      ret = Rotate90(src_image->width, src_image->height, (int*)src_image->data[i], (int*)dst_image->data[i], dir, counter);
      break;
    case IM_FLOAT:
      ret = Rotate90(src_image->width, src_image->height, (float*)src_image->data[i], (float*)dst_image->data[i], dir, counter);
      break;
    case IM_CFLOAT:
      ret = Rotate90(src_image->width, src_image->height, (imcfloat*)src_image->data[i], (imcfloat*)dst_image->data[i], dir, counter);
      break;
    case IM_DOUBLE:
      ret = Rotate90(src_image->width, src_image->height, (double*)src_image->data[i], (double*)dst_image->data[i], dir, counter);
      break;
    case IM_CDOUBLE:
      ret = Rotate90(src_image->width, src_image->height, (imcdouble*)src_image->data[i], (imcdouble*)dst_image->data[i], dir, counter);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

int imProcessRotate180(const imImage* src_image, imImage* dst_image)
{
  int src_depth = src_image->has_alpha && dst_image->has_alpha? src_image->depth+1: src_image->depth;

  int ret = 0;
  int counter = imProcessCounterBegin("Rotate180");
  imCounterTotal(counter, src_depth*src_image->height, "Processing...");

  for (int i = 0; i < src_depth; i++)
  {
    switch(src_image->data_type)
    {
    case IM_BYTE:
      ret = Rotate180(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)dst_image->data[i], counter);
      break;
    case IM_SHORT:
      ret = Rotate180(src_image->width, src_image->height, (short*)src_image->data[i], (short*)dst_image->data[i], counter);
      break;
    case IM_USHORT:
      ret = Rotate180(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)dst_image->data[i], counter);
      break;
    case IM_INT:
      ret = Rotate180(src_image->width, src_image->height, (int*)src_image->data[i], (int*)dst_image->data[i], counter);
      break;
    case IM_FLOAT:
      ret = Rotate180(src_image->width, src_image->height, (float*)src_image->data[i], (float*)dst_image->data[i], counter);
      break;
    case IM_CFLOAT:
      ret = Rotate180(src_image->width, src_image->height, (imcfloat*)src_image->data[i], (imcfloat*)dst_image->data[i], counter);
      break;
    case IM_DOUBLE:
      ret = Rotate180(src_image->width, src_image->height, (double*)src_image->data[i], (double*)dst_image->data[i], counter);
      break;
    case IM_CDOUBLE:
      ret = Rotate180(src_image->width, src_image->height, (imcdouble*)src_image->data[i], (imcdouble*)dst_image->data[i], counter);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

int imProcessLensDistort(const imImage* src_image, imImage* dst_image, double a, double b, double c, int order)
{
  int ret = 0;

  int counter = imProcessCounterBegin("LensDistort");
  int src_depth = src_image->has_alpha && dst_image->has_alpha ? src_image->depth + 1 : src_image->depth;
  imCounterTotal(counter, src_depth*dst_image->height, "Processing...");  /* size of the target image */

  for (int i = 0; i < src_depth; i++)
  {
    switch (src_image->data_type)
    {
    case IM_BYTE:
      ret = LensDistort(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)dst_image->data[i], a, b, c, dst_image->width, dst_image->height, counter, double(0), order);
      break;
    case IM_SHORT:
      ret = LensDistort(src_image->width, src_image->height, (short*)src_image->data[i], (short*)dst_image->data[i], a, b, c, dst_image->width, dst_image->height, counter, double(0), order);
      break;
    case IM_USHORT:
      ret = LensDistort(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)dst_image->data[i], a, b, c, dst_image->width, dst_image->height, counter, double(0), order);
      break;
    case IM_INT:
      ret = LensDistort(src_image->width, src_image->height, (int*)src_image->data[i], (int*)dst_image->data[i], a, b, c, dst_image->width, dst_image->height, counter, double(0), order);
      break;
    case IM_FLOAT:
      ret = LensDistort(src_image->width, src_image->height, (float*)src_image->data[i], (float*)dst_image->data[i], a, b, c, dst_image->width, dst_image->height, counter, double(0), order);
      break;
    case IM_CFLOAT:
      ret = LensDistort(src_image->width, src_image->height, (imcfloat*)src_image->data[i], (imcfloat*)dst_image->data[i], a, b, c, dst_image->width, dst_image->height, counter, imcfloat(0, 0), order);
      break;
    case IM_DOUBLE:
      ret = LensDistort(src_image->width, src_image->height, (double*)src_image->data[i], (double*)dst_image->data[i], a, b, c, dst_image->width, dst_image->height, counter, double(0), order);
      break;
    case IM_CDOUBLE:
      ret = LensDistort(src_image->width, src_image->height, (imcdouble*)src_image->data[i], (imcdouble*)dst_image->data[i], a, b, c, dst_image->width, dst_image->height, counter, imcdouble(0, 0), order);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

int imProcessRadial(const imImage* src_image, imImage* dst_image, double k1, int order)
{
  int ret = 0;

  int counter = imProcessCounterBegin("Radial");
  int src_depth = src_image->has_alpha && dst_image->has_alpha? src_image->depth+1: src_image->depth;
  imCounterTotal(counter, src_depth*dst_image->height, "Processing...");  /* size of the target image */

  for (int i = 0; i < src_depth; i++)
  {
    switch(src_image->data_type)
    {
    case IM_BYTE:
      ret = Radial(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)dst_image->data[i], k1, counter, double(0), order);
      break;
    case IM_SHORT:
      ret = Radial(src_image->width, src_image->height, (short*)src_image->data[i], (short*)dst_image->data[i], k1, counter, double(0), order);
      break;
    case IM_USHORT:
      ret = Radial(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)dst_image->data[i], k1, counter, double(0), order);
      break;
    case IM_INT:
      ret = Radial(src_image->width, src_image->height, (int*)src_image->data[i], (int*)dst_image->data[i], k1, counter, double(0), order);
      break;
    case IM_FLOAT:
      ret = Radial(src_image->width, src_image->height, (float*)src_image->data[i], (float*)dst_image->data[i], k1, counter, double(0), order);
      break;
    case IM_CFLOAT:
      ret = Radial(src_image->width, src_image->height, (imcfloat*)src_image->data[i], (imcfloat*)dst_image->data[i], k1, counter, imcfloat(0,0), order);
      break;
    case IM_DOUBLE:
      ret = Radial(src_image->width, src_image->height, (double*)src_image->data[i], (double*)dst_image->data[i], k1, counter, double(0), order);
      break;
    case IM_CDOUBLE:
      ret = Radial(src_image->width, src_image->height, (imcdouble*)src_image->data[i], (imcdouble*)dst_image->data[i], k1, counter, imcdouble(0, 0), order);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

int imProcessSwirl(const imImage* src_image, imImage* dst_image, double k, int order)
{
  int ret = 0;

  int counter = imProcessCounterBegin("Swirl");
  int src_depth = src_image->has_alpha && dst_image->has_alpha? src_image->depth+1: src_image->depth;
  imCounterTotal(counter, src_depth*dst_image->height, "Processing...");  /* size of the target image */

  for (int i = 0; i < src_depth; i++)
  {
    switch(src_image->data_type)
    {
    case IM_BYTE:
      ret = Swirl(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)dst_image->data[i], k, counter, double(0), order);
      break;
    case IM_SHORT:
      ret = Swirl(src_image->width, src_image->height, (short*)src_image->data[i], (short*)dst_image->data[i], k, counter, double(0), order);
      break;
    case IM_USHORT:
      ret = Swirl(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)dst_image->data[i], k, counter, double(0), order);
      break;
    case IM_INT:
      ret = Swirl(src_image->width, src_image->height, (int*)src_image->data[i], (int*)dst_image->data[i], k, counter, double(0), order);
      break;
    case IM_FLOAT:
      ret = Swirl(src_image->width, src_image->height, (float*)src_image->data[i], (float*)dst_image->data[i], k, counter, double(0), order);
      break;
    case IM_CFLOAT:
      ret = Swirl(src_image->width, src_image->height, (imcfloat*)src_image->data[i], (imcfloat*)dst_image->data[i], k, counter, imcfloat(0,0), order);
      break;
    case IM_DOUBLE:
      ret = Swirl(src_image->width, src_image->height, (double*)src_image->data[i], (double*)dst_image->data[i], k, counter, double(0), order);
      break;
    case IM_CDOUBLE:
      ret = Swirl(src_image->width, src_image->height, (imcdouble*)src_image->data[i], (imcdouble*)dst_image->data[i], k, counter, imcdouble(0, 0), order);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

//*******************************************************************************************
//rotate_transf
//   In this case shift to the origin, rotate, but do NOT shift back
//*******************************************************************************************

static void rotate_transf(double cx, double cy, int x, int y, double *xl, double *yl, double cos0, double sin0)
{
  double xr = x + 0.5 - cx;
  double yr = y + 0.5 - cy;
  *xl = double( xr*cos0 + yr*sin0);
  *yl = double(-xr*sin0 + yr*cos0);
}

void imProcessCalcRotateSize(int width, int height, int *new_width, int *new_height, double cos0, double sin0)
{
  double xl, yl, xmin, xmax, ymin, ymax;
  double wd2 = double(width)/2.0;
  double hd2 = double(height)/2.0;

  rotate_transf(wd2, hd2, 0, 0, &xl, &yl, cos0, sin0);
  xmin = xl; ymin = yl;
  xmax = xl; ymax = yl;

  rotate_transf(wd2, hd2, width-1, height-1, &xl, &yl, cos0, sin0);
  xmin = min_op(xmin, xl); ymin = min_op(ymin, yl);
  xmax = max_op(xmax, xl); ymax = max_op(ymax, yl);

  rotate_transf(wd2, hd2, 0, height-1, &xl, &yl, cos0, sin0);
  xmin = min_op(xmin, xl); ymin = min_op(ymin, yl);
  xmax = max_op(xmax, xl); ymax = max_op(ymax, yl);

  rotate_transf(wd2, hd2, width-1, 0, &xl, &yl, cos0, sin0);
  xmin = min_op(xmin, xl); ymin = min_op(ymin, yl);
  xmax = max_op(xmax, xl); ymax = max_op(ymax, yl);

  *new_width = (int)(xmax - xmin + 2.0);
  *new_height = (int)(ymax - ymin + 2.0);
}

int imProcessRotate(const imImage* src_image, imImage* dst_image, double cos0, double sin0, int order)
{
  int ret = 0;

  int counter = imProcessCounterBegin("Rotate");
  int src_depth = src_image->has_alpha && dst_image->has_alpha? src_image->depth+1: src_image->depth;
  imCounterTotal(counter, src_depth*dst_image->height, "Processing...");  /* size of the target image */

  if (src_image->color_space == IM_MAP)
  {
    ret = RotateCenter(src_image->width, src_image->height, (imbyte*)src_image->data[0], dst_image->width, dst_image->height, (imbyte*)dst_image->data[0], cos0, sin0, counter, double(0), 0);
  }
  else
  {
    for (int i = 0; i < src_depth; i++)
    {
      switch(src_image->data_type)
      {
      case IM_BYTE:
        ret = RotateCenter(src_image->width, src_image->height, (imbyte*)src_image->data[i], dst_image->width, dst_image->height, (imbyte*)dst_image->data[i], cos0, sin0, counter, double(0), order);
        break;
      case IM_SHORT:
        ret = RotateCenter(src_image->width, src_image->height, (short*)src_image->data[i], dst_image->width, dst_image->height, (short*)dst_image->data[i], cos0, sin0, counter, double(0), order);
        break;
      case IM_USHORT:
        ret = RotateCenter(src_image->width, src_image->height, (imushort*)src_image->data[i], dst_image->width, dst_image->height, (imushort*)dst_image->data[i], cos0, sin0, counter, double(0), order);
        break;
      case IM_INT:
        ret = RotateCenter(src_image->width, src_image->height, (int*)src_image->data[i], dst_image->width, dst_image->height, (int*)dst_image->data[i], cos0, sin0, counter, double(0), order);
        break;
      case IM_FLOAT:
        ret = RotateCenter(src_image->width, src_image->height, (float*)src_image->data[i], dst_image->width, dst_image->height, (float*)dst_image->data[i], cos0, sin0, counter, double(0), order);
        break;
      case IM_CFLOAT:
        ret = RotateCenter(src_image->width, src_image->height, (imcfloat*)src_image->data[i], dst_image->width, dst_image->height, (imcfloat*)dst_image->data[i], cos0, sin0, counter, imcfloat(0,0), order);
        break;
      case IM_DOUBLE:
        ret = RotateCenter(src_image->width, src_image->height, (double*)src_image->data[i], dst_image->width, dst_image->height, (double*)dst_image->data[i], cos0, sin0, counter, double(0), order);
        break;
      case IM_CDOUBLE:
        ret = RotateCenter(src_image->width, src_image->height, (imcdouble*)src_image->data[i], dst_image->width, dst_image->height, (imcdouble*)dst_image->data[i], cos0, sin0, counter, imcdouble(0, 0), order);
        break;
      }

      if (!ret)
        break;
    }
   }

  imProcessCounterEnd(counter);

  return ret;
}

int imProcessRotateRef(const imImage* src_image, imImage* dst_image, double cos0, double sin0, int x, int y, int to_origin, int order)
{
  int ret = 0;

  int counter = imProcessCounterBegin("RotateRef");
  int src_depth = src_image->has_alpha && dst_image->has_alpha? src_image->depth+1: src_image->depth;
  imCounterTotal(counter, src_depth*dst_image->height, "Processing...");  /* size of the target image */

  if (src_image->color_space == IM_MAP)
  {
    ret = Rotate(src_image->width, src_image->height, (imbyte*)src_image->data[0], dst_image->width, dst_image->height, (imbyte*)dst_image->data[0], cos0, sin0, x, y, to_origin, counter, double(0), 0);
  }
  else
  {
    for (int i = 0; i < src_depth; i++)
    {
      switch(src_image->data_type)
      {
      case IM_BYTE:
        ret = Rotate(src_image->width, src_image->height, (imbyte*)src_image->data[i], dst_image->width, dst_image->height, (imbyte*)dst_image->data[i], cos0, sin0, x, y, to_origin, counter, double(0), order);
        break;
      case IM_SHORT:
        ret = Rotate(src_image->width, src_image->height, (short*)src_image->data[i], dst_image->width, dst_image->height, (short*)dst_image->data[i], cos0, sin0, x, y, to_origin, counter, double(0), order);
        break;
      case IM_USHORT:
        ret = Rotate(src_image->width, src_image->height, (imushort*)src_image->data[i], dst_image->width, dst_image->height, (imushort*)dst_image->data[i], cos0, sin0, x, y, to_origin, counter, double(0), order);
        break;
      case IM_INT:
        ret = Rotate(src_image->width, src_image->height, (int*)src_image->data[i], dst_image->width, dst_image->height, (int*)dst_image->data[i], cos0, sin0, x, y, to_origin, counter, double(0), order);
        break;
      case IM_FLOAT:
        ret = Rotate(src_image->width, src_image->height, (float*)src_image->data[i], dst_image->width, dst_image->height, (float*)dst_image->data[i], cos0, sin0, x, y, to_origin, counter, double(0), order);
        break;
      case IM_CFLOAT:
        ret = Rotate(src_image->width, src_image->height, (imcfloat*)src_image->data[i], dst_image->width, dst_image->height, (imcfloat*)dst_image->data[i], cos0, sin0, x, y, to_origin, counter, imcfloat(0,0), order);
        break;
      case IM_DOUBLE:
        ret = Rotate(src_image->width, src_image->height, (double*)src_image->data[i], dst_image->width, dst_image->height, (double*)dst_image->data[i], cos0, sin0, x, y, to_origin, counter, double(0), order);
        break;
      case IM_CDOUBLE:
        ret = Rotate(src_image->width, src_image->height, (imcdouble*)src_image->data[i], dst_image->width, dst_image->height, (imcdouble*)dst_image->data[i], cos0, sin0, x, y, to_origin, counter, imcdouble(0, 0), order);
        break;
      }

      if (!ret)
        break;
    }
   }

  imProcessCounterEnd(counter);

  return ret;
}

int imProcessMirror(const imImage* src_image, imImage* dst_image)
{
  int i;
  int src_depth = src_image->has_alpha && dst_image->has_alpha? src_image->depth+1: src_image->depth;

  int ret = 0;
  int counter = imProcessCounterBegin("Mirror");
  imCounterTotal(counter, src_depth*src_image->height, "Processing...");

  for (i = 0; i < src_depth; i++)
  {
    switch(src_image->data_type)
    {
    case IM_BYTE:
      ret = Mirror(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)dst_image->data[i], counter);
      break;
    case IM_SHORT:
      ret = Mirror(src_image->width, src_image->height, (short*)src_image->data[i], (short*)dst_image->data[i], counter);
      break;
    case IM_USHORT:
      ret = Mirror(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)dst_image->data[i], counter);
      break;
    case IM_INT:
      ret = Mirror(src_image->width, src_image->height, (int*)src_image->data[i], (int*)dst_image->data[i], counter);
      break;
    case IM_FLOAT:
      ret = Mirror(src_image->width, src_image->height, (float*)src_image->data[i], (float*)dst_image->data[i], counter);
      break;
    case IM_CFLOAT:
      ret = Mirror(src_image->width, src_image->height, (imcfloat*)src_image->data[i], (imcfloat*)dst_image->data[i], counter);
      break;
    case IM_DOUBLE:
      ret = Mirror(src_image->width, src_image->height, (double*)src_image->data[i], (double*)dst_image->data[i], counter);
      break;
    case IM_CDOUBLE:
      ret = Mirror(src_image->width, src_image->height, (imcdouble*)src_image->data[i], (imcdouble*)dst_image->data[i], counter);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

int imProcessFlip(const imImage* src_image, imImage* dst_image)
{
  int i;
  int src_depth = src_image->has_alpha && dst_image->has_alpha? src_image->depth+1: src_image->depth;

  int ret = 0;
  int counter = imProcessCounterBegin("Flip");
  imCounterTotal(counter, src_depth*src_image->height, "Processing...");

  for (i = 0; i < src_depth; i++)
  {
    switch(src_image->data_type)
    {
    case IM_BYTE:
      ret = Flip(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)dst_image->data[i], counter);
      break;
    case IM_SHORT:
      ret = Flip(src_image->width, src_image->height, (short*)src_image->data[i], (short*)dst_image->data[i], counter);
      break;
    case IM_USHORT:
      ret = Flip(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)dst_image->data[i], counter);
      break;
    case IM_INT:
      ret = Flip(src_image->width, src_image->height, (int*)src_image->data[i], (int*)dst_image->data[i], counter);
      break;
    case IM_FLOAT:
      ret = Flip(src_image->width, src_image->height, (float*)src_image->data[i], (float*)dst_image->data[i], counter);
      break;
    case IM_CFLOAT:
      ret = Flip(src_image->width, src_image->height, (imcfloat*)src_image->data[i], (imcfloat*)dst_image->data[i], counter);
      break;
    case IM_DOUBLE:
      ret = Flip(src_image->width, src_image->height, (double*)src_image->data[i], (double*)dst_image->data[i], counter);
      break;
    case IM_CDOUBLE:
      ret = Flip(src_image->width, src_image->height, (imcdouble*)src_image->data[i], (imcdouble*)dst_image->data[i], counter);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

int imProcessInterlaceSplit(const imImage* src_image, imImage* dst_image1, imImage* dst_image2)
{
  int i;
  int src_depth = src_image->has_alpha && dst_image1->has_alpha && dst_image2->has_alpha ? src_image->depth + 1 : src_image->depth;

  int ret = 0;
  int counter = imProcessCounterBegin("InterlaceSplit");
  imCounterTotal(counter, src_depth*src_image->height, "Processing...");

  for (i = 0; i < src_depth; i++)
  {
    switch(src_image->data_type)
    {
    case IM_BYTE:
      ret = InterlaceSplit(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)dst_image1->data[i], (imbyte*)dst_image2->data[i], counter);
      break;
    case IM_SHORT:
      ret = InterlaceSplit(src_image->width, src_image->height, (short*)src_image->data[i], (short*)dst_image1->data[i], (short*)dst_image2->data[i], counter);
      break;
    case IM_USHORT:
      ret = InterlaceSplit(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)dst_image1->data[i], (imushort*)dst_image2->data[i], counter);
      break;
    case IM_INT:
      ret = InterlaceSplit(src_image->width, src_image->height, (int*)src_image->data[i], (int*)dst_image1->data[i], (int*)dst_image2->data[i], counter);
      break;
    case IM_FLOAT:
      ret = InterlaceSplit(src_image->width, src_image->height, (float*)src_image->data[i], (float*)dst_image1->data[i], (float*)dst_image2->data[i], counter);
      break;
    case IM_CFLOAT:
      ret = InterlaceSplit(src_image->width, src_image->height, (imcfloat*)src_image->data[i], (imcfloat*)dst_image1->data[i], (imcfloat*)dst_image2->data[i], counter);
      break;
    case IM_DOUBLE:
      ret = InterlaceSplit(src_image->width, src_image->height, (double*)src_image->data[i], (double*)dst_image1->data[i], (double*)dst_image2->data[i], counter);
      break;
    case IM_CDOUBLE:
      ret = InterlaceSplit(src_image->width, src_image->height, (imcdouble*)src_image->data[i], (imcdouble*)dst_image1->data[i], (imcdouble*)dst_image2->data[i], counter);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

