/** \file
 * \brief Color Processing Operations
 *
 * See Copyright Notice in im_lib.h
 */

#include <im.h>
#include <im_util.h>
#include <im_color.h>
#include <im_colorhsi.h>
#include <im_palette.h>

#include "im_process_counter.h"
#include "im_process_pnt.h"

#include <stdlib.h>
#include <memory.h>


template <class T>
static void DoPseudoColor(T* src_data, imbyte** dst_data, int count)
{
  imbyte *red = dst_data[0],
    *green = dst_data[1],
    *blue = dst_data[2];
  T min, max;
  unsigned char r, g, b;

  imMinMaxType(src_data, count, min, max);

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    double norm = (double)(src_data[i] - min) / (double)(max - min);  // now 0 <= norm <= 1

    imColorHSI2RGBbyte(norm * 360, 1.0, norm, &r, &g, &b);

    red[i] = r;
    green[i] = g;
    blue[i] = b;
  }
}

void imProcessPseudoColor(const imImage* src_image, imImage* dst_image)
{
  switch (src_image->data_type)
  {
  case IM_BYTE:
    DoPseudoColor((imbyte*)src_image->data[0], (imbyte**)dst_image->data, src_image->count);
    break;
  case IM_SHORT:
    DoPseudoColor((short*)src_image->data[0], (imbyte**)dst_image->data, src_image->count);
    break;
  case IM_USHORT:
    DoPseudoColor((imushort*)src_image->data[0], (imbyte**)dst_image->data, src_image->count);
    break;
  case IM_FLOAT:
    DoPseudoColor((float*)src_image->data[0], (imbyte**)dst_image->data, src_image->count);
    break;
  case IM_DOUBLE:
    DoPseudoColor((double*)src_image->data[0], (imbyte**)dst_image->data, src_image->count);
    break;
  }
}

static void rgb2yrgb(imbyte* r, imbyte* g, imbyte* b, imbyte* y)
{
  int ri,gi,bi;

  *y = imColorRGB2Luma(*r, *g, *b);
  ri = (*r) - (*y) + 128;
  gi = (*g) - (*y) + 128;
  bi = (*b) - (*y) + 128;

  if (ri < 0) ri = 0;
  if (gi < 0) gi = 0;
  if (bi < 0) bi = 0;

  *r = (imbyte)ri;
  *g = (imbyte)gi;
  *b = (imbyte)bi;
}

void imProcessSplitYChroma(const imImage* src_image, imImage* y_image, imImage* chroma_image)
{
  imbyte 
    *red=(imbyte*)src_image->data[0],
    *green=(imbyte*)src_image->data[1],
    *blue=(imbyte*)src_image->data[2],
    *red2=(imbyte*)chroma_image->data[0],
    *green2=(imbyte*)chroma_image->data[1],
    *blue2=(imbyte*)chroma_image->data[2],
    *map1=(imbyte*)y_image->data[0];

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(src_image->count))
#endif
  for (int i = 0; i < src_image->count; i++)
  {
    imbyte R = red[i];
    imbyte G = green[i];
    imbyte B = blue[i];
    imbyte Y;

    rgb2yrgb(&R, &G, &B, &Y);

    map1[i] = Y;

    red2[i] = R;
    green2[i] = G;
    blue2[i] = B;
  }
}

template <class T>
static void DoSplitHSIReal(T** data, T* hue, T* saturation, T* intensity, int count)
{
  T *red = data[0],
      *green=data[1],
       *blue=data[2];
  double h, s, i;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int ii = 0; ii < count; ii++)
  {
    imColorRGB2HSI(red[ii], green[ii], blue[ii], &h, &s, &i);

    hue[ii] = (T)h;
    saturation[ii] = (T)s;
    intensity[ii] = (T)i;
  }
}

template <class T>
static void DoSplitHSIByte(imbyte** data, T* hue, T* saturation, T* intensity, int count)
{
  imbyte *red=data[0],
       *green=data[1],
        *blue=data[2];
  double h, s, i;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int ii = 0; ii < count; ii++)
  {
    imColorRGB2HSIbyte(red[ii], green[ii], blue[ii], &h, &s, &i);

    hue[ii] = (T)h;
    saturation[ii] = (T)s;
    intensity[ii] = (T)i;
  }
}

void imProcessSplitHSI(const imImage* src_image, imImage* dst_image1, imImage* dst_image2, imImage* dst_image3)
{
  switch(src_image->data_type)
  {
  case IM_BYTE:
    if (dst_image1->data_type == IM_FLOAT)
      DoSplitHSIByte((imbyte**)src_image->data, (float*)dst_image1->data[0], (float*)dst_image2->data[0], (float*)dst_image3->data[0], src_image->count);
    else
      DoSplitHSIByte((imbyte**)src_image->data, (double*)dst_image1->data[0], (double*)dst_image2->data[0], (double*)dst_image3->data[0], src_image->count);
    break;
  case IM_FLOAT:                                                                                                                               
    DoSplitHSIReal((float**)src_image->data, (float*)dst_image1->data[0], (float*)dst_image2->data[0], (float*)dst_image3->data[0], src_image->count);
    break;                                                                                
  case IM_DOUBLE:
    DoSplitHSIReal((double**)src_image->data, (double*)dst_image1->data[0], (double*)dst_image2->data[0], (double*)dst_image3->data[0], src_image->count);
    break;
  }

  imImageSetPalette(dst_image1, imPaletteHues(), 256);
}

template <class T>
static void DoMergeHSIReal(T** data, T* hue, T* saturation, T* intensity, int count)
{
  T *red = data[0],
      *green=data[1],
       *blue=data[2];
  double R, G, B;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    imColorHSI2RGB(hue[i], saturation[i], intensity[i], &R, &G, &B);

    red[i] = (T)R;
    green[i] = (T)G;
    blue[i] = (T)B;
  }
}

template <class T>
static void DoMergeHSIByte(imbyte** data, T* hue, T* saturation, T* intensity, int count)
{
  imbyte *red=data[0],
       *green=data[1],
        *blue=data[2];

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    imColorHSI2RGBbyte(hue[i], saturation[i], intensity[i], &red[i], &green[i], &blue[i]);
  }
}

void imProcessMergeHSI(const imImage* src_image1, const imImage* src_image2, const imImage* src_image3, imImage* dst_image)
{
  switch(dst_image->data_type)
  {
  case IM_BYTE:
    if (src_image1->data_type == IM_FLOAT)
      DoMergeHSIByte((imbyte**)dst_image->data, (float*)src_image1->data[0], (float*)src_image2->data[0], (float*)src_image3->data[0], dst_image->count);
    else
      DoMergeHSIByte((imbyte**)dst_image->data, (double*)src_image1->data[0], (double*)src_image2->data[0], (double*)src_image3->data[0], dst_image->count);
    break;                                                                                                                                    
  case IM_FLOAT:                                                                                                                               
    DoMergeHSIReal((float**)dst_image->data, (float*)src_image1->data[0], (float*)src_image2->data[0], (float*)src_image3->data[0], dst_image->count);
    break;                                                                                
  case IM_DOUBLE:
    DoMergeHSIReal((double**)dst_image->data, (double*)src_image1->data[0], (double*)src_image2->data[0], (double*)src_image3->data[0], dst_image->count);
    break;
  }
}

void imProcessSplitComponents(const imImage* src_image, imImage** dst_image)
{
  memcpy(dst_image[0]->data[0], src_image->data[0], src_image->plane_size);
  memcpy(dst_image[1]->data[0], src_image->data[1], src_image->plane_size);
  memcpy(dst_image[2]->data[0], src_image->data[2], src_image->plane_size);
  if (imColorModeDepth(src_image->color_space) == 4 || src_image->has_alpha) 
    memcpy(dst_image[3]->data[0], src_image->data[3], src_image->plane_size);
}

void imProcessMergeComponents(const imImage** src_image, imImage* dst_image)
{
  memcpy(dst_image->data[0], src_image[0]->data[0], dst_image->plane_size);
  memcpy(dst_image->data[1], src_image[1]->data[0], dst_image->plane_size);
  memcpy(dst_image->data[2], src_image[2]->data[0], dst_image->plane_size);
  if (imColorModeDepth(dst_image->color_space) == 4 || dst_image->has_alpha) 
    memcpy(dst_image->data[3], src_image[3]->data[0], dst_image->plane_size);
}

template <class ST, class DT>
static void DoNormalizeComp(ST** src_data, DT** dst_data, int count, int depth)
{
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    int d;

    DT sum = 0;
    for(d = 0; d < depth; d++)
      sum += (DT)(src_data[d][i]);

    for(d = 0; d < depth; d++)
    {
      if (sum == 0)
        dst_data[d][i] = 0;
      else
        dst_data[d][i] = (DT)(src_data[d][i]) / sum;
    }
  }
}

void imProcessNormalizeComponents(const imImage* src_image, imImage* dst_image)
{
  switch(src_image->data_type)
  {
  case IM_BYTE:
    if (dst_image->data_type == IM_FLOAT)
      DoNormalizeComp((imbyte**)src_image->data, (float**)dst_image->data, src_image->count, src_image->depth);
    else
      DoNormalizeComp((imbyte**)src_image->data, (double**)dst_image->data, src_image->count, src_image->depth);
    break;
  case IM_SHORT:                                                                                                                               
    if (dst_image->data_type == IM_FLOAT)
      DoNormalizeComp((short**)src_image->data, (float**)dst_image->data, src_image->count, src_image->depth);
    else
      DoNormalizeComp((short**)src_image->data, (double**)dst_image->data, src_image->count, src_image->depth);
    break;
  case IM_USHORT:                                                                                                                               
    if (dst_image->data_type == IM_FLOAT)
      DoNormalizeComp((imushort**)src_image->data, (float**)dst_image->data, src_image->count, src_image->depth);
    else
      DoNormalizeComp((imushort**)src_image->data, (double**)dst_image->data, src_image->count, src_image->depth);
    break;
  case IM_INT:                                                                                                                               
    if (dst_image->data_type == IM_FLOAT)
      DoNormalizeComp((int**)src_image->data, (float**)dst_image->data, src_image->count, src_image->depth);
    else
      DoNormalizeComp((int**)src_image->data, (double**)dst_image->data, src_image->count, src_image->depth);
    break;
  case IM_FLOAT:                                                                                                                               
    if (dst_image->data_type == IM_FLOAT)
      DoNormalizeComp((float**)src_image->data, (float**)dst_image->data, src_image->count, src_image->depth);
    else
      DoNormalizeComp((float**)src_image->data, (double**)dst_image->data, src_image->count, src_image->depth);
    break;
  case IM_DOUBLE:
    if (dst_image->data_type == IM_FLOAT)
      DoNormalizeComp((double**)src_image->data, (double**)dst_image->data, src_image->count, src_image->depth);
    else
      DoNormalizeComp((double**)src_image->data, (double**)dst_image->data, src_image->count, src_image->depth);
    break;
  }
}

template <class T> 
static void DoReplaceColor(T *src_data, T *dst_data, int count, int depth, double* src_color, double* dst_color)
{
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    int d, equal = 1;
    for (d = 0; d < depth; d++)
    {
      if (src_data[i+d*count] != (T)src_color[d])
      {
        equal = 0;
        break;
      }
    }

    for (d = 0; d < depth; d++)
    {
      if (equal)
        dst_data[i+d*count] = (T)dst_color[d];
      else
        dst_data[i+d*count] = src_data[i+d*count];
    }
  }
}

void imProcessReplaceColor(const imImage* src_image, imImage* dst_image, double* src_color, double* dst_color)
{
  switch(src_image->data_type)
  {
  case IM_BYTE:
    DoReplaceColor((imbyte*)src_image->data[0],   (imbyte*)dst_image->data[0],   src_image->count, src_image->depth, src_color, dst_color);
    break;                                                                                         
  case IM_SHORT:                                                                                   
    DoReplaceColor((short*)src_image->data[0], (short*)dst_image->data[0], src_image->count, src_image->depth, src_color, dst_color);
    break;                                                                                         
  case IM_USHORT:                                                                                   
    DoReplaceColor((imushort*)src_image->data[0], (imushort*)dst_image->data[0], src_image->count, src_image->depth, src_color, dst_color);
    break;                                                                                         
  case IM_INT:                                                                                     
    DoReplaceColor((int*)src_image->data[0],      (int*)dst_image->data[0],      src_image->count, src_image->depth, src_color, dst_color);
    break;                                                                                         
  case IM_FLOAT:                                                                                   
    DoReplaceColor((float*)src_image->data[0],    (float*)dst_image->data[0],    src_image->count, src_image->depth, src_color, dst_color);
    break;                                                                                
  case IM_DOUBLE:
    DoReplaceColor((double*)src_image->data[0], (double*)dst_image->data[0], src_image->count, src_image->depth, src_color, dst_color);
    break;
  }
}

template <class ST, class DT> 
static void DoSetAlphaColor(ST *src_data, DT *dst_data, int count, int depth, double* src_color, double dst_alpha)
{
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    int equal = 1;
    for (int d = 0; d < depth; d++)
    {
      if (src_data[d*count + i] != (ST)src_color[d])
      {
        equal = 0;
        break;
      }
    }

    if (equal)
      dst_data[i] = (DT)dst_alpha;
  }
}

void imProcessSetAlphaColor(const imImage* src_image, imImage* dst_image, double* src_color, double dst_alpha)
{
  int a = 0; // dst_image is a mask to be used as alpha
  if (dst_image->has_alpha)
    a = dst_image->depth; // Index of the alpha channel

  switch(src_image->data_type)
  {
  case IM_BYTE:
    switch(dst_image->data_type)
    {
    case IM_BYTE:
      DoSetAlphaColor((imbyte*)src_image->data[0],   (imbyte*)dst_image->data[a],   src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_SHORT:                                                                                   
      DoSetAlphaColor((imbyte*)src_image->data[0], (short*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_USHORT:                                                                                   
      DoSetAlphaColor((imbyte*)src_image->data[0], (imushort*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_INT:                                                                                    
      DoSetAlphaColor((imbyte*)src_image->data[0],      (int*)dst_image->data[a],      src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_FLOAT:                                                                                   
      DoSetAlphaColor((imbyte*)src_image->data[0],    (float*)dst_image->data[a],    src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                
    case IM_DOUBLE:
      DoSetAlphaColor((imbyte*)src_image->data[0], (double*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    }
    break;                                                                                        
  case IM_SHORT:                                                                                   
    switch(dst_image->data_type)
    {
    case IM_BYTE:
      DoSetAlphaColor((short*)src_image->data[0],   (imbyte*)dst_image->data[a],   src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_SHORT:                                                                                   
      DoSetAlphaColor((short*)src_image->data[0], (short*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_USHORT:                                                                                   
      DoSetAlphaColor((short*)src_image->data[0], (imushort*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_INT:                                                                                    
      DoSetAlphaColor((short*)src_image->data[0],      (int*)dst_image->data[a],      src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_FLOAT:                                                                                   
      DoSetAlphaColor((short*)src_image->data[0],    (float*)dst_image->data[a],    src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                
    case IM_DOUBLE:
      DoSetAlphaColor((short*)src_image->data[0], (double*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    }
    break;                                                                                        
  case IM_USHORT:                                                                                   
    switch(dst_image->data_type)
    {
    case IM_BYTE:
      DoSetAlphaColor((imushort*)src_image->data[0],   (imbyte*)dst_image->data[a],   src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_SHORT:                                                                                   
      DoSetAlphaColor((imushort*)src_image->data[0], (short*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_USHORT:                                                                                   
      DoSetAlphaColor((imushort*)src_image->data[0], (imushort*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_INT:                                                                                    
      DoSetAlphaColor((imushort*)src_image->data[0],      (int*)dst_image->data[a],      src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_FLOAT:                                                                                   
      DoSetAlphaColor((imushort*)src_image->data[0],    (float*)dst_image->data[a],    src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                
    case IM_DOUBLE:
      DoSetAlphaColor((imushort*)src_image->data[0], (double*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    }
    break;                                                                                        
  case IM_INT:                                                                                    
    switch(dst_image->data_type)
    {
    case IM_BYTE:
      DoSetAlphaColor((int*)src_image->data[0],   (imbyte*)dst_image->data[a],   src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_SHORT:                                                                                   
      DoSetAlphaColor((int*)src_image->data[0], (short*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_USHORT:                                                                                   
      DoSetAlphaColor((int*)src_image->data[0], (imushort*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_INT:                                                                                    
      DoSetAlphaColor((int*)src_image->data[0],      (int*)dst_image->data[a],      src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_FLOAT:                                                                                   
      DoSetAlphaColor((int*)src_image->data[0],    (float*)dst_image->data[a],    src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                
    case IM_DOUBLE:
      DoSetAlphaColor((int*)src_image->data[0], (double*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    }
    break;                                                                                        
  case IM_FLOAT:                                                                                   
    switch(dst_image->data_type)
    {
    case IM_BYTE:
      DoSetAlphaColor((float*)src_image->data[0],   (imbyte*)dst_image->data[a],   src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_SHORT:                                                                                   
      DoSetAlphaColor((float*)src_image->data[0], (short*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_USHORT:                                                                                   
      DoSetAlphaColor((float*)src_image->data[0], (imushort*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_INT:                                                                                    
      DoSetAlphaColor((float*)src_image->data[0],      (int*)dst_image->data[a],      src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                        
    case IM_FLOAT:                                                                                   
      DoSetAlphaColor((float*)src_image->data[0],    (float*)dst_image->data[a],    src_image->count, src_image->depth, src_color, dst_alpha);
      break;                                                                                
    case IM_DOUBLE:
      DoSetAlphaColor((double*)src_image->data[0], (double*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    }
    break;                                                                                
  case IM_DOUBLE:
    switch (dst_image->data_type)
    {
    case IM_BYTE:
      DoSetAlphaColor((double*)src_image->data[0], (imbyte*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    case IM_SHORT:
      DoSetAlphaColor((double*)src_image->data[0], (short*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    case IM_USHORT:
      DoSetAlphaColor((double*)src_image->data[0], (imushort*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    case IM_INT:
      DoSetAlphaColor((double*)src_image->data[0], (int*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    case IM_FLOAT:
      DoSetAlphaColor((double*)src_image->data[0], (float*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    case IM_DOUBLE:
      DoSetAlphaColor((double*)src_image->data[0], (double*)dst_image->data[a], src_image->count, src_image->depth, src_color, dst_alpha);
      break;
    }
    break;
  }
}

template <class T>
static void DoFixBGR(T** src_data, T** dst_data, int count)
{
  T* src_b = src_data[0];
  T* src_g = src_data[1];
  T* src_r = src_data[2];
  T* dst_r = dst_data[0];
  T* dst_g = dst_data[1];
  T* dst_b = dst_data[2];

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    // tem values to allow in-place operation
    T b = src_b[i];
    T g = src_g[i];
    T r = src_r[i];

    dst_r[i] = r;
    dst_g[i] = g;
    dst_b[i] = b;
  }
}

void imProcessFixBGR(const imImage* src_image, imImage* dst_image)
{
  switch (src_image->data_type)
  {
  case IM_BYTE:
    DoFixBGR((imbyte**)src_image->data, (imbyte**)dst_image->data, src_image->count);
    break;
  case IM_SHORT:
    DoFixBGR((short**)src_image->data, (short**)dst_image->data, src_image->count);
    break;
  case IM_USHORT:
    DoFixBGR((imushort**)src_image->data, (imushort**)dst_image->data, src_image->count);
    break;
  case IM_INT:
    DoFixBGR((int**)src_image->data, (int**)dst_image->data, src_image->count);
    break;
  case IM_FLOAT:
    DoFixBGR((float**)src_image->data, (float**)dst_image->data, src_image->count);
    break;
  case IM_DOUBLE:
    DoFixBGR((double**)src_image->data, (double**)dst_image->data, src_image->count);
    break;
  }
}

static double iColorNormHue(double H)
{
  while (H < 0.0)
    H += 360;

  if (H > 360.0)
    H = fmod(H, 360.0);

  return H;
}

template <class T>
static void DoSelectHue(T** src_data, T** dst_data, int count, double hue_start, double hue_end)
{
  T *dst_red = dst_data[0],
    *dst_green = dst_data[1],
    *dst_blue = dst_data[2];
  T *src_red = src_data[0],
    *src_green = src_data[1],
    *src_blue = src_data[2];
  T min, max;
  double r, g, b;
  double H, S, I;

  imMinMaxType(src_data[0], count * 3, min, max);

  hue_start = iColorNormHue(hue_start);
  hue_end = iColorNormHue(hue_end);
  int interval360 = 0;
  if (hue_start > hue_end)
    interval360 = 1;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    r = (double)(src_red[i] - min) / (double)(max - min);  // now 0 <= norm <= 1
    g = (double)(src_green[i] - min) / (double)(max - min);  // now 0 <= norm <= 1
    b = (double)(src_blue[i] - min) / (double)(max - min);  // now 0 <= norm <= 1

    imColorRGB2HSI(r, g, b, &H, &S, &I);

    if ((!interval360 && (H < hue_start || H > hue_end)) ||
        (interval360 && (H < hue_start && H > hue_end)))
    {
      dst_red[i] = (T)min;
      dst_green[i] = (T)min;
      dst_blue[i] = (T)min;
    }
    else
    {
      dst_red[i] = src_red[i];
      dst_green[i] = src_green[i];
      dst_blue[i] = src_blue[i];
    }
  }
}

void imProcessSelectHue(const imImage* src_image, imImage* dst_image, double hue_start, double hue_end)
{
  switch (src_image->data_type)
  {
  case IM_BYTE:
    DoSelectHue((imbyte**)src_image->data, (imbyte**)dst_image->data, src_image->count, hue_start, hue_end);
    break;
  case IM_SHORT:
    DoSelectHue((short**)src_image->data, (short**)dst_image->data, src_image->count, hue_start, hue_end);
    break;
  case IM_USHORT:
    DoSelectHue((imushort**)src_image->data, (imushort**)dst_image->data, src_image->count, hue_start, hue_end);
    break;
  case IM_FLOAT:
    DoSelectHue((float**)src_image->data, (float**)dst_image->data, src_image->count, hue_start, hue_end);
    break;
  case IM_DOUBLE:
    DoSelectHue((double**)src_image->data, (double**)dst_image->data, src_image->count, hue_start, hue_end);
    break;
  }
}

template <class T>
static void DoSelectHSI(T** src_data, T** dst_data, int count, double hue_start, double hue_end, double sat_start, double sat_end, double int_start, double int_end)
{
  T *dst_red = dst_data[0],
    *dst_green = dst_data[1],
    *dst_blue = dst_data[2];
  T *src_red = src_data[0],
    *src_green = src_data[1],
    *src_blue = src_data[2];
  T min, max;
  double r, g, b;
  double H, S, I;

  imMinMaxType(src_data[0], count * 3, min, max);

  hue_start = iColorNormHue(hue_start);
  hue_end = iColorNormHue(hue_end);
  int interval360 = 0;
  if (hue_start > hue_end)
    interval360 = 1;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    r = (double)(src_red[i] - min) / (double)(max - min);  // now 0 <= norm <= 1
    g = (double)(src_green[i] - min) / (double)(max - min);  // now 0 <= norm <= 1
    b = (double)(src_blue[i] - min) / (double)(max - min);  // now 0 <= norm <= 1

    imColorRGB2HSI(r, g, b, &H, &S, &I);

    if (((!interval360 && (H < hue_start || H > hue_end)) || (interval360 && (H < hue_start && H > hue_end))) ||
        (S < sat_start || S > sat_end) ||
        (I < int_start || I > int_end))
    {
      dst_red[i] = (T)min;
      dst_green[i] = (T)min;
      dst_blue[i] = (T)min;
    }
    else
    {
      dst_red[i] = src_red[i];
      dst_green[i] = src_green[i];
      dst_blue[i] = src_blue[i];
    }
  }
}

void imProcessSelectHSI(const imImage* src_image, imImage* dst_image, double hue_start, double hue_end, double sat_start, double sat_end, double int_start, double int_end)
{
  switch (src_image->data_type)
  {
  case IM_BYTE:
    DoSelectHSI((imbyte**)src_image->data, (imbyte**)dst_image->data, src_image->count, hue_start, hue_end, sat_start, sat_end, int_start, int_end);
    break;
  case IM_SHORT:
    DoSelectHSI((short**)src_image->data, (short**)dst_image->data, src_image->count, hue_start, hue_end, sat_start, sat_end, int_start, int_end);
    break;
  case IM_USHORT:
    DoSelectHSI((imushort**)src_image->data, (imushort**)dst_image->data, src_image->count, hue_start, hue_end, sat_start, sat_end, int_start, int_end);
    break;
  case IM_FLOAT:
    DoSelectHSI((float**)src_image->data, (float**)dst_image->data, src_image->count, hue_start, hue_end, sat_start, sat_end, int_start, int_end);
    break;
  case IM_DOUBLE:
    DoSelectHSI((double**)src_image->data, (double**)dst_image->data, src_image->count, hue_start, hue_end, sat_start, sat_end, int_start, int_end);
    break;
  }
}

