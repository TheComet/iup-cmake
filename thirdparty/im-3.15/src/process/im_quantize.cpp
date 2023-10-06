/** \file
 * \brief Additional Image Quantization Operations
 *
 * See Copyright Notice in im_lib.h
 */


#include <im.h>
#include <im_util.h>
#include <im_palette.h>
#include <im_math.h>
#include <im_convert.h>

#include "im_process_counter.h"
#include "im_process_pnt.h"
#include "im_process_ana.h"

#include <stdlib.h>
#include <memory.h>


void imProcessQuantizeRGBUniform(const imImage* src_image, imImage* dst_image, int dither)
{
  imbyte *dst_map=(imbyte*)dst_image->data[0], 
         *red_map=(imbyte*)src_image->data[0],
         *green_map=(imbyte*)src_image->data[1],
         *blue_map=(imbyte*)src_image->data[2];

  imImageSetPalette(dst_image, imPaletteUniform(), 256);

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(src_image->height))
#endif
  for (int y = 0; y < src_image->height; y++)
  {
    int line_offset = y*src_image->width;
    for (int x = 0; x < src_image->width; x++)
    {                         
      int offset = line_offset+x;
      if (dither)
        dst_map[offset] = (imbyte)imPaletteUniformIndexHalftoned(imColorEncode(red_map[offset], green_map[offset], blue_map[offset]), x, y);
      else
        dst_map[offset] = (imbyte)imPaletteUniformIndex(imColorEncode(red_map[offset], green_map[offset], blue_map[offset]));
    }
  }
}

void imProcessQuantizeRGBMedianCut(const imImage* image, imImage* NewImage)
{
  imConvertRGB2Map(image->width, image->height, (imbyte*)image->data[0], (imbyte*)image->data[1], (imbyte*)image->data[2], (imbyte*)NewImage->data[0], NewImage->palette, &NewImage->palette_count);
}

void imProcessQuantizeGrayUniform(const imImage* src_image, imImage* dst_image, int grays)
{
  int i;

  imbyte *dst_map=(imbyte*)dst_image->data[0], 
         *src_map=(imbyte*)src_image->data[0];

  imbyte re_map[256];
  memset(re_map, 0, 256);

  double factor = (double)grays / 256.0;
  double factor256 = 256.0 / (double)grays;

  for (i = 0; i < 256; i++)  // for all src values
  {             
    int value = imResampleInt(i, factor);     // from 0-255 to 0-(grays-1)  => this will discart information
    if (dst_image->color_space != IM_MAP)
      value = imResampleInt(value, factor256);  // from 0-(grays-1) back to 0-255
    re_map[i] = (imbyte)IM_BYTECROP(value);
  }

  int total_count = src_image->count*src_image->depth;
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(total_count))
#endif
  for (i = 0; i < total_count; i++)
    dst_map[i] = re_map[src_map[i]];
}

void imProcessQuantizeGrayMedianCut(imImage* src_image, imImage* dst_image, int grays)
{
  int i;
  unsigned long histo[256];
  imbyte *dst_map = (imbyte*)dst_image->data[0],
    *src_map = (imbyte*)src_image->data[0];

  int size = src_image->width * src_image->height;
  imCalcByteHistogram(src_map, size, histo, 1);

  int cut_map[256];
  memset(cut_map, 0, 256 * sizeof(int));

  double offset = double(size) / (grays); // region size
  double cut = offset;                  // first region size
  int map = 0, start = 0, end = 255;

  while (histo[start] == 0) start++;         // ignores zeros at start and end
  while (histo[end - 1] == histo[255]) end--;

  cut_map[0] = start;

  for (i = start + 1; i < end; i++)
  {
    if (histo[i] >(unsigned long)cut)  // outsize the last region
    {
      map++;                            // Change region
      cut_map[map] = i;

      if (map == grays - 1)
      {
        cut_map[grays] = end + 1;          // Last region
        break;
      }

      offset = double(size - histo[i - 1]) / (grays - map);  // recalc region size
      cut = histo[i - 1] + offset;                       // next region size
    }
  }

  imbyte re_map[256];
  memset(re_map, 0, 256);
  imbyte remaped;

  imCalcByteHistogram(src_map, size, histo, 0);

  int j;
  unsigned long sum, N;
  for (i = 0; i < grays; i++)
  {
    N = 0; sum = 0;
    for (j = cut_map[i]; j < cut_map[i + 1]; j++)
    {
      sum += histo[j] * j;
      N += histo[j];
    }
    remaped = imbyte(sum / N);

    for (j = cut_map[i]; j < cut_map[i + 1]; j++)
      re_map[j] = remaped;
  }

  for (i = 0; i < size; i++)
    dst_map[i] = re_map[src_map[i]];
}

