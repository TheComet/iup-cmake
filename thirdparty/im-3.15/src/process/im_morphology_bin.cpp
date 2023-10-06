/** \file
 * \brief Morphology Operations for Binary Images
 *
 * See Copyright Notice in im_lib.h
 */


#include <im.h>
#include <im_util.h>

#include "im_process_counter.h"
#include "im_process_loc.h"
#include "im_process_pnt.h"

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <math.h>


static int DoBinMorphConvolve(imbyte *map, imbyte* new_map, int width, int height, const imImage* kernel, int counter, int hit_value, int miss_value)
{
  int kh2 = kernel->height/2;
  int kw2 = kernel->width/2;

  int* kernel_data = (int*)kernel->data[0];

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINHEIGHT(height))
#endif
  for(int j = 0; j < height; j++)
  {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_BEGIN_PROCESSING;

    int new_offset = j * width;

    for(int i = 0; i < width; i++)
    {
      int hit = 1;
    
      for(int y = -kh2; y <= kh2 && hit; y++)
      {
        int offset;
        int* kernel_line = kernel_data + (y+kh2)*kernel->width;

        if ((j + y < 0) ||       // pass the bottom border
            (j + y >= height))   // pass the top border
          offset = -1;
        else
          offset = (j + y) * width;

        for(int x = -kw2; x <= kw2; x++)
        {
          if ((offset == -1) ||
              (i + x < 0) ||     // pass the left border
              (i + x >= width))  // pass the right border
          {
            if(kernel_line[x+kw2] != -1 && kernel_line[x+kw2] != 0)  // 0 extension beyond borders
              hit = 0;
          }
          else
          {
            if(kernel_line[x+kw2] != -1 && kernel_line[x+kw2] != map[offset + (i + x)])
              hit = 0;
          }
        }
      }

      new_map[new_offset + i] = (imbyte)(hit? hit_value: miss_value);
    }    

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  return processing;
}

int imProcessBinMorphConvolve(const imImage* src_image, imImage* dst_image, const imImage *kernel, int hit_white, int iter)
{
  int j, ret = 0, hit_value, miss_value;
  void *tmp = NULL;
  int counter;

  if (hit_white)
  {
    hit_value = 1;
    miss_value = 0;
  }
  else
  {
    hit_value = 0;
    miss_value = 1;
  }

  if (iter > 1)
  {
    tmp = malloc(src_image->size);
    if (!tmp)
      return IM_ERR_MEM;
  }

  counter = imProcessCounterBegin("BinMorphConvolve");
  imCounterTotal(counter, src_image->height*iter, "Processing...");

  for (j = 0; j < iter; j++)
  {
    if (j == 0)
      ret = DoBinMorphConvolve((imbyte*)src_image->data[0], (imbyte*)dst_image->data[0], src_image->width, src_image->height, kernel, counter, hit_value, miss_value);
    else
    {
      memcpy(tmp, dst_image->data[0], src_image->size);
      ret = DoBinMorphConvolve((imbyte*)tmp, (imbyte*)dst_image->data[0], src_image->width, src_image->height, kernel, counter, hit_value, miss_value);
    }

    if (!ret) 
      break;
  }

  if (tmp) free(tmp);
  imProcessCounterEnd(counter);

  return ret;
}

int imProcessBinMorphErode(const imImage* src_image, imImage* dst_image, int kernel_size, int iter)
{
  imImage* kernel = imImageCreate(kernel_size, kernel_size, IM_GRAY, IM_INT);
  imImageSetAttribute(kernel, "Description", IM_BYTE, -1, (void*)"Erode");

  int* kernel_data = (int*)kernel->data[0];
  for(int i = 0; i < kernel->count; i++)
      kernel_data[i] = 1;

  int ret = imProcessBinMorphConvolve(src_image, dst_image, kernel, 1, iter);
  imImageDestroy(kernel);
  return ret;
}

int imProcessBinMorphDilate(const imImage* src_image, imImage* dst_image, int kernel_size, int iter)
{
  imImage* kernel = imImageCreate(kernel_size, kernel_size, IM_GRAY, IM_INT);
  imImageSetAttribute(kernel, "Description", IM_BYTE, -1, (void*)"Dilate");
  // Kernel is all zeros
  int ret = imProcessBinMorphConvolve(src_image, dst_image, kernel, 0, iter);
  imImageDestroy(kernel);
  return ret;
}

int imProcessBinMorphOpen(const imImage* src_image, imImage* dst_image, int kernel_size, int iter)
{
  imImage*temp = imImageClone(src_image);
  if (!temp)
    return 0;

  if (!imProcessBinMorphErode(src_image, temp, kernel_size, iter)) 
    {imImageDestroy(temp); return 0;}
  if (!imProcessBinMorphDilate(temp, dst_image, kernel_size, iter)) 
    {imImageDestroy(temp); return 0;}

  imImageDestroy(temp);
  return 1;
}

int imProcessBinMorphClose(const imImage* src_image, imImage* dst_image, int kernel_size, int iter)
{
  imImage*temp = imImageClone(src_image);
  if (!temp)
    return 0;

  if (!imProcessBinMorphDilate(src_image, temp, kernel_size, iter)) 
    {imImageDestroy(temp); return 0;}
  if (!imProcessBinMorphErode(temp, dst_image, kernel_size, iter)) 
    {imImageDestroy(temp); return 0;}

  imImageDestroy(temp);
  return 1;
}

int imProcessBinMorphOutline(const imImage* src_image, imImage* dst_image, int kernel_size, int iter)
{
  if (!imProcessBinMorphErode(src_image, dst_image, kernel_size, iter)) 
    return 0;
  imProcessArithmeticOp(src_image, dst_image, dst_image, IM_BIN_DIFF);
  return 1;
}
