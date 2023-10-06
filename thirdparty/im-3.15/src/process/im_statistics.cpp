/** \file
 * \brief Image Statistics Calculations
 *
 * See Copyright Notice in im_lib.h
 */


#include <im.h>
#include <im_util.h>
#include <im_color.h>
#include <im_math_op.h>

#include "im_process_counter.h"
#include "im_process_ana.h"

#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include <stdio.h>


template <class T>
static int DoCalcHisto(T* map, int size, unsigned long* histo, int hcount, int cumulative, int shift, int counter, int width)
{
  memset(histo, 0, hcount * sizeof(unsigned long));

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(size))
#endif
  for (int i = 0; i < size; i++)
  {
    if (i % width == 0)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_BEGIN_PROCESSING;

    int index = map[i] + shift;
#ifdef _OPENMP
#pragma omp atomic
#endif
    histo[index]++;

    if (i % width == 0)
    {
      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_END_PROCESSING;
  }

  if (cumulative)
  {
    /* make cumulative histogram */
    for (int i = 1; i < hcount; i++)
      histo[i] += histo[i-1];
  }

  return processing;
}

void imCalcByteHistogram(const imbyte* map, int size, unsigned long* histo, int cumulative)
{
  DoCalcHisto(map, size, histo, 256, cumulative, 0, -1, size);
}

void imCalcUShortHistogram(const imushort* map, int size, unsigned long* histo, int cumulative)
{
  DoCalcHisto(map, size, histo, 65536, cumulative, 0, -1, size);
}

void imCalcShortHistogram(const short* map, int size, unsigned long* histo, int cumulative)
{
  DoCalcHisto(map, size, histo, 65536, cumulative, 32768, -1, size);
}

int imCalcHistogram(const imImage* src_image, unsigned long* histo, int plane, int cumulative)
{
  int ret = 0;
  int counter = imProcessCounterBegin("Histogram");
  imCounterTotal(counter, src_image->height, "Calculating...");

  switch (src_image->data_type)
  {
  case IM_BYTE:
    ret = DoCalcHisto(  (imbyte*)src_image->data[plane], src_image->count, histo, 256,   cumulative, 0, counter, src_image->width);
    break;
  case IM_SHORT:
    ret = DoCalcHisto(   (short*)src_image->data[plane], src_image->count, histo, 65536, cumulative, 32768, counter, src_image->width);
    break;
  case IM_USHORT:
    ret = DoCalcHisto((imushort*)src_image->data[plane], src_image->count, histo, 65536, cumulative, 0, counter, src_image->width);
    break;
  }

  imProcessCounterEnd(counter);
  return ret;
}

int imCalcGrayHistogram(const imImage* image, unsigned long* histo, int cumulative)
{
  int counter = imProcessCounterBegin("GrayHistogram");

  int hcount = imHistogramCount(image->data_type);

  IM_INT_PROCESSING;

  if (image->color_space == IM_GRAY)
  {
    if (!imCalcHistogram(image, histo, 0, cumulative))
    {
      imProcessCounterEnd(counter);
      return 0;
    }
  }
  else 
  {
    int i;
    memset(histo, 0, hcount * sizeof(unsigned long));
    imCounterTotal(counter, image->height, "Calculating...");

    if (image->color_space == IM_MAP || image->color_space == IM_BINARY)
    {
      imbyte* map = (imbyte*)image->data[0];
      imbyte gray_map[256], r, g, b;

      for (i = 0; i < image->palette_count; i++)
      {
        imColorDecode(&r, &g, &b, image->palette[i]);
        gray_map[i] = imColorRGB2Luma(r, g, b);
      }

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(image->count))
#endif
      for (i = 0; i < image->count; i++)
      {
        if (i % image->width == 0)
        {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
        }
        IM_BEGIN_PROCESSING;

        int index = gray_map[map[i]];
#ifdef _OPENMP
#pragma omp atomic
#endif
        histo[index]++;

        if (i % image->width == 0)
        {
          IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
        }
        IM_END_PROCESSING;
      }
    }
    else   // RGB
    {
      if (image->data_type == IM_USHORT)
      {
        imushort* r = (imushort*)image->data[0];
        imushort* g = (imushort*)image->data[1];
        imushort* b = (imushort*)image->data[2];

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(image->count))
#endif
        for (i = 0; i < image->count; i++)
        {
          if (i % image->width == 0)
          {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
          }
          IM_BEGIN_PROCESSING;

          imushort index = imColorRGB2Luma(*r++, *g++, *b++);
#ifdef _OPENMP
#pragma omp atomic
#endif
          histo[index]++;

          if (i % image->width == 0)
          {
            IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
          }
          IM_END_PROCESSING;
        }
      }
      else if (image->data_type == IM_SHORT)
      {
        short* r = (short*)image->data[0];
        short* g = (short*)image->data[1];
        short* b = (short*)image->data[2];

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(image->count))
#endif
        for (i = 0; i < image->count; i++)
        {
          if (i % image->width == 0)
          {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
          }
          IM_BEGIN_PROCESSING;

          int index = imColorRGB2Luma(*r++, *g++, *b++) + 32768;
#ifdef _OPENMP
#pragma omp atomic
#endif
          histo[index]++;

          if (i % image->width == 0)
          {
            IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
          }
          IM_END_PROCESSING;
        }
      }
      else
      {
        imbyte* r = (imbyte*)image->data[0];
        imbyte* g = (imbyte*)image->data[1];
        imbyte* b = (imbyte*)image->data[2];

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(image->count))
#endif
        for (i = 0; i < image->count; i++)
        {
          if (i % image->width == 0)
          {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
          }
          IM_BEGIN_PROCESSING;

          imbyte index = imColorRGB2Luma(*r++, *g++, *b++);
#ifdef _OPENMP
#pragma omp atomic
#endif
          histo[index]++;

          if (i % image->width == 0)
          {
            IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
          }
          IM_END_PROCESSING;
        }
      }
    }

    if (cumulative)
    {
      /* make cumulative histogram */
      for (i = 1; i < hcount; i++)
        histo[i] += histo[i-1];
    }
  }

  imProcessCounterEnd(counter);
  return processing;
}

static int count_map(const imImage* image, unsigned long *total_count)
{
  int hcount;
  unsigned long* histo = imHistogramNew(image->data_type, &hcount);

  if (!imCalcHistogram(image, histo, 0, 0))
  {
    imHistogramRelease(histo);
    return 0;
  }

  unsigned long numcolor = 0;

  for (int i = 0; i < hcount; i++)
  {             
    if(histo[i] != 0)
      numcolor++;
  }

  imHistogramRelease(histo);

  *total_count = numcolor;
  return 1;
}

static int count_comp(const imImage* image, int counter, int width, unsigned long *total_count)
{
  imbyte *count = (imbyte*)calloc(sizeof(imbyte), 1 << 21); /* (2^24)/8=2^21 ~ 2Mb - using a bit array */
  if (!count)
    return 0;

  imbyte *comp0 = (imbyte*)image->data[0];
  imbyte *comp1 = (imbyte*)image->data[1];
  imbyte *comp2 = (imbyte*)image->data[2];

  int index;
  unsigned long numcolor = 0;

  for(int i = 0; i < image->count; i++)
  {
    index = comp0[i] << 16 | comp1[i] << 8 | comp2[i];

    if(imDataBitGet(count, index) == 0)
      numcolor++;

    imDataBitSet(count, index, 1);

    if (i % width == 0)
    {
      if (!imCounterInc(counter))
        return 0;
    }
  }

  free(count);

  *total_count = numcolor;
  return 1;
}

int imCalcCountColors(const imImage* image, unsigned long* count)
{
  int ret = 0;
  int counter = imProcessCounterBegin("CountColors");

  if (imColorModeDepth(image->color_space) > 1)
  {
    imCounterTotal(counter, image->height, "Calculating...");
    ret = count_comp(image, counter, image->width, count);
  }
  else
    ret = count_map(image, count);

  imProcessCounterEnd(counter);
  return ret;
}

template <class T>
static int DoStats(T* data, int count, imStats* stats, int counter, int width)
{
  memset(stats, 0, sizeof(imStats));

  IM_INT_PROCESSING;

  unsigned long positive = 0;
  unsigned long negative = 0;
  unsigned long zeros = 0;
  double mean = 0;
  double stddev = 0;

  T min, max;
  imMinMax(data, count, min, max);

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count)) \
                         reduction (+:positive, negative, zeros, mean, stddev) 
#endif
  for (int i = 0; i < count; i++)
  {
    if (i % width == 0)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_BEGIN_PROCESSING;

    if (data[i] > 0)
      positive++;

    if (data[i] < 0)
      negative++;

    if (data[i] == 0)
      zeros++;

    mean += (double)data[i];
    stddev += ((double)data[i])*((double)data[i]);

    if (i % width == 0)
    {
      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_END_PROCESSING;
  }

  double dcount = (double)count;
  mean /= dcount;
  stddev = sqrt((stddev - dcount*mean*mean)/(dcount-1.0));

  stats->max = max;
  stats->min = min;
  stats->positive = positive;
  stats->negative = negative;
  stats->zeros = zeros;
  stats->mean = mean;
  stats->stddev = stddev;

  return processing;
}

int imCalcImageStatistics(const imImage* image, imStats* stats)
{
  int ret = 0;
  int counter = imProcessCounterBegin("ImageStatistics");
  imCounterTotal(counter, image->depth*image->height, "Calculating...");

  for (int i = 0; i < image->depth; i++)
  {
    switch(image->data_type)
    {
    case IM_BYTE:
      ret = DoStats((imbyte*)image->data[i], image->count, &stats[i], counter, image->width);
      break;                                                                                
    case IM_SHORT:                                                                           
      ret = DoStats((short*)image->data[i], image->count, &stats[i], counter, image->width);
      break;                                                                                
    case IM_USHORT:                                                                           
      ret = DoStats((imushort*)image->data[i], image->count, &stats[i], counter, image->width);
      break;                                                                                
    case IM_INT:                                                                           
      ret = DoStats((int*)image->data[i], image->count, &stats[i], counter, image->width);
      break;                                                                                
    case IM_FLOAT:                                                                           
      ret = DoStats((float*)image->data[i], image->count, &stats[i], counter, image->width);
      break;                                                                                
    case IM_DOUBLE:
      ret = DoStats((double*)image->data[i], image->count, &stats[i], counter, image->width);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);
  return ret;
}

int imCalcHistogramStatistics(const imImage* image, imStats* stats)
{
  int counter = imProcessCounterBegin("HistogramStatistics");

  int hcount;
  unsigned long* histo = imHistogramNew(image->data_type, &hcount);

  for (int d = 0; d < image->depth; d++)
  {
    if (!imCalcHistogram(image, histo, d, 0))
    {
      imHistogramRelease(histo);
      imProcessCounterEnd(counter);
      return 0;
    }

    DoStats((unsigned long*)histo, hcount, &stats[d], -1, hcount);
  }

  imHistogramRelease(histo);

  imProcessCounterEnd(counter);
  return 1;
}

int imCalcHistoImageStatistics(const imImage* image, int* median, int* mode)
{
  int counter = imProcessCounterBegin("HistoImageStatistics");

  int hcount;
  unsigned long* histo = imHistogramNew(image->data_type, &hcount);

  for (int d = 0; d < image->depth; d++)
  {
    int i;

    if (!imCalcHistogram(image, histo, d, 0))
    {
      imHistogramRelease(histo);
      imProcessCounterEnd(counter);
      return 0;
    }

    unsigned long half = image->count/2;
    unsigned long count = histo[0];
    for (i = 1; i < hcount; i++)
    {
      if (count > half)
      {
        median[d] = i-1;
        median[d] += imHistogramShift(image->data_type);
        break;
      }

      count += histo[i];
    }

    unsigned long max = histo[0];
    for (i = 1; i < hcount; i++)
    {
      if (max < histo[i])
        max = histo[i];
    }

    int found_mode = 0;
    for (i = 0; i < hcount; i++)
    {
      if (histo[i] == max)
      {
        if (found_mode)
        {
          mode[d] = -1;  // must be unique or it is invalid
          break;
        }

        mode[d] = i;
        mode[d] += imHistogramShift(image->data_type);
        found_mode = 1;
      }
    }
  }

  imHistogramRelease(histo);

  imProcessCounterEnd(counter);
  return 1;
}

int imCalcPercentMinMax(const imImage* image, double percent, int ignore_zero, int *min, int *max)
{
  int counter = imProcessCounterBegin("PercentMinMax");

  int zero = -imHistogramShift(image->data_type);
  int hcount;
  unsigned long* histo = imHistogramNew(image->data_type, &hcount);

  if (!imCalcGrayHistogram(image, histo, 0))
  {
    imProcessCounterEnd(counter);
    return 0;
  }

  unsigned long acum, cut = (unsigned long)((image->count * percent) / 100.0f);

  acum = 0;
  for (*min = 0; *min < hcount; (*min)++)
  {  
    if (ignore_zero)
    {
      if (*min != zero)
      {
        acum += histo[*min];
        if (acum > cut)
          break;
      }
    }
    else
    {
      acum += histo[*min];
      if (acum > cut)
        break;
    }
  }

  acum = 0;
  for (*max = hcount-1; *max > 0; (*max)--)
  {  
    acum += histo[*max];
    if (acum > cut)
      break;
  }

  if (*min >= *max)
  {
    *min = 0;
    *max = hcount-1;
  }

  *min += imHistogramShift(image->data_type);
  *max += imHistogramShift(image->data_type);

  imHistogramRelease(histo);

  imProcessCounterEnd(counter);
  return 1;
}

int imCalcSNR(const imImage* image, const imImage* noise_image, double *snr)
{
  int counter = imProcessCounterBegin("SNR");

  imStats stats[4];
  if (!imCalcImageStatistics((imImage*)image, stats))
  {
    imProcessCounterEnd(counter);
    return 0;
  }

  imStats noise_stats[4];
  if (!imCalcImageStatistics((imImage*)noise_image, noise_stats))
  {
    imProcessCounterEnd(counter);
    return 0;
  }

  if (image->color_space == IM_RGB)
  {
    noise_stats[0].stddev += noise_stats[1].stddev;
    noise_stats[0].stddev += noise_stats[2].stddev;
    noise_stats[0].stddev /= 3;
    stats[0].stddev += stats[1].stddev;
    stats[0].stddev += stats[2].stddev;
    stats[0].stddev /= 3;
  }

  if (noise_stats[0].stddev == 0)
    *snr = 0;
  else

    *snr = 20.*log10(stats[0].stddev / noise_stats[0].stddev);

  imProcessCounterEnd(counter);
  return 1;
}

template <class T> 
static int DoRMSOp(T *map1, T *map2, int count, double *rms, int counter, int width)
{
  double rmserror = 0;

  IM_INT_PROCESSING;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:rmserror) if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    if (i % width == 0)
    {
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_BEGIN_PROCESSING;

    double diff = double(map1[i] - map2[i]);
    rmserror += diff * diff;

    if (i % width == 0)
    {
      IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    }
    IM_END_PROCESSING;
  }

  *rms = rmserror;
  return processing;
}
  
int imCalcRMSError(const imImage* image1, const imImage* image2, double *rmserror)
{
  *rmserror = 0;

  int count = image1->count*image1->depth;

  int ret = 0;
  int counter = imProcessCounterBegin("RMSError");
  imCounterTotal(counter, image1->depth*image1->height, "Calculating...");

  switch(image1->data_type)
  {
  case IM_BYTE:
    ret = DoRMSOp((imbyte*)image1->data[0], (imbyte*)image2->data[0], count, rmserror, counter, image1->width);
    break;
  case IM_SHORT:
    ret = DoRMSOp((short*)image1->data[0], (short*)image2->data[0], count, rmserror, counter, image1->width);
    break;
  case IM_USHORT:
    ret = DoRMSOp((imushort*)image1->data[0], (imushort*)image2->data[0], count, rmserror, counter, image1->width);
    break;
  case IM_INT:
    ret = DoRMSOp((int*)image1->data[0], (int*)image2->data[0], count, rmserror, counter, image1->width);
    break;
  case IM_FLOAT:
    ret = DoRMSOp((float*)image1->data[0], (float*)image2->data[0], count, rmserror, counter, image1->width);
    break;
  case IM_CFLOAT:
    ret = DoRMSOp((float*)image1->data[0], (float*)image2->data[0], 2 * count, rmserror, counter, image1->width);
    break;
  case IM_DOUBLE:
    ret = DoRMSOp((double*)image1->data[0], (double*)image2->data[0], count, rmserror, counter, image1->width);
    break;
  case IM_CDOUBLE:
    ret = DoRMSOp((double*)image1->data[0], (double*)image2->data[0], 2 * count, rmserror, counter, image1->width);
    break;
  }

  *rmserror = sqrt((*rmserror) / double((count * image1->depth)));

  imProcessCounterEnd(counter);
  return ret;
}
