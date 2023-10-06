/** \file
 * \brief Fast Fourier Transform using FFTW library
 *
 * See Copyright Notice in im_lib.h
 */

#include <im.h>
#include <im_util.h>
#include <im_complex.h>
#include <im_convert.h>

#include "im_process.h"

#include <stdlib.h>
#include <assert.h>
#include <memory.h>

#ifdef USE_FFTW3
#include "fftw3.h"
#else
#include "fftw.h"
#endif

/* FFTW 2.x can have float or double functions, not both.
   FFTW 3.x can have both. */

template <class T>
static void iCopyCol(imComplex<T> *map1, imComplex<T> *map2, int height, int width1, int width2)
{
  int i;
  for(i = 0; i < height; i++)
  {
    *map1 = *map2;
    map1 += width1;
    map2 += width2;
  }
}

template <class T>
static void iCenterFFT(imComplex<T> *map, int width, int height, int inverse)
{
  imComplex<T> *map1, *map2, *map3, *tmp;
  int i, half1_width, half2_width, half1_height, half2_height;

  if (inverse)
  {
    half1_width = width/2;
    half1_height = height/2;

    half2_width = (width+1)/2;
    half2_height = (height+1)/2;
  }
  else
  {
    half1_width = (width+1)/2;
    half1_height = (height+1)/2;

    half2_width = width/2;
    half2_height = height/2;
  }

  tmp = (imComplex<T>*)malloc(half1_width*sizeof(imComplex<T>));

  map1 = map;
  map2 = map + half1_width;
  map3 = map + half2_width;
  for(i = 0; i < height; i++)
  {
    memcpy(tmp, map1, half1_width*sizeof(imComplex<T>));
    memcpy(map1, map2, half2_width*sizeof(imComplex<T>));
    memcpy(map3, tmp, half1_width*sizeof(imComplex<T>));

    map1 += width;
    map2 += width;
    map3 += width;
  }

  free(tmp);

  tmp = (imComplex<T>*)malloc(half1_height*sizeof(imComplex<T>));

  map1 = map;
  map2 = map + half1_height*width;
  map3 = map + half2_height*width;
  for(i = 0; i < width; i++)
  {
    iCopyCol(tmp, map1, half1_height, 1, width);
    iCopyCol(map1, map2, half2_height, width, width);
    iCopyCol(map3, tmp, half1_height, width, 1);

    map1++;
    map2++;
    map3++;
  }

  free(tmp);
}

template <class T>
static void iNormalize(imComplex<T> *map, int width, int height, int normalize)
{
  T NM = (T)(width * height);
  int count = (int)(2 * NM);

  if (normalize == 1)
    NM = (T)sqrt(NM);

  T *fmap = (T*)map;
  for (int i = 0; i < count; i++)
    *fmap++ /= NM;
}

static void iDoFFT(void *map, int width, int height, int data_type, int inverse, int center, int normalize)
{
  if (inverse && center)
  {
    if (data_type == IM_CFLOAT)
      iCenterFFT((imComplex<float>*)map, width, height, inverse);
    else
      iCenterFFT((imComplex<double>*)map, width, height, inverse);
  }

#ifdef USE_FFTW3
  if (data_type == IM_CFLOAT)
  {
    fftwf_plan plan = fftwf_plan_dft_2d(height, width, 
                                        (fftwf_complex*)map, (fftwf_complex*)map, // in-place transform
                                        inverse?FFTW_BACKWARD:FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
  }
  else
  {
    fftw_plan plan = fftw_plan_dft_2d(height, width, 
                                      (fftw_complex*)map, (fftw_complex*)map, // in-place transform
                                      inverse ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
  }
#else
  if (data_type == IM_CFLOAT)
  {
    fftwnd_plan plan = fftw2d_create_plan(height, width, inverse ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    fftwnd(plan, 1, (FFTW_COMPLEX*)map, 1, 0, 0, 0, 0);
    fftwnd_destroy_plan(plan);
  }
#endif

  if (!inverse && center)
  {
    if (data_type == IM_CFLOAT)
      iCenterFFT((imComplex<float>*)map, width, height, inverse);
    else
      iCenterFFT((imComplex<double>*)map, width, height, inverse);
  }

  if (normalize)
  {
    if (data_type == IM_CFLOAT)
      iNormalize((imComplex<float>*)map, width, height, normalize);
    else
      iNormalize((imComplex<double>*)map, width, height, normalize);
  }
}

void imProcessSwapQuadrants(imImage* image, int inverse)
{
  for (int i = 0; i < image->depth; i++)
  {
    if (image->data_type == IM_CFLOAT)
      iCenterFFT((imComplex<float>*)image->data[i], image->width, image->height, inverse);
    else
      iCenterFFT((imComplex<double>*)image->data[i], image->width, image->height, inverse);
  }
}

void imProcessFFTraw(imImage* image, int inverse, int center, int normalize)
{
  for (int i = 0; i < image->depth; i++)
    iDoFFT(image->data[i], image->width, image->height, image->data_type, inverse, center, normalize);
}

void imProcessFFT(const imImage* src_image, imImage* dst_image)
{
  if (src_image->data_type != IM_CFLOAT && src_image->data_type != IM_CDOUBLE)
    imProcessConvertDataType(src_image, dst_image, 0, 0, 0, 0);
  else
    imImageCopy(src_image, dst_image);

  imProcessFFTraw(dst_image, 0, 1, 0); // forward, centered, unnormalized
}

void imProcessIFFT(const imImage* src_image, imImage* dst_image)
{
  imImageCopy(src_image, dst_image);

  imProcessFFTraw(dst_image, 1, 1, 2); // inverse, uncentered, double normalized
}

void imProcessCrossCorrelation(const imImage* src_image1, const imImage* src_image2, imImage* dst_image)
{
  imImage *tmp_image = imImageCreateBased(src_image2, -1, -1, -1, dst_image->data_type);
  if (!tmp_image) 
    return;

  if (src_image2->data_type != IM_CFLOAT && src_image2->data_type != IM_CDOUBLE)
    imProcessConvertDataType(src_image2, tmp_image, 0, 0, 0, 0);
  else
    imImageCopy(src_image2, tmp_image);

  if (src_image1->data_type != IM_CFLOAT && src_image1->data_type != IM_CDOUBLE)
    imProcessConvertDataType(src_image1, dst_image, 0, 0, 0, 0);
  else
    imImageCopy(src_image1, dst_image);

  imProcessFFTraw(tmp_image, 0, 1, 1);   // forward, centered, normalized
  imProcessFFTraw(dst_image, 0, 1, 1);

  imProcessMultiplyConj(dst_image, tmp_image, dst_image);

  imProcessFFTraw(dst_image, 1, 1, 1);   // inverse, uncentered, normalized
  imProcessSwapQuadrants(dst_image, 0);  // from origin to center

  imImageDestroy(tmp_image);
}

void imProcessAutoCorrelation(const imImage* src_image, imImage* dst_image)
{
  if (src_image->data_type != IM_CFLOAT && src_image->data_type != IM_CDOUBLE)
    imProcessConvertDataType(src_image, dst_image, 0, 0, 0, 0);
  else
    imImageCopy(src_image, dst_image);

  imProcessFFTraw(dst_image, 0, 0, 1);   // forward, at origin, normalized

  imProcessMultiplyConj(dst_image, dst_image, dst_image);

  imProcessFFTraw(dst_image, 1, 0, 1);   // inverse, at origin, normalized
  imProcessSwapQuadrants(dst_image, 0);  // from origin to center
}
