/** \file
 * \brief Binary Arithmetic Operations
 *
 * See Copyright Notice in im_lib.h
 */


#include <im.h>
#include <im_util.h>
#include <im_math.h>
#include <im_complex.h>

#include "im_process_counter.h"
#include "im_process_pnt.h"
#include "im_math_op.h"
#include "im_color.h"

#include <stdlib.h>
#include <memory.h>


template <class T>
static inline T backsub_op(const T& v1, const T& v2, const double& tol, int use_diff)
{
  T diff = diff_op(v1, v2);

  if (diff <= tol)
    return (T)0; // subtract background
  else
  {
    if (use_diff)
      return (T)diff; // show the difference
    else
      return (T)v1;  // show the original image
  }
}

template <class T>
static void DoBackSub(T *map1, T *map2, T *map, int count, double tol, int diff)
{
  for (int i = 0; i < count; i++)
    map[i] = backsub_op(map1[i], map2[i], tol, diff);
}

void imProcessBackSub(const imImage* src_image1, imImage* src_image2, imImage* dst_image, double tol, int diff)
{
  int count = src_image1->count;

  for (int i = 0; i < src_image1->depth; i++)
  {
    switch (src_image1->data_type)
    {
    case IM_BYTE:
      DoBackSub((imbyte*)src_image1->data[i], (imbyte*)src_image2->data[i], (imbyte*)dst_image->data[i], count, tol, diff);
      break;
    case IM_USHORT:
      DoBackSub((imushort*)src_image1->data[i], (imushort*)src_image2->data[i], (imushort*)dst_image->data[i], count, tol, diff);
      break;
    case IM_INT:
      DoBackSub((int*)src_image1->data[i], (int*)src_image2->data[i], (int*)dst_image->data[i], count, tol, diff);
      break;
    case IM_FLOAT:
      DoBackSub((float*)src_image1->data[i], (float*)src_image2->data[i], (float*)dst_image->data[i], count, tol, diff);
      break;
    case IM_CFLOAT:
      DoBackSub((imcfloat*)src_image1->data[i], (imcfloat*)src_image2->data[i], (imcfloat*)dst_image->data[i], count, tol, diff);
      break;
    case IM_DOUBLE:
      DoBackSub((double*)src_image1->data[i], (double*)src_image2->data[i], (double*)dst_image->data[i], count, tol, diff);
      break;
    case IM_CDOUBLE:
      DoBackSub((imcdouble*)src_image1->data[i], (imcdouble*)src_image2->data[i], (imcdouble*)dst_image->data[i], count, tol, diff);
      break;
    }
  }
}

template <class T1, class T2, class T3> 
static void DoBinaryOp(T1 *map1, T2 *map2, T3 *map, int count, int op)
{
  int i;

  switch(op)
  {
  case IM_BIN_ADD:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif 
    for (i = 0; i < count; i++)
      map[i] = add_op((T3)map1[i], (T3)map2[i]);
    break;
  case IM_BIN_SUB:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = sub_op((T3)map1[i], (T3)map2[i]);
    break;
  case IM_BIN_MUL:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = mul_op((T3)map1[i], (T3)map2[i]);
    break;
  case IM_BIN_DIV:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = div_op((T3)map1[i], (T3)map2[i]);
    break;
  case IM_BIN_DIFF:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = diff_op((T3)map1[i], (T3)map2[i]);
    break;
  case IM_BIN_MIN:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = min_op((T3)map1[i], (T3)map2[i]);
    break;
  case IM_BIN_MAX:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = max_op((T3)map1[i], (T3)map2[i]);
    break;
  case IM_BIN_POW:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = pow_op((T3)map1[i], (T3)map2[i]);
    break;
  }
}

static void DoBinaryOpByte(imbyte *map1, imbyte *map2, imbyte *map, int count, int op)
{
  int i;

  switch(op)
  {
  case IM_BIN_ADD:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(add_op((int)map1[i], (int)map2[i]));
    break;
  case IM_BIN_SUB:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(sub_op((int)map1[i], (int)map2[i]));
    break;
  case IM_BIN_MUL:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(mul_op((int)map1[i], (int)map2[i]));
    break;
  case IM_BIN_DIV:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(div_op((int)map1[i], (int)map2[i]));
    break;
  case IM_BIN_DIFF:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(diff_op((int)map1[i], (int)map2[i]));
    break;
  case IM_BIN_MIN:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(min_op((int)map1[i], (int)map2[i]));
    break;
  case IM_BIN_MAX:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(max_op((int)map1[i], (int)map2[i]));
    break;
  case IM_BIN_POW:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(pow_op((int)map1[i], (int)map2[i]));
    break;
  }
}

template <class T>
static void DoBinaryOpCpxReal(imComplex<T> *map1, T *map2, imComplex<T> *map, int count, int op)
{
  int i;

  switch(op)
  {
  case IM_BIN_ADD:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = add_op(map1[i], map2[i]);
    break;
  case IM_BIN_SUB:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = sub_op(map1[i], map2[i]);
    break;
  case IM_BIN_MUL:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = mul_op(map1[i], map2[i]);
    break;
  case IM_BIN_DIV:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = div_op(map1[i], map2[i]);
    break;
  case IM_BIN_DIFF:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = diff_op(map1[i], map2[i]);
    break;
  case IM_BIN_MIN:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = min_op(map1[i], map2[i]);
    break;
  case IM_BIN_MAX:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = max_op(map1[i], map2[i]);
    break;
  case IM_BIN_POW:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = pow_op(map1[i], map2[i]);
    break;
  }
}

void imProcessArithmeticOp(const imImage* src_image1, const imImage* src_image2, imImage* dst_image, int op)
{
  int count = src_image1->count*src_image1->depth;  /* do NOT include alpha here */

  switch(src_image1->data_type)
  {
  case IM_BYTE:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryOp((imbyte*)src_image1->data[0], (imbyte*)src_image2->data[0], (double*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_FLOAT)
      DoBinaryOp((imbyte*)src_image1->data[0], (imbyte*)src_image2->data[0], (float*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_SHORT)
      DoBinaryOp((imbyte*)src_image1->data[0], (imbyte*)src_image2->data[0], (short*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_USHORT)
      DoBinaryOp((imbyte*)src_image1->data[0], (imbyte*)src_image2->data[0], (imushort*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_INT)
      DoBinaryOp((imbyte*)src_image1->data[0], (imbyte*)src_image2->data[0], (int*)dst_image->data[0], count, op);
    else
      DoBinaryOpByte((imbyte*)src_image1->data[0], (imbyte*)src_image2->data[0], (imbyte*)dst_image->data[0], count, op);
    break;
  case IM_SHORT:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryOp((short*)src_image1->data[0], (short*)src_image2->data[0], (double*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_FLOAT)
      DoBinaryOp((short*)src_image1->data[0], (short*)src_image2->data[0], (float*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_INT)
      DoBinaryOp((short*)src_image1->data[0], (short*)src_image2->data[0], (int*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_USHORT)
      DoBinaryOp((short*)src_image1->data[0], (short*)src_image2->data[0], (imushort*)dst_image->data[0], count, op);
    else
      DoBinaryOp((short*)src_image1->data[0], (short*)src_image2->data[0], (short*)dst_image->data[0], count, op);
    break;
  case IM_USHORT:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryOp((imushort*)src_image1->data[0], (imushort*)src_image2->data[0], (double*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_FLOAT)
      DoBinaryOp((imushort*)src_image1->data[0], (imushort*)src_image2->data[0], (float*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_INT)
      DoBinaryOp((imushort*)src_image1->data[0], (imushort*)src_image2->data[0], (int*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_SHORT)
      DoBinaryOp((imushort*)src_image1->data[0], (imushort*)src_image2->data[0], (short*)dst_image->data[0], count, op);
    else
      DoBinaryOp((imushort*)src_image1->data[0], (imushort*)src_image2->data[0], (imushort*)dst_image->data[0], count, op);
    break;
  case IM_INT:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryOp((int*)src_image1->data[0], (int*)src_image2->data[0], (double*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_FLOAT)
      DoBinaryOp((int*)src_image1->data[0], (int*)src_image2->data[0], (float*)dst_image->data[0], count, op);
    else
      DoBinaryOp((int*)src_image1->data[0], (int*)src_image2->data[0], (int*)dst_image->data[0], count, op);
    break;
  case IM_FLOAT:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryOp((float*)src_image1->data[0], (float*)src_image2->data[0], (double*)dst_image->data[0], count, op);
    else
      DoBinaryOp((float*)src_image1->data[0], (float*)src_image2->data[0], (float*)dst_image->data[0], count, op);
    break;
  case IM_CFLOAT:
    if (src_image2->data_type == IM_FLOAT)
      DoBinaryOpCpxReal((imcfloat*)src_image1->data[0], (float*)src_image2->data[0], (imcfloat*)dst_image->data[0], count, op);
    else
      DoBinaryOp((imcfloat*)src_image1->data[0], (imcfloat*)src_image2->data[0], (imcfloat*)dst_image->data[0], count, op);
    break;
  case IM_DOUBLE:
    if (dst_image->data_type == IM_FLOAT)
      DoBinaryOp((double*)src_image1->data[0], (double*)src_image2->data[0], (float*)dst_image->data[0], count, op);
    else
      DoBinaryOp((double*)src_image1->data[0], (double*)src_image2->data[0], (double*)dst_image->data[0], count, op);
    break;
  case IM_CDOUBLE:
    if (src_image2->data_type == IM_DOUBLE)
      DoBinaryOpCpxReal((imcdouble*)src_image1->data[0], (double*)src_image2->data[0], (imcdouble*)dst_image->data[0], count, op);
    else
      DoBinaryOp((imcdouble*)src_image1->data[0], (imcdouble*)src_image2->data[0], (imcdouble*)dst_image->data[0], count, op);
    break;
  }
}

template <class T>
static inline imComplex<T> blend_op(const imComplex<T>& v1, const imComplex<T>& v2, const double& alpha)
{
  return v1*T(alpha) + v2*T(1.0f - alpha);
}

template <class T>
static inline T blend_op(const T& v1, const T& v2, const double& alpha)
{
  return (T)(v1*alpha + v2*(1.0f - alpha));
}

template <class T> 
static void DoBlendConst(T *map1, T *map2, T *map, int count, double alpha)
{
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
    map[i] = blend_op(map1[i], map2[i], alpha);
}

void imProcessBlendConst(const imImage* src_image1, const imImage* src_image2, imImage* dst_image, double alpha)
{
  int count = src_image1->count*src_image1->depth;

  switch(src_image1->data_type)
  {
  case IM_BYTE:
    DoBlendConst((imbyte*)src_image1->data[0], (imbyte*)src_image2->data[0], (imbyte*)dst_image->data[0], count, alpha);
    break;
  case IM_SHORT:
    DoBlendConst((short*)src_image1->data[0], (short*)src_image2->data[0], (short*)dst_image->data[0], count, alpha);
    break;
  case IM_USHORT:
    DoBlendConst((imushort*)src_image1->data[0], (imushort*)src_image2->data[0], (imushort*)dst_image->data[0], count, alpha);
    break;
  case IM_INT:
    DoBlendConst((int*)src_image1->data[0], (int*)src_image2->data[0], (int*)dst_image->data[0], count, alpha);
    break;
  case IM_FLOAT:
    DoBlendConst((float*)src_image1->data[0], (float*)src_image2->data[0], (float*)dst_image->data[0], count, alpha);
    break;
  case IM_CFLOAT:
    DoBlendConst((imcfloat*)src_image1->data[0], (imcfloat*)src_image2->data[0], (imcfloat*)dst_image->data[0], count, alpha);
    break;
  case IM_DOUBLE:
    DoBlendConst((double*)src_image1->data[0], (double*)src_image2->data[0], (double*)dst_image->data[0], count, alpha);
    break;
  case IM_CDOUBLE:
    DoBlendConst((imcdouble*)src_image1->data[0], (imcdouble*)src_image2->data[0], (imcdouble*)dst_image->data[0], count, alpha);
    break;
  }
}

template <class T, class TA> 
static void DoBlend(T *map1, T *map2, TA *alpha, T *map, int count, int alpha_count, double type_max)
{
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
  {
    int a = i % alpha_count;
    map[i] = blend_op(map1[i], map2[i], ((double)alpha[a]) / type_max);
  }
}

void imProcessBlend(const imImage* src_image1, const imImage* src_image2, const imImage* alpha, imImage* dst_image)
{
  int count = src_image1->count*src_image1->depth;
  double type_max = (double)imColorMax(src_image1->data_type);

  switch(src_image1->data_type)
  {
  case IM_BYTE:
    DoBlend((imbyte*)src_image1->data[0], (imbyte*)src_image2->data[0], (imbyte*)alpha->data[0], (imbyte*)dst_image->data[0], count, alpha->count, type_max);
    break;
  case IM_SHORT:
    DoBlend((short*)src_image1->data[0], (short*)src_image2->data[0], (short*)alpha->data[0], (short*)dst_image->data[0], count, alpha->count, type_max);
    break;
  case IM_USHORT:
    DoBlend((imushort*)src_image1->data[0], (imushort*)src_image2->data[0], (imushort*)alpha->data[0], (imushort*)dst_image->data[0], count, alpha->count, type_max);
    break;
  case IM_INT:
    DoBlend((int*)src_image1->data[0], (int*)src_image2->data[0], (int*)alpha->data[0], (int*)dst_image->data[0], count, alpha->count, type_max);
    break;
  case IM_FLOAT:
    DoBlend((float*)src_image1->data[0], (float*)src_image2->data[0], (float*)alpha->data[0], (float*)dst_image->data[0], count, alpha->count, type_max);
    break;
  case IM_CFLOAT:
    DoBlend((imcfloat*)src_image1->data[0], (imcfloat*)src_image2->data[0], (float*)alpha->data[0], (imcfloat*)dst_image->data[0], count, alpha->count, type_max);
    break;
  case IM_DOUBLE:
    DoBlend((double*)src_image1->data[0], (double*)src_image2->data[0], (double*)alpha->data[0], (double*)dst_image->data[0], count, alpha->count, type_max);
    break;
  case IM_CDOUBLE:
    DoBlend((imcdouble*)src_image1->data[0], (imcdouble*)src_image2->data[0], (double*)alpha->data[0], (imcdouble*)dst_image->data[0], count, alpha->count, type_max);
    break;
  }
}

#define COMPOSE_OVER(_SRC, _SRC_ALPHA, _DST, _TMP_MULTI, _TMP_ALPHA) (T)(((_SRC_ALPHA)*(_SRC) + (_TMP_MULTI)*(_DST)) / (_TMP_ALPHA))
#define ALPHA_BLEND(_src,_dst,_alpha) (T)(((_src) * (_alpha) + (_dst) * (max - (_alpha))) / max)

template <class T, class TA>
static inline T compose_op(const T& v1, const T& v2, const T& alpha1, const T& alpha2, const TA& max)
{
  if (alpha1 != max)   /* some transparency */                                                                     
  {                                                                                                                    
    if (alpha1 != 0) /* source not full transparent */                                                             
    {                                                                                                                  
      if (alpha2 == 0) /* target full transparent */                                                            
      {                                                                                                                
        return v1;                                                                                           
      }                                                                                                                
      else if (alpha2 == max) /* target opaque */                                                               
      {                                                                                                                
        return ALPHA_BLEND(v1, v2, alpha1);                                                   
      }                                                                                                                
      else /* (0<alpha2<max && 0<alpha1<max) target and source are semi-transparent */                      
      {                                                                                                                
        /* Closed Compositing SRC over DST  (see smith95a.pdf)        */                                               
        /* Colors NOT Premultiplied by Alpha                          */                                               
        /* DST = SRC * SRC_ALPHA + DST * DST_ALPHA * (1 - SRC_ALPHA)  */                                               
        /* DST_ALPHA = SRC_ALPHA + DST_ALPHA * (1 - SRC_ALPHA)        */                                               
        /* DST /= DST_ALPHA */                                                                                         
        TA _tmp_multi = alpha2 * (max - alpha1);                                                            
        TA _tmp_src_alpha = alpha1*max;                                                                           
        TA _tmp_alpha = _tmp_src_alpha + _tmp_multi;                                                                  
        return COMPOSE_OVER(v1, _tmp_src_alpha, v2, _tmp_multi, _tmp_alpha);                     
      }                                                                                                                
    }                                                                                                                  
    else  /* (alpha1 == 0) source full transparent */                                                              
    {                                                                                                                  
      return v2;                                                                                           
    }                                                                                                                  
  }                                                                                                                    
  else  /* (alpha1 == max) source has no alpha = opaque */                                                         
  {                                                                                                                    
    return v1;                                                                                               
  }                                                                                                                    
}

template <class T, class TA>
static inline T compose_alpha_op(const T& alpha1, const T& alpha2, const TA& max)
{
  if (alpha1 != max)   /* some transparency */                                                                     
  {                                                                                                                    
    if (alpha1 != 0) /* source not full transparent */                                                             
    {                                                                                                                  
      if (alpha2 == 0) /* target full transparent */                                                            
      {                                                                                                                
        return alpha1;                                                                                     
      }                                                                                                                
      else if (alpha2 == max) /* target opaque */                                                               
      {                                                                                                                
        /* alpha2 is not changed */                                                                              
        return alpha2;
      }                                                                                                                
      else /* (0<alpha2<max && 0<alpha1<max) target and source are semi-transparent */                      
      {                                                                                                                
        /* Closed Compositing SRC over DST  (see smith95a.pdf)        */                                               
        /* Colors NOT Premultiplied by Alpha                          */                                               
        /* DST = SRC * SRC_ALPHA + DST * DST_ALPHA * (1 - SRC_ALPHA)  */                                               
        /* DST_ALPHA = SRC_ALPHA + DST_ALPHA * (1 - SRC_ALPHA)        */                                               
        /* DST /= DST_ALPHA */                                                                                         
        TA _tmp_multi = alpha2 * (max - alpha1);                                                            
        TA _tmp_src_alpha = alpha1*max;                                                                           
        TA _tmp_alpha = _tmp_src_alpha + _tmp_multi;                                                                  
        return (T)(_tmp_alpha / max);
      }                                                                                                                
    }                                                                                                                  
    else  /* (alpha1 == 0) source full transparent */                                                              
    {                                                                                                                  
      /* alpha2 is not changed */                                                                                
      return alpha2;
    }                                                                                                                  
  }                                                                                                                    
  else  /* (alpha1 == max) source has no alpha = opaque */                                                         
  {                                                                                                                    
    return (unsigned char)max;   /* set target as opaque */                                                   
  }                                                                                                                    
}

template <class T, class TA> 
static void DoCompose(T *map1, T *map2, T *alpha1, T *alpha2, T *map, int count, TA max)
{
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
    map[i] = compose_op(map1[i], map2[i], alpha1[i], alpha2[i], max);
}

template <class T, class TA> 
static void DoComposeAlpha(T *alpha1, T *alpha2, T *dst_alpha, int count, TA max)
{
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
  for (int i = 0; i < count; i++)
    dst_alpha[i] = compose_alpha_op(alpha1[i], alpha2[i], max);
}

void imProcessCompose(const imImage* src_image1, const imImage* src_image2, imImage* dst_image)
{
  int count = src_image1->count, 
      src_alpha = src_image1->depth;
  int type_max = (int)imColorMax(src_image1->data_type);

  if (!src_image1->has_alpha || !src_image2->has_alpha || !dst_image->has_alpha)
    return;

  for (int i = 0; i < src_image1->depth; i++)
  {
    switch(src_image1->data_type)
    {
    case IM_BYTE:
      DoCompose((imbyte*)src_image1->data[i], (imbyte*)src_image2->data[i], (imbyte*)src_image1->data[src_alpha], (imbyte*)src_image2->data[src_alpha], (imbyte*)dst_image->data[i], count, type_max);
      break;
    case IM_SHORT:
      DoCompose((short*)src_image1->data[i], (short*)src_image2->data[i], (short*)src_image1->data[src_alpha], (short*)src_image2->data[src_alpha], (short*)dst_image->data[i], count, type_max);
      break;
    case IM_USHORT:
      DoCompose((imushort*)src_image1->data[i], (imushort*)src_image2->data[i], (imushort*)src_image1->data[src_alpha], (imushort*)src_image2->data[src_alpha], (imushort*)dst_image->data[i], count, type_max);
      break;
    case IM_INT:
      DoCompose((int*)src_image1->data[i], (int*)src_image2->data[i], (int*)src_image1->data[src_alpha], (int*)src_image2->data[src_alpha], (int*)dst_image->data[i], count, type_max);
      break;
    case IM_FLOAT:
      DoCompose((float*)src_image1->data[i], (float*)src_image2->data[i], (float*)src_image1->data[src_alpha], (float*)src_image2->data[src_alpha], (float*)dst_image->data[i], count, (float)type_max);
      break;
    case IM_DOUBLE:
      DoCompose((double*)src_image1->data[i], (double*)src_image2->data[i], (double*)src_image1->data[src_alpha], (double*)src_image2->data[src_alpha], (double*)dst_image->data[i], count, (double)type_max);
      break;
    }
  }

  /* one more for the alpha channel */
  switch(src_image1->data_type)
  {
  case IM_BYTE:
    DoComposeAlpha((imbyte*)src_image1->data[src_alpha], (imbyte*)src_image2->data[src_alpha], (imbyte*)dst_image->data[src_alpha], count, type_max);
    break;
  case IM_SHORT:
    DoComposeAlpha((short*)src_image1->data[src_alpha], (short*)src_image2->data[src_alpha], (short*)dst_image->data[src_alpha], count, type_max);
    break;
  case IM_USHORT:
    DoComposeAlpha((imushort*)src_image1->data[src_alpha], (imushort*)src_image2->data[src_alpha], (imushort*)dst_image->data[src_alpha], count, type_max);
    break;
  case IM_INT:
    DoComposeAlpha((int*)src_image1->data[src_alpha], (int*)src_image2->data[src_alpha], (int*)dst_image->data[src_alpha], count, type_max);
    break;
  case IM_FLOAT:
    DoComposeAlpha((float*)src_image1->data[src_alpha], (float*)src_image2->data[src_alpha], (float*)dst_image->data[src_alpha], count, (float)type_max);
    break;
  case IM_DOUBLE:
    DoComposeAlpha((double*)src_image1->data[src_alpha], (double*)src_image2->data[src_alpha], (double*)dst_image->data[src_alpha], count, (double)type_max);
    break;
  }
}

template <class T>
static void DoBinaryConstOpCpxReal(imComplex<T> *map1, T value, imComplex<T> *map, int count, int op)
{
  int i;

  switch(op)
  {
  case IM_BIN_ADD:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = add_op(map1[i], value);
    break;
  case IM_BIN_SUB:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = sub_op(map1[i], value);
    break;
  case IM_BIN_MUL:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = mul_op(map1[i], value);
    break;
  case IM_BIN_DIV:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = div_op(map1[i], value);
    break;
  case IM_BIN_DIFF:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = diff_op(map1[i], value);
    break;
  case IM_BIN_MIN:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = min_op(map1[i], value);
    break;
  case IM_BIN_MAX:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = max_op(map1[i], value);
    break;
  case IM_BIN_POW:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = pow_op(map1[i], value);
    break;
  }
}

template <class T1, class T2, class T3> 
static void DoBinaryConstOp(T1 *map1, T2 value, T3 *map, int count, int op)
{
  int i;

  switch(op)
  {
  case IM_BIN_ADD:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (T3)add_op((T2)map1[i], value);
    break;
  case IM_BIN_SUB:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (T3)sub_op((T2)map1[i], value);
    break;
  case IM_BIN_MUL:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (T3)mul_op((T2)map1[i], value);
    break;
  case IM_BIN_DIV:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (T3)div_op((T2)map1[i], value);
    break;
  case IM_BIN_DIFF:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (T3)diff_op((T2)map1[i], value);
    break;
  case IM_BIN_MIN:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (T3)min_op((T2)map1[i], value);
    break;
  case IM_BIN_MAX:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (T3)max_op((T2)map1[i], value);
    break;
  case IM_BIN_POW:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (T3)pow_op((T2)map1[i], value);
    break;
  }
}

template <class T1> 
static void DoBinaryConstOpByte(T1 *map1, int value, imbyte *map, int count, int op)
{
  int i;

  switch(op)
  {
  case IM_BIN_ADD:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(add_op((int)map1[i], value));
    break;
  case IM_BIN_SUB:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(sub_op((int)map1[i], value));
    break;
  case IM_BIN_MUL:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(mul_op((int)map1[i], value));
    break;
  case IM_BIN_DIV:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(div_op((int)map1[i], value));
    break;
  case IM_BIN_DIFF:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(diff_op((int)map1[i], value));
    break;
  case IM_BIN_MIN:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(min_op((int)map1[i], value));
    break;
  case IM_BIN_MAX:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(max_op((int)map1[i], value));
    break;
  case IM_BIN_POW:
#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(count))
#endif
    for (i = 0; i < count; i++)
      map[i] = (imbyte)crop_byte(pow_op((int)map1[i], value));
    break;
  }
}

void imProcessArithmeticConstOp(const imImage* src_image1, double value, imImage* dst_image, int op)
{
  int count = src_image1->count*src_image1->depth;  /* do NOT include alpha here */

  switch(src_image1->data_type)
  {
  case IM_BYTE:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryConstOp((imbyte*)src_image1->data[0], (double)value, (double*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_FLOAT)
      DoBinaryConstOp((imbyte*)src_image1->data[0], (float)value, (float*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_SHORT)
      DoBinaryConstOp((imbyte*)src_image1->data[0], (short)value, (short*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_USHORT)
      DoBinaryConstOp((imbyte*)src_image1->data[0], (imushort)value, (imushort*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_INT)
      DoBinaryConstOp((imbyte*)src_image1->data[0], (int)value, (int*)dst_image->data[0], count, op);
    else
      DoBinaryConstOpByte((imbyte*)src_image1->data[0], (int)value, (imbyte*)dst_image->data[0], count, op);
    break;
  case IM_SHORT:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryConstOp((short*)src_image1->data[0], (double)value, (double*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_FLOAT)
      DoBinaryConstOp((short*)src_image1->data[0], (float)value, (float*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_INT)
      DoBinaryConstOp((short*)src_image1->data[0], (int)value, (int*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_USHORT)
      DoBinaryConstOp((short*)src_image1->data[0], (int)value, (imushort*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_BYTE)
      DoBinaryConstOpByte((short*)src_image1->data[0], (int)value, (imbyte*)dst_image->data[0], count, op);
    else
      DoBinaryConstOp((short*)src_image1->data[0], (int)value, (short*)dst_image->data[0], count, op);
    break;
  case IM_USHORT:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryConstOp((imushort*)src_image1->data[0], (double)value, (double*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_FLOAT)
      DoBinaryConstOp((imushort*)src_image1->data[0], (float)value, (float*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_INT)
      DoBinaryConstOp((imushort*)src_image1->data[0], (int)value, (int*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_SHORT)
      DoBinaryConstOp((imushort*)src_image1->data[0], (int)value, (short*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_BYTE)
      DoBinaryConstOpByte((imushort*)src_image1->data[0], (int)value, (imbyte*)dst_image->data[0], count, op);
    else
      DoBinaryConstOp((imushort*)src_image1->data[0], (int)value, (imushort*)dst_image->data[0], count, op);
    break;
  case IM_INT:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryConstOp((int*)src_image1->data[0], (double)value, (double*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_FLOAT)
      DoBinaryConstOp((int*)src_image1->data[0], (float)value, (float*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_SHORT)
      DoBinaryConstOp((int*)src_image1->data[0], (int)value, (short*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_USHORT)
      DoBinaryConstOp((int*)src_image1->data[0], (int)value, (imushort*)dst_image->data[0], count, op);
    else if (dst_image->data_type == IM_BYTE)
      DoBinaryConstOpByte((int*)src_image1->data[0], (int)value, (imbyte*)dst_image->data[0], count, op);
    else
      DoBinaryConstOp((int*)src_image1->data[0], (int)value, (int*)dst_image->data[0], count, op);
    break;
  case IM_FLOAT:
    if (dst_image->data_type == IM_DOUBLE)
      DoBinaryConstOp((float*)src_image1->data[0], (double)value, (double*)dst_image->data[0], count, op);
    else
      DoBinaryConstOp((float*)src_image1->data[0], (float)value, (float*)dst_image->data[0], count, op);
    break;
  case IM_CFLOAT:
    DoBinaryConstOpCpxReal((imcfloat*)src_image1->data[0], (float)value, (imcfloat*)dst_image->data[0], count, op);
    break;
  case IM_DOUBLE:
    if (dst_image->data_type == IM_FLOAT)
      DoBinaryConstOp((double*)src_image1->data[0], (double)value, (float*)dst_image->data[0], count, op);
    else
      DoBinaryConstOp((double*)src_image1->data[0], (double)value, (double*)dst_image->data[0], count, op);
    break;
  case IM_CDOUBLE:
    DoBinaryConstOpCpxReal((imcdouble*)src_image1->data[0], (double)value, (imcdouble*)dst_image->data[0], count, op);
    break;
  }
}

void imProcessMultipleMean(const imImage** src_image_list, int src_image_count, imImage* dst_image)
{
  const imImage* image1 = src_image_list[0];
  imImage* aux_image = NULL;

  int data_type = image1->data_type;
  if (image1->data_type == IM_BYTE)
    data_type = IM_USHORT;

  imImage *acum_image = imImageCreateBased(image1, -1, -1, -1, data_type);
  if (!acum_image)
    return;

  if (image1->data_type == IM_BYTE)
    aux_image = imImageCreateBased(image1, -1, -1, -1, data_type);

  for(int i = 0; i < src_image_count; i++)
  {
    const imImage *image = src_image_list[i];
    if (aux_image)
    {
      imProcessUnArithmeticOp(image, aux_image, IM_UN_EQL);
      imProcessArithmeticOp(aux_image, acum_image, acum_image, IM_BIN_ADD);  /* acum_image += image */
    }
    else
      imProcessArithmeticOp(image, acum_image, acum_image, IM_BIN_ADD);  /* acum_image += image */
  }

  imProcessArithmeticConstOp(acum_image, double(src_image_count), dst_image, IM_BIN_DIV);

  imImageDestroy(acum_image);
  if (aux_image)
    imImageDestroy(aux_image);
}

void imProcessMultipleStdDev(const imImage** src_image_list, int src_image_count, const imImage *mean_image, imImage* dst_image)
{
  imImage* aux_image = imImageClone(dst_image);
  if (!aux_image)
    return;

  // sdtdev = sqrt( sum(sqr(x - m)) / N)

  // a = sum(sqr(x - m))
  for(int i = 0; i < src_image_count; i++)
  {
    // aux_image = image - mean_image
    imProcessArithmeticOp(src_image_list[i], mean_image, aux_image, IM_BIN_SUB);

    // aux_image = aux_image * aux_image
    imProcessUnArithmeticOp(aux_image, aux_image, IM_UN_SQR);

    // dst_image += aux_image
    imProcessArithmeticOp(aux_image, dst_image, dst_image, IM_BIN_ADD);
  }

  // dst_image = dst_image / src_image_count;
  imProcessArithmeticConstOp(dst_image, double(src_image_count), dst_image, IM_BIN_DIV);

  // dst_image = sqrt(dst_image);
  imProcessUnArithmeticOp(dst_image, dst_image, IM_UN_SQRT);

  imImageDestroy(aux_image);
}

static int compare_imFloat(const void *elem1, const void *elem2)
{
  double* v1 = (double*)elem1;
  double* v2 = (double*)elem2;

  if (*v1 < *v2)
    return -1;

  if (*v1 > *v2)
    return 1;

  return 0;
}

static int imMultiMedianFunc(const double* src_value, double *dst_value, double* params, void* userdata, int x, int y, int d, int src_count)
{
  (void)params;
  (void)userdata;
  (void)x;
  (void)y;
  (void)d;

  qsort((double*)src_value, src_count, sizeof(double), compare_imFloat);
  *dst_value = src_value[src_count / 2];

  return 1;
}

int imProcessMultipleMedian(const imImage** src_image_list, int src_image_count, imImage* dst_image)
{
  return imProcessMultiPointOp(src_image_list, src_image_count, dst_image, imMultiMedianFunc, NULL, NULL, "MultipleMedian");
}

template <class DT>
static double AutoCovCalc(int width, int height, DT *src_map, DT *mean_map, int x, int y, int count)
{
  double value = 0;
  int Ny = height - y;
  int Nx = width - x;
  int offset, offset1, line_offset, line_offset1;

  for (int i = 0; i < Ny; i++)
  {
    line_offset = width*i;
    line_offset1 = width*(i + y);

    for (int j = 0; j < Nx; j++)
    {
      offset = line_offset + j;
      offset1 = line_offset1 + (j + x);
      value += double(src_map[offset] - mean_map[offset]) * 
               double(src_map[offset1] - mean_map[offset1]);
    }
  }

  return value/(double)count;
}

template <class ST, class DT> 
static int doAutoCov(int width, int height, ST *src_map, ST *mean_map, DT *dst_map, int counter)
{
  int count = width*height;
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
      dst_map[line_offset + x] = (DT)AutoCovCalc(width, height, src_map, mean_map, x, y, count);
    }

    IM_COUNT_PROCESSING;
#ifdef _OPENMP
#pragma omp flush (processing)
#endif
    IM_END_PROCESSING;
  }

  return processing;
}

int imProcessAutoCovariance(const imImage* src_image, const imImage* mean_image, imImage* dst_image)
{
  int ret = 0;

  int counter = imProcessCounterBegin("AutoConvariance");
  imCounterTotal(counter, src_image->depth*src_image->height, "Processing...");

  for (int i = 0; i < src_image->depth; i++)
  {
    switch(src_image->data_type)
    {
    case IM_BYTE:
      if (dst_image->data_type == IM_FLOAT)
        ret = doAutoCov(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)mean_image->data[i], (float*)dst_image->data[i], counter);
      else
        ret = doAutoCov(src_image->width, src_image->height, (imbyte*)src_image->data[i], (imbyte*)mean_image->data[i], (double*)dst_image->data[i], counter);
      break;
    case IM_SHORT:
      if (dst_image->data_type == IM_FLOAT)
        ret = doAutoCov(src_image->width, src_image->height, (short*)src_image->data[i], (short*)mean_image->data[i], (float*)dst_image->data[i], counter);
      else
        ret = doAutoCov(src_image->width, src_image->height, (short*)src_image->data[i], (short*)mean_image->data[i], (double*)dst_image->data[i], counter);
      break;
    case IM_USHORT:
      if (dst_image->data_type == IM_FLOAT)
        ret = doAutoCov(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)mean_image->data[i], (float*)dst_image->data[i], counter);
      else
        ret = doAutoCov(src_image->width, src_image->height, (imushort*)src_image->data[i], (imushort*)mean_image->data[i], (double*)dst_image->data[i], counter);
      break;
    case IM_INT:
      if (dst_image->data_type == IM_FLOAT)
        ret = doAutoCov(src_image->width, src_image->height, (int*)src_image->data[i], (int*)mean_image->data[i], (float*)dst_image->data[i], counter);
      else
        ret = doAutoCov(src_image->width, src_image->height, (int*)src_image->data[i], (int*)mean_image->data[i], (double*)dst_image->data[i], counter);
      break;
    case IM_FLOAT:
      ret = doAutoCov(src_image->width, src_image->height, (float*)src_image->data[i], (float*)mean_image->data[i], (float*)dst_image->data[i], counter);
      break;
    case IM_DOUBLE:
      ret = doAutoCov(src_image->width, src_image->height, (double*)src_image->data[i], (double*)mean_image->data[i], (double*)dst_image->data[i], counter);
      break;
    }

    if (!ret)
      break;
  }

  imProcessCounterEnd(counter);

  return ret;
}

void imProcessMultiplyConj(const imImage* src_image1, const imImage* src_image2, imImage* dst_image)
{
  int total_count = src_image1->count*src_image1->depth;

  imcfloat* map = (imcfloat*)dst_image->data[0];
  imcfloat* map1 = (imcfloat*)src_image1->data[0];
  imcfloat* map2 = (imcfloat*)src_image2->data[0];

#ifdef _OPENMP
#pragma omp parallel for if (IM_OMP_MINCOUNT(total_count))
#endif
  for (int i = 0; i < total_count; i++)
  {
    imcfloat tmp; // this will allow an in-place operation
    tmp.real = map1[i].real * map2[i].real + map1[i].imag * map2[i].imag; 
    tmp.imag = map1[i].real * map2[i].imag - map1[i].imag * map2[i].real;
    map[i] = tmp;
  }
}
