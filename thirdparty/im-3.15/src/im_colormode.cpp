/** \file
 * \brief Color Mode Utilities
 *
 * See Copyright Notice in im_lib.h
 */


#include <stdlib.h>
#include <memory.h>
#include <string.h>

#include "im.h"
#include "im_util.h"


const char* imColorModeComponentName(int color_mode, int component)
{
  const char* comp_name[4];
  int color_space = imColorModeSpace(color_mode);

  switch (color_space)
  {
  default:
    comp_name[0] = "Red";
    comp_name[1] = "Green";
    comp_name[2] = "Blue";
    break;
  case IM_CMYK:
    comp_name[0] = "Cyan";
    comp_name[1] = "Magenta";
    comp_name[2] = "Yellow";
    comp_name[3] = "Black";
    break;
  case IM_YCBCR:
    comp_name[0] = "Y";
    comp_name[1] = "Cb";
    comp_name[2] = "Cr";
    break;
  case IM_XYZ:
    comp_name[0] = "X";
    comp_name[1] = "Y";
    comp_name[2] = "Z";
    break;
  case IM_LAB:
    comp_name[0] = "L";
    comp_name[1] = "a";
    comp_name[2] = "b";
    break;
  case IM_LUV:
    comp_name[0] = "L";
    comp_name[1] = "u";
    comp_name[2] = "v";
    break;
  case IM_MAP:    
    return "Index";
  case IM_GRAY:   
    return "Gray";
  case IM_BINARY: 
    return "Binary";
  }

  return comp_name[component];
}

const char* imColorModeSpaceName(int color_mode)
{
  int color_space = imColorModeSpace(color_mode);
  switch (color_space)
  {
  case IM_RGB:    return "RGB";
  case IM_MAP:    return "Map";
  case IM_GRAY:   return "Gray";
  case IM_BINARY: return "Binary";
  case IM_CMYK:   return "CMYK";
  case IM_YCBCR:  return "Y'CbCr";
  case IM_LAB:    return "CIE L*a*b*";
  case IM_LUV:    return "CIE L*u*v*";
  case IM_XYZ:    return "CIE XYZ";
  }

  return NULL;
}

int imColorModeDepth(int color_mode)
{
  int depth = 0;

  int color_space = imColorModeSpace(color_mode);
  switch (color_space)
  {
  case IM_GRAY:
  case IM_BINARY:
  case IM_MAP:   
    depth = 1; 
    break;
  case IM_CMYK:
    depth = 4; 
    break;
  default:
    depth = 3; 
    break;
  }

  if (imColorModeHasAlpha(color_mode))
    depth++;

  return depth;
}

int imColorModeToBitmap(int color_mode)
{
  int color_space = imColorModeSpace(color_mode);
  switch (color_space)
  {
  case IM_BINARY:
  case IM_GRAY:
  case IM_MAP:
    return color_space;
  default:
    return IM_RGB;
  }
}

int imColorModeIsBitmap(int color_mode, int data_type)
{
  if (imColorModeSpace(color_mode) == IM_BINARY || 
      imColorModeSpace(color_mode) == IM_MAP)
    return 1;

  if ((imColorModeSpace(color_mode) == IM_RGB || 
       imColorModeSpace(color_mode) == IM_GRAY) &&
      (data_type == IM_BYTE))
    return 1;

  return 0;
}
