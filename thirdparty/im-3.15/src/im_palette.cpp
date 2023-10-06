/** \file
 * \brief Palette Generators
 * Creates several standard palettes
 *
 * See Copyright Notice in im_lib.h
 */

#include "im.h"
#include "im_util.h"
#include "im_palette.h"
#include "im_colorhsi.h"
#include "im_color.h"

#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>


long* imPaletteNew(int count)
{
  return (long*)malloc(sizeof(long) * count);
}

void imPaletteRelease(long* palette)
{
  free(palette);
}

long* imPaletteDuplicate(const long* palette, int count)
{
  long* new_palette = imPaletteNew(count);

  for (int i = 0; i < count; i++)
    new_palette[i] = palette[i];

  return new_palette;
}

inline int iSqr(int x)
{
  return x*x;
}

inline int iAbs(int x)
{
  return x < 0? -x: x;
}

int imPaletteFindNearest(const long* palette, int palette_count, long color)
{
  assert(palette);
  assert(palette_count);

  int lSqrDiff, lBestDiff = (unsigned int)-1;
  int pIndex = -1;

  imbyte red1, green1, blue1;
  imColorDecode(&red1, &green1, &blue1, color);

  for (int lIndex = 0; lIndex < palette_count; lIndex++, palette++)
  {
    if (color == *palette)
      return lIndex;

    imbyte red2, green2, blue2;
    imColorDecode(&red2, &green2, &blue2, *palette);

    lSqrDiff = iSqr(red1 - red2) +
               iSqr(green1 - green2) +
               iSqr(blue1 - blue2);

    if (lSqrDiff < lBestDiff)
    {
      lBestDiff = lSqrDiff;
      pIndex = lIndex;
    }
  }

  return pIndex;
}

int imPaletteFindColor(const long* palette, int palette_count, long color, unsigned char tol)
{
  assert(palette);
  assert(palette_count);

  /* Divides in two section for faster results when Tolerance is 0.*/
  if (tol == 0)
  {
    for (int lIndex = 0; lIndex < palette_count; lIndex++, palette++)
    {
      if (color == *palette)
        return lIndex;
    }
  }
  else
  {
    imbyte red1, green1, blue1;
    imColorDecode(&red1, &green1, &blue1, color);

    for (int lIndex = 0; lIndex < palette_count; lIndex++, palette++)
    {
      imbyte red2, green2, blue2;
      imColorDecode(&red2, &green2, &blue2, *palette);

      if (iAbs(red1 - red2) < tol &&
          iAbs(green1 - green2) < tol &&
          iAbs(blue1 - blue2) < tol)
      {
        return lIndex;
      }
    }
  }

  return -1;
}

long* imPaletteGray(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;

  for (int lIndex = 0; lIndex < 256; lIndex++)
  {
    /* From (0, 0, 0) to (255, 255, 255)*/
    /* From   Black   to      White     */
    *(ct++) = imColorEncode((imbyte)lIndex, (imbyte)lIndex, (imbyte)lIndex);
  }

  return palette;
}

long* imPaletteRed(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;

  for (int lIndex = 0; lIndex < 256; lIndex++)
  {
    /* From (0, 0, 0) to (255, 0, 0) */
    /* From   Black   to      Red   */
    *(ct++) = imColorEncode((imbyte)lIndex, 0, 0);
  }

  return palette;
}

long* imPaletteGreen(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  
  for (int lIndex = 0; lIndex < 256; lIndex++)
  {
    /* From (0, 0, 0) to (0, 255, 0)*/
    /* From   Black   to    Green   */
    *(ct++) = imColorEncode(0, (imbyte)lIndex, 0);
  }

  return palette;
}

long* imPaletteBlue(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  
  for (int lIndex = 0; lIndex < 256; lIndex++)
  {
    /* From (0, 0, 0) to (0, 0, 255)*/
    /* From   Black   to    Blue    */
    *(ct++) = imColorEncode(0, 0, (imbyte)lIndex);
  }

  return palette;
}

long* imPaletteYellow(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  
  for (int lIndex = 0; lIndex < 256; lIndex++)
  {
    /* From (0, 0, 0) to (255, 255, 0)*/
    /* From   Black   to      Yellow  */
    *(ct++) = imColorEncode((imbyte)lIndex, (imbyte)lIndex, 0);
  }

  return palette;
}

long* imPaletteMagenta(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  
  for (int lIndex = 0; lIndex < 256; lIndex++)
  {
    /* From (0, 0, 0) to (255, 0, 255)*/
    /* From   Black   to    Magenta   */
    *(ct++) = imColorEncode((imbyte)lIndex, 0, (imbyte)lIndex);
  }

  return palette;
}

long* imPaletteCyan(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  
  for (int lIndex = 0; lIndex < 256; lIndex++)
  {
    /* From (0, 0, 0) to (0, 255, 255)*/
    /* From   Black   to     Cyan    */
    *(ct++) = imColorEncode(0, (imbyte)lIndex, (imbyte)lIndex);
  }

  return palette;
}

long* imPaletteHues(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  int i;
  float tone, step1 = 255.0f/41.0f, step2 = 255.0f/42.0f;

  /* 1+42+1+41+1+42+1+41+1+42+1+41+1 = 256 */

  /* red */
  *(ct++) = imColorEncode((imbyte)255, 0, 0);

  for (tone = step2, i = 0; i < 42; i++, tone += step2)
  {
    /* From (255, 0, 0) to (255, 255, 0) */
    /* From      Red    to      Yellow   */
    *(ct++) = imColorEncode((imbyte)255, (imbyte)tone, 0);
  }

  /* yellow */
  *(ct++) = imColorEncode((imbyte)255, (imbyte)255, 0);

  for (tone = step1, i = 0; i < 41; i++, tone += step1)
  {
    /* From (255, 255, 0) to (0, 255, 0)  */
    /* From     Yellow    to    Green    */
    *(ct++) = imColorEncode((imbyte)(255.0f-tone), (imbyte)255, 0);
  }

  /* green */
  *(ct++) = imColorEncode(0, (imbyte)255, 0);;

  for (tone = step2, i = 0; i < 42; i++, tone += step2)
  {
    /* From (0, 255, 0) to (0, 255, 255) */
    /* From    Green    to     Cyan      */
    *(ct++) = imColorEncode(0, (imbyte)255, (imbyte)tone);
  }

  /* cian */
  *(ct++) = imColorEncode(0, (imbyte)255, (imbyte)255);

  for (tone = step1, i = 0; i < 41; i++, tone += step1)
  {
    /* From (0, 255, 255) to (0, 0, 255) */
    /* From     Cyan      to     Blue    */
    *(ct++) = imColorEncode(0, (imbyte)(255.0f-tone), (imbyte)255);
  }

  /* blue */
  *(ct++) = imColorEncode(0, 0, (imbyte)255);

  for (tone = step2, i = 0; i < 42; i++, tone += step2)
  {
    /* From (0, 0, 255) to (255, 0, 255) */
    /* From    Blue     to    Magenta    */
    *(ct++) = imColorEncode((imbyte)tone, 0, (imbyte)255);
  }

  /* magenta */
  *(ct++) = imColorEncode((imbyte)255, 0, (imbyte)255);

  for (tone = step1, i = 0; i < 41; i++, tone += step1)
  {
    /* From (255, 0, 255) to (255, 0, 0) */
    /* From    Magenta    to      Red    */
    *(ct++) = imColorEncode((imbyte)255, 0, (imbyte)(255.0f-tone));
  }

  /* red */
  *(ct++) = imColorEncode((imbyte)255, 0, 0);

  return palette;
}

long* imPaletteRainbow(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  int hue;
  unsigned char r, g, b;
  double h, s, i, factor, H;

  s = 1.0;
  factor = 360.0 / 255.0;

  for (hue = 0; hue < 256; hue++)
  {
    h = hue * factor;
    h = 300.0 - h;
    if (h < 0) h += 360;
    H = h / 57.2957795131;

    i = imColorHSI_ImaxS(H, cos(H), sin(H));

    imColorHSI2RGBbyte(h, s, i, &r, &g, &b);

    *(ct++) = imColorEncode(r, g, b);;
  }

  return palette;
}

long* imPaletteBlueIce(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  
  for (int lIndex = 0; lIndex < 256; lIndex++)
  {
    /* From (0, 0, 255) to (255, 255, 255)*/
    /* From    Blue    to       White     */
    *(ct++) = imColorEncode((imbyte)lIndex, (imbyte)lIndex, 255);
  }

  return palette;
}

long* imPaletteHotIron(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  int lIndex, lSubIndex;

  for (lIndex = 0, lSubIndex = 0; lSubIndex < 128; lSubIndex++, lIndex += 2)
  {
    /* From (0, 0, 0) to (254, 0, 0) */
    /* From   Black   to     ~Red    */
    *(ct++) = imColorEncode((imbyte)lIndex, 0, 0);
  }

  for (lIndex = 0, lSubIndex = 0; lSubIndex < 64; lSubIndex++, lIndex += 2)
  {
    /* From (255, 0, 0) to (255, 126, 0) */
    /* From      Red    to     ~Orange  */
    *(ct++) = imColorEncode(255, (imbyte)lIndex, 0);
  }

  for (lIndex = 0, lSubIndex = 0; lSubIndex < 63; lSubIndex++, lIndex += 2)
  {
    /* From (255, 128, 0) to (255, 252, 252)*/
    /* From     Orange    to     ~White   */
    imbyte red = 255;
    imbyte green = (imbyte)(128 + lIndex);
    imbyte blue = (imbyte)(lIndex * 2 + 4);

    *(ct++) = imColorEncode(red, green, blue);
  }

  *(ct++) = imColorEncode(255, 255, 255);

  return palette;
}

long* imPaletteBlackBody(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  int lIndex, lSubIndex;

  for (lIndex = 0, lSubIndex = 0; lSubIndex < 85; lSubIndex++, lIndex += 3)
  {
    /* From (0, 0, 0) to (252, 0, 0) */
    /* From   Black   to     ~Red   */
    *(ct++) = imColorEncode((imbyte)lIndex, 0, 0);
  }

  for (lIndex = 0, lSubIndex = 0; lSubIndex < 85; lSubIndex++, lIndex += 3)
  {
    /* From (255, 0, 0) to (255, 252, 0)*/
    /* From      Red    to     ~Yellow */
    *(ct++) = imColorEncode(255, (imbyte)lIndex, 0);
  }

  for (lIndex = 0, lSubIndex = 0; lSubIndex < 86; lSubIndex++, lIndex += 3)
  {
    /* From (255, 255, 0) to (255, 255, 255)*/
    /* From     Yellow    to      White  */
    *(ct++) = imColorEncode(255, 255, (imbyte)lIndex);
  }

  return palette;
}

long* imPaletteHighContrast(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  int lIndex;
#define NUM_HC 128

  static struct{ unsigned char r, g, b; } HighContrastColors[NUM_HC] = {
    { 0,0,0 },     

    { 255,0,0 },      { 128,0,0 },      { 64,0,0 },       { 192,0,0 },
    { 0,255,0 },      { 0,128,0 },      { 0,64,0 },       { 0,192,0 },    
    { 0,0,255 },      { 0,0,128 },      { 0,0,64 },       { 0,0,192 },    
    { 255,255,0 },    { 128,128,0 },    { 64,64,0 },      { 192,192,0 },    
    { 255,0,255 },    { 128,0,128 },    { 64,0,64 },      { 192,0,192 },    
    { 0,255,255 },    { 0,128,128 },    { 0,64,64 },      { 0,192,192 },    
    { 255,255,255 },  { 128,128,128 },  { 64,64,64 },     { 192,192,192 },    

    { 255,128,128 },  { 64,255,255 },   { 192,255,255 },   
    { 128,255,128 },  { 255,64,255 },   { 255,192,255 },     
    { 128,128,255 },  { 255,255,64 },   { 255,255,192 },     
    { 255,255,128 },  { 64,64,255 },    { 192,192,255 },     
    { 255,128,255 },  { 64,255,64 },    { 192,255,192 },     
    { 128,255,255 },  { 255,64,64 },    { 255,192,192 },   

    { 128,64,64 },    { 128,192,192 },   
    { 64,128,64 },    { 192,128,192 },   
    { 64,64,128 },    { 192,192,128 },   
    { 128,128,64 },   { 128,128,192 },   
    { 128,64,128 },   { 128,192,128 },   
    { 64,128,128 },   { 192,128,128 },   
    
   { 160,0,0 },      { 32,0,0 },       { 224,0,0 },       { 96,0,0 },    
   { 0,160,0 },      { 0,32,0 },       { 0,224,0 },       { 0,96,0 },    
   { 0,0,160 },      { 0,0,32 },       { 0,0,224 },       { 0,0,96 },    
   { 160,160,0 },    { 32,32,0 },      { 224,224,0 },     { 96,96,0 },   
   { 160,0,160 },    { 32,0,32 },      { 224,0,224 },     { 96,0,96 },   
   { 0,160,160 },    { 0,32,32 },      { 0,224,224 },     { 0,96,96 },   
   { 160,160,160 },  { 32,32,32 },     { 224,224,224 },   { 96,96,96 },  

    { 255,160,160 },  { 32,255,255 },   { 224,255,255 },     { 96,255,255 },  
    { 160,255,160 },  { 255,32,255 },   { 255,224,255 },     { 255,96,255 },   
    { 160,160,255 },  { 255,255,32 },   { 255,255,224 },     { 255,255,96 },   
    { 255,255,160 },  { 32,32,255 },    { 224,224,255 },     { 96,96,255 },   
    { 255,160,255 },  { 32,255,32 },    { 224,255,224 },     { 96,255,96 },   
    { 160,255,255 },  { 255,32,32 },    { 255,224,224 },     { 255,96,96 },  

    { 128,32,32 },    { 64,192,192 },      { 64,224,224 },
    { 32,128,32 },    { 192,64,192 },      { 224,64,224 },
    { 32,32,128 },    { 192,192,64 },      { 224,224,64 },
    { 128,128,32 },   { 64,64,192 },       { 64,64,224 }, 
    { 128,32,128 },   { 64,192,64 },       { 64,224,64 }, 
    { 32,128,128 },   { 192,64,64 },       

 };

  for (lIndex = 0; lIndex < NUM_HC; lIndex++)
  {
    *(ct++) = imColorEncode(HighContrastColors[lIndex].r, 
                            HighContrastColors[lIndex].g, 
                            HighContrastColors[lIndex].b);
  }

  for (; lIndex < 256; lIndex++)
  {
    *(ct++) = imColorEncode(HighContrastColors[lIndex - NUM_HC].r,
                            HighContrastColors[lIndex - NUM_HC].g, 
                            HighContrastColors[lIndex - NUM_HC].b);
  }

  return palette;
}

long* imPaletteLinear(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;
  int lIntensity, lHue;
  unsigned char r, g, b;

  for (lIntensity = 0; lIntensity < 32; lIntensity++)
  {
    double intensity = imColorReconstruct((int)(lIntensity * 8), (int)0, (int)255);  // 8 = 256 / 32         

    for (lHue = 0; lHue < 8; lHue++)
    {
      double hue = imColorReconstruct((int)(lHue * 45), (int)0, (int)360) * 360;  // 45 = 360 / 8

      imColorHSI2RGBbyte(hue, 1.0, intensity, &r, &g, &b);

      *(ct++) = imColorEncode(r, g, b);
    }
  }

  return palette;
}

/* 256 divided in 6 steps results in these steps.*/
static int iSixStepsTable[6] = {0, 51, 102, 153, 204, 255};

long* imPaletteUniform(void)
{
  long* palette = imPaletteNew(256);
  long* ct = palette;

  for (int lRedIndex = 0; lRedIndex < 6; lRedIndex++)
    for (int lGreenIndex = 0; lGreenIndex < 6; lGreenIndex++)
      for (int lBlueIndex = 0; lBlueIndex < 6; lBlueIndex++)
      {
        imbyte red = (imbyte)iSixStepsTable[lRedIndex];
        imbyte green = (imbyte)iSixStepsTable[lGreenIndex];
        imbyte blue = (imbyte)iSixStepsTable[lBlueIndex];

        *(ct++) = imColorEncode(red, green, blue);
      }

  /* We initialized only 216 colors (6x6x6), rest 40 colors.*/
  /* Fill them with a gray scale palette.*/
  for (int lIndex = 6; lIndex < 246; lIndex += 6)
  {
    *(ct++) = imColorEncode((imbyte)lIndex, (imbyte)lIndex, (imbyte)lIndex);
  }

  return palette;
}

/* X divided by 51. Convert to 216 color space. */
static int iDividedBy51Table[256] =
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5
};

/* X multiplied by 36. Shift to red position.*/
static int iTimes36Table[6] = {0, 36, 72, 108, 144, 180};

/* X multiplied by 36. Shift to green position.*/
static int iTimes6Table[6] = {0, 6, 12, 18, 24, 30};

int imPaletteUniformIndex(long color)
{
  imbyte red, green, blue;
  imColorDecode(&red, &green, &blue, color);
  return iTimes36Table[iDividedBy51Table[red]] + iTimes6Table[iDividedBy51Table[green]] + iDividedBy51Table[blue];
}

/* Remainder of X divided by 51. Used to position in the halftone*/
static int iModulo51Table[256] =
{
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
  16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
  32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
  48, 49, 50,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
  29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
  45, 46, 47, 48, 49, 50,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
  26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
  42, 43, 44, 45, 46, 47, 48, 49, 50,  0,  1,  2,  3,  4,  5,  6,
  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
  23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
  39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,  0,  1,  2,  3,
  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
  36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,  0
};

/* Dither matrices for 8 bit to 2.6 bit halftones.*/
static int iHalftone8x8Table[64] =
{
  0, 38,  9, 47,  2, 40, 11, 50,
  25, 12, 35, 22, 27, 15, 37, 24,
  6, 44,  3, 41,  8, 47,  5, 43,
  31, 19, 28, 15, 34, 21, 31, 18,
  1, 39, 11, 49,  0, 39, 10, 48,
  27, 14, 36, 23, 26, 13, 35, 23,
  7, 46,  4, 43,  7, 45,  3, 42,
  33, 20, 30, 17, 32, 19, 29, 16
};

int imPaletteUniformIndexHalftoned(long color, int x, int y)
{
  int lHalf = iHalftone8x8Table[(x % 8) * 8 + y % 8];

  imbyte red, green, blue;
  imColorDecode(&red, &green, &blue, color);

  /* Now, look up each value in the halftone matrix using an 8x8 ordered dither.*/
  int lRed = iDividedBy51Table[red] + (iModulo51Table[red] > lHalf? 1: 0);
  int lGreen = iDividedBy51Table[green] + (iModulo51Table[green] > lHalf? 1: 0);
  int lBlue = iDividedBy51Table[blue] + (iModulo51Table[blue] > lHalf? 1: 0);

  return iTimes36Table[lRed] + iTimes6Table[lGreen] + lBlue;
}
