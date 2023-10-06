/*
  Simple Draw Application

  Shows the same picture on several canvas. Used to quick test the CD library and
  to demonstrate the use of CD library functions.

  This module uses only the CD library,
  there is another module to initialize the Window and its menus.
  */


#include "cd.h"
#include "cdcgm.h"
#include "cddgn.h"
#include "cddxf.h"
#include "cdclipbd.h"
#include "cdemf.h"
#include "cdimage.h"
#include "cdirgb.h"
#include "cdmf.h"
#include "cdprint.h"
#include "cdps.h"
#include "cdpdf.h"
#include "cdpptx.h"
#include "cdsvg.h"
#include "cdwmf.h"
#include "cdiup.h"
#include "cddbuf.h"
#include "cddebug.h"
#include "wd.h"
#include "cdgdiplus.h"
#include "cddirect2d.h"
#include "cdgl.h"
#include <iupdraw_cd.h>

#include "simple.h"

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

/* Global variables */

cdCanvas *winCanvas = NULL;        /* The window drawing canvas */
char* winData = NULL;
cdCanvas *dbCanvas = NULL;        /* The double buffer canvas */
cdCanvas *curCanvas = NULL;        /* The current canvas */

int clipping = CD_CLIPOFF;                     /* Clipping flag, same as the CD */
int write_mode = CD_REPLACE;                   /* Write Mode flag, same as the CD */
int contextplus = 0;
int simple_draw = 0;
int use_transform = 0;
int simulate = 0;
int use_opengl = 0;

enum { DRAW_ALL, DRAW_TEXTFONTS, DRAW_TEXTALIGN, DRAW_TEST };

#define STYLE_SIZE 10 /* A small pattern and stipple size */
long pattern[STYLE_SIZE*STYLE_SIZE];          /* Pattern buffer */
unsigned char stipple[STYLE_SIZE*STYLE_SIZE]; /* Stipple buffer */

#define IMAGE_SIZE 100
unsigned char red[IMAGE_SIZE*IMAGE_SIZE];       /* Red image buffer */
unsigned char green[IMAGE_SIZE*IMAGE_SIZE];     /* Green image buffer */
unsigned char blue[IMAGE_SIZE*IMAGE_SIZE];      /* Blue image buffer */
unsigned char alpha[IMAGE_SIZE*IMAGE_SIZE];     /* Alpha image buffer */


/* Prototype of the function that makes the drawing independent of canvas. */
void SimpleDraw(cdCanvas* canvas);

void SimpleInitAlpha(int width, int height, unsigned char* _alpha)
{
  int c, l;
  /* initialize the alpha image buffer with a degrade from transparent to opaque */
  for (l = 0; l < height; l++)
  for (c = 0; c < width; c++)
    _alpha[l*width + c] = (unsigned char)((c * 255) / (width - 1));
}

void SimpleCreateCanvasWindow(void)
{
  /* creates the canvas based in an existing window */
  if (contextplus) cdUseContextPlus(1);
  winCanvas = cdCreateCanvas(CD_IUP, winData);
  if (contextplus) cdUseContextPlus(0);
  curCanvas = winCanvas;
}

void SimpleCreateCanvas(char* data)
{
  int c, l;

  memset(pattern, 0xFF, STYLE_SIZE*STYLE_SIZE * 4);

  pattern[11] = CD_RED;   /*------------*/
  pattern[21] = CD_RED;   /*  0123456789*/
  pattern[31] = CD_RED;   /*            */
  pattern[41] = CD_RED;   /*9 WWWWWWWWWW*/
  pattern[51] = CD_RED;   /*8 WWWWGGGGGW*/
  pattern[12] = CD_RED;   /*7 WWWGGGGGBW*/
  pattern[22] = CD_RED;   /*6 WWGGGGGBBW*/
  pattern[32] = CD_RED;   /*5 WrrrrrBBBW*/
  pattern[42] = CD_RED;   /*4 WrrrrrBBBW*/
  pattern[52] = CD_RED;   /*3 WrrrrrBBWW*/
  pattern[13] = CD_RED;   /*2 WrrrrrBWWW*/
  pattern[23] = CD_RED;   /*1 WrrrrrWWWW*/
  pattern[33] = CD_RED;   /*0 WWWWWWWWWW*/
  pattern[43] = CD_RED;   /*------------*/
  pattern[53] = CD_RED;
  pattern[14] = CD_RED;   pattern[15] = CD_RED;
  pattern[24] = CD_RED;   pattern[25] = CD_RED;
  pattern[34] = CD_RED;   pattern[35] = CD_RED;
  pattern[44] = CD_RED;   pattern[45] = CD_RED;
  pattern[54] = CD_RED;   pattern[55] = CD_RED;

  pattern[26] = CD_BLUE;  pattern[37] = CD_BLUE;
  pattern[36] = CD_BLUE;  pattern[47] = CD_BLUE;
  pattern[46] = CD_BLUE;  pattern[57] = CD_BLUE;
  pattern[56] = CD_BLUE;  pattern[67] = CD_BLUE;

  pattern[48] = CD_BLUE;  pattern[62] = CD_GREEN;
  pattern[58] = CD_BLUE;  pattern[63] = CD_GREEN;
  pattern[68] = CD_BLUE;  pattern[64] = CD_GREEN;
  pattern[78] = CD_BLUE;  pattern[65] = CD_GREEN;
  pattern[66] = CD_GREEN;

  pattern[73] = CD_GREEN; pattern[84] = CD_GREEN;
  pattern[74] = CD_GREEN; pattern[85] = CD_GREEN;
  pattern[75] = CD_GREEN; pattern[86] = CD_GREEN;
  pattern[76] = CD_GREEN; pattern[87] = CD_GREEN;
  pattern[77] = CD_GREEN; pattern[88] = CD_GREEN;

  /* initialize the stipple buffer with cross pattern */
  for (l = 0; l < STYLE_SIZE; l++)
  for (c = 0; c < STYLE_SIZE; c++)
    stipple[l*STYLE_SIZE + c] = (c % 4) == 0 ? 1 : 0;

  SimpleInitAlpha(IMAGE_SIZE, IMAGE_SIZE, alpha);

  winData = data;
  SimpleCreateCanvasWindow();
  SimpleDrawWindow();
}

int SimpleTransform(void)
{
  use_transform = !use_transform;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleContextPlus(void)
{
#ifdef USE_CONTEXTPLUS
  contextplus = !contextplus;
  SimpleKillCanvas();
  SimpleCreateCanvasWindow();
  SimpleDraw(curCanvas);
#endif
  return 0;
}

void PlayCanvasDriver(cdContext* ctx, char* StrData)
{
  int w, h;
  cdCanvasActivate(curCanvas);
  cdCanvasBackground(curCanvas, CD_WHITE);
  cdCanvasClear(curCanvas);
  cdCanvasGetSize(curCanvas, &w, &h, 0, 0);
  //  cdCanvasPlay(curCanvas, ctx, 100, w-100, 100, h-100, StrData);
  cdCanvasPlay(curCanvas, ctx, 0, w - 1, 0, h - 1, StrData);
  //  cdCanvasPlay(curCanvas, ctx, 0, 0, 0, 0, StrData);
}

int SimplePlayClipboard(void)
{
  PlayCanvasDriver(CD_CLIPBOARD, NULL);
  return 0;
}

int SimplePlayCGMBin(void)
{
  PlayCanvasDriver(CD_CGM, "simple_b.cgm");
  return 0;
}

int SimplePlayCGMText(void)
{
  PlayCanvasDriver(CD_CGM, "simple_t.cgm");
  return 0;
}

int SimplePlayMetafile(void)
{
  PlayCanvasDriver(CD_METAFILE, "simple.mf");
  return 0;
}

int SimplePlayWMF(void)
{
  PlayCanvasDriver(CD_WMF, "simple.wmf");
  return 0;
}

int SimplePlayEMF(void)
{
  PlayCanvasDriver(CD_EMF, "simple.emf");
  return 0;
}

int SimpleRepaint(void)
{
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleDrawWindow(void)
{
  use_opengl = 0;
  curCanvas = winCanvas;
  SimpleDraw(curCanvas);
  return 0;
}

void DrawCanvasDriver(cdContext* ctx, char* StrData)
{
  cdCanvas* tmpCanvas = cdCreateCanvas(ctx, StrData);
  if (tmpCanvas == NULL)
  {
    printf("CreateCanvas(%s) - Failed!\n", StrData);
    return;
  }
  printf("CreateCanvas(%s)\n", StrData);
  SimpleDraw(tmpCanvas);
  cdKillCanvas(tmpCanvas);
  printf("KillCanvas()\n");
}

void DrawCanvasDriverSize(cdContext* ctx, char* name, int pixels, const char* opt)
{
  char StrData[100];
  int w, h;
  double w_mm, h_mm;
  cdCanvasGetSize(curCanvas, &w, &h, &w_mm, &h_mm);
  if (pixels == 1)  /* WMF and EMF */
    sprintf(StrData, "%s %dx%d %s", name, w, h, opt);
  else if (pixels == 2)  /* PDF and PS */
    sprintf(StrData, "%s -w%g -h%g -s%g %s", name, w_mm, h_mm, ((double)w / w_mm)*25.4, opt);
  else  /* others */
    sprintf(StrData, "%s %gx%g %g %s", name, w_mm, h_mm, (double)w / w_mm, opt);
  DrawCanvasDriver(ctx, StrData);
}

void DrawCanvasDriverSizeParam(cdContext* ctx, char* param)
{
  char StrData[100];
  int w, h;
  cdCanvasGetSize(curCanvas, &w, &h, 0, 0);
  sprintf(StrData, "%dx%d %s", w, h, param);
  DrawCanvasDriver(ctx, StrData);
}

int SimpleDrawDebug(void)
{
  DrawCanvasDriverSize(CD_DEBUG, "simple_debug.txt", 0, "");
  return 0;
}

int SimpleDrawCGMText(void)
{
  DrawCanvasDriverSize(CD_CGM, "simple_t.cgm", 0, "-t");
  return 0;
}

int SimpleDrawCGMBin(void)
{
  DrawCanvasDriverSize(CD_CGM, "simple_b.cgm", 0, "");
  return 0;
}

int SimpleDrawDXF(void)
{
  DrawCanvasDriverSize(CD_DXF, "simple.dxf", 0, "-ac2000");  //""
  return 0;
}

int SimpleDrawDGN(void)
{
  DrawCanvasDriverSize(CD_DGN, "simple.dgn", 0, "");
  return 0;
}

int SimpleDrawEMF(void)
{
  if (contextplus) cdUseContextPlus(1);
  DrawCanvasDriverSize(CD_EMF, "simple.emf", 1, "");
  if (contextplus) cdUseContextPlus(0);
  return 0;
}

int SimpleDrawMetafile(void)
{
  DrawCanvasDriverSize(CD_METAFILE, "simple.mf", 0, "");
  return 0;
}

int SimpleDrawPS(void)
{
  DrawCanvasDriverSize(CD_PS, "simple.ps", 2, "-l0 -r0 -t0 -b0");
  return 0;
}

int SimpleDrawSVG(void)
{
  DrawCanvasDriverSize(CD_SVG, "simple.svg", 0, "");
  return 0;
}

int SimpleDrawPDF(void)
{
  DrawCanvasDriverSize(CD_PDF, "simple.pdf", 2, "");
  return 0;
}

int SimpleDrawEPS(void)
{
  DrawCanvasDriverSize(CD_PS, "simple.eps", 2, "-e");
  return 0;
}

int SimpleDrawWMF(void)
{
  DrawCanvasDriverSize(CD_WMF, "simple.wmf", 1, "");
  return 0;
}

int SimpleDrawPPTX(void)
{
  DrawCanvasDriverSize(CD_PPTX, "simple.pptx", 0, "");
  return 0;
}

int SimpleDrawPrint(void)
{
  if (contextplus) cdUseContextPlus(1);
  DrawCanvasDriver(CD_PRINTER, "simple print");
  if (contextplus) cdUseContextPlus(0);
  return 0;
}

int SimpleDrawPrintDialog(void)
{
  if (contextplus) cdUseContextPlus(1);
  DrawCanvasDriver(CD_PRINTER, "simple -d");   /* show dialog */
  if (contextplus) cdUseContextPlus(0);
  return 0;
}

int SimpleDrawClipboardBitmap(void)
{
  if (contextplus) cdUseContextPlus(1);
  DrawCanvasDriverSizeParam(CD_CLIPBOARD, "-b");
  if (contextplus) cdUseContextPlus(0);
  return 0;
}

int SimpleDrawClipboardMetafile(void)
{
  if (contextplus) cdUseContextPlus(1);
  DrawCanvasDriverSizeParam(CD_CLIPBOARD, "-m");
  if (contextplus) cdUseContextPlus(0);
  return 0;
}

int SimpleDrawClipboardEMF(void)
{
  if (contextplus) cdUseContextPlus(1);
  DrawCanvasDriverSizeParam(CD_CLIPBOARD, "");
  if (contextplus) cdUseContextPlus(0);
  return 0;
}

int SimpleReplace(void)
{
  write_mode = CD_REPLACE;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleXor(void)
{
  write_mode = CD_XOR;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleNotXor(void)
{
  write_mode = CD_NOT_XOR;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleClippingOff(void)
{
  clipping = CD_CLIPOFF;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleClippingArea(void)
{
  clipping = CD_CLIPAREA;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleClippingPolygon(void)
{
  clipping = CD_CLIPPOLYGON;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleClippingRegion(void)
{
  clipping = CD_CLIPREGION;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleAll(void)
{
  simple_draw = DRAW_ALL;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleTextAlign(void)
{
  simple_draw = DRAW_TEXTALIGN;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleTextFonts(void)
{
  simple_draw = DRAW_TEXTFONTS;
  SimpleDraw(curCanvas);
  return 0;
}

int SimpleTest(void)
{
  simple_draw = DRAW_TEST;
  SimpleDraw(curCanvas);
  return 0;
}

void* CreateImageRGBA(int w, int h)
{
  void* myImage;
  unsigned char * _alpha = malloc(w * h);
  SimpleInitAlpha(w, h, _alpha);
  cdCanvasSetAttribute(curCanvas, "IMAGEALPHA", (char*)_alpha);
  cdCanvasSetAttribute(curCanvas, "IMAGEFORMAT", "32");    // afetara´ o proximo cdCreateImage
  myImage = cdCanvasCreateImage(curCanvas, w, h);
  cdCanvasSetAttribute(curCanvas, "IMAGEFORMAT", NULL);    // remove o atributo para nao afetar outros cdCreateImage
  return myImage;
}

int SimpleDrawImage(void)
{
  use_opengl = 0;
  if (dbCanvas) cdKillCanvas(dbCanvas);

  if (contextplus) cdUseContextPlus(1);
  dbCanvas = cdCreateCanvas(CD_DBUFFER, winCanvas);
  if (contextplus) cdUseContextPlus(0);

  curCanvas = dbCanvas;
  SimpleDraw(curCanvas);

  return 0;
}

int SimpleDrawImageRGB(void)
{
  use_opengl = 0;
  if (dbCanvas) cdKillCanvas(dbCanvas);

  if (contextplus) cdUseContextPlus(1);
  dbCanvas = cdCreateCanvas(CD_DBUFFERRGB, winCanvas);
  if (contextplus) cdUseContextPlus(0);

  curCanvas = dbCanvas;
  SimpleDraw(curCanvas);

  return 0;
}

#ifdef USE_OPENGL
int SimpleDrawGL(void)
{
  char StrData[100];
  int w, h;
  double w_mm, h_mm;

  if (use_opengl)
    return 0;

  cdCanvasGetSize(curCanvas, &w, &h, &w_mm, &h_mm);

  sprintf(StrData, "%dx%d %g", w, h, ((double)w / w_mm));

  if (dbCanvas) cdKillCanvas(dbCanvas);

  dbCanvas = cdCreateCanvas(CD_GL, StrData);

  curCanvas = dbCanvas;
  use_opengl = 1;
  SimpleDraw(curCanvas);

  return 0;
}
#endif

int SimpleDrawIupDraw(void)
{
  if (dbCanvas) cdKillCanvas(dbCanvas);

  dbCanvas = cdCreateCanvas(CD_IUPDRAW, winData);

  curCanvas = dbCanvas;

  SimpleDraw(curCanvas);

  return 0;
}

#ifdef WIN32
#include <iup.h>
#endif

int SimpleDrawDirect2D(void)
{
#ifdef WIN32
  if (dbCanvas) cdKillCanvas(dbCanvas);

  dbCanvas = cdCreateCanvas(CD_DIRECT2D_NATIVEWINDOW, IupGetAttribute((Ihandle*)winData, "HWND"));

  curCanvas = dbCanvas;

  SimpleDraw(curCanvas);
#endif

  return 0;
}

int SimpleDrawSimulate(void)
{
  simulate = !simulate;

  if (simulate)
    cdCanvasSimulate(curCanvas, CD_SIM_ALL);
  else
    cdCanvasSimulate(curCanvas, CD_SIM_NONE);

  SimpleDraw(curCanvas);

  return 0;
}

void SimpleKillCanvas(void)
{
  if (dbCanvas)
  {
    cdKillCanvas(dbCanvas);
    dbCanvas = NULL;
  }
  if (winCanvas)
  {
    cdKillCanvas(winCanvas);
    winCanvas = NULL;
  }
}

void SimpleDrawTextFonts(cdCanvas* canvas);
void SimpleDrawTextAlign(cdCanvas* canvas);
void SimpleDrawAll(cdCanvas* canvas);
void SimpleDrawTest(cdCanvas* canvas);

void SimpleDraw(cdCanvas* canvas)
{
  if (!canvas)
    return;

#ifdef USE_OPENGL
  if (use_opengl)
    SimpleUpdateSize(canvas);
#endif

  /* refresh CD canvas size, when window size has changed */
  cdCanvasActivate(canvas);

#if 1
  if (cdCanvasGetContext(canvas) == CD_PPTX)
  {
#if 0
    int width, height, dx, dy;

    cdCanvasSetAttribute(canvas, "MASTERSLIDE", "1");

    cdCanvasGetSize(canvas, &width, &height, NULL, NULL);
    dx = width / 10;
    dy = height / 10;

    cdCanvasForeground(canvas, cdEncodeColor(220, 220, 255));
    cdCanvasBox(canvas, dx, width - dx, dy, height - dy);
    cdCanvasLineWidth(canvas, 10);
    cdCanvasForeground(canvas, cdEncodeColor(120, 120, 255));
    cdCanvasRect(canvas, dx, width - dx, dy, height - dy);

    cdCanvasSetAttribute(canvas, "MASTERSLIDE", NULL);
#else
    cdCanvasSetAttribute(canvas, "MASTERSLIDEFILE", "master.pptx");
#endif
}
#endif


  if (simple_draw == DRAW_TEXTFONTS)
    SimpleDrawTextFonts(canvas);
  else if (simple_draw == DRAW_TEXTALIGN)
    SimpleDrawTextAlign(canvas);
  else if (simple_draw == DRAW_TEST)
    SimpleDrawTest(canvas);
  else
    SimpleDrawAll(canvas);

  /* Adds a new page, or
     flushes the file, or
     flushes the screen, or
     swap the double buffer. */
  cdCanvasFlush(canvas);

#ifdef USE_OPENGL
  if (use_opengl)
    SimpleFlush();
#endif
}

void SimpleDrawAll(cdCanvas* canvas)
{
  int w, h;
  cdCanvasGetSize(canvas, &w, &h, NULL, NULL);

  /* Clear the background to be white */
  cdCanvasBackground(canvas, CD_WHITE);
  //  cdBackground(CD_GREEN);
  cdCanvasClear(canvas);

//  cdCanvasSetAttribute(canvas, "ANTIALIAS", "0");

  /* Draw a rectangle and a polygon at the bottom-left area,
     using a thick line with transparency.
     Observe that transparency is only supported in a few drivers,
     and line join is not supported in the IMAGERGB driver. */
  cdCanvasLineWidth(canvas, 3);
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);
  cdCanvasForeground(canvas, cdEncodeAlpha(CD_DARK_MAGENTA, 128));
  cdCanvasRect(canvas, 100, 200, 100, 200);

  cdCanvasBegin(canvas, CD_OPEN_LINES);
  cdCanvasVertex(canvas, 300, 250);
  cdCanvasVertex(canvas, 320, 270);
  cdCanvasVertex(canvas, 350, 260);
  cdCanvasVertex(canvas, 340, 200);
  cdCanvasVertex(canvas, 310, 210);
  cdCanvasEnd(canvas);

  /* Draw the red diagonal line with a custom line style.
     Observe that line styles are not supported in the IMAGERGB driver. */
  cdCanvasForeground(canvas, CD_RED);
  cdCanvasLineWidth(canvas, 3);
  {
    int dashes[] = { 20, 15, 5, 5 };
    cdCanvasLineStyleDashes(canvas, dashes, 4);
  }
  cdCanvasLineStyle(canvas, CD_CUSTOM);
  cdCanvasLine(canvas, 0, 0, w - 1, h - 1);

  /* Draw the blue diagonal line with a pre-defined line style.
     Observe that the pre-defined line style is dependent on the driver. */
  cdCanvasForeground(canvas, CD_BLUE);
  cdCanvasLineWidth(canvas, 10);
  cdCanvasLineStyle(canvas, CD_DOTTED);
  cdCanvasLine(canvas, 0, h - 1, w - 1, 0);

  switch (clipping)
  {
  case CD_CLIPOFF:
    cdCanvasClip(canvas, CD_CLIPOFF);
    break;
  case CD_CLIPAREA:
    /* Defines the clipping area equals the canvas area minus a 100 pixels margin. */
    cdCanvasClipArea(canvas, 100, w - 100, 100, h - 100);
    cdCanvasClip(canvas, CD_CLIPAREA);
    break;
  case CD_CLIPPOLYGON:
    cdCanvasBegin(canvas, CD_CLIP);
    cdCanvasVertex(canvas, 100, 100);
    cdCanvasVertex(canvas, w - 100, 100);
    cdCanvasVertex(canvas, w / 2, h - 100);
    cdCanvasEnd(canvas);
    cdCanvasClip(canvas, CD_CLIPPOLYGON);
    break;
  case CD_CLIPREGION:

    cdCanvasTextAlignment(canvas, CD_CENTER);
    cdCanvasFont(canvas, "Times", CD_BOLD, 50);

    cdCanvasBegin(canvas, CD_REGION);
    cdCanvasRegionCombineMode(canvas, CD_UNION);
    cdCanvasBox(canvas, 100, 200, 100, 200);
    cdCanvasSector(canvas, w / 2 - 50, h / 2 + 50, 150, 150, 0, 360);
    cdCanvasSector(canvas, w / 2 - 50, h / 2 - 50, 150, 150, 0, 360);
    cdCanvasSector(canvas, w / 2 + 50, h / 2 + 50, 150, 150, 0, 360);
    cdCanvasSector(canvas, w / 2 + 50, h / 2 - 50, 150, 150, 0, 360);
    cdCanvasRegionCombineMode(canvas, CD_DIFFERENCE);
    cdCanvasText(canvas, w / 2, h / 2, "TEXT");
    cdCanvasEnd(canvas);
    //    cdCanvasOffsetRegion(canvas, -50, 50);
    cdCanvasClip(canvas, CD_CLIPREGION);

    cdCanvasForeground(canvas, CD_DARK_RED);
    cdCanvasBox(canvas, 0, w, 0, h);
    break;
  }

  switch (write_mode)
  {
  case CD_REPLACE:
    cdCanvasWriteMode(canvas, CD_REPLACE);
    break;
  case CD_XOR:
    cdCanvasWriteMode(canvas, CD_XOR);
    break;
  case CD_NOT_XOR:
    cdCanvasWriteMode(canvas, CD_NOT_XOR);
    break;
  }

  if (use_transform)
  {
    cdCanvasTransform(canvas, NULL);
    cdCanvasTransformTranslate(canvas, w / 2, h / 2);
    cdCanvasTransformRotate(canvas, 30);
    cdCanvasTransformScale(canvas, 0.5, 0.5);
    cdCanvasTransformTranslate(canvas, -w / 2, -h / 2);
  }

  //  cdSetfAttribute("ROTATE", "15 %d %d", w/2, h/2);

  /* Reset line style and width */
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);
  cdCanvasLineWidth(canvas, 1);
  cdCanvasBackOpacity(canvas, CD_OPAQUE);
  cdCanvasBackOpacity(canvas, CD_TRANSPARENT);

  /* Draw an arc at bottom-left, and a sector at bottom-right.
     Notice that counter-clockwise orientation of both. */
  cdCanvasInteriorStyle(canvas, CD_SOLID);
  cdCanvasForeground(canvas, CD_MAGENTA);
  cdCanvasSector(canvas, w - 100, 100, 100, 100, 50, 180);
  cdCanvasForeground(canvas, CD_RED);
  cdCanvasArc(canvas, 100, 100, 100, 100, 50, 180);

  /* Draw a solid filled rectangle at center. */
  cdCanvasForeground(canvas, CD_YELLOW);
  cdCanvasBox(canvas, w / 2 - 100, w / 2 + 100, h / 2 - 100, h / 2 + 100);

  /* Prepare font for text. */
  cdCanvasTextAlignment(canvas, CD_NORTH);
  cdCanvasTextAlignment(canvas, CD_CENTER);
  cdCanvasTextOrientation(canvas, 70);
  cdCanvasFont(canvas, "Times", CD_BOLD, 14);
  //  cdCanvasSetAttribute(canvas, "ADDFONTMAP", "MyWingDings=wingxing");
  //  cdCanvasFont(canvas, "MyWingDings", CD_PLAIN, 24);
  //  cdCanvasFont(canvas, "wingxing", CD_PLAIN, 24); // Current folder
  //  cdCanvasFont(canvas, "WingDings", CD_PLAIN, 24); // Native System - Windows
  //  cdCanvasFont(canvas, "Purisa", CD_PLAIN, 24); // Native System - Linux
  //  cdCanvasFont(canvas, "../../test/simple/wingxing.ttf", CD_PLAIN, 24); // Full path

  //  cdCanvasSetAttribute(canvas, "UTF8MODE", "1");
  /* Draw text at center, with orientation,
     and draw its bounding box.
     Notice that in some drivers the bounding box is not precise. */
  {
    int rect[8];
    //char* mode = cdCanvasGetAttribute(canvas, "UTF8MODE");
    //int utf8mode = mode ? (mode[0] == '1' ? 1 : 0) : 0;
    //if (utf8mode)
    //  cdCanvasGetTextBounds(canvas, w / 2, h / 2, "Simple Draw (p8-Ã§Ã£Ã­)­\nSecond Line", rect);
    //else
    //  cdCanvasGetTextBounds(canvas, w / 2, h / 2, "Simple Draw (p-çãí)\nSecond Line", rect);
    cdCanvasGetTextBounds(canvas, w / 2, h / 2, "cdMin Draw (çãí)", rect);
    cdCanvasForeground(canvas, CD_RED);
    cdCanvasBegin(canvas, CD_CLOSED_LINES);
    cdCanvasVertex(canvas, rect[0], rect[1]);
    cdCanvasVertex(canvas, rect[2], rect[3]);
    cdCanvasVertex(canvas, rect[4], rect[5]);
    cdCanvasVertex(canvas, rect[6], rect[7]);
    cdCanvasEnd(canvas);


    cdCanvasForeground(canvas, CD_BLUE);
    //if (utf8mode)
    //  cdCanvasText(canvas, w / 2, h / 2, "Simple Draw (p8-Ã§Ã£Ã­)\nSecond Line");
    //else
    //  cdCanvasText(canvas, w / 2, h / 2, "Simple Draw (p-çãí)\nSecond Line");
    cdCanvasText(canvas, w / 2, h / 2, "cdMin Draw (çãí)");
    cdCanvasTextOrientation(canvas, 0);
  }

  /* Prepare World Coordinates */
  wdCanvasViewport(canvas, 0, w - 1, 0, h - 1);
  if (w > h)
    wdCanvasWindow(canvas, 0, (double)w / (double)h, 0, 1);
  else
    wdCanvasWindow(canvas, 0, 1, 0, (double)h / (double)w);

  /* Draw a filled blue rectangle in WC */
  wdCanvasBox(canvas, 0.20, 0.30, 0.40, 0.50);
  cdCanvasForeground(canvas, CD_RED);
  /* Draw the diagonal of that rectangle in WC */
  wdCanvasLine(canvas, 0.20, 0.40, 0.30, 0.50);

  //  wdVectorTextDirection(0, 0, 1, 1);
  /* Prepare Vector Text in WC. */
  wdCanvasVectorCharSize(canvas, 0.07);

  //  wdVectorText(0.1, 0.4, "ñç áéíóú àèìòù âêîôû äëïöü");
  //  wdVectorText(0.1, 0.2, "ÑÇ ÁÉÍÓÚ ÀÈÌÒÙ ÂÊÎÔÛ ÄËÏÖÜ");
  //{
  //  int i;
  //  char t[2];
  //  char s[10];
  //  int x = 20;
  //  int y = 0;
  //  t[1] = 0;
  //  for (i = 0; i < 256; i++)
  //  {
  //    int dx = 90;
  //    t[0] = (char)i;
  //    sprintf(s, "%d", i);
  //    cdText(x, y, s);
  //    cdText(x+dx, y, t);
  //    cdVectorText(x+2*dx, y, t);
  //    
  //    x += 3*dx + 2*dx/3;
  //    if ((i+1) % 7 == 0)
  //    {
  //      x = 20;
  //      y += 90;
  //    }

  //  }
  //}

  /* Draw vector text, and draw its bounding box.
     We also use this text to show when we are using a contextplus driver. */
  {
    double rect[8];
    cdCanvasForeground(canvas, CD_RED);
    if (contextplus)
      wdCanvasGetVectorTextBounds(canvas, "WDj-Plus", 0.25, 0.35, rect);
    else
      wdCanvasGetVectorTextBounds(canvas, "WDj", 0.25, 0.35, rect);
    cdCanvasBegin(canvas, CD_CLOSED_LINES);
    wdCanvasVertex(canvas, rect[0], rect[1]);
    wdCanvasVertex(canvas, rect[2], rect[3]);
    wdCanvasVertex(canvas, rect[4], rect[5]);
    wdCanvasVertex(canvas, rect[6], rect[7]);
    cdCanvasEnd(canvas);

    cdCanvasLineWidth(canvas, 2);
    cdCanvasLineStyle(canvas, CD_CONTINUOUS);
    if (contextplus)
      wdCanvasVectorText(canvas, 0.25, 0.35, "WDj-Plus");
    else
      wdCanvasVectorText(canvas, 0.25, 0.35, "WDj");
    cdCanvasLineWidth(canvas, 1);
  }

  /* Draw a filled path at center-right (looks like a weird fish).
     Notice that in PDF the arc is necessarily a circle arc, and not an ellipse. */
  cdCanvasInteriorStyle(canvas, CD_HATCH);
  cdCanvasHatch(canvas, CD_DIAGCROSS);
  cdCanvasForeground(canvas, CD_GREEN);
  cdCanvasBegin(canvas, CD_PATH);
  cdCanvasPathSet(canvas, CD_PATH_MOVETO);
  cdCanvasVertex(canvas, w / 2 + 200, h / 2);
  cdCanvasPathSet(canvas, CD_PATH_LINETO);
  cdCanvasVertex(canvas, w / 2 + 230, h / 2 + 50);
  cdCanvasPathSet(canvas, CD_PATH_LINETO);
  cdCanvasVertex(canvas, w / 2 + 250, h / 2 + 50);
  cdCanvasPathSet(canvas, CD_PATH_CURVETO);
  cdCanvasVertex(canvas, w / 2 + 150 + 150, h / 2 + 200 - 50); /* control point for start */
  cdCanvasVertex(canvas, w / 2 + 150 + 180, h / 2 + 250 - 50); /* control point for end */
  cdCanvasVertex(canvas, w / 2 + 150 + 180, h / 2 + 200 - 50); /* end point */
  cdCanvasPathSet(canvas, CD_PATH_CURVETO);
  cdCanvasVertex(canvas, w / 2 + 150 + 180, h / 2 + 150 - 50);
  cdCanvasVertex(canvas, w / 2 + 150 + 150, h / 2 + 100 - 50);
  cdCanvasVertex(canvas, w / 2 + 150 + 300, h / 2 + 100 - 50);
  cdCanvasPathSet(canvas, CD_PATH_LINETO);
  cdCanvasVertex(canvas, w / 2 + 150 + 300, h / 2 - 50);
  cdCanvasPathSet(canvas, CD_PATH_ARC);
  cdCanvasVertex(canvas, w / 2 + 300, h / 2);  /* center */
  cdCanvasVertex(canvas, 200, 100);  /* width, height */
  cdCanvasVertex(canvas, -30 * 1000, -170 * 1000);  /* start angle, end angle (degrees / 1000) */
  //cdCanvasPathSet(canvas, CD_PATH_CLOSE);
  //cdCanvasPathSet(canvas, CD_PATH_STROKE);
  //cdCanvasPathSet(canvas, CD_PATH_FILL);
  cdCanvasPathSet(canvas, CD_PATH_FILLSTROKE);
  cdCanvasEnd(canvas);

  //cdCanvasForeground(canvas, CD_BLACK);
  //cdCanvasArc(canvas, w / 2 + 300, h / 2, 200, 100, -170, -30);

  //cdCanvasForeground(canvas, CD_DARK_RED);
  //cdCanvasRect(canvas, w / 2 + 300 - 100, w / 2 + 300 + 100, h / 2 - 50, h / 2 + 50);

  /* Draw 3 pixels at center left. */
  cdCanvasPixel(canvas, 10, h / 2 + 0, CD_RED);
  cdCanvasPixel(canvas, 11, h / 2 + 1, CD_GREEN);
  cdCanvasPixel(canvas, 12, h / 2 + 2, CD_BLUE);
  cdCanvasPixel(canvas, 13, h / 2 + 3, CD_RED);
  cdCanvasPixel(canvas, 14, h / 2 + 4, CD_GREEN);
  cdCanvasPixel(canvas, 15, h / 2 + 5, CD_BLUE);

  /* Draw 4 mark types, distributed near each corner.  */
  cdCanvasForeground(canvas, CD_RED);
  cdCanvasMarkSize(canvas, 30);
  cdCanvasMarkType(canvas, CD_PLUS);
  cdCanvasMark(canvas, 200, 200);
  cdCanvasMarkType(canvas, CD_CIRCLE);
  cdCanvasMark(canvas, w - 200, 200);
  cdCanvasMarkType(canvas, CD_HOLLOW_CIRCLE);
  cdCanvasMark(canvas, 200, h - 200);
  cdCanvasMarkType(canvas, CD_DIAMOND);
  cdCanvasMark(canvas, w - 200, h - 200);

  /* Draw all the line style possibilities at bottom.
     Notice that they have some small differences between drivers. */
  cdCanvasLineWidth(canvas, 1);
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);
  cdCanvasLine(canvas, 0, 10, w, 10);
  cdCanvasLineStyle(canvas, CD_DASHED);
  cdCanvasLine(canvas, 0, 20, w, 20);
  cdCanvasLineStyle(canvas, CD_DOTTED);
  cdCanvasLine(canvas, 0, 30, w, 30);
  cdCanvasLineStyle(canvas, CD_DASH_DOT);
  cdCanvasLine(canvas, 0, 40, w, 40);
  cdCanvasLineStyle(canvas, CD_DASH_DOT_DOT);
  cdCanvasLine(canvas, 0, 50, w, 50);

  /* Draw all the hatch style possibilities in the top-left corner.
     Notice that they have some small differences between drivers. */
  cdCanvasHatch(canvas, CD_HORIZONTAL);
  cdCanvasBox(canvas, 0, 50, h - 60, h);
  cdCanvasHatch(canvas, CD_VERTICAL);
  cdCanvasBox(canvas, 50, 100, h - 60, h);
  cdCanvasHatch(canvas, CD_FDIAGONAL);
  cdCanvasBox(canvas, 100, 150, h - 60, h);
  cdCanvasHatch(canvas, CD_BDIAGONAL);
  cdCanvasBox(canvas, 150, 200, h - 60, h);
  cdCanvasHatch(canvas, CD_CROSS);
  cdCanvasBox(canvas, 200, 250, h - 60, h);
  cdCanvasHatch(canvas, CD_DIAGCROSS);
  cdCanvasBox(canvas, 250, 300, h - 60, h);

  /* Draw 4 regions, in diamond shape,
     at top, bottom, left, right,
     using different interior styles. */

  /* At top, not filled polygon, notice that the last line style is used. */
  cdCanvasBegin(canvas, CD_CLOSED_LINES);
  cdCanvasVertex(canvas, w / 2, h - 100);
  cdCanvasVertex(canvas, w / 2 + 50, h - 150);
  cdCanvasVertex(canvas, w / 2, h - 200);
  cdCanvasVertex(canvas, w / 2 - 50, h - 150);
  cdCanvasEnd(canvas);

  /* At left, hatch filled polygon */
  cdCanvasHatch(canvas, CD_DIAGCROSS);
  cdCanvasBegin(canvas, CD_FILL);
  cdCanvasVertex(canvas, 100, h / 2);
  cdCanvasVertex(canvas, 150, h / 2 + 50);
  cdCanvasVertex(canvas, 200, h / 2);
  cdCanvasVertex(canvas, 150, h / 2 - 50);
  cdCanvasEnd(canvas);

  /* At right, pattern filled polygon */
  cdCanvasPattern(canvas, STYLE_SIZE, STYLE_SIZE, pattern);
  cdCanvasBegin(canvas, CD_FILL);
  cdCanvasVertex(canvas, w - 100, h / 2);
  cdCanvasVertex(canvas, w - 150, h / 2 + 50);
  cdCanvasVertex(canvas, w - 200, h / 2);
  cdCanvasVertex(canvas, w - 150, h / 2 - 50);
  cdCanvasEnd(canvas);

  /* At bottom, stipple filled polygon */
  cdCanvasStipple(canvas, STYLE_SIZE, STYLE_SIZE, stipple);
  cdCanvasBegin(canvas, CD_FILL);
  cdCanvasVertex(canvas, w / 2, 100);
  cdCanvasVertex(canvas, w / 2 + 50, 150);
  cdCanvasVertex(canvas, w / 2, 200);
  cdCanvasVertex(canvas, w / 2 - 50, 150);
  cdCanvasEnd(canvas);

  /* Draw two beziers at bottom-left */
  cdCanvasBegin(canvas, CD_BEZIER);
  cdCanvasVertex(canvas, 100, 100);
  cdCanvasVertex(canvas, 150, 200);
  cdCanvasVertex(canvas, 180, 250);
  cdCanvasVertex(canvas, 180, 200);
  cdCanvasVertex(canvas, 180, 150);
  cdCanvasVertex(canvas, 150, 100);
  cdCanvasVertex(canvas, 300, 100);
  cdCanvasEnd(canvas);

  /* Initialize the image buffer contents */
  //#define IMAGE_SIZE 16
  memset(red, 0xFF, IMAGE_SIZE*IMAGE_SIZE / 2);
  memset(green, 0x5F, IMAGE_SIZE*IMAGE_SIZE / 2);
  memset(blue, 0x5F, IMAGE_SIZE*IMAGE_SIZE / 2);
  memset(red + IMAGE_SIZE*IMAGE_SIZE / 2, 0x5F, IMAGE_SIZE*IMAGE_SIZE / 2);
  memset(green + IMAGE_SIZE*IMAGE_SIZE / 2, 0x8F, IMAGE_SIZE*IMAGE_SIZE / 2);
  memset(blue + IMAGE_SIZE*IMAGE_SIZE / 2, 0x5F, IMAGE_SIZE*IMAGE_SIZE / 2);
  memset(red + IMAGE_SIZE*(IMAGE_SIZE - 1), 0, IMAGE_SIZE);
  memset(green + IMAGE_SIZE*(IMAGE_SIZE - 1), 0, IMAGE_SIZE);
  memset(blue + IMAGE_SIZE*(IMAGE_SIZE - 1), 0, IMAGE_SIZE);
  memset(red, 0, IMAGE_SIZE);
  memset(green, 0, IMAGE_SIZE);
  memset(blue, 0, IMAGE_SIZE);
  {
    int i, offset;
    for (i = 0; i < IMAGE_SIZE; i++)
    {
      offset = i*IMAGE_SIZE;
      red[offset] = 0;
      green[offset] = 0;
      blue[offset] = 0;
      red[offset + IMAGE_SIZE - 1] = 0;
      green[offset + IMAGE_SIZE - 1] = 0;
      blue[offset + IMAGE_SIZE - 1] = 0;
    }
  }

//    cdCanvasPutImageRectRGB(canvas, IMAGE_SIZE, IMAGE_SIZE, red, green, blue, w - 400, h - 310, IMAGE_SIZE, IMAGE_SIZE, 0, 0, 0, 0);
  //  cdCanvasPutImageRectRGBA(canvas, IMAGE_SIZE, IMAGE_SIZE, red, green, blue, alpha, w - 400, h - 310, IMAGE_SIZE, IMAGE_SIZE, 0, 0, 0, 0);
  //  cdCanvasPutImageRectRGB(canvas, IMAGE_SIZE, IMAGE_SIZE, red, green, blue,         w - 400, h - 310, 3*IMAGE_SIZE, 3*IMAGE_SIZE, 0, 0, 0, 0);
  cdCanvasPutImageRectRGBA(canvas, IMAGE_SIZE, IMAGE_SIZE, red, green, blue, alpha, w - 400, h - 310, 3*IMAGE_SIZE, 3*IMAGE_SIZE, 0, 0, 0, 0);

  if (contextplus || cdCanvasGetAttribute(canvas, "LINEARGRADIENT") != NULL)
//  if (0)
  {
    cdCanvasForeground(canvas, CD_YELLOW);
    cdCanvasBackground(canvas, CD_GREEN);
    cdCanvasSetfAttribute(canvas, "LINEARGRADIENT", "%d 200 %d 250", w - 150, w - 100);  /* x1 y1 x2 y2 */
    cdCanvasBox(canvas, w - 150, w - 100, 200, 250);

    cdCanvasSetfAttribute(canvas, "RADIALGRADIENT", "%d %d 50", w - 125, 325);
    cdCanvasSector(canvas, w - 125, 325, 50, 50, 0, 360);
  }

  cdCanvasSetAttribute(canvas, "ROTATE", NULL);
  //cdCanvasSetAttribute(canvas, "ANTIALIAS", "0");
  if (use_transform)
    cdCanvasTransform(canvas, NULL);
  cdCanvasClip(canvas, CD_CLIPOFF);
}

void DrawVectorTextBox(cdCanvas* canvas, int x, int y, char* text)
{
  int rect[8], draw_box;

  cdCanvasLineWidth(canvas, 1);
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);

  draw_box = 0;
  if (draw_box)
  {
    int xmin, xmax, ymin, ymax;
    cdCanvasGetVectorTextBox(canvas, x, y, text, &xmin, &xmax, &ymin, &ymax);
    cdCanvasForeground(canvas, CD_GREEN);
    cdCanvasRect(canvas, xmin, xmax, ymin, ymax);

    if (cdCanvasTextOrientation(canvas, CD_QUERY) == 0)
    {
      cdCanvasForeground(canvas, CD_RED);
      cdCanvasLine(canvas, xmin, y, xmax, y);
    }
  }
  else
  {
    /* bounding box */
    cdCanvasGetVectorTextBounds(canvas, text, x, y, rect);
    cdCanvasForeground(canvas, CD_GREEN);
    cdCanvasBegin(canvas, CD_CLOSED_LINES);
    cdCanvasVertex(canvas, rect[0], rect[1]);
    cdCanvasVertex(canvas, rect[2], rect[3]);
    cdCanvasVertex(canvas, rect[4], rect[5]);
    cdCanvasVertex(canvas, rect[6], rect[7]);
    cdCanvasEnd(canvas);
  }

  /* reference point */
  cdCanvasForeground(canvas, CD_BLUE);
  cdCanvasMarkType(canvas, CD_PLUS);
  cdCanvasMarkSize(canvas, 30);
  cdCanvasMark(canvas, x, y);

  cdCanvasForeground(canvas, CD_BLACK);
  cdCanvasVectorText(canvas, x, y, text);
}

void DrawTextBox(cdCanvas* canvas, int x, int y, char* text)
{
  int rect[8], draw_box;

  cdCanvasLineWidth(canvas, 1);
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);

  draw_box = 0;
  if (draw_box)
  {
    int xmin, xmax, ymin, ymax;
    cdCanvasGetTextBox(canvas, x, y, text, &xmin, &xmax, &ymin, &ymax);
    cdCanvasRect(canvas, xmin, xmax, ymin, ymax);

    if (cdCanvasTextOrientation(canvas, CD_QUERY) == 0)
    {
      cdCanvasForeground(canvas, CD_RED);
      cdCanvasLine(canvas, xmin, y, xmax, y);
    }
  }
  else
  {
    /* bounding box */
    cdCanvasGetTextBounds(canvas, x, y, text, rect);
    cdCanvasForeground(canvas, CD_GREEN);
    cdCanvasBegin(canvas, CD_CLOSED_LINES);
    cdCanvasVertex(canvas, rect[0], rect[1]);
    cdCanvasVertex(canvas, rect[2], rect[3]);
    cdCanvasVertex(canvas, rect[4], rect[5]);
    cdCanvasVertex(canvas, rect[6], rect[7]);
    cdCanvasEnd(canvas);
  }

  /* reference point */
  cdCanvasForeground(canvas, CD_BLUE);
  cdCanvasMarkType(canvas, CD_PLUS);
  cdCanvasMarkSize(canvas, 30);
  cdCanvasMark(canvas, x, y);

  cdCanvasForeground(canvas, CD_BLACK);
  cdCanvasText(canvas, x, y, text);
}

void SimpleDrawTextAlign(cdCanvas* canvas)
{
  int w, h, i, xoff, yoff, use_vector;

  int text_aligment[] = {
    CD_NORTH,
    CD_SOUTH,
    CD_EAST,
    CD_WEST,
    CD_NORTH_EAST,
    CD_NORTH_WEST,
    CD_SOUTH_EAST,
    CD_SOUTH_WEST,
    CD_CENTER,
    CD_BASE_CENTER,
    CD_BASE_RIGHT,
    CD_BASE_LEFT
  };

#if 1
  char* text_aligment_str[] = {
    "North (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "South (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "East (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "West (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "North East (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "North West (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "South East (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "South West (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "Center (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "Base Center (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "Base Right (Ãyj)\nSecond Line (Ãyj)\nThird Line",
    "Base Left (Ãyj)\nSecond Line (Ãyj)\nThird Line"
  };
#else
  char* text_aligment_str[] = {
    "North (Ãyj)",
    "South (Ãyj)",
    "East (Ãyj)",
    "West (Ãyj)",
    "North East (Ãyj)",
    "North West (Ãyj)",
    "South East (Ãyj)",
    "South West (Ãyj)",
    "Center (Ãyj)",
    "Base Center (Ãyj)",
    "Base Right (Ãyj)",
    "Base Left (Ãyj)"
  };
#endif

  cdCanvasGetSize(canvas, &w, &h, 0, 0);

  cdCanvasBackground(canvas, CD_WHITE);
  cdCanvasClear(canvas);

  use_vector = 0;

#if 0
  if (use_vector)
    cdCanvasVectorTextDirection(canvas, 0, 0, 1, 1);
  else
    cdCanvasTextOrientation(canvas, 45);
#endif

  xoff = w / 4;
  yoff = h / 7;

  if (use_vector)
    cdCanvasVectorCharSize(canvas, 30);
  else
  {
    cdCanvasFont(canvas, "Times", CD_PLAIN, 14);
    //cdCanvasFont(canvas, "Helvetica", CD_PLAIN, 24);
  }

  for (i = 0; i < 12; i++)
  {
    cdCanvasTextAlignment(canvas, text_aligment[i]);
    if (i < 6)
    {
      if (use_vector)
        DrawVectorTextBox(canvas, xoff, yoff*(i + 1), text_aligment_str[i]);
      else
        DrawTextBox(canvas, xoff, yoff*(i + 1), text_aligment_str[i]);
    }
    else
    {
      if (use_vector)
        DrawVectorTextBox(canvas, 3 * xoff, yoff*(i - 5), text_aligment_str[i]);
      else
        DrawTextBox(canvas, 3 * xoff, yoff*(i - 5), text_aligment_str[i]);
    }
  }
}

void DrawTextFont(cdCanvas* canvas, const char* font, int size, int xoff, int yoff, char* text)
{
  cdCanvasFont(canvas, font, CD_PLAIN, size);
  DrawTextBox(canvas, xoff, yoff, text);

  cdCanvasFont(canvas, font, CD_BOLD, size);
  DrawTextBox(canvas, 2 * xoff, yoff, text);

  cdCanvasFont(canvas, font, CD_ITALIC, size);
  DrawTextBox(canvas, 3 * xoff, yoff, text);

  cdCanvasFont(canvas, font, CD_BOLD_ITALIC, size);
  DrawTextBox(canvas, 4 * xoff, yoff, text);
}

void SimpleDrawTextFonts(cdCanvas* canvas)
{
  int xoff, yoff, size;

  cdCanvasBackground(canvas, CD_WHITE);
  cdCanvasClear(canvas);

  xoff = 470;
  yoff = 150;
  //  size = -30;
  xoff = 270;
  size = 10;

  cdCanvasTextAlignment(canvas, CD_BASE_CENTER);

  //cdCanvasTextOrientation(canvas, 45);

  cdCanvasTextOrientation(canvas, 0);
  DrawTextFont(canvas, "Courier", size, xoff, yoff, "Courier");
  DrawTextFont(canvas, "Times", size, xoff, 2 * yoff, "Times");
  cdCanvasTextOrientation(canvas, 45);
  DrawTextFont(canvas, "Helvetica", size, xoff, 3 * yoff, "Helvetica");
  DrawTextFont(canvas, "Calibri", size, xoff, 4 * yoff, "Calibri");
  DrawTextFont(canvas, "Segoe UI", size, xoff, 5 * yoff, "Segoe UI");

  {
    //    static char native[50] = "Tecmedia, -60";
    //    static char native[50] = "-*-helvetica-medium-r-*-*-8-*";
    //    static char native[50] = "Edwardian Script ITC, 24";
    //    cdSetAttribute("ADDFONTMAP","Edwardian Script ITC=ITCEDSCR");

    //    char native[50] = "Book Antiqua, 24";
    //    cdSetAttribute("ADDFONTMAP", "Book Antiqua=BKANT");

    //    cdNativeFont("-d");
    //    cdNativeFont(native);
    //    DrawTextBox(xoff, yoff, native);
    //    DrawTextBox(xoff, yoff, "The quick brown fox.");
  }

  //cdNativeFont("Tecmedia, 36");

  //cdSetAttribute("ADDFONTMAP", "WingDings=WingDing");
  //cdNativeFont("WingDings, 36");

  //cdText(500, 50, "X");
  //cdText(500, 50, "abcdefghijklmnopqrstuvwxyz");
  //cdText(500, 150, "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  //cdText(500, 250, "1234567890");
  //cdText(500, 350, "'\"!@#$%¨&*()_+-=[]^/;.,");

  //cdFont(CD_COURIER, 0, 22);
  //cdText(10, 60, "abcdefghijklmnopqrstuvwxyz");
  //cdText(10, 160, "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  //cdText(10, 260, "1234567890");
  //cdText(500, 360, "'\"!@#$%¨&*()_+-=[]^/;.,");
}

//void SimpleDrawTest(cdCanvas* canvas)
void SimpleDrawTest1(cdCanvas* canvas)
{
  long pattern_p[16];  /* 4x4 pattern */
  int w, h;
  int xmin, xmax, ymin, ymax;

  /* notice that if we are not using world coordinates
     it is harder to position all the objetcs we want. */
  cdCanvasGetSize(canvas, &w, &h, 0, 0);

  cdCanvasBackground(canvas, CD_WHITE);
  cdCanvasClear(canvas);

  /* pattern_p initialization */
  pattern_p[0] = CD_RED;    pattern_p[1] = CD_RED;    /* first line */
  pattern_p[2] = CD_YELLOW; pattern_p[3] = CD_YELLOW;
  pattern_p[4] = CD_RED;    pattern_p[5] = CD_RED;    /* second line */
  pattern_p[6] = CD_YELLOW; pattern_p[7] = CD_YELLOW;
  pattern_p[8] = CD_YELLOW; pattern_p[9] = CD_YELLOW; /* third line */
  pattern_p[10] = CD_YELLOW; pattern_p[11] = CD_YELLOW;
  pattern_p[12] = CD_YELLOW; pattern_p[13] = CD_YELLOW; /* fourth line */
  pattern_p[14] = CD_YELLOW; pattern_p[15] = CD_YELLOW;

  /* set the line attributes */
  cdCanvasLineWidth(canvas, 4);
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);

  /* in the center draw a pattern_p pizza
     with a slice missing */
  cdCanvasPattern(canvas, 4, 4, pattern_p);
  cdCanvasSector(canvas, w / 2, h / 2, w / 2, h / 2, 45, 0);
  /* draws a dark red border */
  cdCanvasForeground(canvas, CD_DARK_RED);
  cdCanvasInteriorStyle(canvas, CD_HOLLOW);
  cdCanvasSector(canvas, w/2, h/2, w/2, h/2, 45, 0);

  /* on the left a red hash diamond */
  /* notice the the default back opacity is transparent
     and the pattern of the sector will still be visible
     inside the hatch where the two objects intersect */
  cdCanvasForeground(canvas, CD_RED);
  cdCanvasHatch(canvas, CD_DIAGCROSS); 
  cdCanvasBegin(canvas, CD_FILL);
  cdCanvasVertex(canvas, w/4, h/4); 
  cdCanvasVertex(canvas, w/2-w/8, h/2); 
  cdCanvasVertex(canvas, w/4, 3*h/4); 
  cdCanvasVertex(canvas, w/8, h/2); 
  cdCanvasEnd(canvas);

  /* draws a blue roof.*/
  cdCanvasForeground(canvas, CD_BLUE);
  cdCanvasLine(canvas, w/8, h/2, w/4, 3*h/4);
  cdCanvasLine(canvas, w/4, 3*h/4, w/2-w/8, h/2);

  /* draws a dashed ribbon on the right 
     with a custom color */
  cdCanvasForeground(canvas, cdEncodeColor(100, 25, 200));
  cdCanvasLineStyle(canvas, CD_DASH_DOT);
  cdCanvasBegin(canvas, CD_BEZIER);
  cdCanvasVertex(canvas, 3*w/4-20, h/2-50); 
  cdCanvasVertex(canvas, 3*w/4+150, 3*h/4-50); 
  cdCanvasVertex(canvas, 3*w/4-150, 3*h/4-50); 
  cdCanvasVertex(canvas, 3*w/4+20, h/2-50); 
  cdCanvasEnd(canvas);

  cdCanvasFont(canvas, "Helvetica", CD_BOLD, 40);
  cdCanvasTextAlignment(canvas, CD_CENTER);
  cdCanvasText(canvas, w/2, h/4-50, "Canvas Draw");
  cdCanvasGetTextBox(canvas, w/2, h/4-50, "Canvas Draw", &xmin, &xmax, &ymin, &ymax);
  cdCanvasRect(canvas, xmin, xmax, ymin, ymax);
}

#if 0
#include <cd_old.h>

void draw_wd(void)
{
  char* text;
  double x, y;
  double rect[8];

  cdBackground(CD_WHITE);
  cdClear();
  cdLineStyle(CD_CONTINUOUS);
  cdLineWidth(1);

  //  wdVectorTextDirection(0, 0, 1, 1);
  cdTextAlignment(CD_NORTH_WEST);

  //  text = "Vector Text";
  text = "Vector Text\nSecond Line\nThird Line";
  x = 0.25;
  y = 0.40;

  cdForeground(CD_BLACK);
  wdLine(0, 0, 1, 1);
  wdLine(0, 1, 1, 0);

  cdForeground(CD_GREEN);
  cdMarkType(CD_STAR);
  wdMark(x, y);

  cdForeground(CD_BLUE);
  wdVectorCharSize(0.1);
  wdVectorText(x, y, text);

  cdForeground(CD_RED);
  wdGetVectorTextBounds(text, x, y, rect);
  cdBegin(CD_CLOSED_LINES);
  wdVertex(rect[0], rect[1]);
  wdVertex(rect[2], rect[3]);
  wdVertex(rect[4], rect[5]);
  wdVertex(rect[6], rect[7]);
  cdEnd();
}

void SimpleDrawTest(cdCanvas* canvas)
//void SimpleDrawTestHardCopy(cdCanvas* canvas)
{
  int w, h;

  cdActivate(canvas);

  cdGetCanvasSize(&w, &h, 0, 0);

  wdViewport(0,w-1,0,h-1);
  if (w>h)
    wdWindow(0,(double)w/(double)h,0,1);
  else
    wdWindow(0,1,0,(double)h/(double)w);

  draw_wd();

  //wdHardcopy(CD_CLIPBOARD, "800x600", cdActiveCanvas(), draw_wd );
}

//void SimpleDrawTest(cdCanvas* canvas)
void SimpleDrawTestImageRGB1(cdCanvas* canvas)
{
  int size = 2048*2048;
  unsigned char *red_p, *green_p, *blue_p;
  cdCanvas* canvas_p = cdCreateCanvas(CD_IMAGERGB, "2048x2048");
  cdCanvasActivate(canvas_p);

  red_p = calloc(size, 1);
  green_p = calloc(size, 1);
  blue_p = calloc(size, 1);

  cdCanvasPutImageRectRGB(canvas_p, 2048, 2048, red_p, green_p, blue_p, 0, 3, 2048, 2017, 0, 2047, 3, 2020);

  free(red_p);
  free(green_p);
  free(blue_p);

  cdKillCanvas(canvas_p);
}

//void SimpleDrawTest(cdCanvas* canvas)
void SimpleDrawVectorFont(cdCanvas* canvas)
{
  cdActivate(canvas);

  cdBackground(CD_WHITE);
  cdClear();
  cdLineStyle(CD_CONTINUOUS);
  cdLineWidth(1);

  //  wdVectorText(0.1, 0.4, "ãõñç áéíóú àèìòù âêîôû äëïöü");
  //  wdVectorText(0.1, 0.2, "ÃÕÑÇ ÁÉÍÓÚ ÀÈÌÒÙ ÂÊÎÔÛ ÄËÏÖÜ ");
  cdVectorFont("../etc/vectorfont26.txt"); /* original Simplex II */
  {
    int i;
    char t[2];
    char s[10];
    int x = 10;
    int y = 600;
    t[1] = 0;
    cdFont(CD_COURIER, CD_BOLD, 14);
    cdVectorCharSize(25);
    for (i = 128; i < 256; i++)
    {
      int dx = 30;
      t[0] = (char)i;
      sprintf(s, "%3d", i);
      cdForeground(CD_DARK_RED);
      cdText(x, y, s);
      //      cdText(x+dx, y, t);
      cdForeground(CD_BLACK);
      cdVectorText(x+2*dx-10, y, t);

      x += 3*dx;
      if ((i+1) % 8 == 0)
      {
        x = 10;
        y -= 30;
      }
    }
    //cdFont(CD_TIMES_ROMAN, CD_PLAIN, 24);
    //cdVectorCharSize(24);
    //  for (i = 192; i < 256; i++)
    //  {
    //    int dx = 92;
    //    t[0] = (char)i;
    //    sprintf(s, "%d", i);
    //    cdText(x, y, s);
    //    cdText(x+dx, y, t);
    //    cdVectorText(x+2*dx, y, t);
    //    
    //    x += 3*dx + 2*dx/3;
    //    if ((i+1) % 4 == 0)
    //    {
    //      x = 30;
    //      y += 52;
    //    }
    //  }
  }
}

//void SimpleDrawTest(cdCanvas* canvas)
void SimpleDrawPoly(cdCanvas* canvas)
{
  int w, h;

  cdActivate(canvas);

  cdGetCanvasSize(&w, &h, 0, 0);

  cdBackground(CD_WHITE);
  cdClear();

  //cdSetAttribute("ANTIALIAS", "0");
  cdForeground(cdEncodeAlpha(cdEncodeColor(255, 0, 0), 100));

  cdfCanvasArc(cdActiveCanvas(), 255, 255, 100, 100, 0, 360);

  cdLine(0, 0, 200, 200);

  cdBegin(CD_BEZIER);
  cdVertex(100, 100);
  cdVertex(150, 200);
  cdVertex(180, 250);
  cdVertex(180, 200);
  cdVertex(180, 150);
  cdVertex(150, 100);
  cdVertex(300, 100);
  cdEnd();


  cdEnd();
}
#endif

//void SimpleDrawTest(cdCanvas* canvas)
void SimpleDrawTestImageRGB(cdCanvas* canvas1)
{
  //cdCanvas* canvas = cdCreateCanvas(CD_IMAGERGB, "570x569");
  cdCanvas* canvas = canvas1;

  if (canvas == canvas1)
  {
    cdCanvasClear(canvas);
    cdCanvasInteriorStyle(canvas, CD_SOLID);
  }

  cdCanvasBegin(canvas, CD_FILL);
  cdCanvasVertex(canvas, 279, 81);
  cdCanvasVertex(canvas, 280, 81);
  cdCanvasVertex(canvas, 281, 81);
  cdCanvasVertex(canvas, 282, 81);
  cdCanvasVertex(canvas, 283, 81);
  cdCanvasVertex(canvas, 284, 82);
  cdCanvasVertex(canvas, 284, 81);
  cdCanvasVertex(canvas, 286, 81);
  cdCanvasVertex(canvas, 287, 81);
  cdCanvasVertex(canvas, 288, 80);
  cdCanvasVertex(canvas, 288, 79);
  cdCanvasVertex(canvas, 288, 77);
  cdCanvasVertex(canvas, 288, 76);
  cdCanvasVertex(canvas, 286, 76);
  cdCanvasVertex(canvas, 285, 76);
  cdCanvasVertex(canvas, 284, 76);
  cdCanvasVertex(canvas, 283, 76);
  cdCanvasVertex(canvas, 282, 77);
  cdCanvasVertex(canvas, 281, 77);
  cdCanvasVertex(canvas, 279, 78);
  cdCanvasVertex(canvas, 279, 79);
  cdCanvasVertex(canvas, 280, 80);
  cdCanvasVertex(canvas, 280, 81);
  cdCanvasVertex(canvas, 279, 81);
  cdCanvasEnd(canvas);

  if (canvas != canvas1)
    cdKillCanvas(canvas);
}

void SimpleDrawTest2(cdCanvas* canvas)
{
  cdCanvasBackground(canvas, CD_WHITE);
  cdCanvasClear(canvas);

  cdCanvasInteriorStyle(canvas, CD_SOLID);
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);
  cdCanvasLineWidth(canvas, 1);

  cdCanvasForeground(canvas, CD_RED);
  cdCanvasRect(canvas, 100, 200, 100, 200);
  cdCanvasForeground(canvas, CD_GREEN);
  cdCanvasBox(canvas, 100, 200, 100, 200);

  cdCanvasForeground(canvas, CD_GREEN);
  cdCanvasBox(canvas, 300, 400, 100, 200);
  cdCanvasForeground(canvas, CD_RED);
  cdCanvasRect(canvas, 300, 400, 100, 200);
}

void cdCanvasSectorXY(cdCanvas* canvas, int x1, int y1, int x2, int y2, double angle1, double angle2)
{
  int xc = (x2 + x1) / 2;
  int yc = (y2 + y1) / 2;
  int w = x2 - x1 + 1;
  int h = y2 - y1 + 1;
  cdCanvasSector(canvas, xc, yc, w, h, angle1, angle2);
}

//void SimpleDrawTest(cdCanvas* canvas)
void SimpleDrawTestPrecision(cdCanvas* canvas)
{
  /* white background */
  cdCanvasBackground(canvas, CD_WHITE);
  cdCanvasClear(canvas);

  cdCanvasInteriorStyle(canvas, CD_SOLID);
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);
  cdCanvasLineWidth(canvas, 1);

  /* Guide Lines */
  cdCanvasForeground(canvas, cdEncodeColor(255, 0, 0));
  cdCanvasLine(canvas, 10, 5, 10, 19);
  cdCanvasLine(canvas, 14, 5, 14, 19);
  cdCanvasLine(canvas, 5, 10, 19, 10);
  cdCanvasLine(canvas, 5, 14, 19, 14);

  /* Stroke Rectangle, must cover guide lines */
  cdCanvasForeground(canvas, cdEncodeColor(0, 0, 0));
  cdCanvasRect(canvas, 10, 14, 10, 14);

  /* Guide Lines */
  cdCanvasForeground(canvas, cdEncodeColor(255, 0, 0));
  cdCanvasLine(canvas, 10, 5 + 30, 10, 19 + 30);
  cdCanvasLine(canvas, 14, 5 + 30, 14, 19 + 30);
  cdCanvasLine(canvas, 5, 10 + 30, 19, 10 + 30);
  cdCanvasLine(canvas, 5, 14 + 30, 19, 14 + 30);

  /* Fill Rectangle, must cover guide lines */
  cdCanvasForeground(canvas, cdEncodeColor(0, 0, 0));
  cdCanvasBox(canvas, 10, 14, 10 + 30, 14 + 30);

  cdCanvasForeground(canvas, cdEncodeColor(255, 0, 0));
  cdCanvasBox(canvas, 30, 50, 10, 30);

  cdCanvasForeground(canvas, cdEncodeColor(0, 0, 0));
  cdCanvasSectorXY(canvas, 30, 10, 50, 30, 0, 360);

  cdCanvasForeground(canvas, cdEncodeColor(255, 0, 0));
  cdCanvasBox(canvas, 60, 80, 10, 30);

  cdCanvasForeground(canvas, cdEncodeColor(0, 0, 0));
  cdCanvasSectorXY(canvas, 60, 10, 80, 30, 0, 360);


  cdCanvasForeground(canvas, cdEncodeColor(255, 0, 0));
  cdCanvasBox(canvas, 30, 50, 10 + 30, 30 + 30);

  cdCanvasForeground(canvas, cdEncodeColor(0, 0, 0));
  cdCanvasSectorXY(canvas, 30, 10 + 30, 50, 30 + 30, 45, 135);

  cdCanvasForeground(canvas, cdEncodeColor(255, 0, 0));
  cdCanvasBox(canvas, 60, 80, 10 + 30, 30 + 30);

  cdCanvasForeground(canvas, cdEncodeColor(0, 0, 0));
  cdCanvasSectorXY(canvas, 60, 10 + 30, 80, 30 + 30, 45, 135);

  cdCanvasLine(canvas, 20, 70 - 2, 20, 70 + 2);
  cdCanvasLine(canvas, 20 - 2, 70, 20 + 2, 70);

  cdCanvasTextAlignment(canvas, CD_SOUTH_WEST);
  cdCanvasFont(canvas, "Helvetica", CD_PLAIN, -80);
  cdCanvasText(canvas, 20, 70, "Text");
}

void SimpleDrawTest3(cdCanvas* canvas)
{
  cdCanvas *pdf;                                                                                                               
  int canvas_height;
  int canvas_width;
  int margin;
  char *text;
  int height;
  int width;

  pdf = cdCreateCanvasf(CD_PDF, "\"%s\" -p%d", "tmp.pdf", CD_A4);

  cdCanvasSetAttribute(pdf, "UTF8MODE", "1");
  cdCanvasActivate(pdf);

  margin = 5;
  cdCanvasGetSize(pdf, &canvas_width, &canvas_height, NULL, NULL);

  text = "F = ((1 / 4 * π * ε0) * (Q1 * Q2) / R²";

  cdCanvasFont(pdf, "Courier", CD_PLAIN, 12);

  cdCanvasGetTextSize(pdf, text, &width, &height);
  cdCanvasText(pdf, margin, canvas_height - margin - height, text);

  cdCanvasDeactivate(pdf);
  cdKillCanvas(pdf);
}

void SimpleDrawTest4(cdCanvas* canvas)
{
  int w, h;

  cdCanvasClear(canvas);

  cdCanvasForeground(canvas, CD_BLACK);
  cdCanvasLineStyle(canvas, CD_CONTINUOUS);
  cdCanvasLineWidth(canvas, 1);
  cdCanvasFont(canvas, "Tahoma", CD_PLAIN, 20);

  cdCanvasTextAlignment(canvas, CD_BASE_LEFT);

  cdCanvasText(canvas, 200, 400, "It's a test");

  cdCanvasForeground(canvas, CD_RED);
  cdCanvasGetTextSize(canvas, " ", &w, &h);
  cdCanvasGetTextSize(canvas, "It's a test", &w, &h);

  cdCanvasRect(canvas, 200, 200 + w, 400, 400 + h);
}

#if 1

#include <Windows.h>

void SimpleDrawTest(cdCanvas* canvas)
{
  cdCanvas * cn;
  int w, h, pixel_size;
  HDC hDC;
  SIZE sz;
  RECT rect;
  POINT pnt;
  HPEN pen;
  HFONT font;

  //cn = canvas;
  cn = cdCreateCanvas(CD_PRINTER, "%s -d");
  cdCanvasForeground(cn, CD_WHITE);

  /* Generic initiation */

  cdCanvasActivate(cn);
  cdCanvasBackground(cn, CD_WHITE);
  cdCanvasClear(cn);

  cdCanvasForeground(cn, CD_BLACK);
  cdCanvasLineStyle(cn, CD_CONTINUOUS);
  cdCanvasLineWidth(cn, 1);

  /* Getting native hDC */
  hDC = (HDC)cdCanvasGetAttribute(cn, "HDC");

  /* Creating Font with identical properties to the CD canvas one */
  pixel_size  = -MulDiv(20, GetDeviceCaps(hDC, LOGPIXELSY), 72);    /* -27 */
  font = CreateFontA(pixel_size, 0, 0, 0, FW_NORMAL, 0, 0, 0, 
                     RUSSIAN_CHARSET, OUT_TT_ONLY_PRECIS, CLIP_DEFAULT_PRECIS, ANTIALIASED_QUALITY, DEFAULT_PITCH, 
                     "Tahoma");

  /*
  hFont = CreateFontW(-27, 0, 0, 0, 0, 0, 0, 0,
                     DEFAULT_CHARSET, OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, FF_DONTCARE | DEFAULT_PITCH,
                     L"Tahoma");
  */
  SelectObject(hDC, font);

  /* Drawing test string temporarily ignoring different Y-axizs direction */
  /* In my work I internally use ANSI and only input/output with UTF8 so don't mind I'm using A API-functions */

  rect.left = 500;
  rect.bottom = 300;
  rect.right = 900;
  rect.top = 200;

  SetTextAlign(hDC, TA_LEFT|TA_TOP);

  ExtTextOutA(hDC, 500, 300, 0, &rect, "It's a test", 11, NULL);

  /* Getting string image extent */

  GetTextExtentPoint32A(hDC, "It's a test", 11, &sz);

  /* Creating thin red pen for visualization of text extent */

  pen = CreatePen(PS_SOLID, 1, 0x000000FF);
  SelectObject(hDC, pen);

  /* Drawing visualization of text extent itself */

  MoveToEx(hDC, 500, 300, &pnt);
  LineTo(hDC, 500 + sz.cx, 300);
  LineTo(hDC, 500 + sz.cx, 300 + sz.cy);
  LineTo(hDC, 500, 300 + sz.cy);
  LineTo(hDC, 500, 300);

  /* Now doing the same sequence but using CD canvas.
  Due to different Y-axis direction results will be in lower part of the page */

  cdCanvasFont(cn, "Tahoma", CD_PLAIN, 20);
  // cdCanvasNativeFont(cn, "-d");
  cdCanvasTextAlignment(cn, CD_SOUTH_WEST);

  cdCanvasText(cn, 200, 400, "It's a test");

  cdCanvasForeground(cn, CD_RED);
  cdCanvasGetTextSize(cn, "It's a test", &w, &h);

  cdCanvasRect(cn, 200, 200 + w, 400, 400 + h);

  if (cn != canvas)
    cdKillCanvas(cn);
}

#endif
