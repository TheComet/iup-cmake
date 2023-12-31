/** \file
 * \brief Cairo Base Driver
 *
 * See Copyright Notice in cd.h
 */

#ifndef __CDCAIROCTX_H
#define __CDCAIROCTX_H

#include <cairo.h>
#include <pango/pango.h>

#include "cd.h"
#include "cd_private.h"

#ifndef __GTK_PRINT_UNIX_DIALOG_H__
typedef struct _GtkPrintJob  GtkPrintJob;
#endif

#define CAIRO_VERSION_110 CAIRO_VERSION_ENCODE(1, 10, 0)

struct _cdCtxImage {
  unsigned int w, h;
  double w_mm, h_mm;   /* size in mm                                  */                  
  double xres, yres;   /* resolution in pixels/mm                     */     
  int bpp;
  cairo_t* cr;
};

struct _cdCtxCanvas
{
  cdCanvas* canvas;

  cairo_t* cr;

  /* text attributes */
  PangoContext *fontcontext;
  PangoFontDescription *fontdesc;
  PangoLayout *fontlayout;

  /* fill attributes */
  cairo_pattern_t *pattern, *solid;
  int last_source;
  int hatchboxsize;

  /* custom attributes */

  int img_format;

  double rotate_angle;
  int rotate_center_x;
  int rotate_center_y;

  int poly_holes[500];
  int holes;

#if CAIRO_VERSION >= CAIRO_VERSION_110
  cairo_region_t* new_rgn;
#endif

  void* window;    /* used in NativeWindow in GDK */

#ifdef WIN32
  void* hWnd;        /* used in NativeWindow in Win32 */
  void* hDC;  
  int isOwnedDC;
#else
  void* dpy;         /* used in NativeWindow in X11 */
  unsigned long wnd;          
#endif

  int user_image;   /* used in ImageRGB */
  unsigned char *rgb;

  int eps;  /* used in PS */

  cdImage* image_dbuffer;       /* Used in double buffer driver */
  cdCanvas* canvas_dbuffer;
  int kill_dbuffer;

  GtkPrintJob* job;    /* used in Printer (GTK) */
  char* printername;   /* used in Printer (Win32) */
  int scale_points;
  double scale;

  char* utf8_buffer;
  int utf8mode;
};

#define cdCairoGetRed(_)   (((double)cdRed(_))/255.)
#define cdCairoGetGreen(_) (((double)cdGreen(_))/255.)
#define cdCairoGetBlue(_)  (((double)cdBlue(_))/255.)
#define cdCairoGetAlpha(_)  (((double)cdAlpha(_))/255.)

cdCtxCanvas *cdcairoCreateCanvas(cdCanvas* canvas, cairo_t* cr);
void cdcairoInitTable(cdCanvas* canvas);
void cdcairoKillCanvas(cdCtxCanvas *ctxcanvas);

#endif
