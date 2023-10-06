#include <stdlib.h>
#include <stdio.h>

#include <iup.h>
#include <iupgl.h>
#include <cd.h>
#include "cddirect2d.h"

#include "simple.h"


int cmdExit(void)
{
  return IUP_CLOSE;
}

void simple_loadled (void);

#ifdef USE_OPENGL
/* USE_OPENGL - add to linker:
cdgl
iupgl
ftgl
glu32
opengl32
*/

void SimpleUpdateSize(cdCanvas* cnv)
{
  Ihandle* canvas = IupGetHandle("SimpleCanvas");
  IupGLMakeCurrent(canvas);

  if (cnv)
  {
    int w, h;
    double res = IupGetDouble(NULL, "SCREENDPI") / 25.4;
    IupGetIntInt(canvas, "DRAWSIZE", &w, &h);

    cdCanvasSetfAttribute(cnv, "SIZE", "%dx%d %g", w, h, res);  /* no need to update resolution */
  }
}

void SimpleFlush(void)
{
  IupGLSwapBuffers(IupGetHandle("SimpleCanvas"));
}
#endif

int main(int argc, char** argv)
{
  IupOpen(&argc, &argv);                        

#ifdef USE_CONTEXTPLUS
  cdInitContextPlus();
#endif
#ifdef USE_OPENGL
  IupGLCanvasOpen();
#endif
#ifdef WIN32
  cdInitDirect2D();
#endif

  simple_loadled();
#ifdef USE_OPENGL
  {
    Ihandle* dialog = IupGetHandle("SimpleDialog");
    Ihandle* canvas = IupGetHandle("SimpleCanvas");
    IupDestroy(canvas);
    canvas = IupGLCanvas("SimpleRepaint");
    IupSetAttribute(canvas, "BUFFER", "DOUBLE");
    IupSetHandle("SimpleCanvas", canvas);
    IupAppend(dialog, canvas);
  }
#endif

  IupSetAttribute(IupGetHandle("SimpleDialog"), "SIZE", "HALFxHALF");
  IupSetAttribute(IupGetHandle("SimpleDialog"), "PLACEMENT", "MAXIMIZED");
  IupShow(IupGetHandle("SimpleDialog"));
  IupSetAttribute(IupGetHandle("SimpleDialog"), "USERSIZE", NULL);

  SimpleCreateCanvas((char*)IupGetHandle("SimpleCanvas"));

  IupSetFunction("cmdExit", (Icallback) cmdExit);

  IupSetFunction("SimplePlayClipboard", (Icallback) SimplePlayClipboard);
  IupSetFunction("SimplePlayCGMText", (Icallback) SimplePlayCGMText);
  IupSetFunction("SimplePlayCGMBin", (Icallback) SimplePlayCGMBin);
  IupSetFunction("SimplePlayMetafile", (Icallback) SimplePlayMetafile);
  IupSetFunction("SimplePlayWMF", (Icallback) SimplePlayWMF);
  IupSetFunction("SimplePlayEMF", (Icallback) SimplePlayEMF);

  IupSetFunction("SimpleDrawDebug", (Icallback) SimpleDrawDebug);
  IupSetFunction("SimpleDrawWindow", (Icallback) SimpleDrawWindow);
  IupSetFunction("SimpleDrawCGMText", (Icallback) SimpleDrawCGMText);
  IupSetFunction("SimpleDrawCGMBin", (Icallback) SimpleDrawCGMBin);
  IupSetFunction("SimpleDrawDXF", (Icallback) SimpleDrawDXF);
  IupSetFunction("SimpleDrawDGN", (Icallback) SimpleDrawDGN);
  IupSetFunction("SimpleDrawEMF", (Icallback) SimpleDrawEMF);
  IupSetFunction("SimpleDrawMetafile", (Icallback) SimpleDrawMetafile);
  IupSetFunction("SimpleDrawPDF", (Icallback) SimpleDrawPDF);
  IupSetFunction("SimpleDrawPS", (Icallback) SimpleDrawPS);
  IupSetFunction("SimpleDrawEPS", (Icallback) SimpleDrawEPS);
  IupSetFunction("SimpleDrawSVG", (Icallback) SimpleDrawSVG);
  IupSetFunction("SimpleDrawWMF", (Icallback)SimpleDrawWMF);
  IupSetFunction("SimpleDrawPPTX", (Icallback)SimpleDrawPPTX);
  IupSetFunction("SimpleDrawPrint", (Icallback)SimpleDrawPrint);
  IupSetFunction("SimpleDrawPrintDialog", (Icallback) SimpleDrawPrintDialog);
  IupSetFunction("SimpleDrawClipboardBitmap", (Icallback) SimpleDrawClipboardBitmap);
  IupSetFunction("SimpleDrawClipboardMetafile", (Icallback) SimpleDrawClipboardMetafile);
  IupSetFunction("SimpleDrawClipboardEMF", (Icallback) SimpleDrawClipboardEMF);
  IupSetFunction("SimpleDrawImage", (Icallback) SimpleDrawImage);
  IupSetFunction("SimpleDrawImageRGB", (Icallback) SimpleDrawImageRGB);
  IupSetFunction("SimpleDrawSimulate", (Icallback) SimpleDrawSimulate);

#ifdef USE_OPENGL
  IupSetFunction("SimpleDrawGL", (Icallback) SimpleDrawGL);
#endif
  IupSetFunction("SimpleDrawIupDraw", (Icallback)SimpleDrawIupDraw);
  IupSetFunction("SimpleDrawDirect2D", (Icallback)SimpleDrawDirect2D);

  IupSetFunction("SimpleNotXor", (Icallback) SimpleNotXor);
  IupSetFunction("SimpleXor", (Icallback) SimpleXor);
  IupSetFunction("SimpleReplace", (Icallback) SimpleReplace);
  IupSetFunction("SimpleClippingOff", (Icallback) SimpleClippingOff);
  IupSetFunction("SimpleClippingArea", (Icallback) SimpleClippingArea);
  IupSetFunction("SimpleClippingPolygon", (Icallback) SimpleClippingPolygon);
  IupSetFunction("SimpleClippingRegion", (Icallback) SimpleClippingRegion);
  IupSetFunction("SimpleContextPlus", (Icallback) SimpleContextPlus);
  IupSetFunction("SimpleTransform", (Icallback) SimpleTransform);

  IupSetFunction("SimpleAll", (Icallback) SimpleAll);
  IupSetFunction("SimpleTextAlign", (Icallback)SimpleTextAlign);
  IupSetFunction("SimpleTextFonts", (Icallback)SimpleTextFonts);
  IupSetFunction("SimpleTest", (Icallback) SimpleTest);

  IupSetFunction("SimpleRepaint", (Icallback) SimpleRepaint);

  SimpleDrawWindow();

#ifdef USE_OPENGL
  SimpleUpdateSize(NULL);
  IupUpdate(IupGetHandle("SimpleCanvas"));
#endif

  IupMainLoop();

  SimpleKillCanvas();

#ifdef USE_CONTEXTPLUS
  cdFinishContextPlus();
#endif
#ifdef WIN32
  cdFinishDirect2D();
#endif

  IupClose();

  return EXIT_SUCCESS;
}
