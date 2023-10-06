/* IM 3 sample that shows an image.

  Needs "im.lib", "iup.lib", "cd.lib" and "iupcd.lib".

  Usage: im_view <file_name>

    Example: im_view test.tif

    
  Click on image to open another file.
*/

#include <im_plus.h>
#include <cd_plus.h>
#include <iup_plus.h>
#include <iup_class_cbs.hpp>

#include <im_format_jp2.h>
#include <im_format_avi.h>
#include <im_format_wmv.h>

#include <stdio.h>
#include <string.h>


static void PrintError(int error)
{
  switch (error)
  {
  case IM_ERR_OPEN:
    Iup::Message("IM", "Error Opening File.");
    break;
  case IM_ERR_MEM:
    Iup::Message("IM", "Insufficient memory.");
    break;
  case IM_ERR_ACCESS:
    Iup::Message("IM", "Error Accessing File.");
    break;
  case IM_ERR_DATA:
    Iup::Message("IM", "Image type not Supported.");
    break;
  case IM_ERR_FORMAT:
    Iup::Message("IM", "Invalid Format.");
    break;
  case IM_ERR_COMPRESS:
    Iup::Message("IM", "Invalid or unsupported compression.");
    break;
  default:
    Iup::Message("IM", "Unknown Error.");
  }
}

class imView
{
  cd::Canvas* canvas_draw;
  im::Image* image;
  Iup::Canvas canvas;

  int disable_repaint;  /* used to optimize repaint, while opening a new file */

public:
  imView()
    :canvas_draw(NULL), image(NULL), disable_repaint(0), canvas()
  {
    Iup::Dialog iup_dialog(canvas);
    iup_dialog.SetAttribute("SIZE", "HALFxHALF");  /* initial size */

    // 1) Register "this" object as a callback receiver (need only once)
    IUP_PLUS_INITCALLBACK(iup_dialog, imView);

    // 2) Associate the callback with the button
    IUP_PLUS_SETCALLBACK(canvas, "BUTTON_CB", CanvasButton);
    IUP_PLUS_SETCALLBACK(canvas, "ACTION", CanvasRepaint);
    IUP_PLUS_SETCALLBACK(canvas, "MAP_CB", CanvasMap);

    iup_dialog.Show();

    iup_dialog.SetAttribute("SIZE", NULL); /* remove initial size limitation */
  };

  ~imView()
  {
    if (canvas_draw) delete canvas_draw;
    if (image) delete image;
  }

  void ShowImage(const char* file_name);

protected:
  // 3) Declare the callback as a member function
  IUP_CLASS_DECLARECALLBACK_IFn(imView, CanvasRepaint);
  IUP_CLASS_DECLARECALLBACK_IFnii(imView, CanvasButton);
  IUP_CLASS_DECLARECALLBACK_IFn(imView, CanvasMap);
  IUP_CLASS_DECLARECALLBACK_IFn(imView, DialogClose);
};

int imView::CanvasRepaint(Ihandle*)
{
  if (!canvas_draw || disable_repaint)
    return IUP_DEFAULT;

  canvas_draw->Activate();
  canvas_draw->Clear();

  if (!image)
    return IUP_DEFAULT;

  canvas_draw->PutImage(*image, 0, 0, image->Width(), image->Height());
  
  canvas_draw->Flush();
  
  return IUP_DEFAULT;
}

void imView::ShowImage(const char* file_name)
{
  int error = 0;
  im::Image* new_image = new im::Image(file_name, 0, error, true);
  if (error) PrintError(error);
  if (!new_image) return;

  if (image) delete image;
  image = new_image;

  Iup::Dialog dialog(canvas.GetParent());
  dialog.SetString("TITLE", file_name);

  canvas.Update();
}

int imView::CanvasButton(Ihandle*, int but, int pressed)
{
  char file_name[200] = "*.*";

  if (but != IUP_BUTTON1 || !pressed)
    return IUP_DEFAULT;
  
  disable_repaint = 1;
  if (Iup::GetFile(file_name) != 0)
  {
    disable_repaint = 0;
    return IUP_DEFAULT;
  }

  disable_repaint = 0;
  ShowImage(file_name);
  
  return IUP_DEFAULT;
}

int imView::CanvasMap(Ihandle*)
{
  canvas_draw = new cd::CanvasIup(canvas);
  return IUP_DEFAULT;
}

int main(int argc, char* argv[])
{
  //  imFormatRegisterJP2();
  //  imFormatRegisterAVI();
  //  imFormatRegisterWMV();   

  Iup::Open(argc, argv);

  imView view;
 
  /* Try to get a file name from the command line. */
  if (argc > 1)
    view.ShowImage(argv[1]);
  else   
  {
    char file_name[1024] = "*.*";
    if (IupGetFile(file_name) == 0)
      view.ShowImage(file_name);
  }
                                   
  Iup::MainLoop();
  
  Iup::Close();

  return 0;
}
