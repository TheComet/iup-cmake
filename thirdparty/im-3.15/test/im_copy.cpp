/* IM 3 sample that copies an image from one file to another. 
   It is good to test the file formats read and write.
   If the target does not supports the input image it aborts and returns an error.

  Needs "im.lib".

  Usage: im_copy <input_file_name> <output_file_name> [<output_format> [<output_compression]]

    Example: im_copy test.tif test_proc.tif
*/

#include <im_plus.h>

#ifdef WIN32
#include <im_format_avi.h>
#include <im_format_wmv.h>
#endif

#include <stdio.h>
#include <stdlib.h>


void PrintError(int error)
{
  switch (error)
  {
  case IM_ERR_OPEN:
    printf("Error Opening File.\n");
    break;
  case IM_ERR_MEM:
    printf("Insufficient memory.\n");
    break;
  case IM_ERR_ACCESS:
    printf("Error Accessing File.\n");
    break;
  case IM_ERR_DATA:
    printf("Image type not Supported.\n");
    break;
  case IM_ERR_FORMAT:
    printf("Invalid Format.\n");
    break;
  case IM_ERR_COMPRESS:
    printf("Invalid or unsupported compression.\n");
    break;
  default:
    printf("Unknown Error.\n");
  }
}

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    printf("Invalid number of arguments.\n");
    return 1;
  }

#ifdef WIN32
//  imFormatRegisterAVI();
//  imFormatRegisterWMV();
#endif

  int error;
  im::File ifile(argv[1], error);
  if (ifile.Failed())
  {
    PrintError(error);
    return 1;
  }

  char format[10];
  char compression[20];
  int image_count;
  ifile.GetInfo(format, compression, image_count);

  im::File ofile(argv[2], (argc < 3) ? format : argv[3], error);
  if (ofile.Failed())
  {
    PrintError(error);
    return 1;
  }

  if (argc < 4)
    ofile.SetInfo(compression);
  else
    ofile.SetInfo(argv[4]);

  for (int i = 0; i < image_count; i++)
  {
    im::Image image = ifile.LoadImage(i, error);
    if (error != IM_ERR_NONE)
    {
      PrintError(error);
      return 1;
    }

    error = ofile.SaveImage(image);
    if (error != IM_ERR_NONE)
    {
      PrintError(error);
      return 1;
    }

    printf(".");
  }
  printf("done");

  return 0;
}
