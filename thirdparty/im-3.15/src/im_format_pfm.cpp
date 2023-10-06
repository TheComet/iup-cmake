/** \file
 * \brief PFM - Portable FloatMap Image Format
 *
 * See Copyright Notice in im_lib.h
 */

#include "im_format.h"
#include "im_format_all.h"
#include "im_util.h"
#include "im_counter.h"

#include "im_binfile.h"

#include <stdlib.h>
#include <string.h>
#include <memory.h>


static const char* iPFMCompTable[1] = 
{
  "NONE",
};

class imFileFormatPFM: public imFileFormatBase
{
  imBinFile* handle;          /* the binary file handle */
  unsigned char image_type;

public:
  imFileFormatPFM(const imFormat* _iformat): imFileFormatBase(_iformat) {}
  ~imFileFormatPFM() {}

  int Open(const char* file_name);
  int New(const char* file_name);
  void Close();
  void* Handle(int index);
  int ReadImageInfo(int index);
  int ReadImageData(void* data);
  int WriteImageInfo();
  int WriteImageData(void* data);
};

class imFormatPFM: public imFormat
{
public:
  imFormatPFM()
    :imFormat("PFM", 
              "Portable FloatMap Image Format", 
              "*.pfm;", 
              iPFMCompTable, 
              1, 
              0)
    {}
  ~imFormatPFM() {}

  imFileFormatBase* Create(void) const { return new imFileFormatPFM(this); }
  int CanWrite(const char* compression, int color_mode, int data_type) const;
};


void imFormatRegisterPFM(void)
{
  imFormatRegister(new imFormatPFM());
}

int imFileFormatPFM::Open(const char* file_name)
{
  unsigned char sig[2];

  /* opens the binary file for reading */
  handle = imBinFileOpen(file_name);
  if (!handle)
    return IM_ERR_OPEN;

  /* reads the PFM format identifier */
  imBinFileRead(handle, sig, 2, 1);
  if (imBinFileError(handle))
  {
    imBinFileClose(handle);
    return IM_ERR_ACCESS;
  }

  if (sig[0] != 'P' || (sig[1] != 'f' && sig[1] != 'F'))
  {
    imBinFileClose(handle);
    return IM_ERR_FORMAT;
  }

  this->image_type = sig[1];  // 'F' means color, 'f' means grayscale
  this->image_count = 1;

  strcpy(this->compression, "NONE");

  return IM_ERR_NONE;
}

int imFileFormatPFM::New(const char* file_name)
{
  /* opens the binary file for writing */
  handle = imBinFileNew(file_name);
  if (!handle)
    return IM_ERR_OPEN;

  this->image_count = 1;  

  return IM_ERR_NONE;
}

void imFileFormatPFM::Close()
{
  imBinFileClose(handle);
}

void* imFileFormatPFM::Handle(int index)
{
  if (index == 0)
    return (void*)this->handle;
  else
    return NULL;
}

int imFileFormatPFM::ReadImageInfo(int index)
{
  (void)index;

  if (this->image_type == 'f')
    this->file_color_mode = IM_GRAY;
  else
    this->file_color_mode = IM_RGB | IM_PACKED;

  if (!imBinFileReadInteger(handle, &this->width))
    return IM_ERR_ACCESS;

  if (!imBinFileReadInteger(handle, &this->height))
    return IM_ERR_ACCESS;

  if (this->height <= 0 || this->width <= 0)
    return IM_ERR_DATA;

  double byte_order = 1.0;
  if (!imBinFileReadReal(handle, &byte_order))
    return IM_ERR_ACCESS;

  if (byte_order > 0)
    imBinFileByteOrder(handle, IM_BIGENDIAN);
  else
    imBinFileByteOrder(handle, IM_LITTLEENDIAN);

  this->file_data_type = IM_FLOAT;

  return IM_ERR_NONE;
}

int imFileFormatPFM::WriteImageInfo()
{
  this->file_data_type = this->user_data_type;
  this->file_color_mode = imColorModeSpace(this->user_color_mode);

  if (this->file_color_mode == IM_GRAY)
    this->image_type = 'f';
  else
  {
    this->image_type = 'F';
    this->file_color_mode |= IM_PACKED;
  }

  imBinFilePrintf(handle, "P%c\n", (int)this->image_type);

  if (imBinFileError(handle))
    return IM_ERR_ACCESS;

  imBinFilePrintf(handle, "%d ", this->width);
  imBinFilePrintf(handle, "%d\n", this->height);

  if (imBinCPUByteOrder() == IM_BIGENDIAN)
    imBinFilePrintf(handle, "1.0\n");
  else
    imBinFilePrintf(handle, "-1.0\n");

  /* tests if everything was ok */
  if (imBinFileError(handle))
    return IM_ERR_ACCESS;

  return IM_ERR_NONE;
}

int imFileFormatPFM::ReadImageData(void* data)
{
  imCounterTotal(this->counter, this->height, "Reading PFM...");

  int line_raw_size = imImageLineSize(this->width, this->file_color_mode, this->file_data_type);

  for (int lin = 0; lin < this->height; lin++)
  {
    imBinFileRead(handle, this->line_buffer, line_raw_size, 1);

    if (imBinFileError(handle))
      return IM_ERR_ACCESS;     

    imFileLineBufferRead(this, data, lin, 0);

    if (!imCounterInc(this->counter))
      return IM_ERR_COUNTER;
  }

  return IM_ERR_NONE;
}

int imFileFormatPFM::WriteImageData(void* data)
{
  imCounterTotal(this->counter, this->height, "Writing PFM...");

  int line_raw_size = imImageLineSize(this->width, this->file_color_mode, this->file_data_type);

  for (int lin = 0; lin < this->height; lin++)
  {
    imFileLineBufferWrite(this, data, lin, 0);

    imBinFileWrite(handle, this->line_buffer, line_raw_size, 1);

    if (imBinFileError(handle))
      return IM_ERR_ACCESS;     

    if (!imCounterInc(this->counter))
      return IM_ERR_COUNTER;
  }

  return IM_ERR_NONE;
}

int imFormatPFM::CanWrite(const char* compression, int color_mode, int data_type) const
{
  int color_space = imColorModeSpace(color_mode);

  if (color_space != IM_GRAY && color_space != IM_RGB)
    return IM_ERR_DATA;                       
                                              
  if (data_type != IM_FLOAT)
    return IM_ERR_DATA;

  if (!compression || compression[0] == 0)
    return IM_ERR_NONE;

  if (!imStrEqual(compression, "NONE"))
    return IM_ERR_COMPRESS;

  return IM_ERR_NONE;
}
