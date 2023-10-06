/** \file
 * \brief KRN - IM Kernel File Format
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
#include <math.h>


static const char* iKRNCompTable[1] = 
{
  "NONE"
};

class imFileFormatKRN: public imFileFormatBase
{
  imBinFile* handle;          /* the binary file handle */

public:
  imFileFormatKRN(const imFormat* _iformat) : imFileFormatBase(_iformat), handle(0) {}
  ~imFileFormatKRN() {}

  int Open(const char* file_name);
  int New(const char* file_name);
  void Close();
  void* Handle(int index);
  int ReadImageInfo(int index);
  int ReadImageData(void* data);
  int WriteImageInfo();
  int WriteImageData(void* data);
};

class imFormatKRN: public imFormat
{
public:
  imFormatKRN()
    :imFormat("KRN", 
              "IM Kernel File Format", 
              "*.krn;", 
              iKRNCompTable, 
              1, 
              0)
    {}
  ~imFormatKRN() {}

  imFileFormatBase* Create(void) const { return new imFileFormatKRN(this); }
  int CanWrite(const char* compression, int color_mode, int data_type) const;
};

void imFormatRegisterKRN(void)
{
  imFormatRegister(new imFormatKRN());
}

int imFileFormatKRN::Open(const char* file_name)
{
  char sig[9];

  /* opens the binary file for reading */
  handle = imBinFileOpen(file_name);
  if (!handle)
    return IM_ERR_OPEN;

  /* reads the KRN format identifier */
  imBinFileRead(handle, sig, 8, 1);
  if (imBinFileError(handle))
  {
    imBinFileClose(handle);
    return IM_ERR_ACCESS;
  }

  sig[8] = 0;

  if (!imStrEqual(sig, "IMKERNEL"))
  {
    imBinFileClose(handle);
    return IM_ERR_FORMAT;
  }

  imBinFileSkipLine(handle);

  this->image_count = 1;
  strcpy(this->compression, "NONE");

  return IM_ERR_NONE;
}

int imFileFormatKRN::New(const char* file_name)
{
  /* opens the binary file for writing */
  handle = imBinFileNew(file_name);
  if (!handle)
    return IM_ERR_OPEN;

  this->image_count = 1;  
  if (!imBinFileWrite(handle, (void*)"IMKERNEL\n", 9, 1))
  {
    imBinFileClose(handle);
    return IM_ERR_ACCESS;
  }

  return IM_ERR_NONE;
}

void imFileFormatKRN::Close()
{
  imBinFileClose(handle);
}

void* imFileFormatKRN::Handle(int index)
{
  if (index == 0)
    return (void*)this->handle;
  else
    return NULL;
}

int imFileFormatKRN::ReadImageInfo(int index)
{
  (void)index;
  this->file_color_mode = IM_GRAY|IM_TOPDOWN;

  char desc[512];
  int desc_size = 512;
  if (!imBinFileReadLine(handle, desc, &desc_size))
    return IM_ERR_ACCESS;

  imAttribTable* attrib_table = AttribTable();
  if (desc_size)
    attrib_table->Set("Description", IM_BYTE, desc_size, desc);

  if (!imBinFileReadInteger(handle, &this->width))
    return IM_ERR_ACCESS;

  if (!imBinFileReadInteger(handle, &this->height))
    return IM_ERR_ACCESS;

  int type;
  if (!imBinFileReadInteger(handle, &type))
    return IM_ERR_ACCESS;

  if (type == 0)
    this->file_data_type = IM_INT;
  else
    this->file_data_type = IM_FLOAT;

  return IM_ERR_NONE;
}

int imFileFormatKRN::WriteImageInfo()
{
  this->file_data_type = this->user_data_type;
  this->file_color_mode = IM_GRAY|IM_TOPDOWN;

  imAttribTable* attrib_table = AttribTable();

  int attrib_size;
  const void* attrib_data = attrib_table->Get("Description", NULL, &attrib_size);
  if (attrib_data)
  {
    char* desc = (char*)attrib_data;
    int size = 0;
    while(size < (attrib_size-1) && (desc[size] != '\r' && desc[size] != '\n'))
      size++;

    imBinFileWrite(handle, desc, size, 1);
  }
  imBinFileWrite(handle, (void*)"\n", 1, 1);

  imBinFilePrintf(handle, "%d\n", this->width);
  imBinFilePrintf(handle, "%d\n", this->height);

  if (this->file_data_type == IM_INT)
    imBinFileWrite(handle, (void*)"0\n", 1, 1);
  else
    imBinFileWrite(handle, (void*)"1\n", 1, 1);
  
  /* tests if everything was ok */
  if (imBinFileError(handle))
    return IM_ERR_ACCESS;

  return IM_ERR_NONE;
}

int imFileFormatKRN::ReadImageData(void* data)
{
  imCounterTotal(this->counter, this->height, "Reading KRN...");

  for (int lin = 0; lin < this->height; lin++)
  {
    for (int col = 0; col < this->width; col++)
    {
      if (this->file_data_type == IM_INT)
      {
        int value;
        if (!imBinFileReadInteger(handle, &value))
          return IM_ERR_ACCESS;

        ((int*)this->line_buffer)[col] = value;
      }
      else
      {
        double value;
        if (!imBinFileReadReal(handle, &value))
          return IM_ERR_ACCESS;

        ((float*)this->line_buffer)[col] = (float)value;
      }
    }

    imFileLineBufferRead(this, data, lin, 0);

    if (!imCounterInc(this->counter))
      return IM_ERR_COUNTER;
  }

  return IM_ERR_NONE;
}

int imFileFormatKRN::WriteImageData(void* data)
{
  imCounterTotal(this->counter, this->height, "Writing KRN...");

  for (int lin = 0; lin < this->height; lin++)
  {
    imFileLineBufferWrite(this, data, lin, 0);

    for (int col = 0; col < this->width; col++)
    {
      if (this->file_data_type == IM_INT)
      {
        int value = ((int*)this->line_buffer)[col];

        if (!imBinFilePrintf(handle, "%d ", value))
          return IM_ERR_ACCESS;
      }
      else
      {
        float value = ((float*)this->line_buffer)[col];

        if (!imBinFilePrintf(handle, "%.9f ", value))
          return IM_ERR_ACCESS;
      }
    }

    imBinFileWrite(handle, (void*)"\n", 1, 1);

    if (imBinFileError(handle))
      return IM_ERR_ACCESS;     

    if (!imCounterInc(this->counter))
      return IM_ERR_COUNTER;
  }

  return IM_ERR_NONE;
}

int imFormatKRN::CanWrite(const char* compression, int color_mode, int data_type) const
{
  int color_space = imColorModeSpace(color_mode);

  if (color_space != IM_GRAY)
    return IM_ERR_DATA;                       
                                              
  if (data_type != IM_INT && data_type != IM_FLOAT)
    return IM_ERR_DATA;

  if (!compression || compression[0] == 0)
    return IM_ERR_NONE;

  if (!imStrEqual(compression, "NONE"))
    return IM_ERR_COMPRESS;

  return IM_ERR_NONE;
}
