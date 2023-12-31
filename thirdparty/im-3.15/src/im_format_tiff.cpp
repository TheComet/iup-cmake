/** \file
 * \brief TIFF - Tagged Image File Format
 *
 * See Copyright Notice in im_lib.h
 * See libTIFF Copyright Notice in tiff.h
 */

#include "im_format.h"
#include "im_util.h"
#include "im_format_all.h"
#include "im_counter.h"
#include "im_binfile.h"

#include "tiffiop.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>


static tmsize_t iTIFFReadProc(thandle_t fd, void* buf, tmsize_t size)
{
  imBinFile* file_bin = (imBinFile*)fd;
  return imBinFileRead(file_bin, buf, (unsigned long)size, 1);
}

static tmsize_t iTIFFWriteProc(thandle_t fd, void* buf, tmsize_t size)
{
  imBinFile* file_bin = (imBinFile*)fd;
  return imBinFileWrite(file_bin, buf, (unsigned long)size, 1);
}

static toff_t iTIFFSeekProc(thandle_t fd, toff_t off, int whence)
{
  imBinFile* file_bin = (imBinFile*)fd;
  switch (whence)
  {
  case SEEK_SET:
    imBinFileSeekTo(file_bin, (unsigned long)off);
    break;
  case SEEK_CUR:
    imBinFileSeekOffset(file_bin, (unsigned long)off);
    break;
  case SEEK_END:
    imBinFileSeekFrom(file_bin, (unsigned long)off);
    break;
  }

  return imBinFileTell(file_bin);
}

static int iTIFFCloseProc(thandle_t fd)
{
  imBinFile* file_bin = (imBinFile*)fd;
  imBinFileClose(file_bin);
  return 0;
}

static toff_t iTIFFSizeProc(thandle_t fd)
{
  imBinFile* file_bin = (imBinFile*)fd;
  return imBinFileSize(file_bin);
}

static int iTIFFMapProc(thandle_t fd, void** pbase, toff_t* psize)
{
  (void)fd; (void)pbase; (void)psize;
  return (0);
}

static void iTIFFUnmapProc(thandle_t fd, void* base, toff_t size)
{
  (void)fd; (void)base; (void)size;
}

static TIFF* iTIFFOpen(const char* name, const char* mode)
{
  imBinFile* bin_file;
  TIFF* tiff;

  if (mode[0] == 'r')
    bin_file = imBinFileOpen(name);
  else
    bin_file = imBinFileNew(name);

  if (!bin_file)
    return NULL;

  tiff = TIFFClientOpen(name, mode, (thandle_t)bin_file, iTIFFReadProc, iTIFFWriteProc,
                        iTIFFSeekProc, iTIFFCloseProc,
                        iTIFFSizeProc, iTIFFMapProc,
                        iTIFFUnmapProc);
  if (!tiff)
    imBinFileClose(bin_file);

  return tiff;
}

//Used to debug TIFF loading and decoding
//#define IM_TIFF_DEBUG_RGBA 1

#define TIFFTAG_GEOPIXELSCALE        33550
#define TIFFTAG_INTERGRAPH_MATRIX    33920
#define TIFFTAG_GEOTIEPOINTS         33922
#define TIFFTAG_GEOTRANSMATRIX       34264
#define TIFFTAG_GEOKEYDIRECTORY      34735
#define TIFFTAG_GEODOUBLEPARAMS      34736
#define TIFFTAG_GEOASCIIPARAMS       34737

#define TIFFTAG_CFAREPEATPATTERNDIM	33421	/* dimensions of CFA pattern */
#define TIFFTAG_CFAPATTERN		33422	/* color filter array pattern */
#define	    PHOTOMETRIC_CFA		32803	/* color filter array */
#define	    PHOTOMETRIC_LINEARRAW		34892

static const TIFFFieldInfo iTiffFieldInfo[] = 
{
  /* Patch from Dave Coffin (Used for DNG) */
  { TIFFTAG_WHITELEVEL,	-2, -1,	TIFF_LONG,	FIELD_CUSTOM,  0,	1,	"WhiteLevel" },
  { TIFFTAG_WHITELEVEL,	-2, -1,	TIFF_SHORT,	FIELD_CUSTOM,  0,	1,	"WhiteLevel" },
  { TIFFTAG_CFAREPEATPATTERNDIM, 2, 2, TIFF_SHORT,	FIELD_CUSTOM, 0,	0,	"CFARepeatPatternDim" },
  { TIFFTAG_CFAPATTERN,	-1, -1,	TIFF_BYTE,	FIELD_CUSTOM, 0,	1,	"CFAPattern" },

  /* GeoTIFF Tags */
  { TIFFTAG_GEOPIXELSCALE,	-1,-1, TIFF_DOUBLE,	FIELD_CUSTOM, TRUE,	TRUE,	"GeoPixelScale" },
  { TIFFTAG_INTERGRAPH_MATRIX,-1,-1, TIFF_DOUBLE,	FIELD_CUSTOM, TRUE,	TRUE,	"Intergraph TransformationMatrix" },
  { TIFFTAG_GEOTIEPOINTS,	-1,-1, TIFF_DOUBLE,	FIELD_CUSTOM, TRUE,	TRUE,	"GeoTiePoints" },
  { TIFFTAG_GEOTRANSMATRIX,	-1,-1, TIFF_DOUBLE,	FIELD_CUSTOM, TRUE,	TRUE,	"GeoTransformationMatrix" },
  { TIFFTAG_GEOKEYDIRECTORY,-1,-1, TIFF_SHORT,	FIELD_CUSTOM, TRUE,	TRUE,	"GeoKeyDirectory" },
  { TIFFTAG_GEODOUBLEPARAMS,	-1,-1, TIFF_DOUBLE,	FIELD_CUSTOM, TRUE,	TRUE,	"GeoDoubleParams" },
  { TIFFTAG_GEOASCIIPARAMS,	-1,-1, TIFF_ASCII,	FIELD_CUSTOM, TRUE,	FALSE, "GeoASCIIParams" }
};

#define IMTIFF_NUMCOMP 15

/* this list must be sorted because of bsearch */
static uint16 iTIFFCompIdTable [IMTIFF_NUMCOMP] = 
{
  COMPRESSION_NONE,                       
  COMPRESSION_CCITTRLE,                   
  COMPRESSION_CCITTFAX3,                  
  COMPRESSION_CCITTFAX4,                  
  COMPRESSION_LZW,                        
  COMPRESSION_JPEG,                       
  COMPRESSION_ADOBE_DEFLATE,              
  COMPRESSION_NEXT,                       
  COMPRESSION_CCITTRLEW,                  
  COMPRESSION_PACKBITS,                   
  COMPRESSION_THUNDERSCAN,                
  COMPRESSION_PIXARLOG,                   
  COMPRESSION_DEFLATE,                    
  COMPRESSION_SGILOG,
  COMPRESSION_SGILOG24
};

static int iTIFFCompareCompID(const void *elem1, const void *elem2)
{
  const uint16 *tiff_comp_elem1 = (const uint16 *)elem1;
  const uint16 *tiff_comp_elem2 = (const uint16 *)elem2;

  if (*tiff_comp_elem1 > *tiff_comp_elem2)
    return 1;

  if (*tiff_comp_elem1 < *tiff_comp_elem2)
    return -1;

  return 0;
}

static int iTIFFGetCompIndex(uint16 Compression)
{
  if (Compression == COMPRESSION_OJPEG)
    Compression = COMPRESSION_JPEG;

  uint16* comp_result = (uint16 *)bsearch(&Compression, iTIFFCompIdTable, sizeof(iTIFFCompIdTable)/sizeof(uint16), sizeof(uint16), iTIFFCompareCompID);

  if (comp_result == NULL)
  {
    return -1;
  }

  return (int)(comp_result - iTIFFCompIdTable);
}

/* this list must follow iTIFFCompIdTable order */
static const char* iTIFFCompTable[IMTIFF_NUMCOMP] = 
{
  "NONE",
  "CCITTRLE",
  "CCITTFAX3",
  "CCITTFAX4",
  "LZW",
  "JPEG",
  "ADOBEDEFLATE",
  "NEXT",
  "CCITTRLEW",
  "RLE",
  "THUNDERSCAN",
  "PIXARLOG",    
  "DEFLATE",
  "SGILOG",
  "SGILOG24"
};

static uint16 iTIFFCompFind(const char* compression)
{
  for(int i = 0; i < IMTIFF_NUMCOMP; i++)
  {
    if (imStrEqual(compression, iTIFFCompTable[i]))
      return iTIFFCompIdTable[i];
  }

  return (uint16)-1;
}

static uint16 iTIFFCompDefault(int color_space, int data_type)
{
  if (color_space == IM_BINARY)
    return COMPRESSION_CCITTRLE;

  if (color_space == IM_MAP)
    return COMPRESSION_PACKBITS;

  if (color_space == IM_YCBCR && data_type == IM_BYTE)
    return COMPRESSION_JPEG;

  if (color_space == IM_XYZ)
    return COMPRESSION_SGILOG;

  if (data_type >= IM_FLOAT)
    return COMPRESSION_NONE;

  return COMPRESSION_LZW;
}

static uint16 iTIFFCompCalc(const char* compression, int color_mode, int data_type)
{
  uint16 Compression;
  if (compression[0] == 0)
    Compression = iTIFFCompDefault(imColorModeSpace(color_mode), data_type);
  else
    Compression = iTIFFCompFind(compression);

  return Compression;
}

static int iTIFFGetDataType(TIFFDataType field_type)
{
  switch(field_type)
  {
  case TIFF_UNDEFINED:
  case TIFF_ASCII: 
  case TIFF_BYTE:
  case TIFF_SBYTE:
    return IM_BYTE;
  case TIFF_SSHORT:
    return IM_SHORT;
  case TIFF_SHORT:
    return IM_USHORT;
  case TIFF_LONG:
  case TIFF_SLONG:
    return IM_INT;
  case TIFF_RATIONAL:
  case TIFF_SRATIONAL:
  case TIFF_FLOAT:
    return IM_FLOAT;
  case TIFF_DOUBLE:
    return IM_DOUBLE;
  default:  /* TIFF_NOTYPE, TIFF_IFD, TIFF_LONG8, TIFF_SLONG8, TIFF_IFD8 */
    return -1;
  }
}

static int iTIFFWriteTag(TIFF* tiff, int index, const char* name, int data_type, int count, const void* data)
{
  const TIFFField *fld = TIFFFieldWithName(tiff, name);
  if (fld)
  {
    /* ignored tags */
    if (fld->field_tag == TIFFTAG_EXIFIFD ||         /* offset */
        fld->field_tag == TIFFTAG_GPSIFD ||          
        fld->field_tag == TIFFTAG_INTEROPERABILITYIFD ||   
	      fld->field_tag == TIFFTAG_SUBIFD ||          
	      fld->field_tag == TIFFTAG_COLORMAP ||        /* handled elsewhere */
	      fld->field_tag == TIFFTAG_EXTRASAMPLES ||
	      fld->field_tag == TIFFTAG_TRANSFERFUNCTION ||
	      fld->field_tag == TIFFTAG_RESOLUTIONUNIT ||
	      fld->field_tag == TIFFTAG_XRESOLUTION ||
	      fld->field_tag == TIFFTAG_YRESOLUTION ||
        fld->field_tag == TIFFTAG_INKNAMES)
      return 1;

    /* test if tag type is the same as the attribute type */
    if (iTIFFGetDataType(fld->field_type) != data_type)
      return 1;

    if (fld->field_passcount)
    {
			if (fld->field_writecount == TIFF_VARIABLE2)
      {
        uint32 value_count = (uint32)count;
        if (TIFFSetField(tiff, fld->field_tag, value_count, data) != 1)
          return 1;
      }
      else
      {
        uint16 value_count = (uint16)count;
        if (TIFFSetField(tiff, fld->field_tag, value_count, data) != 1)
          return 1;
      }
    } 
    else
    {
      if (fld->field_tag == TIFFTAG_PAGENUMBER ||
			    fld->field_tag == TIFFTAG_HALFTONEHINTS ||
			    fld->field_tag == TIFFTAG_YCBCRSUBSAMPLING ||
          fld->field_tag == TIFFTAG_DOTRANGE)
      {
        // there are 2 separated ushort values
        uint16* ushort_value = (uint16*)data;
        TIFFSetField(tiff, fld->field_tag, ushort_value[0], ushort_value[1]);
        return 1;
      }

      if (count > 1 || fld->field_type == TIFF_ASCII)
        TIFFSetField(tiff, fld->field_tag, data);
      else
      {
        switch(data_type)
        {
        case IM_BYTE:
          {
            imbyte* byte_data = (imbyte*)data;
            TIFFSetField(tiff, fld->field_tag, *byte_data);
          }
          break;
        case IM_SHORT:
          {
            short* short_data = (short*)data;
            TIFFSetField(tiff, fld->field_tag, *short_data);
          }
          break;
        case IM_USHORT:
          {
            imushort* short_data = (imushort*)data;
            TIFFSetField(tiff, fld->field_tag, *short_data);
          }
          break;
        case IM_INT:
          {
            int* long_data = (int*)data;
            TIFFSetField(tiff, fld->field_tag, *long_data);
          }
          break;
        case IM_FLOAT:
          {
            float* float_data = (float*)data;
            TIFFSetField(tiff, fld->field_tag, *float_data);
          }
          break;
        case IM_DOUBLE:
          {
            double* double_data = (double*)data;
            TIFFSetField(tiff, fld->field_tag, *double_data);
          }
          break;
        default:
          break;
        }
      }
    } 
  }

  (void)index;
  return 1;
}

static void iTIFFWriteCustomTags(TIFF* tiff, imAttribTable* attrib_table)
{
  attrib_table->ForEach(tiff, (imAttribTableCallback)iTIFFWriteTag);
}

static void iTIFFReadCustomTags(TIFF* tiff, imAttribTable* attrib_table)
{
  int  i;
  short tag_count;

  tag_count = (short) TIFFGetTagListCount(tiff);
  for( i = 0; i < tag_count; i++ )
  {
    ttag_t tag = TIFFGetTagListEntry(tiff, i);
    const TIFFField *fld;

    fld = TIFFFieldWithTag(tiff, tag);
    if (fld == NULL)
      continue;

    /* offsets */
    if (fld->field_tag == TIFFTAG_EXIFIFD ||         
        fld->field_tag == TIFFTAG_GPSIFD ||          
        fld->field_tag == TIFFTAG_INTEROPERABILITYIFD ||   
	      fld->field_tag == TIFFTAG_SUBIFD ||          
	      fld->field_tag == TIFFTAG_COLORMAP)
      continue;
        
    /* handled elsewhere */
	  if (fld->field_tag == TIFFTAG_EXTRASAMPLES ||
	      fld->field_tag == TIFFTAG_TRANSFERFUNCTION ||
	      fld->field_tag == TIFFTAG_RESOLUTIONUNIT ||
	      fld->field_tag == TIFFTAG_XRESOLUTION ||
	      fld->field_tag == TIFFTAG_YRESOLUTION ||
        fld->field_tag == TIFFTAG_INKNAMES)
      continue;

      if (fld->field_tag == TIFFTAG_BLACKLEVEL ||
          fld->field_tag == TIFFTAG_DEFAULTCROPSIZE ||
          fld->field_tag == TIFFTAG_DEFAULTCROPORIGIN)
      {
        //TODO re-check this
        /* libTIFF issue. When reading custom tags there is an incorrect interpretation of the tag
        that leads to return always type=RATIONAL for these tags. */
        continue;
      }

    int data_type = iTIFFGetDataType(fld->field_type);
    if (data_type == -1)
          continue;

    int data_count = -1;
    void* data = NULL;

    if (fld->field_passcount)
    {
			if (fld->field_readcount == TIFF_VARIABLE2)
      {
        uint32 value_count;
        if (TIFFGetField(tiff, tag, &value_count, &data) != 1)
          continue;
        data_count = value_count;
      }
      else
      {
        uint16 value_count;
        if (TIFFGetField(tiff, tag, &value_count, &data) != 1)
          continue;
        data_count = value_count;
      }

      if (data && data_count > 0)
        attrib_table->Set(fld->field_name, data_type, data_count, data);
    } 
    else
    {
      data_count = fld->field_readcount;

      if (fld->field_tag == TIFFTAG_PAGENUMBER ||
			    fld->field_tag == TIFFTAG_HALFTONEHINTS ||
			    fld->field_tag == TIFFTAG_YCBCRSUBSAMPLING ||
          fld->field_tag == TIFFTAG_DOTRANGE)
      {
        // there are 2 separated ushort values
        uint16 ushort_value[2];
        if (TIFFGetField(tiff, fld->field_tag, &ushort_value[0], &ushort_value[1]))
          attrib_table->Set(fld->field_name, IM_USHORT, 2, ushort_value);
        continue;
      }

		  if (fld->field_type == TIFF_ASCII ||
		      fld->field_readcount == TIFF_VARIABLE ||
		      fld->field_readcount == TIFF_VARIABLE2 ||
		      fld->field_readcount == TIFF_SPP ||
		      data_count > 1) 
      {
        if (TIFFGetField(tiff, tag, &data) != 1)
          continue;

        if (data)
        {
          if (fld->field_type == TIFF_ASCII && data_count == -1)
            data_count = (int)strlen((char*)data)+1;

          if (data_count > 0)
          {
            char* newstr = NULL;

            if (fld->field_type == TIFF_ASCII && ((char*)data)[data_count-1] != 0)
            {
              int j = data_count-1;
              char* p = (char*)data;
              while (j > 0 && p[j] != 0)
                j--;
              if (j == 0)
              {
                if (fld->field_tag == TIFFTAG_DATETIME ||
			              fld->field_tag == EXIFTAG_DATETIMEORIGINAL ||
                    fld->field_tag == EXIFTAG_DATETIMEDIGITIZED)
                {
                  /* sometimes theses tags get non standard strings,
                      libTIIF does not returns the actual number os bytes read,
                      it returns the standard value of 20.
                      so we will try to find the actual string size, but we risk in a memory invalid access. */
                  j = data_count;
                  while (j < data_count+6 && p[j] != 0)
                    j++;
                  if (j < data_count+6)
                    data_count = j+1;
                }
                else
                {
                  newstr = (char*)malloc(data_count+1);
                  memcpy(newstr, data, data_count);
                  newstr[data_count] = 0;
                  data_count++;
                }
              }
              else
                data_count = j;
            }

            attrib_table->Set(fld->field_name, data_type, data_count, data);

            if (newstr) free(newstr);
          }
        }
      }
      else if (data_count == 1)
      {
        int size = imDataTypeSize(data_type);
        data = malloc(size);
        if (TIFFGetField(tiff, tag, data) == 1)
          attrib_table->Set(fld->field_name, data_type, data_count, data);
        free(data);
        data = NULL;
      }
    } 
  }
}

static void iTIFFReadAttributes(TIFF* tiff, imAttribTable* attrib_table)
{
  uint16 ResolutionUnit = RESUNIT_NONE;
  TIFFGetField(tiff, TIFFTAG_RESOLUTIONUNIT, &ResolutionUnit);
  if (ResolutionUnit != RESUNIT_NONE)
  {
    float xres = 0, yres = 0;

    TIFFGetField(tiff, TIFFTAG_XRESOLUTION, &xres);
    TIFFGetField(tiff, TIFFTAG_YRESOLUTION, &yres);

    if (xres != 0 && yres != 0)
    {
      if (ResolutionUnit == RESUNIT_INCH)
        attrib_table->Set("ResolutionUnit", IM_BYTE, -1, "DPI");
      else
        attrib_table->Set("ResolutionUnit", IM_BYTE, -1, "DPC");

      attrib_table->Set("XResolution", IM_FLOAT, 1, (void*)&xres);
      attrib_table->Set("YResolution", IM_FLOAT, 1, (void*)&yres);
    }
  }

  uint16 *transferfunction[3]; 
  if (TIFFGetField(tiff, TIFFTAG_TRANSFERFUNCTION, &transferfunction[0], &transferfunction[1], &transferfunction[2]))
  {
    uint16 SamplesPerPixel = 1, BitsPerSample = 1, ExtraSamples = 0, *SampleInfo;
    TIFFGetFieldDefaulted(tiff, TIFFTAG_BITSPERSAMPLE, &BitsPerSample);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_EXTRASAMPLES, &ExtraSamples, &SampleInfo);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL, &SamplesPerPixel);

    int num = (SamplesPerPixel - ExtraSamples) > 1 ? 3 : 1;
    int count = 1L<<BitsPerSample;
    if (num == 1)
      attrib_table->Set("TransferFunction0", IM_USHORT, count, transferfunction[0]);
    else
    {
      attrib_table->Set("TransferFunction0", IM_USHORT, count, transferfunction[0]);
      attrib_table->Set("TransferFunction1", IM_USHORT, count, transferfunction[1]);
      attrib_table->Set("TransferFunction2", IM_USHORT, count, transferfunction[2]);
    }
  }

  char *inknames;
  if (TIFFGetField(tiff, TIFFTAG_INKNAMES, &inknames))
  {
    // Ink names are separated by '0', so strlen will measure only the first string
    uint16 numinks;
    TIFFGetField(tiff, TIFFTAG_NUMBEROFINKS, &numinks);
    int inknameslen = 0;
    for (int k = 0; k < (int)numinks; k++)
      inknameslen += (int)strlen(inknames + inknameslen) + 1;
    attrib_table->Set("InkNames", IM_BYTE, inknameslen, inknames);
  }

  iTIFFReadCustomTags(tiff, attrib_table);

  uint64 offset;
  if (TIFFGetField(tiff, TIFFTAG_EXIFIFD, &offset))
  {
    tdir_t cur_dir = TIFFCurrentDirectory(tiff);

    if (!TIFFReadEXIFDirectory(tiff, offset))
    {
      TIFFSetDirectory(tiff, cur_dir);
      return;
    }

    iTIFFReadCustomTags(tiff, attrib_table);
    TIFFSetDirectory(tiff, cur_dir);
  }
}

static void iTIFFWriteAttributes(TIFF* tiff, imAttribTable* attrib_table)
{
  char* res_unit = (char*)attrib_table->Get("ResolutionUnit");
  if (res_unit)
  {
    float* xres = (float*)attrib_table->Get("XResolution");
    float* yres = (float*)attrib_table->Get("YResolution");

    if (xres && yres)
    {
      uint16 tiff_res_unit = RESUNIT_CENTIMETER;
      if (imStrEqual(res_unit, "DPI"))
        tiff_res_unit = RESUNIT_INCH;

      TIFFSetField(tiff, TIFFTAG_RESOLUTIONUNIT, tiff_res_unit);
      TIFFSetField(tiff, TIFFTAG_XRESOLUTION, *xres);
      TIFFSetField(tiff, TIFFTAG_YRESOLUTION, *yres);
    }
  }

  uint16 *transferfunction0 = (uint16*)attrib_table->Get("TransferFunction0"); 
  if (transferfunction0)
  {
    uint16 SamplesPerPixel = 1, ExtraSamples = 0, *SampleInfo;
    TIFFGetFieldDefaulted(tiff, TIFFTAG_EXTRASAMPLES, &ExtraSamples, &SampleInfo);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL, &SamplesPerPixel);

    int num = (SamplesPerPixel - ExtraSamples) > 1 ? 3 : 1;
    if (num == 1)
      TIFFSetField(tiff, TIFFTAG_TRANSFERFUNCTION, transferfunction0);
    else
    {
      uint16 *transferfunction1 = (uint16*)attrib_table->Get("TransferFunction1"); 
      uint16 *transferfunction2 = (uint16*)attrib_table->Get("TransferFunction2"); 

      if (transferfunction1 && transferfunction2)
        TIFFSetField(tiff, TIFFTAG_TRANSFERFUNCTION, transferfunction0, transferfunction1, transferfunction2);
    }
  }

  char* inknames = (char*)attrib_table->Get("InkNames");
  if (inknames)
    TIFFSetField(tiff, TIFFTAG_INKNAMES, inknames);

  int proflength;
  const void* profdata = attrib_table->Get("ICCProfile", (int*)NULL, &proflength);
  if (profdata)
    TIFFSetField(tiff, TIFFTAG_ICCPROFILE, proflength, profdata);

  iTIFFWriteCustomTags(tiff, attrib_table);
}

class imFileFormatTIFF: public imFileFormatBase
{
  TIFF* tiff;

  int invert,      // must invert black and white reference
      lab_fix;     // convert CIE Lab to unsigned

  int cpx_int,     // original data is a complex integer (when reading)
      extra_sample_size, // extra samples fix, eliminate if more than one (when reading)
      sample_size_no_extra, // extra samples fix (when reading)
      h_subsample, v_subsample, // (when reading)
      start_plane; // first band to read in a multiband image (when reading)

  void** tile_buf;
  int tile_buf_count, tile_width, tile_height, 
      tile_start_lin, tile_line_size, tile_line_raw_size;

  int ReadTileline(void* line_buffer, int lin, int plane);
  void InvertBits(void* line_buffer, int size);

public:
  imFileFormatTIFF(const imFormat* _iformat): imFileFormatBase(_iformat) {}
  ~imFileFormatTIFF() {}

  int Open(const char* file_name);
  int New(const char* file_name);
  void Close();
  void* Handle(int index);
  int ReadImageInfo(int index);
  int ReadImageData(void* data);
  int WriteImageInfo();
  int WriteImageData(void* data);
};

class imFormatTIFF: public imFormat
{
public:
  imFormatTIFF()
    :imFormat("TIFF", 
              "Tagged Image File Format", 
              "*.tif;*.tiff;", 
              iTIFFCompTable, 
              IMTIFF_NUMCOMP, 
              1)
    { extra = "LIBTIFF Version 4.0.0"; }
  ~imFormatTIFF() {}

  imFileFormatBase* Create(void) const { return new imFileFormatTIFF(this); }
  int CanWrite(const char* compression, int color_mode, int data_type) const;
};

static void iTIFFDefaultDirectory(TIFF *tiff)
{
  /* Install the IM Tag field info */
  TIFFMergeFieldInfo(tiff, iTiffFieldInfo, TIFFArrayCount(iTiffFieldInfo));
}

void imFormatRegisterTIFF(void)
{
  TIFFSetTagExtender(iTIFFDefaultDirectory);
  imFormatRegister(new imFormatTIFF());
}

int imFileFormatTIFF::Open(const char* file_name)
{
  this->tiff = iTIFFOpen(file_name, "r");
  if (this->tiff == NULL)
    return IM_ERR_FORMAT;

  // Return the compression of the first image in the file.
  uint16 Compression = COMPRESSION_NONE;
  TIFFGetField(this->tiff, TIFFTAG_COMPRESSION, &Compression);
  int comp_index = iTIFFGetCompIndex(Compression);
  if (comp_index == -1) return IM_ERR_COMPRESS;
  strcpy(this->compression, iTIFFCompTable[comp_index]);

  this->image_count = TIFFNumberOfDirectories(this->tiff);
  this->tile_buf = NULL;
  this->start_plane = 0;

  return IM_ERR_NONE;
}

int imFileFormatTIFF::New(const char* file_name)
{
  this->tiff = iTIFFOpen(file_name, "w");
  if (this->tiff == NULL)
    return IM_ERR_OPEN;

  this->tile_buf = NULL;

  return IM_ERR_NONE;
}

void imFileFormatTIFF::Close()
{
  if (this->tile_buf)
  {
    for (int i = 0; i < this->tile_buf_count; i++)
      free(this->tile_buf[i]);
    free(this->tile_buf);
  }

  TIFFClose(this->tiff);
}

void* imFileFormatTIFF::Handle(int index)
{
  if (index == 0)
    return (void*)this->tiff->tif_fd;
  else if (index == 1)
    return (void*)this->tiff;
  else
    return NULL;
}
                
int imFileFormatTIFF::ReadImageInfo(int index)
{
  int sub_ifd = -1;

  // Defaults that can be different for each image
  this->cpx_int = 0;
  this->invert = 0;
  this->lab_fix = 0;
  this->extra_sample_size = 0;
  this->sample_size_no_extra = 0;
  this->h_subsample = 1;
  this->v_subsample = 1;

  if (!TIFFSetDirectory(this->tiff, (tdir_t)index))
    return IM_ERR_ACCESS;

  imAttribTable* attrib_table = AttribTable();

  uint16* attrib_start_plane = (uint16*)attrib_table->Get("MultiBandSelect");
  if (attrib_start_plane)
    this->start_plane = *attrib_start_plane;
  else
    this->start_plane = 0;

  uint16* sub_ifd_atrib = (uint16*)attrib_table->Get("SubIFDSelect");
  if (sub_ifd_atrib) sub_ifd = *sub_ifd_atrib;

  /* must clear the attribute list, because it can have multiple images and 
     has many attributes that may exists only for specific images. */
  attrib_table->RemoveAll();
  imFileSetBaseAttributes(this);

  void* data = NULL;
  if (TIFFGetField(this->tiff, TIFFTAG_DNGVERSION, &data) == 1 && data)
  {
    uint32 SubFileType = 0;
    TIFFGetField(this->tiff, TIFFTAG_SUBFILETYPE, &SubFileType);

    uint16 SubIFDsCount = 0;
    uint64* SubIFDs = NULL;
    TIFFGetField(this->tiff, TIFFTAG_SUBIFD, &SubIFDsCount, &SubIFDs);
    attrib_table->Set("SubIFDCount", IM_USHORT, 1, (void*)&SubIFDsCount);

    /* If it is a DNG file and has SubIFDs, 
       then ignore the thumbnail and position at the image SubIFD. */
    if (SubFileType == FILETYPE_REDUCEDIMAGE && SubIFDsCount != 0)
    {
      if (sub_ifd < 0 || sub_ifd >= SubIFDsCount) sub_ifd = SubIFDsCount-1;
      uint64 SubIFDOffset = SubIFDs[sub_ifd];

      /* Load the main image attributes, the SubIFD contains only a few attributes. */
      iTIFFReadAttributes(this->tiff, attrib_table);

      TIFFSetSubDirectory(this->tiff, SubIFDOffset);
    }
  }

  uint16 Compression = COMPRESSION_NONE;
  TIFFGetField(this->tiff, TIFFTAG_COMPRESSION, &Compression);
  int comp_index = iTIFFGetCompIndex(Compression);
  if (comp_index == -1) return IM_ERR_COMPRESS;
  strcpy(this->compression, iTIFFCompTable[comp_index]);

  if (Compression == COMPRESSION_JPEG)
    TIFFSetField(this->tiff, TIFFTAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB);

  uint32 Width;
  if (!TIFFGetField(this->tiff, TIFFTAG_IMAGEWIDTH, &Width))       
    return IM_ERR_FORMAT;
  this->width = Width;

  uint32 Height;
  if (!TIFFGetField(this->tiff, TIFFTAG_IMAGELENGTH, &Height))
    return IM_ERR_FORMAT;
  this->height = Height;

  uint16 Photometric;
  if (!TIFFGetField(this->tiff, TIFFTAG_PHOTOMETRIC, &Photometric))
    return IM_ERR_FORMAT;
  attrib_table->Set("Photometric", IM_USHORT, 1, (void*)&Photometric);

  switch(Photometric)
  {
  case PHOTOMETRIC_MINISWHITE:
    this->invert = 1;
  case PHOTOMETRIC_LINEARRAW:
  case PHOTOMETRIC_CFA:
  case PHOTOMETRIC_LOGL:
  case PHOTOMETRIC_MASK:
  case PHOTOMETRIC_MINISBLACK:
    this->file_color_mode = IM_GRAY;
    break;
  case PHOTOMETRIC_PALETTE:
    this->file_color_mode = IM_MAP;
    break;
  case PHOTOMETRIC_RGB:
    this->file_color_mode = IM_RGB;
    break;
  case PHOTOMETRIC_SEPARATED:
    this->file_color_mode = IM_CMYK;
    break;
  case PHOTOMETRIC_YCBCR:
    if (Compression == COMPRESSION_JPEG)
      this->file_color_mode = IM_RGB;
    else
      this->file_color_mode = IM_YCBCR;
    break;
  case PHOTOMETRIC_CIELAB:
    this->lab_fix = 1;
  case PHOTOMETRIC_ITULAB:
  case PHOTOMETRIC_ICCLAB:
    this->file_color_mode = IM_LAB;
    break;
  case PHOTOMETRIC_LOGLUV:
    this->file_color_mode = IM_XYZ;
    break;
  default: 
    return IM_ERR_DATA;
  }

  if (Photometric == PHOTOMETRIC_LOGLUV || Photometric == PHOTOMETRIC_LOGL)
    TIFFSetField(this->tiff, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);

  uint16 SamplesPerPixel = 1, BitsPerSample = 1;
  TIFFGetFieldDefaulted(this->tiff, TIFFTAG_BITSPERSAMPLE, &BitsPerSample);
  TIFFGetFieldDefaulted(this->tiff, TIFFTAG_SAMPLESPERPIXEL, &SamplesPerPixel);

  if (BitsPerSample == 1 && this->file_color_mode == IM_GRAY)
    this->file_color_mode = IM_BINARY;

  /* consistency checks */
  if (Photometric == PHOTOMETRIC_PALETTE && (SamplesPerPixel != 1 || BitsPerSample > 8))
    return IM_ERR_DATA;

  if (Photometric == PHOTOMETRIC_MASK && (SamplesPerPixel != 1 || BitsPerSample != 1))
    return IM_ERR_DATA;

  if ((Photometric == PHOTOMETRIC_CFA || Photometric == PHOTOMETRIC_LINEARRAW) && SamplesPerPixel == 3)  /* when there are 3 sensors */
    this->file_color_mode = IM_RGB;

  if ((Photometric == PHOTOMETRIC_CFA || Photometric == PHOTOMETRIC_LINEARRAW) && BitsPerSample == 12)
    this->convert_bpp = 12;

  if (Photometric == PHOTOMETRIC_YCBCR && this->file_color_mode == IM_YCBCR)
  {
    uint16 ycbcrsubsampling[2];
    TIFFGetFieldDefaulted(this->tiff, TIFFTAG_YCBCRSUBSAMPLING, &ycbcrsubsampling[0], &ycbcrsubsampling[1]);
    this->h_subsample = ycbcrsubsampling[0];
    this->v_subsample = ycbcrsubsampling[1];

    /* add space for the line buffer (this is more than necessary) */
    this->line_buffer_extra = (int)TIFFScanlineSize(this->tiff);
  }

  uint16 PlanarConfig = PLANARCONFIG_CONTIG;
  TIFFGetFieldDefaulted(this->tiff, TIFFTAG_PLANARCONFIG, &PlanarConfig);

  if (PlanarConfig == PLANARCONFIG_CONTIG && SamplesPerPixel > 1)
    this->file_color_mode |= IM_PACKED;

  uint16 ExtraSamples = 0, *SampleInfo;
  TIFFGetFieldDefaulted(this->tiff, TIFFTAG_EXTRASAMPLES, &ExtraSamples, &SampleInfo);
  if (ExtraSamples == 1)
  {
    switch (SampleInfo[0]) 
    {
    case EXTRASAMPLE_UNSPECIFIED: /* !unspecified data */
    case EXTRASAMPLE_ASSOCALPHA:  /* data is pre-multiplied */
    case EXTRASAMPLE_UNASSALPHA:  /* data is not pre-multiplied */
      this->file_color_mode |= IM_ALPHA;
      break;
    }
    attrib_table->Set("ExtraSampleInfo", IM_USHORT, 1, (void*)&SampleInfo[0]);
  }
  else if ((ExtraSamples > 1) && (PlanarConfig == PLANARCONFIG_CONTIG))
  {
    /* usually a multi band image, we read only one band */
    this->sample_size_no_extra = (BitsPerSample*(SamplesPerPixel-ExtraSamples) + 7)/8; 
    this->extra_sample_size = (BitsPerSample*SamplesPerPixel + 7)/8;

    /* add space for the line buffer (this is more than necessary) */
    this->line_buffer_extra = (int)TIFFScanlineSize(this->tiff);
  }

  uint16 SampleFormat = SAMPLEFORMAT_UINT;
  TIFFGetField(this->tiff, TIFFTAG_SAMPLEFORMAT, &SampleFormat);
  switch(SampleFormat)
  {
  case SAMPLEFORMAT_VOID:
  case SAMPLEFORMAT_UINT:
    if (BitsPerSample < 8)
    {
      if (BitsPerSample != 1 && BitsPerSample != 2 && BitsPerSample != 4)
        return IM_ERR_DATA;

      this->file_data_type = IM_BYTE;
      this->convert_bpp = BitsPerSample;
    }
    else if (BitsPerSample == 8)
      this->file_data_type = IM_BYTE;
    else if (BitsPerSample <= 16)
      this->file_data_type = IM_USHORT;
    else if (BitsPerSample <= 32)
    {
      this->switch_type = 1;             // switch unsigned to signed
      this->file_data_type = IM_INT;
    }
    else
      return IM_ERR_DATA;
    break;
  case SAMPLEFORMAT_INT:
    if (BitsPerSample <= 8)
    {
      this->switch_type = 1;             // switch signed to unsigned
      this->file_data_type = IM_BYTE;
    }
    else if (BitsPerSample <= 16)
      this->file_data_type = IM_SHORT;
    else if (BitsPerSample <= 32)
      this->file_data_type = IM_INT;
    else
      return IM_ERR_DATA;
    break;
  case SAMPLEFORMAT_IEEEFP:
    if (BitsPerSample == 32)
      this->file_data_type = IM_FLOAT;      
    else if (BitsPerSample == 64)
      this->file_data_type = IM_DOUBLE;   
    else
      return IM_ERR_DATA;
    break;
  case SAMPLEFORMAT_COMPLEXINT:
    if (BitsPerSample == 32)
    {
      this->cpx_int = 1;
      this->file_data_type = IM_CFLOAT;  // convert short to float
    }
    else if (BitsPerSample == 64)
    {
      this->cpx_int = 2;
      this->file_data_type = IM_CFLOAT;  // convert int to float     
    }
    else
      return IM_ERR_DATA;
    break;
  case SAMPLEFORMAT_COMPLEXIEEEFP:
    if (BitsPerSample == 64)
      this->file_data_type = IM_CFLOAT;      
    else if (BitsPerSample == 128)
      this->file_data_type = IM_CDOUBLE;
    else
      return IM_ERR_DATA;
    break;
  default:
    return IM_ERR_DATA;
  }

  uint16 *rmap, *gmap, *bmap; 
  if (TIFFGetField(this->tiff, TIFFTAG_COLORMAP, &rmap, &gmap, &bmap))
  {
    long pal[256];
    int pal_count = 1 << BitsPerSample;

    for (int c = 0; c < pal_count; c++)
    {
      pal[c] = imColorEncode((unsigned char)(rmap[c] >> 8),
                             (unsigned char)(gmap[c] >> 8),
                             (unsigned char)(bmap[c] >> 8));
    }

    imFileSetPalette(this, pal, pal_count);
  }

  if (TIFFIsTiled(this->tiff))
  {
    if (this->tile_buf)
    {
      for (int i = 0; i < this->tile_buf_count; i++)
        free(this->tile_buf[i]);
      free(this->tile_buf);
    }

    uint32 tileWidth, tileLength;
    TIFFGetField(this->tiff, TIFFTAG_TILEWIDTH, &tileWidth);
    TIFFGetField(this->tiff, TIFFTAG_TILELENGTH, &tileLength);
    this->tile_width = (int)tileWidth;
    this->tile_height = (int)tileLength;

    this->tile_buf_count = (Width + tileWidth-1) / tileWidth;
    if (PlanarConfig == PLANARCONFIG_SEPARATE)
      this->tile_buf_count *= SamplesPerPixel;
    this->tile_line_size = (int)TIFFTileRowSize(this->tiff);
    this->tile_line_raw_size = (int)TIFFScanlineSize(this->tiff);
    this->tile_start_lin = 0;

    this->tile_buf = (void**)malloc(sizeof(void*)*this->tile_buf_count);
    size_t tile_size = TIFFTileSize(this->tiff);
    for (int t = 0; t < this->tile_buf_count; t++)
      this->tile_buf[t] = malloc(tile_size);
  }

  if (SamplesPerPixel < imColorModeDepth(this->file_color_mode))
    return IM_ERR_DATA;

  if (SamplesPerPixel > 1 && imColorModeSpace(this->file_color_mode) == IM_GRAY)
  {
    /* multiband data, we read only one band */
    attrib_table->Set("MultiBandCount", IM_USHORT, 1, (void*)&SamplesPerPixel);
  }

  uint16 Orientation;
  TIFFGetFieldDefaulted(this->tiff, TIFFTAG_ORIENTATION, &Orientation);
  switch (Orientation) 
  {
  case ORIENTATION_TOPRIGHT:
  case ORIENTATION_RIGHTTOP:  
  case ORIENTATION_LEFTTOP:  
  case ORIENTATION_TOPLEFT:
    this->file_color_mode |= IM_TOPDOWN;
    break;
  }
  attrib_table->Set("Orientation", IM_USHORT, 1, (void*)&Orientation);

  iTIFFReadAttributes(this->tiff, attrib_table);

#ifdef IM_TIFF_DEBUG_RGBA
  this->file_color_mode = IM_RGB;
  this->file_data_type = IM_BYTE;
#endif

  return IM_ERR_NONE;
}

int imFileFormatTIFF::WriteImageInfo()
{
  this->file_color_mode = this->user_color_mode;
  this->file_data_type = this->user_data_type;

  this->lab_fix = 0;
  this->invert = 0;

  uint16 Compression = iTIFFCompCalc(this->compression, this->file_color_mode, this->file_data_type);
  if (Compression == (uint16)-1)
    return IM_ERR_COMPRESS;

  int comp_index = iTIFFGetCompIndex(Compression);
  strcpy(this->compression, iTIFFCompTable[comp_index]);

  TIFFSetField(this->tiff, TIFFTAG_COMPRESSION, Compression);

  uint32 Width = this->width;
  TIFFSetField(this->tiff, TIFFTAG_IMAGEWIDTH, Width);

  uint32 Height = this->height;
  TIFFSetField(this->tiff, TIFFTAG_IMAGELENGTH, Height);

  static uint16 colorspace2photometric [] =
  {
    PHOTOMETRIC_RGB,    
    PHOTOMETRIC_PALETTE,    
    PHOTOMETRIC_MINISBLACK,   
    PHOTOMETRIC_MINISBLACK, 
    PHOTOMETRIC_SEPARATED,   
    PHOTOMETRIC_YCBCR,  
    PHOTOMETRIC_CIELAB,
    (uint16)-1,          // Pure Luv not supported
    PHOTOMETRIC_LOGLUV   // LogLuv Saved as XYZ
  };

  uint16 Photometric = colorspace2photometric[imColorModeSpace(this->file_color_mode)];

  // Correction for sgi LogL
  if (Compression == COMPRESSION_SGILOG && Photometric == PHOTOMETRIC_MINISBLACK)
    Photometric = PHOTOMETRIC_LOGL;

  // Corrections for JPEG, automatic convert from RGB to YCbCr when writing
  if (Compression == COMPRESSION_JPEG && Photometric == PHOTOMETRIC_RGB)
    Photometric = PHOTOMETRIC_YCBCR;

  /* The �normal� PhotometricInterpretation for bilevel CCITT compressed data is WhiteIsZero.
     Although TIFF readers should process PhotometricInterpretation in BlackIsZero as well,
     some applications assume WhiteIsZero. So if these compressions are used switch to WhiteIsZero. */
  if (imColorModeSpace(this->file_color_mode)==IM_BINARY && 
      (Compression >= COMPRESSION_CCITTRLE && Compression <= COMPRESSION_CCITT_T6))
    Photometric = PHOTOMETRIC_MINISWHITE;

  imAttribTable* attrib_table = AttribTable();

  uint16* photometric_attrib = (uint16*)attrib_table->Get("Photometric");
  if (photometric_attrib)
  {
    if (*photometric_attrib == PHOTOMETRIC_MASK && Photometric == PHOTOMETRIC_MINISBLACK)
      Photometric = PHOTOMETRIC_MASK;
    else if (*photometric_attrib == PHOTOMETRIC_MINISWHITE && Photometric == PHOTOMETRIC_MINISBLACK)
      Photometric = PHOTOMETRIC_MINISWHITE;
    else if (*photometric_attrib == PHOTOMETRIC_ICCLAB && Photometric == PHOTOMETRIC_CIELAB)
      Photometric = PHOTOMETRIC_ICCLAB;
    else if (*photometric_attrib == PHOTOMETRIC_ITULAB && Photometric == PHOTOMETRIC_CIELAB)
      Photometric = PHOTOMETRIC_ITULAB;
  }

  if (Photometric == PHOTOMETRIC_CIELAB)
    this->lab_fix = 1;

  if (Photometric == PHOTOMETRIC_MINISWHITE)
    this->invert = 1;

  TIFFSetField(this->tiff, TIFFTAG_PHOTOMETRIC, Photometric);

  // This is the default, and many software assume/handle only this, so we force it.
  uint16 PlanarConfig = PLANARCONFIG_CONTIG;
  TIFFSetField(this->tiff, TIFFTAG_PLANARCONFIG, PlanarConfig);
  if (imColorModeDepth(this->file_color_mode) > 1)
    this->file_color_mode |= IM_PACKED;

  // Corrections for JPEG, must be set after Photometric and PlanarConfig
  if (Compression == COMPRESSION_JPEG && imColorModeSpace(this->file_color_mode) == IM_RGB)
    TIFFSetField(this->tiff, TIFFTAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB);

  // Compression options
  int* zip_quality = (int*)attrib_table->Get("ZIPQuality");
  if (zip_quality && (Compression == COMPRESSION_DEFLATE || Compression == COMPRESSION_ADOBE_DEFLATE))
    TIFFSetField(this->tiff, TIFFTAG_ZIPQUALITY, *zip_quality);

  if (Compression == COMPRESSION_JPEG)
  {
    int* jpeg_quality = (int*)attrib_table->Get("JPEGQuality");
    if (jpeg_quality)
      TIFFSetField(this->tiff, TIFFTAG_JPEGQUALITY, *jpeg_quality);
  }

  // This is the default, and many software assume/handle only this, so we force it.
  uint16 Orientation = ORIENTATION_TOPLEFT; 
  TIFFSetField(this->tiff, TIFFTAG_ORIENTATION, Orientation);
  this->file_color_mode |= IM_TOPDOWN;

  static uint16 datatype2format[] =
  {
    SAMPLEFORMAT_UINT,    
    SAMPLEFORMAT_INT,     
    SAMPLEFORMAT_UINT,  
    SAMPLEFORMAT_INT,     
    SAMPLEFORMAT_IEEEFP,  
    SAMPLEFORMAT_IEEEFP,
    SAMPLEFORMAT_COMPLEXIEEEFP,
    SAMPLEFORMAT_COMPLEXIEEEFP
  };
  uint16 SampleFormat = datatype2format[this->file_data_type];
  TIFFSetField(this->tiff, TIFFTAG_SAMPLEFORMAT, SampleFormat);

  uint16 BitsPerSample = (uint16)(imDataTypeSize(this->file_data_type)*8);
  if (imColorModeSpace(this->file_color_mode) == IM_BINARY) 
  {
    BitsPerSample = 1;
    this->convert_bpp = 1;
  }
  TIFFSetField(this->tiff, TIFFTAG_BITSPERSAMPLE, BitsPerSample);

  // Correction for Luv, this will change BitsperSample and SampleFormat
  if (Photometric == PHOTOMETRIC_LOGLUV || Photometric == PHOTOMETRIC_LOGL)
    TIFFSetField(this->tiff, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);

  uint16 SamplesPerPixel = (uint16)imColorModeDepth(this->file_color_mode);
  TIFFSetField(this->tiff, TIFFTAG_SAMPLESPERPIXEL, SamplesPerPixel);

  if (imColorModeHasAlpha(this->file_color_mode))
  {
    uint16 ExtraSamples = 1, SampleInfo[1] = {EXTRASAMPLE_UNASSALPHA};
    uint16* sample_info = (uint16*)attrib_table->Get("ExtraSampleInfo");
    if (sample_info) SampleInfo[0] = *sample_info;
    TIFFSetField(this->tiff, TIFFTAG_EXTRASAMPLES, ExtraSamples, SampleInfo);
  }

  if (imColorModeSpace(this->file_color_mode) == IM_MAP)
  {
    uint16 rmap[256], gmap[256], bmap[256];
    memset(rmap, 0, 256 * sizeof(uint16));
    memset(gmap, 0, 256 * sizeof(uint16));
    memset(bmap, 0, 256 * sizeof(uint16));

    unsigned char r, g, b;
    for (int c = 0; c < this->palette_count; c++)
    {
      imColorDecode(&r, &g, &b, this->palette[c]);
      rmap[c] = r;
      gmap[c] = g;
      bmap[c] = b;
      rmap[c] <<= 8;
      gmap[c] <<= 8;
      bmap[c] <<= 8;
    }

    TIFFSetField(this->tiff, TIFFTAG_COLORMAP, rmap, gmap, bmap);
  }

  // Force libTIFF to calculate best RowsPerStrip
  uint32 RowsPerStrip = (uint32)-1; 
  RowsPerStrip = TIFFDefaultStripSize(this->tiff, RowsPerStrip);
  TIFFSetField(this->tiff, TIFFTAG_ROWSPERSTRIP, RowsPerStrip);

  iTIFFWriteAttributes(this->tiff, attrib_table);

  return IM_ERR_NONE;
}

static void iTIFFExpandComplexInt(void* line_buffer, int count, int cpx_int)
{
  count *= 2;

  // conversion will be done in-place

  if (cpx_int == 1)
  {
    // convert short to float, expanding from 16 to 32 bits
    short* short_buffer = (short*)line_buffer;
    float* float_buffer = (float*)line_buffer;

    float_buffer += count-1; // from end to start
    short_buffer += count-1;

    for (int i = 0; i < count; i++)
      *float_buffer-- = (float)(*short_buffer--);
  }
  else
  {
    // convert int to float, same size not expanding    
    int*   int_buffer   = (int*)line_buffer;
    float* float_buffer = (float*)line_buffer;

    for (int i = 0; i < count; i++)
      *float_buffer++ = (float)(*int_buffer++);
  }
}

static void iTIFFExtraSamplesFix(unsigned char* line_buffer, int width, int sample_size_no_extra, int extra_sample_size, int plane)
{
  /* ignore all the other extra samples, here the samples are packed */
  for (int i = 1; i < width; i++)
  {
    memcpy(line_buffer + i*sample_size_no_extra, line_buffer + i*extra_sample_size + plane, sample_size_no_extra);
  }
}

static void iTIFFExpandSubSamplePacked(const unsigned char* src_line_buffer, unsigned char* dst_line_buffer, int width, int lin, int h_subsample, int v_subsample)
{
  //TODO: Check other data types.
  // (2,1) YYCbCr...
  // (2,2) Y0Y0Y1Y1CbCr...
  // (4,1) YYYYCbCr...
  // (4,2) Y0Y0Y0Y0Y1Y1Y1Y1CbCr...
  // (4,4) Y0Y0Y0Y0Y1Y1Y1Y1Y2Y2Y2Y2Y3Y3Y3Y3CbCr...  
  int yy = h_subsample*v_subsample + 2;
  int l = (lin%v_subsample)*v_subsample;
  for (int i=0; i <width; i++)
  {
    int dst = i*3;
    int j = i/h_subsample, r = i%h_subsample;
    int src = j*yy;
    dst_line_buffer[dst+0] = src_line_buffer[src+r+l];    // Y
    dst_line_buffer[dst+1] = src_line_buffer[src+yy-2];   // Cb
    dst_line_buffer[dst+2] = src_line_buffer[src+yy-1];   // Cr
  }
}

static void iTIFFExpandSubSamplePlanar(const unsigned char* src_line_buffer, unsigned char* dst_line_buffer, int width, int plane, int h_subsample)
{
  if (plane != 0)
  {
    for (int i=0; i <width; i++)
    {
      int j = i/h_subsample;
      dst_line_buffer[i] = src_line_buffer[j];
    }
  }
}

/*
For CIELab (PhotometricInterpretation = 8), the L* component is encoded in 8 bits as an unsigned integer
range [0,255], and encoded in 16 bits as an unsigned integer range [0,65535]. The a* and b* components
are encoded in 8 bits as signed integers range [-128,127], and encoded in 16 bits as signed integers range [-
32768,32767]. The 8 bit chrominance values are exactly equal to the 1976 CIE a* and b* values, while the
16 bit values are equal to 256 times the 1976 CIE a* and b* values.

For ICCLab (PhotometricInterpretation = 9), the L* component is encoded in 8 bits as an unsigned integer
range [0,255], and encoded in 16 bits as an unsigned integer range [0,65280]. The a* and b* components
are encoded in 8 bits as unsigned integers range [0,255], and encoded in 16 bits as unsigned integers range
[0,65535]. The 8 bit chrominance values are exactly equal to the 1976 CIE a* and b* values plus 128,
while the 16 bit values are equal to 256 times the 1976 CIE a* and b* values plus 32768 (this is also 256
times the 8 bit encoding). PhotometricInterpretation 9 is designed to match the encoding used by the ICC
profile specification.
*/

static void iTIFFLabFix(void* line_buffer, int width, int data_type, int is_new)
{
  if (data_type == IM_BYTE)
  {
    imbyte* byte_buffer = (imbyte*)line_buffer;

    int offset = 128;
    if (is_new) offset = -128;

    for (int i = 0; i < width; i++)
    {
      *(byte_buffer+1) = (imbyte)(*((char*)byte_buffer+1) + offset);
      *(byte_buffer+2) = (imbyte)(*((char*)byte_buffer+2) + offset);

      byte_buffer += 3;
    }
  }
  else if (data_type == IM_USHORT)
  {
    imushort* ushort_buffer = (imushort*)line_buffer;

    int offset = 32768;
    if (is_new) offset = -32768;

    for (int i = 0; i < width; i++)
    {
      *(ushort_buffer+1) = (imushort)(*((short*)ushort_buffer+1) + offset);
      *(ushort_buffer+2) = (imushort)(*((short*)ushort_buffer+2) + offset);

      ushort_buffer += 3;
    }
  }
  //TODO: do NOT know how it is encoded for other data types.
}

int imFileFormatTIFF::ReadTileline(void* buffer, int lin, int plane)
{
  int t;

  if (lin == 0)
    this->tile_start_lin = 0;

  if (lin == this->tile_start_lin + this->tile_height)
    this->tile_start_lin = lin;

  // load a line of tiles
  if (lin == this->tile_start_lin)
  {
    int x = 0;
    for (t = 0; t < this->tile_buf_count; t++)
    {
      if (TIFFReadTile(this->tiff, this->tile_buf[t], x, tile_start_lin, 0, (tsample_t)plane) <= 0)
        return -1;

      x += this->tile_width;
    }
  }

  int line_size = this->tile_line_size;
  int tile_line = lin - this->tile_start_lin;

  for (t = 0; t < this->tile_buf_count; t++)
  {
    if (t == this->tile_buf_count-1)
    {
      // At the last tile, compute the correct size
      int extra = this->tile_line_size*this->tile_buf_count - this->tile_line_raw_size;
      line_size -= extra;
    }

    memcpy(buffer, (imbyte*)(this->tile_buf[t]) + tile_line*this->tile_line_size, line_size);
    buffer = (imbyte*)(buffer) + line_size;
  }

  return 1;
}

#ifdef IM_TIFF_DEBUG_RGBA
static void iTIFFReadRGBA(TIFF* tif, int w, int h, imbyte* data)
{
  uint32* raster = (uint32*)malloc(w*h*4);
  if (TIFFReadRGBAImage(tif, w, h, raster, 0))
  {
    int count = w*h;
    imbyte* red = data;
    imbyte* green = data + count;
    imbyte* blue = data + 2*count;

    for (int i=0; i<count; i++)
    {
      // ABGR
      red[i] = imbyte((raster[i] >> 0) & 0xFF);
      green[i] = imbyte((raster[i] >> 8) & 0xFF);
      blue[i] = imbyte((raster[i] >> 16) & 0xFF);
    }
  }
  free(raster);
}
#endif

static void iTIFFInvertBits(void* line_buffer, int size)
{
  unsigned char* buf = (unsigned char*)line_buffer;
  for (int b = 0; b < size; b++)
  {
    *buf = ~(*buf);
    buf++;
  }
}

int imFileFormatTIFF::ReadImageData(void* data)
{
  int count = imFileLineBufferCount(this);

  imCounterTotal(this->counter, count, "Reading TIFF...");

#ifdef IM_TIFF_DEBUG_RGBA
  iTIFFReadRGBA(this->tiff, this->width, this->height, (imbyte*)data);
#else
  int lin = 0, plane = this->start_plane;
  for (int i = 0; i < count; i++)
  {
    if (TIFFIsTiled(this->tiff))
    {
      if (this->h_subsample != 1 || this->v_subsample != 1)
      {
        if (i%this->v_subsample==0)
        {
          if (ReadTileline((imbyte*)this->line_buffer+line_buffer_size, lin/this->v_subsample, (tsample_t)plane) <= 0)
            return IM_ERR_ACCESS;
        }

        if (this->file_color_mode & IM_PACKED)
          iTIFFExpandSubSamplePacked((imbyte*)this->line_buffer+line_buffer_size, (imbyte*)this->line_buffer, this->width, lin, this->h_subsample, this->v_subsample);
        else
          iTIFFExpandSubSamplePlanar((imbyte*)this->line_buffer+line_buffer_size, (imbyte*)this->line_buffer, this->width, plane, this->h_subsample);
      }
      else
      {
        if (ReadTileline(this->line_buffer, lin, (tsample_t)plane) <= 0)
          return IM_ERR_ACCESS;
      }
    }
    else
    {
      if (this->h_subsample != 1 || this->v_subsample != 1)
      {
        if (i%this->v_subsample==0)
        {
          if (TIFFReadScanline(this->tiff, (imbyte*)this->line_buffer+line_buffer_size, lin/this->v_subsample, (tsample_t)plane) <= 0)
            return IM_ERR_ACCESS;
        }

        if (this->file_color_mode & IM_PACKED)
          iTIFFExpandSubSamplePacked((imbyte*)this->line_buffer+line_buffer_size, (imbyte*)this->line_buffer, this->width, lin, this->h_subsample, this->v_subsample);
        else
          iTIFFExpandSubSamplePlanar((imbyte*)this->line_buffer+line_buffer_size, (imbyte*)this->line_buffer, this->width, plane, this->h_subsample);
      }
      else
      {
        if (TIFFReadScanline(this->tiff, this->line_buffer, lin, (tsample_t)plane) <= 0)
          return IM_ERR_ACCESS;
      }
    }

    if (this->invert && this->file_data_type == IM_BYTE)
      iTIFFInvertBits(this->line_buffer, this->line_buffer_size);

    if (this->cpx_int)
    {
      int line_count = imImageLineCount(this->width, this->user_color_mode);
      iTIFFExpandComplexInt(this->line_buffer, line_count, this->cpx_int);
    }

    if (this->lab_fix)
      iTIFFLabFix(this->line_buffer, this->width, this->file_data_type, 0);

    if (this->extra_sample_size)
      iTIFFExtraSamplesFix((imbyte*)this->line_buffer, this->width, this->sample_size_no_extra, this->extra_sample_size, plane);

    imFileLineBufferRead(this, data, lin, plane);

    if (!imCounterInc(this->counter))
      return IM_ERR_COUNTER;

    imFileLineBufferInc(this, &lin, &plane);
  }
#endif

  return IM_ERR_NONE;
}

int imFileFormatTIFF::WriteImageData(void* data)
{
  int count = imFileLineBufferCount(this);

  imCounterTotal(this->counter, count, "Writing TIFF...");

  int lin = 0, plane = 0;
  for (int i = 0; i < count; i++)
  {
    imFileLineBufferWrite(this, data, lin, plane);

    if (this->invert && this->file_data_type == IM_BYTE)
      iTIFFInvertBits(this->line_buffer, this->line_buffer_size);

    if (this->lab_fix)
      iTIFFLabFix(this->line_buffer, this->width, this->file_data_type, 1);

    if (TIFFWriteScanline(this->tiff, this->line_buffer, lin, (tsample_t)plane) <= 0)
      return IM_ERR_ACCESS;

    if (!imCounterInc(this->counter))
      return IM_ERR_COUNTER;

    imFileLineBufferInc(this, &lin, &plane);
  }

  this->image_count++;

  if (!TIFFWriteDirectory(this->tiff))
    return IM_ERR_ACCESS;

   return IM_ERR_NONE;
}

int imFormatTIFF::CanWrite(const char* compression, int color_mode, int data_type) const
{
  if (!compression)
    return IM_ERR_NONE;

  if (imColorModeSpace(color_mode) == IM_LUV)
    return IM_ERR_DATA;

  uint16 Compression = iTIFFCompCalc(compression, color_mode, data_type);
  if (Compression == (uint16)-1)
    return IM_ERR_COMPRESS;

  /* no support for 2 bpp or 4 bpp */
  if (Compression == COMPRESSION_THUNDERSCAN || Compression == COMPRESSION_NEXT)
    return IM_ERR_COMPRESS;

  /* Binary compression restrictions */
  if ((Compression == COMPRESSION_CCITTRLE || Compression == COMPRESSION_CCITTRLEW ||
       Compression == COMPRESSION_CCITTFAX3  || Compression == COMPRESSION_CCITTFAX4) &&
      imColorModeSpace(color_mode) != IM_BINARY)
    return IM_ERR_COMPRESS;

  /* JPEG compression restrictions */
  if (Compression == COMPRESSION_JPEG && 
      (data_type != IM_BYTE || 
       imColorModeSpace(color_mode) == IM_MAP || imColorModeSpace(color_mode) == IM_BINARY))
    return IM_ERR_COMPRESS;

  /* Pixar log accepts only 3 types */
  if (Compression == COMPRESSION_PIXARLOG && 
      data_type != IM_BYTE && data_type != IM_USHORT && data_type != IM_FLOAT)
    return IM_ERR_COMPRESS;

  /* SGI Luv compression restrictions */
  if ((Compression == COMPRESSION_SGILOG || Compression == COMPRESSION_SGILOG24) &&
      (imColorModeSpace(color_mode) != IM_XYZ || data_type != IM_FLOAT))
    return IM_ERR_COMPRESS;

  return IM_ERR_NONE;
}
