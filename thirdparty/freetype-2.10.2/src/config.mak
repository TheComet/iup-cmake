PROJNAME = freetype
LIBNAME = freetype
OPT = YES   

# Changes to freetype source code => search for "CDLIB"

SRC  := \
  autofit/autofit.c bdf/bdf.c bzip2/ftbzip2.c cache/ftcache.c cff/cff.c cid/type1cid.c \
  gxvalid/gxvalid.c gzip/ftgzip.c lzw/ftlzw.c otvalid/otvalid.c pcf/pcf.c \
  pfr/pfr.c psaux/psaux.c pshinter/pshinter.c psnames/psnames.c raster/raster.c \
  sfnt/sfnt.c smooth/smooth.c truetype/truetype.c type1/type1.c \
  type42/type42.c winfonts/winfnt.c \
  base/ftapi.c base/ftbase.c base/ftbbox.c base/ftbdf.c base/ftbitmap.c base/ftcid.c base/ftdebug.c \
  base/ftfstype.c base/ftgasp.c base/ftglyph.c base/ftgxval.c base/ftinit.c \
  base/ftmm.c base/ftotval.c base/ftpatent.c base/ftpfr.c base/ftstroke.c \
  base/ftsynth.c base/ftsystem.c base/fttype1.c base/ftwinfnt.c

DEFINES += FT2_BUILD_LIBRARY FT_CONFIG_OPTION_SYSTEM_ZLIB
INCLUDES = ../include
LDIR = ../lib/$(TEC_UNAME)

USE_ZLIB = Yes

ifneq ($(findstring dll, $(TEC_UNAME)), )
  DEFINES += DLL_EXPORT
  SRC += base/ftver.rc
  DEF_FILE = freetype.def
endif

ifneq ($(findstring Win, $(TEC_SYSNAME)), )
  # To be compatible with the existing DLLs of gnuwin32
  LIBNAME = freetype6
endif
ifneq ($(findstring cygw, $(TEC_UNAME)), )
  # To be compatible with the existing DLLs of cygwin
  #LIBNAME = freetype-6
endif

ifneq ($(findstring bc5, $(TEC_UNAME)), )
  FLAGS = -w-8004
endif

ifneq ($(findstring MacOS, $(TEC_UNAME)), )
  ifneq ($(TEC_SYSMINOR), 4)
    BUILD_DYLIB=Yes
  endif
  DEFINES += DARWIN_NO_CARBON
endif
