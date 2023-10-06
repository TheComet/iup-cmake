PROJNAME = im
LIBNAME = im_fftw3
OPT = YES

DEF_FILE = im_fftw3.def
DEPENDDIR = dep

SRC = process/im_fft.cpp

INCLUDES = ../include

DEFINES = USE_FFTW3

USE_IM = Yes
IM = ..
LIBS = im_process

ifneq ($(findstring Win, $(TEC_SYSNAME)), )
  # Windows use local
  FFTW = $(TECTOOLS_HOME)/fftw3
  ifneq ($(findstring _64, $(TEC_UNAME)), )
    LDIR = $(FFTW)/lib/Win64
  else
    LDIR = $(FFTW)/lib/Win32
  endif
  INCLUDES += $(FFTW)/include
  LIBS += libfftw3f-3 libfftw3-3
else  
  # Linux use system
  LIBS += fftw3f fftw3
endif


ifneq ($(findstring Win, $(TEC_SYSNAME)), )
  ifneq ($(findstring ow, $(TEC_UNAME)), )
    DEFINES += IM_DEFMATHFLOAT
  endif         
  ifneq ($(findstring bc, $(TEC_UNAME)), )
    DEFINES += IM_DEFMATHFLOAT
  endif
else
  ifneq ($(findstring MacOS, $(TEC_UNAME)), )
    ifneq ($(TEC_SYSMINOR), 4)
      BUILD_DYLIB=Yes
    endif
  endif
  ifneq ($(findstring AIX, $(TEC_UNAME)), )
    DEFINES += IM_DEFMATHFLOAT
  endif
  ifneq ($(findstring SunOS, $(TEC_UNAME)), )
    DEFINES += IM_DEFMATHFLOAT
  endif
endif

ifneq ($(findstring MacOS, $(TEC_UNAME)), )
  BUILD_DYLIB=Yes
endif

