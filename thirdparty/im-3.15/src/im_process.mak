PROJNAME = im
LIBNAME = im_process
OPT = YES

SRC = \
    im_arithmetic_bin.cpp  im_morphology_gray.cpp  im_quantize.cpp   \
    im_arithmetic_un.cpp   im_geometric.cpp        im_render.cpp     \
    im_color.cpp           im_histogram.cpp        im_resize.cpp     \
    im_convolve.cpp        im_houghline.cpp        im_statistics.cpp \
    im_convolve_rank.cpp   im_logic.cpp            im_threshold.cpp  \
    im_effects.cpp         im_morphology_bin.cpp   im_tonegamut.cpp  \
    im_canny.cpp           im_distance.cpp         im_analyze.cpp    \
    im_kernel.cpp          im_remotesens.cpp       im_point.cpp      \
    im_process_counter.cpp im_binary.cpp
SRC := $(addprefix process/, $(SRC))

SRC += im_convertbitmap.cpp im_convertcolor.cpp im_converttype.cpp
                                       
USE_IM = Yes
IM = ..
DEFINES += IM_PROCESS
DEPENDDIR = dep

ifdef USE_OPENMP
  DEF_FILE := $(LIBNAME).def
  LIBNAME := $(LIBNAME)_omp
endif

ifneq ($(findstring Win, $(TEC_SYSNAME)), )
    ifneq ($(findstring ow, $(TEC_UNAME)), )
      DEFINES += IM_DEFMATHFLOAT
    endif  
    ifneq ($(findstring bc, $(TEC_UNAME)), )
      DEFINES += IM_DEFMATHFLOAT
    endif  
else
  ifneq ($(findstring AIX, $(TEC_UNAME)), )
    DEFINES += IM_DEFMATHFLOAT 
  endif
  ifneq ($(findstring SunOS, $(TEC_UNAME)), )
    DEFINES += IM_DEFMATHFLOAT
  endif
  ifneq ($(findstring HP-UX, $(TEC_UNAME)), )
    DEFINES += IM_DEFMATHFLOAT
  endif
  ifneq ($(findstring MacOS, $(TEC_UNAME)), )
    ifneq ($(TEC_SYSMINOR), 4)
      BUILD_DYLIB=Yes
    endif
  endif
endif
