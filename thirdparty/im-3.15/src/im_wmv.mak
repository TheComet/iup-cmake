PROJNAME = im
LIBNAME = im_wmv
OPT = YES

SRC = im_format_wmv.cpp

ifneq ($(findstring vc9, $(TEC_UNAME)), )
  USE_WMFSDK = Yes
endif
ifneq ($(findstring vc10, $(TEC_UNAME)), )
  USE_WMFSDK = Yes
endif
ifneq ($(findstring dll9, $(TEC_UNAME)), )
  USE_WMFSDK = Yes
endif
ifneq ($(findstring dll10, $(TEC_UNAME)), )
  USE_WMFSDK = Yes
endif

ifndef USE_WMFSDK
  #vc9-vc10 needs an external SDK
  ifneq ($(findstring _64, $(TEC_UNAME)), )
    WMFSDK = d:/lng/wmfsdk95
  else
  #  WMFSDK = d:/lng/wmfsdk11
    WMFSDK = d:/lng/wmfsdk9
  endif
  INCLUDES = $(WMFSDK)/include
  LDIR = $(WMFSDK)/lib
else
  # > vc11, wmf sdk is inside Windows SDK
endif
  
DEFINES = _CRT_NON_CONFORMING_SWPRINTFS                                     
DEPENDDIR = dep

LIBS = wmvcore

USE_IM = Yes
IM = ..

ifneq ($(findstring gcc, $(TEC_UNAME)), )
  $(error No support for WMFSDK in Cygwin)
endif
ifneq ($(findstring mingw, $(TEC_UNAME)), )
  $(error No support for WMFSDK in MingW)
endif
ifneq ($(findstring cygw, $(TEC_UNAME)), )
  $(error No support for WMFSDK in Cygwin)
endif
ifneq ($(findstring dllw, $(TEC_UNAME)), )
  $(error No support for WMFSDK in MingW)
endif
ifneq ($(findstring dllg, $(TEC_UNAME)), )
  $(error No support for WMFSDK in Cygwin)
endif
ifneq ($(findstring owc, $(TEC_UNAME)), )
  $(error No support for WMFSDK in OpenWatcom)
endif
ifneq ($(findstring bc, $(TEC_UNAME)), )
  $(error No support for WMFSDK in BorlandC)
endif
