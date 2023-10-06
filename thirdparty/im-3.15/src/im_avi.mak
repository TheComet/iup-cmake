PROJNAME = im
LIBNAME = im_avi
OPT = YES

SRC = im_format_avi.cpp
DEPENDDIR = dep
                      
LIBS = vfw32

USE_IM=Yes
IM = ..

ifneq ($(findstring cygw, $(TEC_UNAME)), )
  $(error No support for AVI in Cygwin when using Posix mode)
endif
