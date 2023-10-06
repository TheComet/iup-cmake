APPNAME = im_copy
APPTYPE = console

SRC = im_copy.cpp

USE_IM = Yes

IM = ..

USE_STATIC = Yes

ifneq ($(findstring Win, $(TEC_SYSNAME)), )
  LIBS = im_wmv im_avi vfw32 wmvcore
endif
