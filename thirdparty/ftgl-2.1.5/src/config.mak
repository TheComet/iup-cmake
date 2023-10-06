PROJNAME = ftgl
LIBNAME = ftgl
OPT = YES

ifdef DBG
  DEFINES += IUP_ASSERT
  ifndef DBG_DIR
    ifneq ($(findstring Win, $(TEC_SYSNAME)), )
      LIBNAME := $(LIBNAME)_debug
    endif
  endif  
endif  

DEF_FILE = ftgl.def

ftglyph_sources = \
    FTGlyph/FTGlyph.cpp \
    FTGlyph/FTGlyphGlue.cpp \
    FTGlyph/FTBitmapGlyph.cpp \
    FTGlyph/FTBufferGlyph.cpp \
    FTGlyph/FTExtrudeGlyph.cpp \
    FTGlyph/FTOutlineGlyph.cpp \
    FTGlyph/FTPixmapGlyph.cpp \
    FTGlyph/FTPolygonGlyph.cpp \
    FTGlyph/FTTextureGlyph.cpp

ftfont_sources = \
    FTFont/FTFont.cpp \
    FTFont/FTFontGlue.cpp \
    FTFont/FTBitmapFont.cpp \
    FTFont/FTBufferFont.cpp \
    FTFont/FTExtrudeFont.cpp \
    FTFont/FTOutlineFont.cpp \
    FTFont/FTPixmapFont.cpp \
    FTFont/FTPolygonFont.cpp \
    FTFont/FTTextureFont.cpp

ftlayout_sources = \
    FTLayout/FTLayout.cpp \
    FTLayout/FTLayoutGlue.cpp \
    FTLayout/FTSimpleLayout.cpp
    
SRC = \
    FTBuffer.cpp \
    FTCharmap.cpp \
    FTContour.cpp \
    FTFace.cpp \
    FTGlyphContainer.cpp \
    FTLibrary.cpp \
    FTPoint.cpp \
    FTSize.cpp \
    FTVectoriser.cpp \
    $(ftglyph_sources) \
    $(ftfont_sources) \
    $(ftlayout_sources)
    
INCLUDES := ../include .
LDIR = ../lib/$(TEC_UNAME)

DEFINES = FTGL_LIBRARY_STATIC
ifneq ($(findstring dll, $(TEC_UNAME)), )
  DEFINES = FTGL_LIBRARY
endif

USE_FREETYPE = Yes
USE_OPENGL = Yes

ifneq ($(findstring MacOS, $(TEC_UNAME)), )
  ifneq ($(TEC_SYSMINOR), 4)
    BUILD_DYLIB=Yes
  endif
endif
