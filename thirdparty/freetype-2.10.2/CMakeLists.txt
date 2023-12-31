cmake_minimum_required (VERSION 3.21)

project (freetype
    VERSION 2.10.2
    LANGUAGES C)

option (FREETYPE_PIC "Compile with position independent code" ${FTGL_PIC})

add_library (freetype
    "src/autofit/autofit.c"
    "src/bdf/bdf.c"
    "src/bzip2/ftbzip2.c"
    "src/cache/ftcache.c"
    "src/cff/cff.c"
    "src/cid/type1cid.c"
    "src/gxvalid/gxvalid.c"
    "src/gzip/ftgzip.c"
    "src/lzw/ftlzw.c"
    "src/otvalid/otvalid.c"
    "src/pcf/pcf.c"
    "src/pfr/pfr.c"
    "src/psaux/psaux.c"
    "src/pshinter/pshinter.c"
    "src/psnames/psnames.c"
    "src/raster/raster.c"
    "src/sfnt/sfnt.c"
    "src/smooth/smooth.c"
    "src/truetype/truetype.c"
    "src/type1/type1.c"
    "src/type42/type42.c"
    "src/winfonts/winfnt.c"
    "src/base/ftapi.c"
    "src/base/ftbase.c"
    "src/base/ftbbox.c"
    "src/base/ftbdf.c"
    "src/base/ftbitmap.c"
    "src/base/ftcid.c"
    "src/base/ftdebug.c"
    "src/base/ftfstype.c"
    "src/base/ftgasp.c"
    "src/base/ftglyph.c"
    "src/base/ftgxval.c"
    "src/base/ftinit.c"
    "src/base/ftmm.c"
    "src/base/ftotval.c"
    "src/base/ftpatent.c"
    "src/base/ftpfr.c"
    "src/base/ftstroke.c"
    "src/base/ftsynth.c"
    "src/base/ftsystem.c"
    "src/base/fttype1.c"
    "src/base/ftwinfnt.c")
target_include_directories (freetype
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_compile_definitions (freetype
    PRIVATE
        FT2_BUILD_LIBRARY
        FT_CONFIG_OPTION_SYSTEM_ZLIB)
target_link_libraries (freetype
    PRIVATE
        ZLIB::ZLIB)
set_target_properties (freetype PROPERTIES
    MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>
    POSITION_INDEPENDENT_CODE ${FREETYPE_PIC})
