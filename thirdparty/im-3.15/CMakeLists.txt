project (im
    LANGUAGES C CXX
    VERSION 3.15)

option (IM_PIC "Compile with position independent code" OFF)

include (TestBigEndian)
test_big_endian (WORDS_BIGENDIAN)

set (SRCTIFF
    src/libtiff/tif_aux.c       src/libtiff/tif_dirwrite.c   src/libtiff/tif_jpeg.c      src/libtiff/tif_print.c
    src/libtiff/tif_close.c     src/libtiff/tif_dumpmode.c   src/libtiff/tif_luv.c       src/libtiff/tif_read.c
    src/libtiff/tif_codec.c     src/libtiff/tif_error.c      src/libtiff/tif_lzw.c       src/libtiff/tif_strip.c
    src/libtiff/tif_color.c     src/libtiff/tif_extension.c  src/libtiff/tif_next.c      src/libtiff/tif_swab.c
    src/libtiff/tif_compress.c  src/libtiff/tif_fax3.c       src/libtiff/tif_open.c      src/libtiff/tif_thunder.c
    src/libtiff/tif_dir.c       src/libtiff/tif_fax3sm.c     src/libtiff/tif_packbits.c  src/libtiff/tif_tile.c
    src/libtiff/tif_dirinfo.c   src/libtiff/tif_flush.c      src/libtiff/tif_pixarlog.c  src/libtiff/tif_zip.c
    src/libtiff/tif_dirread.c   src/libtiff/tif_getimage.c   src/libtiff/tif_predict.c   src/libtiff/tif_version.c
    src/libtiff/tif_write.c     src/libtiff/tif_warning.c    src/libtiff/tif_ojpeg.c     src/libtiff/tif_lzma.c
    src/libtiff/tif_jbig.c)
set (SRCPNG
    src/libpng/png.c       src/libpng/pngget.c    src/libpng/pngread.c   src/libpng/pngrutil.c  src/libpng/pngwtran.c
    src/libpng/pngerror.c  src/libpng/pngmem.c    src/libpng/pngrio.c    src/libpng/pngset.c    src/libpng/pngwio.c
    src/libpng/pngpread.c  src/libpng/pngrtran.c  src/libpng/pngtrans.c  src/libpng/pngwrite.c  src/libpng/pngwutil.c)
set (SRCJPEG
    src/libjpeg/jcapimin.c  src/libjpeg/jcmarker.c  src/libjpeg/jdapimin.c  src/libjpeg/jdinput.c   src/libjpeg/jdtrans.c
    src/libjpeg/jcapistd.c  src/libjpeg/jcmaster.c  src/libjpeg/jdapistd.c  src/libjpeg/jdmainct.c  src/libjpeg/jerror.c   src/libjpeg/jmemmgr.c
    src/libjpeg/jccoefct.c  src/libjpeg/jcomapi.c   src/libjpeg/jdatadst.c  src/libjpeg/jdmarker.c  src/libjpeg/jfdctflt.c src/libjpeg/jmemnobs.c
    src/libjpeg/jccolor.c   src/libjpeg/jcparam.c   src/libjpeg/jdatasrc.c  src/libjpeg/jdmaster.c  src/libjpeg/jfdctfst.c src/libjpeg/jquant1.c
    src/libjpeg/jcdctmgr.c  src/libjpeg/jdcoefct.c  src/libjpeg/jdmerge.c   src/libjpeg/jfdctint.c  src/libjpeg/jquant2.c
    src/libjpeg/jchuff.c    src/libjpeg/jcprepct.c  src/libjpeg/jdcolor.c   src/libjpeg/jidctflt.c  src/libjpeg/jutils.c   src/libjpeg/jdarith.c
    src/libjpeg/jcinit.c    src/libjpeg/jcsample.c  src/libjpeg/jddctmgr.c  src/libjpeg/jdpostct.c  src/libjpeg/jidctfst.c src/libjpeg/jaricom.c
    src/libjpeg/jcmainct.c  src/libjpeg/jctrans.c   src/libjpeg/jdhuff.c    src/libjpeg/jdsample.c  src/libjpeg/jidctint.c src/libjpeg/jcarith.c)
set (SRCEXIF
    src/libexif/fuji/exif-mnote-data-fuji.c  src/libexif/fuji/mnote-fuji-entry.c  src/libexif/fuji/mnote-fuji-tag.c
    src/libexif/canon/exif-mnote-data-canon.c  src/libexif/canon/mnote-canon-entry.c  src/libexif/canon/mnote-canon-tag.c
    src/libexif/olympus/exif-mnote-data-olympus.c  src/libexif/olympus/mnote-olympus-entry.c  src/libexif/olympus/mnote-olympus-tag.c
    src/libexif/pentax/exif-mnote-data-pentax.c  src/libexif/pentax/mnote-pentax-entry.c  src/libexif/pentax/mnote-pentax-tag.c
    src/libexif/exif-byte-order.c  src/libexif/exif-entry.c  src/libexif/exif-utils.c    src/libexif/exif-format.c  src/libexif/exif-mnote-data.c
    src/libexif/exif-content.c  src/libexif/exif-ifd.c  src/libexif/exif-tag.c src/libexif/exif-data.c  src/libexif/exif-loader.c src/libexif/exif-log.c src/libexif/exif-mem.c)
set (SRCLZF
    src/liblzf/lzf_c.c src/liblzf/lzf_d.c)
set (SRCLZ4
    src/lz4/lz4.c)

add_library (im
    src/im_oldcolor.c         src/im_oldresize.c      src/im_converttype.cpp
    src/im_attrib.cpp         src/im_format.cpp       src/im_format_tga.cpp    src/im_filebuffer.cpp
    src/im_bin.cpp            src/im_format_all.cpp   src/im_format_raw.cpp    src/im_convertopengl.cpp
    src/im_binfile.cpp        src/im_format_sgi.cpp   src/im_datatype.cpp      src/im_format_pcx.cpp
    src/im_colorhsi.cpp       src/im_format_bmp.cpp   src/im_image.cpp         src/im_rgb2map.cpp
    src/im_colormode.cpp      src/im_format_gif.cpp   src/im_lib.cpp           src/im_format_pnm.cpp
    src/im_colorutil.cpp      src/im_format_ico.cpp   src/im_palette.cpp       src/im_format_ras.cpp
    src/im_convertbitmap.cpp  src/im_format_led.cpp   src/im_counter.cpp       src/im_str.cpp
    src/im_convertcolor.cpp   src/im_fileraw.cpp      src/im_format_krn.cpp    src/im_compress.cpp
    src/im_file.cpp           src/im_old.cpp          src/im_format_pfm.cpp
    src/im_format_tiff.cpp    src/im_format_png.cpp   src/im_format_jpeg.cpp
    ${SRCTIFF} src/tiff_binfile.c  # NOT necessary if using another distribution
    ${SRCLZF}
    ${SRCLZ4}
    ${SRCJPEG})
target_compile_definitions (im
    PRIVATE
        $<$<BOOL:${WORDS_BIGENDIAN}>:WORDS_BIGENDIAN>)
target_include_directories (im
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/lz4>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/liblzf>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/libtiff>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/libjpeg>)

# tif requires zlib -- higher level project should include zlib before im
target_link_libraries (im PRIVATE ZLIB::ZLIB)

if (CMAKE_SYSTEM_NAME MATCHES Windows)
    find_package (PNG)
    if (NOT PNG_FOUND)
        target_sources (im PRIVATE
            ${SRCPNG})
        target_include_directories (im PRIVATE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/libpng>)
    endif ()

    target_sources (im PRIVATE
        src/im_sysfile_win32.cpp
        src/im_dib.cpp
        src/im_dibxbitmap.cpp
        src/im.rc)
else ()
    # On Linux, use the installed files in the system package libpng-dev
    # If using GTK, then must use the same libpng they use
    find_package (PNG REQUIRED)
    target_link_libraries (im PRIVATE PNG::PNG)
endif ()

find_package (EXIF)
if (NOT EXIF_FOUND)
    target_sources (im PRIVATE
        ${SRCEXIF})
    target_compile_definitions (im PRIVATE USE_EXIF)
    target_include_directories (im PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/libexif>)
endif ()

set_target_properties (im PROPERTIES
    DEBUG_POSTFIX "d"
    MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>
    POSITION_INDEPENDENT_CODE ${IM_PIC})
