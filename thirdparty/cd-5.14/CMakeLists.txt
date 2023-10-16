project (cd
    VERSION 5.14
    LANGUAGES C)

option (CD_PIC "Compile with position independent code" OFF)

add_library (cd
    # common
    "src/cd.c"
    "src/wd.c"
    "src/wdhdcpy.c"
    "src/rgb2map.c"
    "src/cd_vectortext.c"
    "src/cd_active.c"
    "src/cd_attributes.c"
    "src/cd_bitmap.c"
    "src/cd_image.c"
    "src/cd_primitives.c"
    "src/cd_text.c"
    "src/cd_util.c"

    # svg
    "src/svg/base64.c"
    "src/svg/lodepng.c"
    "src/svg/cdsvg.c"

    #intcgm
    "src/intcgm/cd_intcgm.c"
    "src/intcgm/cgm_bin_get.c"
    "src/intcgm/cgm_bin_parse.c"
    "src/intcgm/cgm_list.c"
    "src/intcgm/cgm_play.c"
    "src/intcgm/cgm_sism.c"
    "src/intcgm/cgm_txt_get.c"
    "src/intcgm/cgm_txt_parse.c"

    # drv
    "src/drv/cddgn.c"
    "src/drv/cdcgm.c"
    "src/drv/cgm.c"
    "src/drv/cddxf.c"
    "src/drv/cdirgb.c"
    "src/drv/cdmf.c"
    "src/drv/cdps.c"
    "src/drv/cdpicture.c"
    "src/drv/cddebug.c"
    "src/drv/cdpptx.c"
    "src/drv/pptx.c"

    # minizip
    "src/minizip/ioapi.c"
    "src/minizip/minizip.c"
    "src/minizip/zip.c"
    "src/minizip/miniunzip.c"
    "src/minizip/unzip.c"

    # sim
    "src/sim/cdfontex.c"
    "src/sim/sim.c"
    "src/sim/cd_truetype.c"
    "src/sim/sim_primitives.c"
    "src/sim/sim_text.c"
    "src/sim/sim_linepolyfill.c")
target_include_directories (cd
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/drv>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/x11>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/win32>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/intcgm>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/sim>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/cairo>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/svg>)
target_link_libraries (cd
    PRIVATE
        freetype
        ZLIB::ZLIB)
set_target_properties (cd PROPERTIES
    MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>
    POSITION_INDEPENDENT_CODE ${CD_PIC})

if (IUP_TOOLKIT MATCHES "GTK2" OR IUP_TOOLKIT MATCHES "GTK3")
    find_package (Cairo REQUIRED)
    find_package (Pango 1.0 REQUIRED)
    find_package (HarfBuzz REQUIRED)

    target_sources (cd PRIVATE
        "src/cairo/cdcairodbuf.c"
        "src/cairo/cdcairopdf.c"
        "src/cairo/cdcairosvg.c"
        "src/cairo/cdcairo.c"
        "src/cairo/cdcairoimg.c"
        "src/cairo/cdcairoirgb.c"
        "src/cairo/cdcairops.c"
        "src/cairo/cdcaironative_gdk.c")
    target_include_directories (cd
        PRIVATE
            ${HARFBUZZ_INCLUDE_DIR})
    target_link_libraries (cd
        PRIVATE
            Cairo::Cairo
            pango-1.0
            freetype
            ${HARFBUZZ_LIBRARY})

    if (IUP_TOOLKIT MATCHES "GTK3")
        find_package (GTK3 REQUIRED)
        target_sources (cd PRIVATE
            "src/gdk/cdgdkclp.c")
        target_compile_definitions (cd PRIVATE
            GDK_DISABLE_DEPRECATED
            GSEAL_ENABLE
            USE_GTK3)
        target_include_directories (cd PRIVATE
            ${GTK3_INCLUDE_DIRS})
    else ()
        find_package (GTK2 REQUIRED)
        target_sources (cd PRIVATE
            "src/gdk/cdgdk.c"
            "src/gdk/cdgdkclp.c"
            "src/gdk/cdgdkdbuf.c"
            "src/gdk/cdgdkimg.c"
            "src/gdk/cdgdknative.c"
            "src/drv/cd0prn.c"
            "src/drv/cd0emf.c")
        target_link_libraries (cd
            PRIVATE
                GTK2::gtk
                GTK2::gdk_pixbuf)
    endif ()
elseif (IUP_TOOLKIT MATCHES "Win32")
    target_sources (cd PRIVATE
        "src/win32/cdwclp.c"
        "src/win32/cdwemf.c"
        "src/win32/cdwimg.c"
        "src/win32/cdwin.c"
        "src/win32/cdwnative.c"
        "src/win32/cdwprn.c"
        "src/win32/cdwwmf.c"
        "src/win32/wmf_emf.c"
        "src/win32/cdwdbuf.c"
        "src/win32/cdwdib.c")
    target_compile_definitions (cd PRIVATE UNICODE WIN32)
endif ()
