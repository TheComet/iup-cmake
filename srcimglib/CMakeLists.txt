project (iupimglib
    LANGUAGES C
    VERSION ${IUP_VERSION})

add_library (iupimglib
    "iup_image_library.c"
    "iup_imglib_circleprogress.c")
target_include_directories (iupimglib
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../src>)
target_compile_definitions (iupimglib
    PRIVATE
        $<$<PLATFORM_ID:Windows>:WIN32>
        $<$<CONFIG:Debug>:IUP_ASSERT>
        IUPIMGLIB_BUILD_LIBRARY
        IUP_IMGLIB_LARGE_ICON
        $<$<BOOL:${IUP_IMGLIB_LARGE}>:IUP_IMGLIB_LARGE>)

if (IUP_TOOLKIT MATCHES "Motif")
    target_sources (iupimglib PRIVATE
        "iup_imglib_basemot16x16_8bpp.c"
        "iup_imglib_logos48x48_8bpp.c"
        "iup_imglib_logos32x32_8bpp.c")
    target_compile_definitions (iupimglib
        PRIVATE
            MOTIF)
elseif (IUP_TOOLKIT MATCHES "Win32")
    target_sources (iupimglib PRIVATE
        "iup_imglib_logos48x48.c"
        "iup_imglib_logos32x32.c")
elseif (IUP_TOOLKIT MATCHES "GTK2")
    target_sources (iupimglib PRIVATE
        "iup_imglib_logos48x48.c"
        "iup_imglib_logos32x32.c"
            "iup_imglib_basegtk24x24.c"
            "iup_imglib_iconsgtk48x48.c")
elseif (IUP_TOOLKIT MATCHES "GTK3")
    target_sources (iupimglib PRIVATE
        "iup_imglib_logos48x48.c"
        "iup_imglib_logos32x32.c"
        "iup_imglib_basegtk324x24.c"
        "iup_imglib_iconsgtk348x48.c")
    target_compile_definitions (iupimglib PRIVATE GTK3)
else ()
    target_sources (iupimglib PRIVATE
        "iup_imglib_basemot16x16_8bpp.c"
        "iup_imglib_logos48x48_8bpp.c"
        "iup_imglib_logos32x32_8bpp.c")
endif ()
