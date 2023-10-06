project (iupimglib
    LANGUAGES C
    VERSION ${IUP_VERSION})

add_library (iupimglib ${IUP_LIB_TYPE}
    "iup_image_library.c"
    "iup_imglib_circleprogress.c")
target_include_directories (iupimglib
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../src>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_compile_definitions (iupimglib
    PRIVATE
        $<$<PLATFORM_ID:Windows>:WIN32>
        $<$<CONFIG:Debug>:IUP_ASSERT>
        IUPIMGLIB_BUILD_LIBRARY
        IUP_IMGLIB_LARGE_ICON
        $<$<BOOL:${IUP_IMGLIB_LARGE}>:IUP_IMGLIB_LARGE>)
target_link_libraries (iupimglib
    PRIVATE
        iup)
set_target_properties (iupimglib
    PROPERTIES
        DEBUG_POSTFIX "d"
        MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>
        RUNTIME_OUTPUT_DIRECTORY ${IUP_BUILD_BINDIR}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG ${IUP_BUILD_BINDIR}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE ${IUP_BUILD_BINDIR}
        LIBRARY_OUTPUT_DIRECTORY ${IUP_BUILD_LIBDIR}
        LIBRARY_OUTPUT_DIRECTORY_DEBUG ${IUP_BUILD_LIBDIR}
        LIBRARY_OUTPUT_DIRECTORY_RELEASE ${IUP_BUILD_LIBDIR}
        ARCHIVE_OUTPUT_DIRECTORY ${IUP_BUILD_LIBDIR}
        ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${IUP_BUILD_LIBDIR}
        ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${IUP_BUILD_LIBDIR})

if (CMAKE_SYSTEM_NAME MATCHES "Windows" AND IUP_LIB_TYPE MATCHES "SHARED")
    target_sources (iupimglib PRIVATE "iupimglib.def")
endif ()

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
        "iup_imglib_basewin32x32.c"
        "iup_imglib_iconswin48x48.c"
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

add_library (IUP::imglib ALIAS iupimglib)

install (
    TARGETS iupimglib
    EXPORT IUPTargets
    INCLUDES DESTINATION ${IUP_INSTALL_INCLUDEDIR}
    ARCHIVE DESTINATION ${IUP_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${IUP_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${IUP_INSTALL_BINDIR})
