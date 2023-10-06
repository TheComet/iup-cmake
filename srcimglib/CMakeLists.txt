project (iupimglib
    LANGUAGES C
    VERSION ${IUP_VERSION})

option (USE_IUP_IMGLIB_LARGE "" ON)

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
        $<$<BOOL:${USE_IUP_IMGLIB_LARGE}>:IUP_IMGLIB_LARGE>)

if (USE_WIN32)
    target_sources (iupimglib PRIVATE
        "iup_imglib_basewin16x16.c"
        "iup_imglib_basewin32x32.c"
        "iup_imglib_logos32x32.c"
        "iup_imglib_logos48x48.c"
        "iup_imglib_iconswin48x48.c")
endif ()
