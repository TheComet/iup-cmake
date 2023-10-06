project (iupcd
    LANGUAGES C
    VERSION ${IUP_VERSION})

add_library (iupcd ${IUP_LIB_TYPE}
    "iup_cd.c"
    "iup_cdutil.c"
    "iup_draw_cd.c")
target_include_directories (iupcd
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../src>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_compile_definitions (iupcd
    PRIVATE
        CD_NO_OLD_INTERFACE)
target_link_libraries (iupcd
    PRIVATE
        freetype
        ftgl
        cd
        iup)
set_target_properties (iupcd
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
    target_sources (iupcd PRIVATE "iupcd.def")
endif ()

add_library (IUP::cd ALIAS iupcd)

install (
    TARGETS iupcd
    EXPORT IUPTargets
    INCLUDES DESTINATION ${IUP_INSTALL_INCLUDEDIR}
    ARCHIVE DESTINATION ${IUP_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${IUP_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${IUP_INSTALL_BINDIR})
