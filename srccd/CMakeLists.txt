project (iupcd
    LANGUAGES C
    VERSION ${IUP_VERSION})

add_library (iupcd
    "iup_cd.c"
    "iup_cdutil.c"
    "iup_draw_cd.c")
target_include_directories (iupcd
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../src>)
target_compile_definitions (iupcd
    PRIVATE
        CD_NO_OLD_INTERFACE)

target_link_libraries (iupcd
    PRIVATE
        freetype
        ftgl
    PUBLIC
        cd)
