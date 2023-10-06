project (iupcd
    LANGUAGES C
    VERSION ${IUP_VERSION})

add_library (iupcd
    "iup_cd.c"
    "iup_cdutil.c"
    "iup_draw_cd.c")
target_include_directories (iupcd
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_compile_definitions (iupcd
    PRIVATE
        CD_NO_OLD_INTERFACE)

target_link_libraries (iupcd
    PRIVATE
        iup
        freetype
        ftgl
    PUBLIC
        cd)
