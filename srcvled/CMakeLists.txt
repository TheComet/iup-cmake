project (iupvled
    VERSION ${IUP_VERSION}
    LANGUAGES C)

add_executable (iupvled
    "iup_vled.c"
    "iup_vled_imgs.c")
target_link_libraries (iupvled PRIVATE
    iup
    iupimglib
    iupcd)

if (IUP_IM)
    target_link_libraries (iupvled PRIVATE
        iupim
        im
        cd)
    target_compile_definitions (iupvled PRIVATE
        USE_IM)
endif ()

