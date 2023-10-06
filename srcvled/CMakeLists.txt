project (iupvled
    VERSION ${IUP_VERSION}
    LANGUAGES C)

add_executable (iupvled
    "iup_vled.c"
    "iup_vled_imgs.c")
target_link_libraries (iupvled PRIVATE
    iup
    iupim
    iupimglib
    iupcd)
