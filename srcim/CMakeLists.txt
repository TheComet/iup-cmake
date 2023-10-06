project (iupim
    VERSION ${IUP_VERSION}
    LANGUAGES C)

add_library (iupim
    "iup_im.c")
target_compile_definitions (iupim
    PRIVATE
        $<$<CONFIG:Debug>:IUP_ASSERT>)
target_link_libraries (iupim
    PRIVATE
        iup
        im)
