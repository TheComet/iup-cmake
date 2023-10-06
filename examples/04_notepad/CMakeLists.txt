project ("04_notepad"
    LANGUAGES C)
    
add_executable (04_notepad
    "main.c")
target_link_libraries (04_notepad
    PRIVATE
        iup
        iupimglib)
