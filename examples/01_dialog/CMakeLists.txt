project ("01_dialog"
    LANGUAGES C)
    
add_executable (01_dialog
    "main.c")
target_link_libraries (01_dialog
    PRIVATE
        iup)
