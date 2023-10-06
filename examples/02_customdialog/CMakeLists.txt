project ("02_customdialog"
    LANGUAGES C)
    
add_executable (02_customdialog
    "main.c")
target_link_libraries (02_customdialog
    PRIVATE
        iup)
