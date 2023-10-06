project ("03_okbutton"
    LANGUAGES C)
    
add_executable (03_okbutton
    "main.c")
target_link_libraries (03_okbutton
    PRIVATE
        iup)
