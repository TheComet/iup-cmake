project ("iup"
    LANGUAGES C
    VERSION ${IUP_VERSION})

# Windows 7 minimum
set (WIN32VER "0x0601")
# OLD Windows XP Minimum
# set (WIN32VER "0x0501")

if (IUP_TOOLKIT MATCHES "Haiku")
    enable_language (CXX)
endif ()

add_library (iup
    "iup_array.h"
    "iup_assert.h"
    "iup_attrib.h"
    "iup_box.h"
    "iup_button.h"
    "iup_canvas.h"
    "iup_childtree.h"
    "iup_class.h"
    "iup_classbase.h"
    "iup_colorhsi.h"
    "iup_dialog.h"
    "iup_dlglist.h"
    "iup_draw.h"
    "iup_drv.h"
    "iup_drvdraw.h"
    "iup_drvfont.h"
    "iup_drvinfo.h"
    "iup_flatscrollbar.h"
    "iup_flattree.h"
    "iup_focus.h"
    "iup_frame.h"
    "iup_func.h"
    "iup_globalattrib.h"
    "iup_image.h"
    "iup_key.h"
    "iup_label.h"
    "iup_layout.h"
    "iup_ledlex.h"
    "iup_linefile.h"
    "iup_list.h"
    "iup_lng_czech_utf8.h"
    "iup_lng_english.h"
    "iup_lng_portuguese.h"
    "iup_lng_portuguese_utf8.h"
    "iup_lng_russian_utf8.h"
    "iup_lng_spanish.h"
    "iup_lng_spanish_utf8.h"
    "iup_loop.h"
    "iup_mask.h"
    "iup_maskmatch.h"
    "iup_maskparse.h"
    "iup_menu.h"
    "iup_names.h"
    "iup_normalizer.h"
    "iup_object.h"
    "iup_predialogs.h"
    "iup_progressbar.h"
    "iup_register.h"
    "iup_stdcontrols.h"
    "iup_str.h"
    "iup_strmessage.h"
    "iup_table.h"
    "iup_tabs.h"
    "iup_text.h"
    "iup_timer.h"
    "iup_toggle.h"
    "iup_tree.h"
    "iup_val.h"

    "iup_array.c"
    "iup_callback.c"
    "iup_dlglist.c"
    "iup_attrib.c"
    "iup_focus.c"
    "iup_font.c"
    "iup_globalattrib.c"
    "iup_object.c"
    "iup_key.c"
    "iup_layout.c"
    "iup_ledlex.c"
    "iup_names.c"
    "iup_ledparse.c"
    "iup_predialogs.c"
    "iup_register.c"
    "iup_scanf.c"
    "iup_show.c"
    "iup_str.c"
    "iup_func.c"
    "iup_childtree.c"
    "iup.c"
    "iup_classattrib.c"
    "iup_dialog.c"
    "iup_assert.c"
    "iup_messagedlg.c"
    "iup_timer.c"
    "iup_image.c"
    "iup_label.c"
    "iup_fill.c"
    "iup_zbox.c"
    "iup_colordlg.c"
    "iup_fontdlg.c"
    "iup_filedlg.c"
    "iup_strmessage.c"
    "iup_menu.c"
    "iup_frame.c"
    "iup_user.c"
    "iup_button.c"
    "iup_radio.c"
    "iup_toggle.c"
    "iup_progressbar.c"
    "iup_text.c"
    "iup_val.c"
    "iup_box.c"
    "iup_hbox.c"
    "iup_vbox.c"
    "iup_cbox.c"
    "iup_class.c"
    "iup_classbase.c"
    "iup_maskmatch.c"
    "iup_mask.c"
    "iup_maskparse.c"
    "iup_tabs.c"
    "iup_spin.c"
    "iup_list.c"
    "iup_getparam.c"
    "iup_link.c"
    "iup_sbox.c"
    "iup_scrollbox.c"
    "iup_normalizer.c"
    "iup_tree.c"
    "iup_split.c"
    "iup_layoutdlg.c"
    "iup_recplay.c"
    "iup_progressdlg.c"
    "iup_expander.c"
    "iup_open.c"
    "iup_table.c"
    "iup_canvas.c"
    "iup_gridbox.c"
    "iup_detachbox.c"
    "iup_backgroundbox.c"
    "iup_linefile.c"
    "iup_config.c"
    "iup_flatbutton.c"
    "iup_animatedlabel.c"
    "iup_draw.c"
    "iup_flatframe.c"
    "iup_flattabs.c"
    "iup_flatscrollbar.c"
    "iup_flatscrollbox.c"
    "iup_gauge.c"
    "iup_dial.c"
    "iup_colorbar.c"
    "iup_colorbrowser.c"
    "iup_colorhsi.c"
    "iup_flatlabel.c"
    "iup_dropbutton.c"
    "iup_flattoggle.c"
    "iup_flatseparator.c"
    "iup_space.c"
    "iup_multibox.c"
    "iup_flatlist.c"
    "iup_globalsdlg.c"
    "iup_flatval.c"
    "iup_loop.c"
    "iup_thread.c"
    "iup_classinfo.c"
    "iup_elempropdlg.c"
    "iup_export.c"
    "iup_flattree.c")
target_include_directories (iup
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)
target_compile_definitions (iup
    PRIVATE
        $<$<PLATFORM_ID:Windows>:WIN32>
        $<$<CONFIG:Debug>:IUP_ASSERT>)

if (IUP_TOOLKIT MATCHES "Haiku")
    target_sources (iup PRIVATE
        "haiku/iuphaiku_button.cpp"
        "haiku/iuphaiku_canvas.cpp"
        "haiku/iuphaiku_clipboard.c"
        "haiku/iuphaiku_common.cpp"
        "haiku/iuphaiku_dialog.cpp"
        "haiku/iuphaiku_dragdrop.c"
        "haiku/iuphaiku_draw.c"
        "haiku/iuphaiku_filedlg.cpp"
        "haiku/iuphaiku_focus.cpp"
        "haiku/iuphaiku_font.cpp"
        "haiku/iuphaiku_fontdlg.c"
        "haiku/iuphaiku_frame.cpp"
        "haiku/iuphaiku_globalattrib.c"
        "haiku/iuphaiku_help.c"
        "haiku/iuphaiku_image.cpp"
        "haiku/iuphaiku_info.cpp"
        "haiku/iuphaiku_label.cpp"
        "haiku/iuphaiku_list.cpp"
        "haiku/iuphaiku_loop.cpp"
        "haiku/iuphaiku_menu.cpp"
        "haiku/iuphaiku_messagedlg.c"
        "haiku/iuphaiku_open.cpp"
        "haiku/iuphaiku_progressbar.cpp"
        "haiku/iuphaiku_tabs.cpp"
        "haiku/iuphaiku_text.cpp"
        "haiku/iuphaiku_timer.cpp"
        "haiku/iuphaiku_tips.c"
        "haiku/iuphaiku_toggle.cpp"
        "haiku/iuphaiku_tree.cpp"
        "haiku/iuphaiku_val.cpp"
        "gtk/iupgtk_draw_cairo.c"
        "gtk/iupgtk_key.c")
elseif (IUP_TOOLKIT MATCHES "GTK2" OR IUP_TOOLKIT MATCHES "GTK3")
    target_sources (iup PRIVATE
        "gtk/iupgtk_focus.c"
        "gtk/iupgtk_clipboard.c"
        "gtk/iupgtk_val.c"
        "gtk/iupgtk_globalattrib.c"
        "gtk/iupgtk_key.c"
        "gtk/iupgtk_tips.c"
        "gtk/iupgtk_loop.c"
        "gtk/iupgtk_open.c"
        "gtk/iupgtk_messagedlg.c"
        "gtk/iupgtk_timer.c"
        "gtk/iupgtk_label.c"
        "gtk/iupgtk_font.c"
        "gtk/iupgtk_fontdlg.c"
        "gtk/iupgtk_filedlg.c"
        "gtk/iupgtk_button.c"
        "gtk/iupgtk_toggle.c"
        "gtk/iupgtk_dragdrop.c"
        "gtk/iupgtk_text.c"
        "gtk/iupgtk_frame.c"
        "gtk/iupgtk_progressbar.c"
        "gtk/iupgtk_tabs.c"
        "gtk/iupgtk_list.c"
        "gtk/iupgtk_tree.c"
        "gtk/iupgtk_canvas.c"
        "gtk/iupgtk_image.c"
        "gtk/iupgtk_dialog.c"
        "gtk/iupgtk_common.c"
        "gtk/iupgtk_str.c"
        "gtk/iupgtk_menu.c"
        "gtk/iupgtk_calendar.c"
        "iup_datepick.c")
    target_compile_definitions (iup
        PRIVATE
            GTK_DISABLE_DEPRECATED)

    if (IUP_TOOLKIT MATCHES "GTK3")
        find_package (GTK3 REQUIRED)
        target_include_directories (iup PRIVATE ${GTK3_INCLUDE_DIRS})
        target_link_libraries (iup PRIVATE ${GTK3_LIBRARIES})
        target_compile_definitions (iup
            PRIVATE
                GDK_DISABLE_DEPRECATED
                GSEAL_ENABLE)
        if (USE_NEW_DRAW)
            target_sources (iup PRIVATE
                "gtk/iupgtk_draw_cairo.c")
        endif ()
    else ()
        find_package (GTK2 REQUIRED)
        target_link_libraries (iup PRIVATE GTK2::gtk)

        if (IUP_NEW_DRAW)
            target_sources (iup PRIVATE
                "gtk/iupgtk_draw_cairo.c")
        else ()
            target_sources (iup PRIVATE
                "gtk/iupgtk_draw_gdk.c")
        endif ()
    endif ()

    if (CMAKE_PLATFORM_ID MATCHES Windows)
        target_compile_definitions (iup PRIVATE
            _WIN32_WINNT=${WIN32VER} _WIN32_IE=0x0600 WINVER=${WIN32VER} NOTREEVIEW)
    else ()
        if (CMAKE_PLATFORM_ID MATCHES Darwin)
            target_sources (iup PRIVATE
                "gtk/iupmac_help.c"
                "gtk/iupmac_info.c")
            target_compile_definitions (iup PRIVATE GTK_MAC)
        else ()
            if (GDK_NULL)
                #to completely remove X11 calls
                # works only for GTK3
                target_compile_definitions (iup PRIVATE GDK_NULL)
            else ()
                find_package (X11 REQUIRED)
                target_link_libraries (iup
                    PRIVATE
                        X11::X11
                        X11::Xt
                        m)
            endif ()
            target_sources (iup PRIVATE
                "gtk/iupgtk_help.c"
                "mot/iupunix_info.c")

            if (IUP_TOOLKIT MATCHES "GTK3")
                target_sources (iup PRIVATE
                    "gtk/iupgtk_info.c")
            else ()
                # Because of iupdrvGetScreenSize limitation
                if (GDK_NULL)
                    target_sources (iup PRIVATE
                        "gtk/iupgtk_info.c")
                else ()
                    target_sources (iup PRIVATE
                        "mot/iupx11_info.c")
                endif ()
            endif ()
        endif ()
    endif ()
elseif (IUP_TOOLKIT MATCHES "Motif")
    find_package (Motif REQUIRED)
    find_package (X11 REQUIRED)

    target_sources (iup PRIVATE
        "mot/iupmot_common.c"
        "mot/iupmot_color.c"
        "mot/iupmot_focus.c"
        "mot/iupmot_font.c"
        "mot/iupmot_key.c"
        "mot/iupmot_loop.c"
        "mot/iupmot_open.c"
        "mot/iupmot_tips.c"
        "mot/iupmot_globalattrib.c"
        "mot/iupmot_dialog.c"
        "mot/iupmot_messagedlg.c"
        "mot/iupmot_draw_x11.c"
        "mot/iupmot_timer.c"
        "mot/iupmot_image.c"
        "mot/iupmot_label.c"
        "mot/iupmot_canvas.c"
        "mot/iupmot_fontdlg.c"
        "mot/iupmot_filedlg.c"
        "mot/iupmot_frame.c"
        "mot/iupmot_button.c"
        "mot/iupmot_toggle.c"
        "mot/iupmot_progressbar.c"
        "mot/iupmot_clipboard.c"
        "mot/iupmot_text.c"
        "mot/iupmot_val.c"
        "mot/iupmot_tabs.c"
        "mot/iupmot_menu.c"
        "mot/iupmot_list.c"
        "mot/iupmot_tree.c"
        "mot/iupmot_dragdrop.c"
        "mot/iupmot_str.c"
        "mot/iupmot_calendar.c"
        "iup_datepick.c"

        "mot/iupunix_help.c"
        "mot/iupunix_info.c"
        "mot/iupx11_info.c")
   target_include_directories (iup
        PRIVATE
            $<BUILD_INTERFACE:${MOTIF_INCLUDE_DIR}>
            $<BUILD_INTERFACE:${X11_INCLUDE_DIR}>)
    target_link_libraries (iup
        PRIVATE
            ${MOTIF_LIBRARIES}
            X11::X11
            X11::Xt
            m)
elseif (IUP_TOOLKIT MATCHES "Win32")
    target_sources (iup PRIVATE
        "win/iupwin_common.c"
        "win/iupwin_brush.c"
        "win/iupwin_focus.c"
        "win/iupwin_font.c"
        "win/iupwin_globalattrib.c"
        "win/iupwin_handle.c"
        "win/iupwin_key.c"
        "win/iupwin_str.c"
        "win/iupwin_loop.c"
        "win/iupwin_open.c"
        "win/iupwin_tips.c"
        "win/iupwin_info.c"
        "win/iupwin_dialog.c"
        "win/iupwin_messagedlg.c"
        "win/iupwin_timer.c"
        "win/iupwin_image.c"
        "win/iupwin_label.c"
        "win/iupwin_canvas.c"
        "win/iupwin_frame.c"
        "win/iupwin_fontdlg.c"
        "win/iupwin_filedlg.c"
        "win/iupwin_dragdrop.c"
        "win/iupwin_button.c"
        "win/iupwin_draw.c"
        "win/iupwin_toggle.c"
        "win/iupwin_clipboard.c"
        "win/iupwin_progressbar.c"
        "win/iupwin_text.c"
        "win/iupwin_val.c"
        "win/iupwin_touch.c"
        "win/iupwin_tabs.c"
        "win/iupwin_menu.c"
        "win/iupwin_list.c"
        "win/iupwin_tree.c"
        "win/iupwin_calendar.c"
        "win/iupwin_datepick.c"

        "win/iupwindows_main.c"
        "win/iupwindows_help.c"
        "win/iupwindows_info.c")
    target_compile_definitions (iup
        PRIVATE
            _WIN32_WINNT=${WIN32VER} _WIN32_IE=0x0900 WINVER=${WIN32VER} NOTREEVIEW UNICODE)
    target_link_libraries (iup
        PRIVATE
            Comctl32)

    if (IUP_NEW_DRAW)
        target_sources (iup PRIVATE
            "win/iupwin_image_wdl.c"
            "win/iupwin_draw_wdl.c"
            "win/iupwin_draw_gdi.c"
            "win/wdl/backend-d2d.c"
            "win/wdl/backend-dwrite.c"
            "win/wdl/backend-gdix.c"
            "win/wdl/backend-wic.c"
            "win/wdl/bitblt.c"
            "win/wdl/brush.c"
            "win/wdl/cachedimage.c"
            "win/wdl/canvas.c"
            "win/wdl/draw.c"
            "win/wdl/fill.c"
            "win/wdl/font.c"
            "win/wdl/image.c"
            "win/wdl/init.c"
            "win/wdl/memstream.c"
            "win/wdl/misc.c"
            "win/wdl/path.c"
            "win/wdl/string.c"
            "win/wdl/strokestyle.c")
        target_include_directories (iup
            PRIVATE
                $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/win/wdl>)
        target_compile_definitions (iup
            PRIVATE
                COBJMACROS _UNICODE USE_NEW_DRAW)
    else ()
        target_sources (iup PRIVATE
            "win/iupwin_draw_gdi.c")
    endif ()
endif ()
