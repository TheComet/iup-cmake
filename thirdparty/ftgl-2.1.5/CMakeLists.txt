cmake_minimum_required (VERSION 3.21)

project (ftgl
    VERSION 2.1.5
    LANGUAGES CXX)

option (FTGL_PIC "Compile with position independent code" OFF)

add_library (ftgl
    "src/FTBuffer.cpp"
    "src/FTCharmap.cpp"
    "src/FTContour.cpp"
    "src/FTFace.cpp"
    "src/FTGlyphContainer.cpp"
    "src/FTLibrary.cpp"
    "src/FTPoint.cpp"
    "src/FTSize.cpp"
    "src/FTVectoriser.cpp"

    # ftglyph sources
    "src/FTGlyph/FTGlyph.cpp"
    "src/FTGlyph/FTGlyphGlue.cpp"
    "src/FTGlyph/FTBitmapGlyph.cpp"
    "src/FTGlyph/FTBufferGlyph.cpp"
    "src/FTGlyph/FTExtrudeGlyph.cpp"
    "src/FTGlyph/FTOutlineGlyph.cpp"
    "src/FTGlyph/FTPixmapGlyph.cpp"
    "src/FTGlyph/FTPolygonGlyph.cpp"
    "src/FTGlyph/FTTextureGlyph.cpp"

    # ftfont sources
    "src/FTFont/FTFont.cpp"
    "src/FTFont/FTFontGlue.cpp"
    "src/FTFont/FTBitmapFont.cpp"
    "src/FTFont/FTBufferFont.cpp"
    "src/FTFont/FTExtrudeFont.cpp"
    "src/FTFont/FTOutlineFont.cpp"
    "src/FTFont/FTPixmapFont.cpp"
    "src/FTFont/FTPolygonFont.cpp"
    "src/FTFont/FTTextureFont.cpp"

    # ftlayout sources
    "src/FTLayout/FTLayout.cpp"
    "src/FTLayout/FTLayoutGlue.cpp"
    "src/FTLayout/FTSimpleLayout.cpp")

target_include_directories (ftgl
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
target_compile_definitions (ftgl
    PRIVATE
        $<$<CONFIG:Debug>:IUP_ASSERT>)
target_link_libraries (ftgl
    PRIVATE
        freetype)
set_target_properties (ftgl PROPERTIES
    MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>
    POSITION_INDEPENDENT_CODE ${FTGL_PIC})
