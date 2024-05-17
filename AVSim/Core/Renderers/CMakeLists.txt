set(RENDERERS_SOURCES Renderer.cpp GuiBackInterconnection.cpp $<$<BOOL:${USE_MAGNUM_GUI}>:GuiApplication.cpp> PARENT_SCOPE)
set(RENDERERS_HEADERS Renderer.h GuiBackInterconnection.h GuiApplication.h PARENT_SCOPE)

