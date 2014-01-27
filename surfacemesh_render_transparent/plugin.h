#pragma once
#include "SurfaceMeshPlugins.h"

class plugin : public SurfaceMeshRenderPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "surfacemesh_render_transparent.plugin.starlab")
    Q_INTERFACES(RenderPlugin)
  
    QString name() { return "Transparent"; }
    QIcon icon(){ return QIcon(":/icons/transparent_mesh.png"); }
    Renderer* instance();
};
