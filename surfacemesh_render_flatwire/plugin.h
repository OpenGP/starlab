#pragma once
#include "SurfaceMeshPlugins.h"

class plugin : public SurfaceMeshRenderPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "surfacemesh_render_flatwire.plugin.starlab")
    Q_INTERFACES(RenderPlugin)

    QString name() { return "Flat Wire"; }
    QIcon icon(){ return QIcon(":/icons/flatwire.png"); }
    Renderer* instance();
};
