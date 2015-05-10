#pragma once
#include "SurfaceMeshPlugins.h"
class surfacemesh_filter_intersection_mesh : public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "surfacemesh_filter_intersection_mesh.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Intersection of meshes"; }
    QString description() { return "Find intersection between two meshes"; }
    void applyFilter(RichParameterSet*);
};
