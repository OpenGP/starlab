#pragma once
#include "SurfaceMeshPlugins.h"
class surfacemesh_filter_simplification : public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "surfacemesh_filter_simplification.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Simplification | Quadric based"; }
    QString description() { return "Mesh decimation using error quadrics."; }

    void initParameters(RichParameterSet *pars);
    void applyFilter(RichParameterSet*);
};
