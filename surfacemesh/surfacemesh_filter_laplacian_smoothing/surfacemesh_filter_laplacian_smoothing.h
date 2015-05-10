#pragma once
#include "SurfaceMeshPlugins.h"

class surfacemesh_filter_laplacian_smoothing : public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "surfacemesh_filter_laplacian_smoothing.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Smoothing | Laplacian"; }
    QString description() { return "A new position for each vertex is chosen based on local information."; }

    void initParameters(RichParameterSet *pars);
    void applyFilter(RichParameterSet*);
};
