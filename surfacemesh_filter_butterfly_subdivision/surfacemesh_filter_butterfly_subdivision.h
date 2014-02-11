#pragma once
#include "SurfaceMeshPlugins.h"

class surfacemesh_filter_butterfly_subdivision : public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "surfacemesh_filter_butterfly_subdivision.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Subdivision | Modified Butterfly"; }
    QString description() { return "A new position for each vertex is chosen based on local information."; }

    void initParameters(RichParameterSet *pars);
    void applyFilter(RichParameterSet*);
};
