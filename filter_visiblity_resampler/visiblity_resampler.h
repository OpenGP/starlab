#pragma once

#include "SurfaceMeshPlugins.h"
#include "SurfaceMeshHelper.h"
#include "RichParameterSet.h"

class visiblity_resampler : public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "visiblity_resampler.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Visibility Resampler"; }
    QString description() { return "Resample a mesh by considering visible points only"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);

	QMap<QString,QVariant> my;
};
