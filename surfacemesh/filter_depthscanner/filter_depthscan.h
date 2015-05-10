#pragma once
#include "FilterPlugin.h"
class filter_depthscan : public FilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "filter_depthscan.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Virtual Depth Scan"; }
    QString description() { return "Take a depth screenshot"; }
    bool isApplicable(StarlabModel*){ return true; }
    void applyFilter(RichParameterSet*);
    void initParameters(RichParameterSet*);
};
