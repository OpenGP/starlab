#include "SurfaceMeshHelper.h"
#include "surfacemesh_filter_laplacian_smoothing.h"

using namespace SurfaceMesh;

void surfacemesh_filter_laplacian_smoothing::initParameters(RichParameterSet *pars){
    pars->addParam(new RichInt("Iterations", 1, "Iterations"));
}

void surfacemesh_filter_laplacian_smoothing::applyFilter(RichParameterSet* pars){
    int iterations = 1;
    if(pars) iterations = pars->getInt("Iterations");
    SurfaceMeshHelper(mesh()).smoothVertexProperty<Vector3>(VPOINT, iterations);
}
