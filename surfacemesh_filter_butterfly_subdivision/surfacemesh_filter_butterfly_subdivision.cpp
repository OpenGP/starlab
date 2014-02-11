#include "SurfaceMeshHelper.h"
#include "surfacemesh_filter_butterfly_subdivision.h"
#include "ModifiedButterflySubdivision.h"

using namespace SurfaceMesh;

void surfacemesh_filter_butterfly_subdivision::initParameters(RichParameterSet *pars){
    pars->addParam(new RichInt("Iterations", 1, "Iterations"));
}

void surfacemesh_filter_butterfly_subdivision::applyFilter(RichParameterSet* pars){
    mesh()->isVisible = false;

    int iterations = 1;
    if(pars) iterations = pars->getInt("Iterations");

    ModifiedButterfly butterfly;
    butterfly.subdivide(*mesh(), iterations);

    foreach(Vertex v, mesh()->vertices()) if(mesh()->is_isolated(v)) mesh()->remove_vertex(v);
    mesh()->garbage_collection();
    mesh()->update_face_normals();
    mesh()->update_vertex_normals();
    mesh()->updateBoundingBox();

    mesh()->isVisible = true;
}
