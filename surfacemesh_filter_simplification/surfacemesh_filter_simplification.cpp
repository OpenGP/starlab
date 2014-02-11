#include "surfacemesh_filter_simplification.h"
#include "decimater.h"

using namespace SurfaceMesh;

void surfacemesh_filter_simplification::initParameters(RichParameterSet *pars){
    pars->addParam(new RichFloat("Percentage", 0.9f, "Percentage"));

    /// We are interested in seeing the mesh edges
    if(drawArea()) drawArea()->setRenderer(mesh(),"Flat Wire");
}

void surfacemesh_filter_simplification::applyFilter(RichParameterSet* pars){
    mesh()->isVisible = false;

    double percentage = 0.9;
    if(pars) percentage = pars->getFloat("Percentage");

    int numFacesBefore = mesh()->n_faces();

    Decimater::simplify(mesh(), percentage);

    foreach(Vertex v, mesh()->vertices()) if(mesh()->is_isolated(v)) mesh()->remove_vertex(v);
    mesh()->garbage_collection();
    mesh()->update_face_normals();
    mesh()->update_vertex_normals();
    mesh()->updateBoundingBox();

    int numFacesAfter = mesh()->n_faces();

    // Log information
    QString report = QString("Faces before (%1), after(%2)").arg(numFacesBefore).arg(numFacesAfter);
    qDebug() << report;
    mainWindow()->setStatusBarMessage(report);

    mesh()->isVisible = true;
}

