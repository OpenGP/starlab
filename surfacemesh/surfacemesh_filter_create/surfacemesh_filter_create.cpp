#include "surfacemesh_filter_create.h"
#include "StarlabDrawArea.h"

using namespace SurfaceMesh;

void surfacemesh_filter_create::initParameters(RichParameterSet *pars){
    pars->addParam(new RichBool("sphere", true, "Sphere", "Sphere"));
    pars->addParam(new RichBool("cube", false, "Cube", "Cube"));
    pars->addParam(new RichBool("cylinder", false, "Cylinder", "Cylinder"));
    pars->addParam(new RichBool("plane", false, "Plane", "Plane"));
    pars->addParam(new RichInt("resolution", 3, "Resolution", "Resolution"));
}

void surfacemesh_filter_create::applyFilter(RichParameterSet *pars){
    SurfaceMeshModel * m = nullptr;

    if(pars->getBool("sphere")){
        m = new SurfaceMeshModel("sphere.off", "sphere");
        makeSphere(*m, pars->getInt("resolution"));
    }
    else if(pars->getBool("cube")){
        m = new SurfaceMeshModel("cube.off", "cube");
        makeCube(*m, pars->getInt("resolution"));
    }
    else if(pars->getBool("cylinder")){
        m = new SurfaceMeshModel("cylinder.off", "cylinder");
        makeCylinder(*m, pars->getInt("resolution"));
    }
    else if(pars->getBool("plane")){
        m = new SurfaceMeshModel("plane.off", "plane");
        makePlane(*m, pars->getInt("resolution"));
    }

    if(m){
        m->updateBoundingBox();
        m->update_face_normals();
        m->update_vertex_normals();

        document()->addModel(m);
        document()->setSelectedModel(m);
    }
}
