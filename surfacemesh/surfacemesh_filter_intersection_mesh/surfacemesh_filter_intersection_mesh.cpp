#include "surfacemesh_filter_intersection_mesh.h"
#include "StarlabDrawArea.h"

using namespace SurfaceMesh;

#define PQP_SUPPORT_SURFACEMESH
#include "PQPLib.h"

void surfacemesh_filter_intersection_mesh::applyFilter(RichParameterSet*)
{
    auto models = document()->models();

    if( models.size() < 2 ){
        QMessageBox::information(mainWindow(), "Mesh intersection", "Please load at least two models.");
        return;
    }

    QVector<SurfaceMeshModel*> meshes;
    meshes << (SurfaceMeshModel*)models.at(0);
    meshes << (SurfaceMeshModel*)models.at(1);

    PQP::Manager manager( meshes.size() );
    manager.addModel( makeModelPQP( meshes[0] ) );
    manager.addModel( makeModelPQP( meshes[1] ) );
    auto isects_pairs = manager.testIntersection();

    // Prepare for better visualization
    for(auto & m : meshes) drawArea()->setRenderer(m, "Wireframe");
    drawArea()->clear();

    // Check if intersected
    if(isects_pairs.size() > 1)
    {
        for(auto & isect : isects_pairs)
        {
            for(int i = 0; i < 2; i++){
                QVector<Vector3> vp;
                for(auto v : meshes[i]->vertices( Face( i == 0 ? isect.id1 : isect.id2 ) ))
                    vp << meshes[i]->vertex_coordinates()[v];
                drawArea()->drawTriangle(vp[0], vp[1], vp[2], i == 0 ? Qt::blue : Qt::green);
            }

            Vector3 p(isect.p[0],isect.p[1],isect.p[2]);
            Vector3 q(isect.q[0],isect.q[1],isect.q[2]);
        }
    }
    else
    {
        PQP::IntersectResult isect = isects_pairs.front();

        Vector3 p(isect.p[0],isect.p[1],isect.p[2]);
        Vector3 q(isect.q[0],isect.q[1],isect.q[2]);

        drawArea()->drawSegment(p, q, 3);

        for(int i = 0; i < 2; i++){
            QVector<Vector3> vp;
            for(auto v : meshes[i]->vertices( Face( i == 0 ? isect.id1 : isect.id2 ) ))
                vp << meshes[i]->vertex_coordinates()[v];
            drawArea()->drawTriangle(vp[0], vp[1], vp[2], i == 0 ? Qt::blue : Qt::green);
        }
    }

    qDebug() << QString("Distance between meshes = %1").arg(isects_pairs.front().distance);
}

