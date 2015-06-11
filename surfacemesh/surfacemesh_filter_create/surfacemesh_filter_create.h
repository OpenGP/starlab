#pragma once
#include "SurfaceMeshPlugins.h"
class surfacemesh_filter_create : public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "surfacemesh_filter_create.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Create new mesh"; }
    QString description() { return "Creates a new mesh"; }
    bool isApplicable(Starlab::Model*) { return true; }
    void initParameters(RichParameterSet *pars);
    void applyFilter(RichParameterSet *pars);
};

/// Helper functions:
static inline void makeSphere(SurfaceMesh::SurfaceMeshModel & geometry, int recursionLevel)
{
    // Cache
    std::map<size_t, int> middlePointIndexCache;

    // Vertex index
    int index = 0;

    auto getMiddlePoint = [&index,&geometry,&middlePointIndexCache](int p1, int p2)
    {
        // first check if we have it already
        bool firstIsSmaller = p1 < p2;
        size_t smallerIndex = firstIsSmaller ? p1 : p2;
        size_t greaterIndex = firstIsSmaller ? p2 : p1;
        size_t key = (smallerIndex << 32) + greaterIndex;

        if(middlePointIndexCache.find(key) != middlePointIndexCache.end())
            return middlePointIndexCache[key];

        // not in cache, calculate it
        Vector3 v1 = geometry.vertex_property<Vector3>("v:point")[Surface_mesh::Vertex(p1)];
        Vector3 v2 = geometry.vertex_property<Vector3>("v:point")[Surface_mesh::Vertex(p2)];
        Vector3 middle = ( v1 + v2 ) * 0.5;

        // add vertex makes sure point is on unit sphere
        geometry.add_vertex( middle.normalized() );
        int i = index++;

        // store it, return index
        middlePointIndexCache[key] = i;
        return i;
    };

    auto addVertex = [&index,&geometry](const Vector3 &p){
        geometry.add_vertex( p.normalized() );
        return index++;
    };

    struct TriangleIndices{
        int v1,v2,v3;
        TriangleIndices(int v1, int v2, int v3) : v1(v1), v2(v2), v3(v3){}
    };

    // create 12 vertices of a icosahedron
    auto t = (1.0 + std::sqrt(5.0)) / 2.0;

    addVertex(Vector3(-1,  t,  0));
    addVertex(Vector3( 1,  t,  0));
    addVertex(Vector3(-1, -t,  0));
    addVertex(Vector3( 1, -t,  0));

    addVertex(Vector3( 0, -1,  t));
    addVertex(Vector3( 0,  1,  t));
    addVertex(Vector3( 0, -1, -t));
    addVertex(Vector3( 0,  1, -t));

    addVertex(Vector3( t,  0, -1));
    addVertex(Vector3( t,  0,  1));
    addVertex(Vector3(-t,  0, -1));
    addVertex(Vector3(-t,  0,  1));

    // create 20 triangles of the icosahedron
    std::list<TriangleIndices> faces;

    // 5 faces around point 0
    faces.push_back(TriangleIndices(0, 11, 5));
    faces.push_back(TriangleIndices(0, 5, 1));
    faces.push_back(TriangleIndices(0, 1, 7));
    faces.push_back(TriangleIndices(0, 7, 10));
    faces.push_back(TriangleIndices(0, 10, 11));

    // 5 adjacent faces
    faces.push_back(TriangleIndices(1, 5, 9));
    faces.push_back(TriangleIndices(5, 11, 4));
    faces.push_back(TriangleIndices(11, 10, 2));
    faces.push_back(TriangleIndices(10, 7, 6));
    faces.push_back(TriangleIndices(7, 1, 8));

    // 5 faces around point 3
    faces.push_back(TriangleIndices(3, 9, 4));
    faces.push_back(TriangleIndices(3, 4, 2));
    faces.push_back(TriangleIndices(3, 2, 6));
    faces.push_back(TriangleIndices(3, 6, 8));
    faces.push_back(TriangleIndices(3, 8, 9));

    // 5 adjacent faces
    faces.push_back(TriangleIndices(4, 9, 5));
    faces.push_back(TriangleIndices(2, 4, 11));
    faces.push_back(TriangleIndices(6, 2, 10));
    faces.push_back(TriangleIndices(8, 6, 7));
    faces.push_back(TriangleIndices(9, 8, 1));

    // refine triangles
    for (int i = 0; i < recursionLevel; i++)
    {
        std::list<TriangleIndices> faces2;
        for (auto tri : faces)
        {
            // replace triangle by 4 triangles
            int a = getMiddlePoint(tri.v1, tri.v2);
            int b = getMiddlePoint(tri.v2, tri.v3);
            int c = getMiddlePoint(tri.v3, tri.v1);

            faces2.push_back(TriangleIndices(tri.v1, a, c));
            faces2.push_back(TriangleIndices(tri.v2, b, a));
            faces2.push_back(TriangleIndices(tri.v3, c, b));
            faces2.push_back(TriangleIndices(a, b, c));
        }
        faces = faces2;
    }

    // done, now add triangles to mesh
    for (auto tri : faces)
    {
        geometry.add_triangle( Surface_mesh::Vertex(tri.v1),
                               Surface_mesh::Vertex(tri.v2),
                               Surface_mesh::Vertex(tri.v3) );
    }
}

// Based on trimesh2
static inline void makePlane(SurfaceMesh::SurfaceMeshModel & mesh, int resolution)
{
    auto mkpoint = [&](SurfaceMesh::SurfaceMeshModel & mesh, float x, float y, float z){
        mesh.add_vertex(SurfaceMesh::Vector3(x,y,z));
    };

    auto mkquad = [&](SurfaceMesh::SurfaceMeshModel & mesh, int _ll, int _lr, int _ul, int _ur){
        SurfaceMesh::Vertex ll(_ll), lr(_lr), ul(_ul), ur(_ur);
        mesh.add_triangle(ll, lr, ur);
        mesh.add_triangle(ll, ur, ul);
    };

    int tess_x = resolution, tess_y = resolution;

    for (int j = 0; j < tess_y+1; j++) {
        float y = -1.0f + 2.0f * j / tess_y;
        for (int i = 0; i < tess_x+1; i++) {
            float x = -1.0f + 2.0f * i / tess_x;
            mkpoint(mesh, x, y, 0);
        }
    }

    for (int j = 0; j < tess_y; j++) {
        for (int i = 0; i < tess_x; i++) {
            int ind = i + j * (tess_x+1);
            mkquad(mesh, ind, ind+1, ind+tess_x+1, ind+tess_x+2);
        }
    }
}

static inline void makeCube(SurfaceMesh::SurfaceMeshModel & mesh, int resolution)
{
    auto mkpoint = [&](SurfaceMesh::SurfaceMeshModel & mesh, float x, float y, float z){
        mesh.add_vertex(SurfaceMesh::Vector3(x,y,z));
    };

    auto mkquad = [&](SurfaceMesh::SurfaceMeshModel & mesh, int _ll, int _lr, int _ul, int _ur){
        SurfaceMesh::Vertex ll(_ll), lr(_lr), ul(_ul), ur(_ur);
        mesh.add_triangle(ll, lr, ur);
        mesh.add_triangle(ll, ur, ul);
    };

    int tess = resolution;

    for (int j = 0; j < tess+1; j++) {
        float y = 1.0f - 2.0f * j / tess;
        for (int i = 0; i < tess+1; i++) {
            float x = 1.0f - 2.0f * i / tess;
            mkpoint(mesh, x, y, -1);
        }
    }
    for (int j = 1; j < tess; j++) {
        float z = -1.0f + 2.0f * j / tess;
        for (int i = 0; i < tess; i++) {
            float x = -1.0f + 2.0f * i / tess;
            mkpoint(mesh, x, -1, z);
        }
        for (int i = 0; i < tess; i++) {
            float y = -1.0f + 2.0f * i / tess;
            mkpoint(mesh, 1, y, z);
        }
        for (int i = 0; i < tess; i++) {
            float x = 1.0f - 2.0f * i / tess;
            mkpoint(mesh, x, 1, z);
        }
        for (int i = 0; i < tess; i++) {
            float y = 1.0f - 2.0f * i / tess;
            mkpoint(mesh, -1, y, z);
        }
    }
    for (int j = 0; j < tess+1; j++) {
        float y = -1.0f + 2.0f * j / tess;
        for (int i = 0; i < tess+1; i++) {
            float x = -1.0f + 2.0f * i / tess;
            mkpoint(mesh, x, y, 1);
        }
    }

    for (int j = 0; j < tess; j++) {
        for (int i = 0; i < tess; i++) {
            int ind = i + j * (tess+1);
            mkquad(mesh, ind, ind+tess+1, ind+1, ind+tess+2);
        }
    }

    auto sqr = [&](int x){
        return x*x;
    };

    int topstart = sqr(tess+1) + 4*tess*(tess-1);
    for (int j = 0; j < tess; j++) {
        int next = sqr(tess+1) + 4*tess*(j-1);
        for (int i = 0; i < tess; i++) {
            int ll = next++;
            int lr = ll + 1;
            int ul = ll + 4*tess;
            int ur = ul + 1;
            if (j == 0) {
                ll = sqr(tess+1)-1 - i;
                lr = ll - 1;
            }
            mkquad(mesh, ll, lr, ul, ur);
        }
        for (int i = 0; i < tess; i++) {
            int ll = next++;
            int lr = ll + 1;
            int ul = ll + 4*tess;
            int ur = ul + 1;
            if (j == 0) {
                ll = tess*(tess+1) - i*(tess+1);
                lr = ll - (tess+1);
            }
            if (j == tess-1) {
                ul = topstart + tess + i*(tess+1);
                ur = ul + (tess+1);
            }
            mkquad(mesh, ll, lr, ul, ur);
        }
        for (int i = 0; i < tess; i++) {
            int ll = next++;
            int lr = ll + 1;
            int ul = ll + 4*tess;
            int ur = ul + 1;
            if (j == 0) {
                ll = i;
                lr = i + 1;
            }
            if (j == tess-1) {
                ul = topstart + sqr(tess+1)-1 - i;
                ur = ul - 1;
            }
            mkquad(mesh, ll, lr, ul, ur);
        }
        for (int i = 0; i < tess; i++) {
            int ll = next++;
            int lr = ll + 1;
            int ul = ll + 4*tess;
            int ur = ul + 1;
            if (j == 0) {
                ll = tess + i*(tess+1);
                lr = ll + (tess+1);
            }
            if (j == tess-1) {
                ul = topstart + tess*(tess+1) - i*(tess+1);
                ur = ul - (tess+1);
            }
            if (i == tess-1) {
                if (j != 0)
                    lr -= 4*tess;
                if (j != tess-1)
                    ur -= 4*tess;
            }
            mkquad(mesh, ll, lr, ul, ur);
        }
    }
    for (int j = 0; j < tess; j++) {
        for (int i = 0; i < tess; i++) {
            int ind = topstart + i + j * (tess+1);
            mkquad(mesh, ind, ind+1, ind+tess+1, ind+tess+2);
        }
    }
}

static inline void makeCylinder(SurfaceMesh::SurfaceMeshModel & mesh, int resolution)
{
    auto mkpoint = [&](SurfaceMesh::SurfaceMeshModel & mesh, float x, float y, float z){
        mesh.add_vertex(SurfaceMesh::Vector3(x,y,z));
    };

    auto mkface = [&](SurfaceMesh::SurfaceMeshModel & mesh, int v1, int v2, int v3){
         mesh.add_triangle(SurfaceMesh::Vertex(v1), SurfaceMesh::Vertex(v2), SurfaceMesh::Vertex(v3));
    };

    auto mkquad = [&](SurfaceMesh::SurfaceMeshModel & mesh, int _ll, int _lr, int _ul, int _ur){
        SurfaceMesh::Vertex ll(_ll), lr(_lr), ul(_ul), ur(_ur);
        mesh.add_triangle(ll, lr, ur);
        mesh.add_triangle(ll, ur, ul);
    };

    float r = 0.25;

    int tess_th = 3 * resolution;
    int tess_h = 1 * resolution;

#ifndef M_TWOPIf
# define M_TWOPIf 6.2831855f
#endif

    mkpoint(mesh, 0, 0, -1);
    for (int j = 1; j <= tess_h; j++) {
        float rr = r * j / tess_h;
        for (int i = 0; i < tess_th; i++) {
            float th = M_TWOPIf * i / tess_th;
            mkpoint(mesh, rr*cos(th), rr*sin(th), -1);
        }
    }
    int side_start = mesh.n_vertices();
    for (int j = 1; j < tess_h; j++) {
        float z = -1.0f + 2.0f * j / tess_h;
        for (int i = 0; i < tess_th; i++) {
            float th = M_TWOPIf * i / tess_th;
            mkpoint(mesh, r*cos(th), r*sin(th), z);
        }
    }
    int top_start = mesh.n_vertices();
    for (int j = tess_h; j > 0; j--) {
        float rr = r * j / tess_h;
        for (int i = 0; i < tess_th; i++) {
            float th = M_TWOPIf * i / tess_th;
            mkpoint(mesh, rr*cos(th), rr*sin(th), 1);
        }
    }
    mkpoint(mesh, 0, 0, 1);

    for (int i = 0; i < tess_th; i++)
        mkface(mesh, 0, ((i+1)%tess_th)+1, i+1);
    for (int j = 1; j < tess_h; j++) {
        int base = 1 + (j-1) * tess_th;
        for (int i = 0; i < tess_th; i++) {
            int i1 = (i+1)%tess_th;
            mkquad(mesh, base+tess_th+i1, base+tess_th+i,base+i1, base+i);
        }
    }

    for (int j = 0; j < tess_h; j++) {
        int base = side_start - tess_th + j * tess_th;
        for (int i = 0; i < tess_th; i++) {
            int i1 = (i+1)%tess_th;
            mkquad(mesh, base + i, base + i1,base+tess_th+i, base+tess_th+i1);
        }
    }

    for (int j = 0; j < tess_h-1; j++) {
        int base = top_start + j * tess_th;
        for (int i = 0; i < tess_th; i++) {
            int i1 = (i+1)%tess_th;
            mkquad(mesh, base+tess_th+i1, base+tess_th+i,base+i1, base+i);
        }
    }
    int base = top_start + (tess_h-1)*tess_th;
    for (int i = 0; i < tess_th; i++)
        mkface(mesh, base+i, base+((i+1)%tess_th), base+tess_th);
}
