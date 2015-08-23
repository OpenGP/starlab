#pragma once
#include "PQP.h"

#include <vector>
typedef std::vector<double> PQPPointType;
typedef std::vector< PQPPointType > PQPTriangleType;
typedef std::vector< PQPTriangleType > PQPTrianglesType;

#ifdef PQP_SUPPORT_SURFACEMESH
static inline PQPTrianglesType makeModelPQP( SurfaceMesh::SurfaceMeshModel * m ){
    PQPTrianglesType triangles;
    for(auto & f : m->faces()){
        PQPTriangleType tri;
        for(auto & v : m->vertices(f)){
            PQPPointType point;
            for(int i = 0; i < 3; i++) point.push_back( m->vertex_coordinates()[v][i] );
            tri.push_back(point);
        }
        triangles.push_back(tri);
    }
    return triangles;
}
#endif

namespace PQP{

struct IntersectResult{
    PQP_REAL distance;
    size_t id1, id2;
    IntersectResult(size_t id1, size_t id2, PQP_REAL p1[3], PQP_REAL p2[3], PQP_REAL distance = 0) :
        distance(distance), id1(id1), id2(id2){
        for(int i = 0; i < 3; i++){
            p[i] = p1[i];
            q[i] = p2[i];
        }
    }
    PQP_REAL p[3], q[3];
};

struct Manager{
    Manager(int num_models = 2){ models.reserve(num_models); }

    void addModel( const PQPTrianglesType & triangles )
    {
        models.push_back(PQP_Model());
        PQP_Model & m = models.back();

        int fid = 0;
        m.BeginModel();
        for(auto & tri : triangles) m.AddTri(&tri[0][0], &tri[1][0], &tri[2][0], fid++);
        m.EndModel();
    }

    std::vector<IntersectResult> testIntersection( size_t model_id1 = 0, size_t model_id2 = 1, PQP_REAL threshold = 1e-12 )
    {
        std::vector<IntersectResult> results;
        if(models.size() < 2) return results;

        PQP_Model & m1 = models[model_id1];
        PQP_Model & m2 = models[model_id2];

        // Leave models as-is
        PQP_REAL R1[3][3], R2[3][3];
        PQP_REAL T1[3], T2[3];
        makeIdentity(R1,T1);
        makeIdentity(R2,T2);

        // Perform collision detection
        PQP_Checker checker;
        PQP_CollideResult collisions;
        checker.PQP_Collide(&collisions, R1, T1, &m1, R2, T2, &m2);

        Tri_Processor isect_tri;
        PQP_REAL p[3],q[3];
        PQP_REAL t1[3][3], t2[3][3];

        for(int i = 0; i < collisions.num_pairs; i++)
        {
            auto & pair = collisions.pairs[i];
            results.push_back( IntersectResult( pair.id1, pair.id2, p, q, 0.0 ) );
        }

        // If no collisions detected return closest points
        if (results.empty())
        {
            PQP_DistanceResult distance;
            checker.PQP_Distance(&distance, R1, T1, &m1, R2, T2, &m2, threshold, threshold);

            for (int i = 0; i < 3; i++){
                p[i] = distance.p1[i];
                q[i] = distance.p2[i];
            }

            results.push_back(IntersectResult( distance.p1ID, distance.p2ID, p, q, distance.Distance() ));
        }

        return results;
    }

    std::vector<PQP_Model> models;

    static inline void makeIdentity( PQP_REAL R[3][3], PQP_REAL T[3] ){
        T[0] = T[1] = T[2] = 0;
        for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) R[i][j] = (i == j) ? 1 : 0;
    }

    static inline void getTriangle( Tri* t, PQP_REAL tri[3][3]){
        tri[0][0] = t->p1[0]; tri[0][1] = t->p1[1]; tri[0][2] = t->p1[2];
        tri[1][0] = t->p2[0]; tri[1][1] = t->p2[1]; tri[1][2] = t->p2[2];
        tri[2][0] = t->p3[0]; tri[2][1] = t->p3[1]; tri[2][2] = t->p3[2];
    }
};
}
