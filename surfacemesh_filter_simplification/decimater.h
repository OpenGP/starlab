// Adapted from OpenFlipper
#pragma once
#include "SurfaceMeshModel.h"

typedef Eigen::Matrix<double,4,4,(Eigen::DontAlign|Eigen::RowMajor)> QuadricMatrix;

class Decimater{
private:
    size_t target_num_faces;
    Surface_mesh * mesh;
    Surface_mesh::Vertex_property<Point> points;
    Surface_mesh::Vertex_property<Normal> vnormal;

    Surface_mesh::Vertex_property<QuadricMatrix> quadrics;
    Surface_mesh::Face_property< std::vector<double> > fplane;
    Surface_mesh::Edge_property< double > errors;

public:
    Decimater(Surface_mesh* mesh, double percent = 0.75)
    {
        this->mesh = mesh;
        this->target_num_faces = mesh->n_faces() * percent;

        points = mesh->vertex_property<Point>("v:point");
        vnormal = mesh->vertex_property<Point>("v:normal");
    }

private:
    void prepare()
    {
        quadrics = mesh->vertex_property<QuadricMatrix>("v:quadrics");
        errors = mesh->edge_property< double >("e:errors");
        fplane = mesh->face_property< std::vector<double> >("f:planes");

        initial_quadrics();
        calculate_errors();
    }

    void cleanUp()
    {
        mesh->remove_vertex_property(quadrics);
        mesh->remove_face_property(fplane);
        mesh->remove_edge_property(errors);
    }

    QuadricMatrix plane(std::vector<double> & plane)
    {
        QuadricMatrix m;
        double a = plane[0];
        double b = plane[1];
        double c = plane[2];
        double d = plane[3];
        m(0) = a*a;  m(1) = a*b;  m(2) = a*c;  m(3) = a*d;
        m(4) = a*b;  m(5) = b*b;  m(6) = b*c;  m(7) = b*d;
        m(8) = a*c;  m(9) = b*c; m(10) = c*c; m(11) = c*d;
        m(12)= a*d; m(13) = b*d; m(14) = c*d; m(15) = d*d;
        return m;
    }

    void initial_quadrics()
    {
        Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();
        for(vit = mesh->vertices_begin(); vit != vend; ++vit)
            quadrics[vit] = QuadricMatrix::Zero();

        /* compute initial quadric */
        Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
        for(fit = mesh->faces_begin(); fit != fend; ++fit)
        {
            fplane[fit] = planeFace(fit);

            /* faces are triangles */
            Surface_mesh::Vertex_around_face_circulator fvit, fvend;
            fvit = fvend = mesh->vertices(fit);

            do{ quadrics[ fvit ] += plane(fplane[fit]); } while (++fvit != fvend);
        }
    }


    /* calculate determinant of a 3*3 Matrix
     * a(mn)  - index of Matrix element
     */
    double det3(const QuadricMatrix & m,
                int a11, int a12, int a13,
                int a21, int a22, int a23,
                int a31, int a32, int a33)
    {
        double det =  m(a11)*m(a22)*m(a33) + m(a13)*m(a21)*m(a32) + m(a12)*m(a23)*m(a31)
                    - m(a13)*m(a22)*m(a31) - m(a11)*m(a23)*m(a32) - m(a12)*m(a21)*m(a33);
        return det;
    }

    double calculate_edge_error(Surface_mesh::Edge edge, Point & p)
    {
        double min_error;
        QuadricMatrix q_bar;
        QuadricMatrix q_delta;

        Surface_mesh::Vertex v1 = mesh->vertex(edge, 0);
        Surface_mesh::Vertex v2 = mesh->vertex(edge, 1);

        /* computer quadric of virtual vertex vf */
        q_bar = quadrics[v1] + quadrics[v2];

        /* test if q_bar is symmetric */
        if (q_bar(1) != q_bar(4) || q_bar(2) != q_bar(8) || q_bar(6) != q_bar(9) ||
            q_bar(3) != q_bar(12) || q_bar(7) != q_bar(13) || q_bar(11) != q_bar(14)){
            return std::numeric_limits<double>::max();
        }

        q_delta <<  q_bar(0), q_bar(1),  q_bar(2),  q_bar(3),
                    q_bar(4), q_bar(5),  q_bar(6),  q_bar(7),
                    q_bar(8), q_bar(9), q_bar(10), q_bar(11),
                           0,        0,         0,        1;

        /* if q_delta is invertible */
        if ( double det = det3(q_delta,0, 1, 2, 4, 5, 6, 8, 9, 10) )   /* note that det(q_delta) equals to M44 */
        {
            p.x() = -1/det*(det3(q_delta,1, 2, 3, 5, 6, 7, 9, 10, 11));
            p.y() =  1/det*(det3(q_delta,0, 2, 3, 4, 6, 7, 8, 10, 11));
            p.z() = -1/det*(det3(q_delta,0, 1, 3, 4, 5, 7, 8, 9, 11));
        }
        /*
        * if q_delta is NOT invertible, select
        * vertex from v1, v2, and (v1+v2)/2
        */
        else
        {
            Point p1 = points[v1], p2 = points[v2];
            Point p3 = (p1 + p2) / 2.0;

            double error1 = vertex_error(q_bar, p1);
            double error2 = vertex_error(q_bar, p2);
            double error3 = vertex_error(q_bar, p3);

            min_error = std::min(error1, std::min(error2, error3));
            if (error1 == min_error) { p = p1; }
            if (error2 == min_error) { p = p2; }
            if (error3 == min_error) { p = p3; }
        }

        min_error = vertex_error(q_bar, p);

        return min_error;
    }

    void calculate_errors()
    {
        // Populate errors on edges
        Surface_mesh::Edge_iterator edgesEnd = mesh->edges_end();
        for (  Surface_mesh::Edge_iterator eit = mesh->edges_begin(); eit != edgesEnd; ++eit)
        {
            Point p(0,0,0);
            errors[eit] = calculate_edge_error(eit,p);
        }
    }

    inline double vertex_error(const QuadricMatrix& q, Point p)
    {
        double x = p.x(), y = p.y(), z = p.z();
        return q(0)*x*x + 2*q(1)*x*y + 2*q(2)*x*z + 2*q(3)*x + q(5)*y*y
            + 2*q(6)*y*z + 2*q(7)*y + q(10)*z*z + 2*q(11)*z + q(15);
    }

    std::vector<double> planeFace(Surface_mesh::Face f)
    {
        std::vector<double> f_plane(4);
        std::vector<Point> v;
        double a,b,c,M;

        // Get face points
        Surface_mesh::Vertex_around_face_circulator fvit, fvend;
        fvit = fvend = mesh->vertices(f);
        do{ v.push_back(points[fvit]);	} while (++fvit != fvend);

        // Compute plane coefficient for face
        a = (v[1].y()-v[0].y())*(v[2].z()-v[0].z()) - (v[1].z()-v[0].z())*(v[2].y()-v[0].y());   /* a1*b2 - a2*b1; */
        b = (v[1].z()-v[0].z())*(v[2].x()-v[0].x()) - (v[1].x()-v[0].x())*(v[2].z()-v[0].z());   /* a2*b0 - a0*b2; */
        c = (v[1].x()-v[0].x())*(v[2].y()-v[0].y()) - (v[1].y()-v[0].y())*(v[2].x()-v[0].x());   /* a0*b1 - a1*b0; */
        M = sqrt(a*a + b*b + c*c);
        a = a/M; b = b/M; c = c/M;
        f_plane[0] = a;	f_plane[1] = b;	f_plane[2] = c;
        f_plane[3] = -1*(a*v[0].x() + b*v[0].y() + c*v[0].z());

        return f_plane;
    }

    void doSimplify()
    {
        prepare();

        int originalNumFaces = mesh->n_faces();

        while (mesh->n_faces() > target_num_faces)
        {
            /* find least-error edge */
            Surface_mesh::Edge minEdge = Surface_mesh::Edge(mesh->edges_begin());

            Surface_mesh::Edge_iterator edgesEnd = mesh->edges_end();
            for (Surface_mesh::Edge_iterator eit = mesh->edges_begin(); eit != edgesEnd; ++eit)
                if(errors[eit] < errors[minEdge])
                    minEdge = eit;

            Surface_mesh::Halfedge minHEdge = mesh->halfedge(minEdge, 0);

            /* update coordinate for modified v1 */
            Surface_mesh::Vertex v1 = mesh->to_vertex(minHEdge), v2 = mesh->from_vertex(minHEdge);
            Point newPos(0,0,0); calculate_edge_error(minEdge, newPos);
            points[v1] = newPos;

            /* update quadric of v1 */
            quadrics[v1] = quadrics[v1] + quadrics[v2];

            /* merge pairs of v2 to v1 */
            mesh->collapse(minHEdge);

            /* update error of pairs involving v1 */
            Surface_mesh::Halfedge_around_vertex_circulator hend = mesh->halfedges(v1);
            for (Surface_mesh::Halfedge_around_vertex_circulator h = mesh->halfedges(v1); ; )
            {
                Surface_mesh::Edge e = mesh->edge(h);

                Point pi(0,0,0);
                errors[e] = calculate_edge_error(e, pi);

                ++h; if(h == hend) break;
            }

            // report percent
            int progress = (1.0 - (double(mesh->n_faces() - target_num_faces) / (originalNumFaces - target_num_faces))) * 100;

            printf("decimation progress: %d %% \r", progress);
        }

        cleanUp();
    }

public:
    static void simplify(Surface_mesh *mesh, double percent = 0.75)
    {
        Decimater d(mesh, percent);
        d.doSimplify();
    }
};
