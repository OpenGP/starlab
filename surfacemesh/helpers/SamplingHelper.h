#pragma once
#include "SurfaceMeshHelper.h"

// Helper structures
struct SamplePoint{
    Eigen::Vector3d pos, n;
    int findex; // index of sampled face
    double u,v;
    double weight;
    int flag;
    SamplePoint(const Eigen::Vector3d& position = Eigen::Vector3d(), const Eigen::Vector3d& normal = Eigen::Vector3d(),
        int face_index = -1.0, double U = 0.0, double V = 0.0, double Weight = 0.0, int flags = 0){
        pos = position;
        n = normal;
        weight = Weight;
        findex = face_index;
        u = U;
        v = V;
        flag = flags;
    }
};

class SamplingHelper : public virtual SurfaceMeshHelper{
protected:

public:
	SamplingHelper(SurfaceMeshModel* mesh) : SurfaceMeshHelper(mesh){}
	
	/// Stratified sampling:
	class SimilarSampler
	{
		// Helpers
		static inline Scalar deg_to_rad(const Scalar& _angle){ return M_PI*(_angle/180.0); }
        static inline Vector3 barycentric(Vector3 p, Vector3 a, Vector3 b, Vector3 c){
            Vector3 v0 = b - a, v1 = c - a, v2 = p - a;
			double d00 = dot(v0, v0); double d01 = dot(v0, v1);
			double d11 = dot(v1, v1); double d20 = dot(v2, v0);
			double d21 = dot(v2, v1); double denom = d00 * d11 - d01 * d01;
			double v = (d11 * d20 - d01 * d21) / denom;
			double w = (d00 * d21 - d01 * d20) / denom;
			double u = 1.0 - v - w;
            return Vector3(u,v,w);
		}
        static inline bool isValidBaryCoord(Vector3 coord){
			if(	coord[0] < 0 || coord[1] < 0 || coord[2] < 0 ||
				coord[0] > 1 || coord[1] > 1 || coord[2] > 1) return false;
			return true;
		}

        // Sort by QMap second value
        template<class F, class S>
        static inline bool sortByFirst(const QPair<F,S>& e1, const QPair<F,S>& e2) {
            return e1.first < e2.first;
        }

        template<class F, class S>
        static inline QList< QPair<S, F> > sortQMapByValue(const QMap<F,S> & map){
            QList< QPair<S, F> > result;

            QMapIterator<F, S> i(map);
            while (i.hasNext()) {
                i.next();
                result.push_back(qMakePair(i.value(), i.key()));
            }

            // Sort that list
            qSort(result.begin(), result.end(), sortByFirst<S,F>);
            return result;
        }

    public:
        static QVector<Vector3> FaceSamples(SurfaceMeshModel * m, int sampleNum,
					QVector<Vector3> & samplesNormals, double * avgSpacing, QVector<SamplePoint> * fullSamples = NULL)
		{
			QVector<Vector3> samples;

			// Compute face areas
			SurfaceMeshHelper h( m );
			h.computeFaceAreas();
			ScalarFaceProperty farea = m->face_property<Scalar>(FAREA);
			   
			Scalar area = 0;
			foreach(Face f, m->faces()) area += farea[f];

			// Compute edge lengths
			ScalarEdgeProperty elength = h.computeEdgeLengths();

			m->update_face_normals();
			Vector3FaceProperty fnormals = m->get_face_property<Vector3>(FNORMAL);

			Scalar samplePerAreaUnit = sampleNum / area;

			// Mesh points
			Vector3VertexProperty points = h.getVector3VertexProperty(VPOINT);

			// Similar Triangles sampling
			foreach(SurfaceMeshModel::Face f, m->faces())
			{
				// Collect vector of triangle points, and map to vertices
                std::vector<Vector3, Eigen::aligned_allocator<Vector3> > triangle, virtualTri;
                QMap<SurfaceMesh::Vertex, size_t> verts;
				Surface_mesh::Vertex_around_face_circulator vit = m->vertices(f),vend=vit;
				do{ verts[vit] = triangle.size(); triangle.push_back(points[vit]);  } while(++vit != vend);
				virtualTri = triangle;

				// Force virtual triangle to be isosceles
				{
					QMap<SurfaceMesh::Halfedge,double> edgeMap;

					// Classify edges by their lengths
					Surface_mesh::Halfedge_around_face_circulator hj(m, f), hend = hj;
					do{ edgeMap[hj] = elength[m->edge(hj)]; } while (++hj != hend);
					QList< QPair<double,SurfaceMesh::Halfedge> > edges = sortQMapByValue(edgeMap);
					SurfaceMesh::Halfedge S = edges.at(0).second, M = edges.at(1).second, L = edges.at(2).second;

					SurfaceMesh::Vertex vP = m->to_vertex(m->next_halfedge(L));
					SurfaceMesh::Vertex v0 = (vP == m->to_vertex(S)) ? m->from_vertex(S) : m->to_vertex(S);
					SurfaceMesh::Vertex vM = (vP == m->to_vertex(M)) ? m->from_vertex(M) : m->to_vertex(M);
                    Vector3 deltaS = (points[vP] - points[v0]).normalized();
                    Vector3 deltaL = (points[vM] - points[v0]).normalized();

					// Push vertex towards triangle with two equal edges
					virtualTri[ verts[vP] ] = points[v0] + (deltaS * elength[m->edge(M)]);

					// Shrink triangle to avoid sampling edges
					triangle[ verts[vP] ] = points[vP] + (-deltaS * elength[m->edge(S)] * 0.001);
					triangle[ verts[vM] ] = points[vM] + (-deltaL * elength[m->edge(L)] * 0.001);
				}

				double varea = 0.5 * cross(Vector3(virtualTri[1] - virtualTri[0]), Vector3(virtualTri[2] - virtualTri[0])).norm();
				
				// compute # samples in the current face
				int n_samples = (int) (0.5 * varea * samplePerAreaUnit);

				// Minimum number of samples per face
				n_samples = qMax(n_samples, 2);

				if(n_samples > 1)
				{
					int n_samples_per_edge = (int)((sqrt(1.0 + 8.0 * (Scalar)n_samples) + 5.0) / 2.0);

					Scalar segmentNum = n_samples_per_edge - 1 ;
					Scalar segmentLen = 1.0 / segmentNum;

					// face sampling
					for(int i = 1; i < (n_samples_per_edge - 1); i++)
					{
						for(int j = 1; j < (n_samples_per_edge - 1) - i; j++)
						{
							Scalar uvw[] = {i*segmentLen, j*segmentLen, 1.0 - (i*segmentLen + j*segmentLen)};

							// Get point from current barycentric coordinate
                            Vector3 p = Vector3::Zero();
							for(int vi = 0; vi < 3; vi++) 
								p = p + (virtualTri[vi] * uvw[vi]);

                            Vector3 coord = barycentric(p, triangle[0], triangle[1], triangle[2]);
							if( !isValidBaryCoord (coord) ) continue;

							samples.push_back( p );
							samplesNormals.push_back( fnormals[f] );

							if(avgSpacing && j > 1 && samples.size() > 1)
							{
								double dist = (samples.back() - samples[samples.size()-2]).norm();
								if(dist != 0.0)	*avgSpacing = dist;
							}

							// Detailed record of samples
							if( fullSamples )
								fullSamples->push_back(SamplePoint(p, fnormals[f], f.idx(), coord[0], coord[1]));
						}
					}
				}
			}

			return samples;
		}

        static QVector<Vector3> EdgeUniform(SurfaceMeshModel * m, int sampleNum, QVector<Vector3> & samplesNormals)
		{
			SurfaceMeshHelper h( m );

			// First loop compute total edge length;
			Scalar edgeSum = 0;
			if(!m->has_edge_property<Scalar>(ELENGTH)) h.computeEdgeLengths();
			ScalarEdgeProperty elength = m->edge_property<Scalar>(ELENGTH);
			foreach(Edge e, m->edges())	edgeSum += elength[e];

			Scalar sampleLen = edgeSum / sampleNum;

			return EdgeUniformFixed(m,samplesNormals,sampleLen);
		}

        static QVector<Vector3> EdgeUniformFixed( SurfaceMeshModel * m, QVector<Vector3> & samplesNormals, double sampleLen )
		{ 
			QVector<Vector3> samples;

			SurfaceMeshHelper h( m );

			// Mesh points
			Vector3VertexProperty points = h.getVector3VertexProperty(VPOINT);

			// Face normals
			m->update_face_normals();
			Vector3FaceProperty fnormals = m->get_face_property<Vector3>(FNORMAL);

			if(!m->has_edge_property<Scalar>(ELENGTH)) SurfaceMeshHelper(m).computeEdgeLengths();
			ScalarEdgeProperty elength = m->edge_property<Scalar>(ELENGTH);

			Scalar rest = 0;

			// This is a check for zero "sampleLen".. better solution is to fix it up in Face sampling
			Scalar sumEdgeLengths = 0;
			foreach(Edge ei, m->edges()) sumEdgeLengths += elength[ei];
			if(sampleLen == 0.0) sampleLen = (sumEdgeLengths / m->n_edges()) / 2;

			foreach(Edge ei, m->edges()){
				Scalar len = elength[ei];
				Scalar samplePerEdge = floor((len+rest)/sampleLen);
				rest = (len+rest) - samplePerEdge * sampleLen;
				Scalar step = 1.0 / (samplePerEdge + 1);
				for(int i = 0; i < samplePerEdge; ++i)
				{
					Scalar alpha = step*(i+1);
					Scalar beta = 1.0 - step*(i+1);
					samples.push_back( (alpha * points[m->vertex(ei,0)]) + (beta * points[m->vertex(ei,1)]) );

					// Normal = average of adj faces
					{
                        Vector3 normal(0,0,0);
						Face f1 = m->face(m->halfedge(ei,0)),f2 = m->face(m->halfedge(ei,1));
						if(f1.is_valid()) normal += fnormals[f1];
						if(f2.is_valid()) normal += fnormals[f1];
						if(f1.is_valid() && f2.is_valid()) normal /= 2.0;

						samplesNormals.push_back( normal );
					}
				}
			}

			return samples;
		}

        static QVector<Vector3> Vertices(SurfaceMeshModel * m, QVector<Vector3> & samplesNormals)
        {
            QVector<Vector3> samples;

            // Mesh points
            Vector3VertexProperty points = m->vertex_property<Vector3>(VPOINT);

            // Normals
            m->update_face_normals();
            m->update_vertex_normals();
            Vector3VertexProperty normals = m->vertex_property<Vector3>(VNORMAL);

            foreach(Vertex v, m->vertices())
            {
                samples.push_back(points[v]);
                samplesNormals.push_back(normals[v]);
            }

            return samples;
        }

        static QVector<Vector3> SimilarSampler::All( SurfaceMeshModel * m, int sampleNum, QVector<Vector3> & samplesNormals )
		{
			QVector<Vector3> samples;

			double sampleSpacing = 0;
			samples += SimilarSampler::FaceSamples( m, sampleNum, samplesNormals, &sampleSpacing );
			samples += SimilarSampler::EdgeUniformFixed( m, samplesNormals, sampleSpacing );
			samples += SimilarSampler::Vertices( m, samplesNormals );

			return samples;
		}
	};
	
    QVector<Vector3> similarSampling(int number_samples, QVector<Vector3> & normals ){
        return SimilarSampler::All(mesh, number_samples, normals);
	}
};
