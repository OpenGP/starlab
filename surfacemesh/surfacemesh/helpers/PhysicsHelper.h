#pragma once
#include "SurfaceMeshHelper.h"
#include "StarlabDrawArea.h"

class PhysicsHelper : public SurfaceMeshHelper{
public:
	PhysicsHelper(SurfaceMeshModel* mesh) : SurfaceMeshHelper(mesh){}

	inline Scalar surfaceArea(){
		Scalar area = 0;
		ScalarFaceProperty fareas = SurfaceMeshHelper(mesh).computeFaceAreas();
		for(auto f : mesh->faces()) area += fareas[f];
		return area;
	}

	inline Scalar volume(){
		Scalar volume = 0;
		Vector3VertexProperty points = mesh->vertex_coordinates();

		for(auto f: mesh->faces()){
			std::vector<Vector3> p;
			for(auto v : mesh->vertices(f)) p.push_back(points[v]);
			Scalar signedVolumeTri = p[0].dot(p[1].cross(p[2])) / 6.0;
			volume += signedVolumeTri;
		}

		return volume;
	}

	SurfaceMesh::Vector3 centerOfMass( Eigen::Matrix3d & inertiaTensor = Eigen::Matrix3d(), 
			Scalar density = 1.0, Scalar * meshMass = NULL )
	{
		Scalar T0;
		Vector3 T1, T2, TP, r;
		int X = 0, Y = 1, Z = 2;

		compVolumeIntegrals(T0, T1, T2, TP);

		Scalar mass = density * T0;

		/* compute center of mass */
		r[X] = T1[X] / T0;
		r[Y] = T1[Y] / T0;
		r[Z] = T1[Z] / T0;

		/* compute inertia tensor */
		Eigen::Matrix3d J;
		J(X,X) = density * (T2[Y] + T2[Z]);
		J(Y,Y) = density * (T2[Z] + T2[X]);
		J(Z,Z) = density * (T2[X] + T2[Y]);
		J(X,Y) = J(Y,X) = - density * TP[X];
		J(Y,Z) = J(Z,Y) = - density * TP[Y];
		J(Z,X) = J(X,Z) = - density * TP[Z];

		/* translate inertia tensor to center of mass */
		J(X,X) -= mass * (r[Y]*r[Y] + r[Z]*r[Z]);
		J(Y,Y) -= mass * (r[Z]*r[Z] + r[X]*r[X]);
		J(Z,Z) -= mass * (r[X]*r[X] + r[Y]*r[Y]);
		J(X,Y) = J(Y,X) += mass * r[X] * r[Y]; 
		J(Y,Z) = J(Z,Y) += mass * r[Y] * r[Z]; 
		J(Z,X) = J(X,Z) += mass * r[Z] * r[X]; 
		inertiaTensor = J;

		if(meshMass) *meshMass = mass;

		return r;
	}

private:

	/* Face wrapper */
	typedef struct {
		Scalar norm[3];
		Scalar w;
		std::vector<SurfaceMesh::Vertex> verts;
	} FACE;

	/*	Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties," 
		journal of graphics tools, volume 1, number 1, 1996. */
	void compVolumeIntegrals(Scalar & T0, Vector3 & T1, Vector3 & T2, Vector3 & TP)
	{
		int X = 0, Y = 1, Z = 2;

		Scalar nx, ny, nz;

		int A, B, C;

		// Clear values
		T0 = T1[X] = T1[Y] = T1[Z] = T2[X] = T2[Y] = T2[Z] = TP[X] = TP[Y] = TP[Z] = 0;

		Scalar P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;
		Scalar Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

		#define SQR(x) ((x)*(x))
		#define CUBE(x) ((x)*(x)*(x))

		Vector3VertexProperty pverts = mesh->vertex_coordinates();
		Vector3FaceProperty fnormals = mesh->face_normals(true);

		for(Face face : mesh->faces())
		{
			// Fill container with face data
			FACE faceContainer;
			FACE * f = &faceContainer;
			{
				for(Vertex v : mesh->vertices(face)) f->verts.push_back(v);

				f->norm[0] = fnormals[face][0];
				f->norm[1] = fnormals[face][1];
				f->norm[2] = fnormals[face][2];
				f->w =	- f->norm[X] * pverts[f->verts[0]][X]
						- f->norm[Y] * pverts[f->verts[0]][Y]
						- f->norm[Z] * pverts[f->verts[0]][Z];
			}
			
			nx = fabs(f->norm[X]);
			ny = fabs(f->norm[Y]);
			nz = fabs(f->norm[Z]);

			if (nx > ny && nx > nz) C = X;
			else C = (ny > nz) ? Y : Z;
			A = (C + 1) % 3;
			B = (A + 1) % 3;

			//void compFaceIntegrals(FACE *f)
			{
				Scalar *n, w;
				Scalar k1, k2, k3, k4;

				//void compProjectionIntegrals(FACE *f)
				{
					Scalar a0, a1, da;
					Scalar b0, b1, db;
					Scalar a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
					Scalar a1_2, a1_3, b1_2, b1_3;
					Scalar C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
					Scalar Cab, Kab, Caab, Kaab, Cabb, Kabb;

					P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

					for (size_t i = 0; i < f->verts.size(); i++) 
					{
						a0 = pverts[f->verts[i]][A];
						b0 = pverts[f->verts[i]][B];
						a1 = pverts[f->verts[(i+1) % f->verts.size()]][A];
						b1 = pverts[f->verts[(i+1) % f->verts.size()]][B];

						da = a1 - a0;
						db = b1 - b0;
						a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
						b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
						a1_2 = a1 * a1; a1_3 = a1_2 * a1; 
						b1_2 = b1 * b1; b1_3 = b1_2 * b1;

						C1 = a1 + a0;
						Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
						Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
						Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
						Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
						Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
						Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

						P1 += db*C1;
						Pa += db*Ca;
						Paa += db*Caa;
						Paaa += db*Caaa;
						Pb += da*Cb;
						Pbb += da*Cbb;
						Pbbb += da*Cbbb;
						Pab += db*(b1*Cab + b0*Kab);
						Paab += db*(b1*Caab + b0*Kaab);
						Pabb += da*(a1*Cabb + a0*Kabb);
					}

					P1 /= 2.0;
					Pa /= 6.0;
					Paa /= 12.0;
					Paaa /= 20.0;
					Pb /= -6.0;
					Pbb /= -12.0;
					Pbbb /= -20.0;
					Pab /= 24.0;
					Paab /= 60.0;
					Pabb /= -60.0;
				}

				w = f->w;
				n = f->norm;
				k1 = 1 / n[C]; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

				Fa = k1 * Pa;
				Fb = k1 * Pb;
				Fc = -k2 * (n[A]*Pa + n[B]*Pb + w*P1);

				Faa = k1 * Paa;
				Fbb = k1 * Pbb;
				Fcc = k3 * (SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb
					+ w*(2*(n[A]*Pa + n[B]*Pb) + w*P1));

				Faaa = k1 * Paaa;
				Fbbb = k1 * Pbbb;
				Fccc = -k4 * (CUBE(n[A])*Paaa + 3*SQR(n[A])*n[B]*Paab 
					+ 3*n[A]*SQR(n[B])*Pabb + CUBE(n[B])*Pbbb
					+ 3*w*(SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb)
					+ w*w*(3*(n[A]*Pa + n[B]*Pb) + w*P1));

				Faab = k1 * Paab;
				Fbbc = -k2 * (n[A]*Pabb + n[B]*Pbbb + w*Pbb);
				Fcca = k3 * (SQR(n[A])*Paaa + 2*n[A]*n[B]*Paab + SQR(n[B])*Pabb
					+ w*(2*(n[A]*Paa + n[B]*Pab) + w*Pa));
			}

			T0 += f->norm[X] * ((A == X) ? Fa : ((B == X) ? Fb : Fc));

			T1[A] += f->norm[A] * Faa;
			T1[B] += f->norm[B] * Fbb;
			T1[C] += f->norm[C] * Fcc;
			T2[A] += f->norm[A] * Faaa;
			T2[B] += f->norm[B] * Fbbb;
			T2[C] += f->norm[C] * Fccc;
			TP[A] += f->norm[A] * Faab;
			TP[B] += f->norm[B] * Fbbc;
			TP[C] += f->norm[C] * Fcca;
		}

		T1[X] /= 2; T1[Y] /= 2; T1[Z] /= 2;
		T2[X] /= 3; T2[Y] /= 3; T2[Z] /= 3;
		TP[X] /= 2; TP[Y] /= 2; TP[Z] /= 2;
	}
};
