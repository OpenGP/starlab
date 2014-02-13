#include "visiblity_resampler.h"

#include "Octree.h"

#include "SamplingHelper.h"

#include "StarlabDrawArea.h"
#include "RenderObjectExt.h"

// OpenMP
#include <omp.h>

void traverseOctree( Octree & octree, starlab::CubeSoup & cs )
{
    if( octree.children.empty() && !octree.triangleData.empty() ){
        cs.addCube( starlab::QVector3( octree.boundingBox.Center() ), octree.boundingBox.yExtent * 2 );
		return;
	}

    foreach(Octree t, octree.children){
		traverseOctree(t,cs);
	}
}

void visiblity_resampler::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichBool("uniformSampling",true,"Uniform sampling"));
    pars->addParam(new RichBool("saveFileXYZ",false,"Save points to XYZ file"));
	pars->addParam(new RichBool("addModel",false,"Add as layer"));
	pars->addParam(new RichBool("viz",true,"Visualize"));
	pars->addParam(new RichInt("randSamples",1e4,"Number of random samples"));
	pars->addParam(new RichInt("sphereSamples",5,"Samples on unit sphere"));
    pars->addParam(new RichBool("triSoup", false, "Load as triangle soup"));
    pars->addParam(new RichBool("exportAsSoup", false, "Export visible soup"));
	pars->addParam(new RichBool("faceCenters", false, "Sample only face centers"));
}

std::vector<Vector3> randomSampleSphere( int numSamples, bool isRegular )
{
	double theta,rho,phi;
	double x,y,z;

    std::vector<Vector3> samples;

	int iAQuantize = 200;
	int iBQuantize = 50;
	double fAQuantize = (double)(iAQuantize-1);
	double fBQuantize = (double)(iBQuantize-1);

	for (int i = 0; i < numSamples; i++)
	{
		// choose a theta & rho
		if (isRegular)
		{
			theta = ((double)(rand() % iAQuantize)) / fAQuantize;
			rho   = ((double)(rand() % iBQuantize)) / fBQuantize;
		}
		else
		{
			theta = ((double)rand()) / (double)RAND_MAX;
			rho   = ((double)rand()) / (double)RAND_MAX;
		}

		// theta and rho now between 0 and 1
		theta *= 2.0 * 3.14159265359;
		rho   *= 2.0;

		// convert rho to phi;
		phi   = rho * 3.14159265359f * 0.5f;

		// spherical to Cartesian
		x = sin(phi) * cos(theta);
		y = cos(phi);
		z = sin(phi) * sin(theta);

        samples.push_back(Vector3(x,y,z));
	}

	return samples;
}

std::vector<Vector3> uniformSampleSphere( int numSamples = 5 )
{
	double x,y,z;

    std::vector<Vector3> samples;

	double d_theta = (M_PI) / numSamples;
	double d_psi = (M_PI) / numSamples;

	for(double theta = 0; theta <= M_PI; theta += d_theta)
	{
		for(double psi = 0; psi <= 2 * M_PI; psi += d_psi)
		{
			x = sin(theta) * cos(psi);
			y = sin(theta) * sin(psi);
			z = cos(theta);

            samples.push_back(Vector3(x,y,z));
		}
	}

	return samples;
}

void visiblity_resampler::applyFilter(RichParameterSet *pars)
{		
    SurfaceMeshModel * meshUsed = NULL;

    if( !pars->getBool("triSoup") )
    {
        meshUsed = mesh();
    }
    else
    {
        // Tri soup only supported for '.obj' files
        QFile file(mesh()->path);
        if (mesh()->path.endsWith("obj") && file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            meshUsed = new SurfaceMeshModel("tri_soup.obj", "tri_soup");
            QVector<Vector3> realPoints;

            QTextStream inF(&file);
            while( !inF.atEnd() ){
                QString line = inF.readLine();
                if(!line.size()) continue;

                if(line.at(0).toLatin1() == 'v' && (line.at(1).toLatin1() == ' ' || line.at(1).toLatin1() == '\t'))
                {
                    QStringList v = line.simplified().split(" ", QString::SkipEmptyParts);
                    realPoints.push_back( Vector3( v[1].toDouble(), v[2].toDouble(), v[3].toDouble() ) );
                }

                if(line.at(0).toLatin1() == 'f')
                {
                    QStringList f = line.simplified().split(" ", QString::SkipEmptyParts);
                    std::vector<SurfaceMesh::Vertex> verts;
                    verts.push_back(meshUsed->add_vertex(realPoints[f[1].split("/", QString::SkipEmptyParts).front().toInt() - 1]));
                    verts.push_back(meshUsed->add_vertex(realPoints[f[2].split("/", QString::SkipEmptyParts).front().toInt() - 1]));
                    verts.push_back(meshUsed->add_vertex(realPoints[f[3].split("/", QString::SkipEmptyParts).front().toInt() - 1]));
                    meshUsed->add_face(verts);
                }
            }

			meshUsed->updateBoundingBox();
			meshUsed->update_face_normals();
			meshUsed->update_vertex_normals();
			document()->addModel( meshUsed );
        }
        else
        {
            meshUsed = mesh();
        }
		file.close();
    }

    SurfaceMeshHelper h(meshUsed);
    Vector3VertexProperty points = meshUsed->vertex_property<Vector3>(VPOINT);

	double surfaceOffset = 1e-6;

    Octree octree(meshUsed);

    std::vector<Vector3> sphere = uniformSampleSphere( pars->getInt("sphereSamples") );
    int rayCount = (int) sphere.size();

	std::vector<SamplePoint> all_samples;

	if( !pars->getBool("faceCenters") )
	{
		if( pars->getBool("uniformSampling") )
		{
			// Similar sampling
            QVector<Vector3> points, normals;
            points = SamplingHelper(mesh()).similarSampling(pars->getInt("randSamples"), normals);
			for(int i = 0; i < (int)points.size(); i++)
				all_samples.push_back(SamplePoint(points[i],normals[i]));
		}
		else
		{
			// Random sampling
            //Sampler sampler(meshUsed);
            //all_samples = sampler.getSamples( pars->getInt("randSamples") );
		}
	}
	else
	{
		Vector3FaceProperty fcenters = h.vector3FaceProperty("f:center", Vector3(0,0,0));
		Vector3FaceProperty fnormals = meshUsed->face_normals(true);
		foreach( Face f, meshUsed->faces() )
		{
			foreach(Vertex v, meshUsed->vertices(f)) fcenters[f] += points[v];
			fcenters[f] /= meshUsed->valence(f);

			all_samples.push_back(SamplePoint( fcenters[f], fnormals[f], 1, f.idx() ));
		}
	}

	int N = (int) all_samples.size();
	std::vector<bool> isUse(N,false);

	#pragma omp parallel for
	for(int i = 0; i < N; i++)
	{
		const SamplePoint & sp = all_samples[i];

		for(int r = 0; r < rayCount; r++)
		{
			const Eigen::Vector3d & d = sphere[r];

			Eigen::Vector3d rayStart = sp.pos + (d * surfaceOffset);
			Ray ray( rayStart, d );
			
            if(octree.intersectRay(ray,0,true).empty())
			{
				isUse[i] = true;
				break;
			}
		}
	}

	document()->pushBusy();

    starlab::PointSoup * ps = NULL;
	SurfaceMesh::Model * m = NULL;

    QString newMeshName = QString("%1_sampled").arg(meshUsed->name);

	if(pars->getBool("addModel")) m = new SurfaceMesh::Model( newMeshName+".obj", newMeshName );
    if(pars->getBool("viz")) ps = new starlab::PointSoup;

	QTextStream * out = NULL;
	QString filename = "";
	if(pars->getBool("saveFileXYZ")) filename = QFileDialog::getSaveFileName(NULL,"Save XYZ","./", "XYZ file (*.xyz)");
	QFile file(filename); 
	if(pars->getBool("saveFileXYZ")){
		file.open(QFile::WriteOnly | QFile::Text);
		out = new QTextStream(&file);
	}

	std::vector<SamplePoint> used_samples;
    std::vector<Vector3> used_samples_points;

	for(int i = 0; i < N; i++){
		if(isUse[i]) 
		{
			used_samples.push_back(all_samples[i]);
			used_samples_points.push_back(all_samples[i].pos);
		}
	}

    /// Remove duplicates
    //std::vector<size_t> corner_xrefs;
    //weld(used_samples_points, corner_xrefs, std::hash_Vector3d(), std::equal_to<Vector3d>());

	QSet<int> visibleFaces;

	// Use samples
	for(int j = 0; j < (int)used_samples.size(); j++)
	{
        //int i = corner_xrefs[j];
        int i = j;

		const SamplePoint& sp = used_samples[i];

        if(pars->getBool("viz")) ps->addPointNormal( Vector3(sp.pos), Vector3(sp.n) );
		if(pars->getBool("addModel")) m->add_vertex( sp.pos );

		if(pars->getBool("saveFileXYZ"))
		{
            Vector3 p = sp.pos;
            Vector3 n = sp.n;
			(*out) << QString("%1 %2 %3 %4 %5 %6\n").arg(p[0]).arg(p[1]).arg(p[2]).arg(n[0]).arg(n[1]).arg(n[2]);
		}

		if(pars->getBool("exportAsSoup"))
		{
			visibleFaces.insert( sp.findex );
		}
	}

	file.close();

	if(pars->getBool("exportAsSoup"))
	{
		SurfaceMeshModel * visibleSoup = new SurfaceMeshModel("visible_soup.obj", "visible_soup");

		foreach(int fid, visibleFaces)
		{
			SurfaceMesh::Face f(fid);
			if(!meshUsed->is_valid(f)) continue;

			std::vector<Vector3> origPoints;
			foreach(Vertex v, meshUsed->vertices(f)) 
				origPoints.push_back( points[v] );

			std::vector<SurfaceMesh::Vertex> verts;
			verts.push_back(visibleSoup->add_vertex( origPoints[0] ));
			verts.push_back(visibleSoup->add_vertex( origPoints[1] ));
			verts.push_back(visibleSoup->add_vertex( origPoints[2] ));
			visibleSoup->add_face(verts);
		}

		document()->addModel( visibleSoup );
	}

	if(pars->getBool("addModel"))
	{
		m->updateBoundingBox();
		document()->addModel(m);
	}

	drawArea()->deleteAllRenderObjects();
	if(pars->getBool("viz")) 
	{
		drawArea()->addRenderObject(ps);
	}

	document()->popBusy();

	drawArea()->updateGL();
}
