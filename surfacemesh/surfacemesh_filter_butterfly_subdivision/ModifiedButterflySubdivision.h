#pragma once
#include "SurfaceMeshModel.h"

#define get_point(x) (points[x])

// Adapted from OpenMesh
/** Modified Butterfly subdivision algorithm
*
* Implementation of the modified butterfly scheme of Denis Zorin, Peter Schröder and Wim Sweldens, 
* 'Interpolating subdivision for meshes with arbitrary topology' in Proceedings
* of SIGGRAPH 1996, ACM SIGGRAPH, 1996, pp. 189-192.
*/
class ModifiedButterfly{
public:

	typedef double									real_t;
	typedef Surface_mesh							mesh_t;

	typedef std::vector< std::vector<real_t> >      weights_t;
	typedef std::vector<real_t>                     weight_t;

public:

	ModifiedButterfly()
	{ init_weights(); }

public:

	/// Pre-compute weights
	void init_weights(size_t _max_valence=50)
	{
		weights.resize(_max_valence);

		//special case: K==3, K==4
		weights[3].resize(4);
		weights[3][0] = real_t(5.0)/12;
		weights[3][1] = real_t(-1.0)/12;
		weights[3][2] = real_t(-1.0)/12;
		weights[3][3] = real_t(3.0)/4;

		weights[4].resize(5);
		weights[4][0] = real_t(3.0)/8;
		weights[4][1] = 0;
		weights[4][2] = real_t(-1.0)/8;
		weights[4][3] = 0;
		weights[4][4] = real_t(3.0)/4;

		for(unsigned int K = 5; K<_max_valence; ++K)
		{
			weights[K].resize(K+1);
			// s(j) = ( 1/4 + cos(2*pi*j/K) + 1/2 * cos(4*pi*j/K) )/K
			real_t   invK  = 1.0/real_t(K);
			real_t sum = 0;
			for(unsigned int j=0; j<K; ++j)
			{
				weights[K][j] = (0.25 + cos(2.0*M_PI*j*invK) + 0.5*cos(4.0*M_PI*j*invK))*invK;
				sum += weights[K][j];
			}
			weights[K][K] = (real_t)1.0 - sum;
		}
	}


protected:


	bool prepare( mesh_t& _m )
	{
		points = _m.vertex_property<Point>("v:point");

		vp_pos_ = _m.vertex_property<Point>("v:point2");
		ep_pos_ = _m.edge_property<Point>("e:point");
		return true;
	}


	bool cleanup( mesh_t& _m )
	{
		_m.remove_vertex_property(vp_pos_);
		_m.remove_edge_property(ep_pos_);

		return true;
	}

public:
    bool subdivide( Surface_mesh& _m, size_t _n)
	{
		prepare(_m);

		///TODO:Implement fixed positions
		mesh_t::Face_iterator   fit, f_end;
		mesh_t::Edge_iterator   eit, e_end;
		mesh_t::Vertex_iterator vit;

		// Do _n subdivisions
		for (size_t i=0; i < _n; ++i)
		{

			// This is an interpolating scheme, old vertices remain the same.
			mesh_t::Vertex_iterator initialVerticesEnd = _m.vertices_end();
			for ( vit  = _m.vertices_begin(); vit != initialVerticesEnd; ++vit)
				vp_pos_[vit] = get_point(vit);

			// Compute position for new vertices and store them in the edge property
			for (eit=_m.edges_begin(); eit != _m.edges_end(); ++eit)
				compute_midpoint( _m, eit );
			
			// Split each edge at midpoint and store precomputed positions (stored in
			// edge property ep_pos_) in the vertex property vp_pos_;

			// Attention! Creating new edges, hence make sure the loop ends correctly.
			e_end = _m.edges_end();
			for (eit=_m.edges_begin(); eit != e_end; ++eit)
				split_edge(_m, eit );


			// Commit changes in topology and reconstitute consistency

			// Attention! Creating new faces, hence make sure the loop ends correctly.
			f_end   = _m.faces_end();
			for (fit = _m.faces_begin(); fit != f_end; ++fit)
				split_face(_m, fit );


			// Commit changes in geometry
			for ( vit  = /*initialVerticesEnd;*/_m.vertices_begin();
				vit != _m.vertices_end(); ++vit)
				points[vit] = vp_pos_[vit];
		}

		return cleanup(_m);
	}

private: // topological modifiers

	void split_face(mesh_t& _m, const mesh_t::Face& _fh)
	{
		mesh_t::Halfedge
			heh1(_m.halfedge(_fh)),
			heh2(_m.next_halfedge(_m.next_halfedge(heh1))),
			heh3(_m.next_halfedge(_m.next_halfedge(heh2)));

		// Cutting off every corner of the 6_gon
		corner_cutting( _m, heh1 );
		corner_cutting( _m, heh2 );
		corner_cutting( _m, heh3 );
	}


	void corner_cutting(mesh_t& _m, const mesh_t::Halfedge& _he)
	{
		// Define Halfedge Handles
		mesh_t::Halfedge
			heh1(_he),
			heh5(heh1),
			heh6(_m.next_halfedge(heh1));

		// Cycle around the polygon to find correct Halfedge
		for (; _m.next_halfedge(_m.next_halfedge(heh5)) != heh1;
			heh5 = _m.next_halfedge(heh5)){}

		mesh_t::Vertex
			vh1 = _m.to_vertex(heh1),
			vh2 = _m.to_vertex(heh5);

		mesh_t::Halfedge
			heh2(_m.next_halfedge(heh5)),
			heh3(_m.new_edge( vh1, vh2)),
			heh4(_m.opposite_halfedge(heh3));

		/* Intermediate result
		*
		*            *
		*         5 /|\
		*          /_  \
		*    vh2> *     *
		*        /|\3   |\
		*       /_  \|4   \
		*      *----\*----\*
		*          1 ^   6
		*            vh1 (adjust_outgoing halfedge!)
		*/

		// Old and new Face
		mesh_t::Face     fh_old(_m.face(heh6));
		mesh_t::Face     fh_new(_m.new_face());


		// Re-Set Handles around old Face
		_m.set_next_halfedge(heh4, heh6);
		_m.set_next_halfedge(heh5, heh4);

		_m.set_face(heh4, fh_old);
		_m.set_face(heh5, fh_old);
		_m.set_face(heh6, fh_old);
		_m.set_halfedge(fh_old, heh4);

		// Re-Set Handles around new Face
		_m.set_next_halfedge(heh1, heh3);
		_m.set_next_halfedge(heh3, heh2);

		_m.set_face(heh1, fh_new);
		_m.set_face(heh2, fh_new);
		_m.set_face(heh3, fh_new);

		_m.set_halfedge(fh_new, heh1);
	}


	void split_edge(mesh_t& _m, const mesh_t::Edge& _eh)
	{
		mesh_t::Halfedge heh     = _m.halfedge(_eh, 0), opp_heh = _m.halfedge(_eh, 1);

		mesh_t::Halfedge new_heh, opp_new_heh, t_heh;
		mesh_t::Vertex   vh;
		mesh_t::Vertex   vh1(_m.to_vertex(heh));
		Point    zero(0,0,0);

		// new vertex
		vh                = _m.add_vertex( zero );

		// memorize position, will be set later
		vp_pos_[vh] = ep_pos_[_eh];

		// Re-link mesh entities
		if (_m.is_boundary(_eh))
		{
			for (t_heh = heh;
				_m.next_halfedge(t_heh) != opp_heh;
				t_heh = _m.opposite_halfedge(_m.next_halfedge(t_heh))){}
		}
		else
		{
			for (t_heh = _m.next_halfedge(opp_heh);
				_m.next_halfedge(t_heh) != opp_heh;
				t_heh = _m.next_halfedge(t_heh) ){}
		}

		new_heh     = _m.new_edge(vh, vh1);
		opp_new_heh = _m.opposite_halfedge(new_heh);
		_m.set_vertex( heh, vh );

		_m.set_next_halfedge(t_heh, opp_new_heh);
		_m.set_next_halfedge(new_heh, _m.next_halfedge(heh));
		_m.set_next_halfedge(heh, new_heh);
		_m.set_next_halfedge(opp_new_heh, opp_heh);

		if (_m.face(opp_heh).is_valid())
		{
			_m.set_face(opp_new_heh, _m.face(opp_heh));
			_m.set_halfedge(_m.face(opp_new_heh), opp_new_heh);
		}

		_m.set_face( new_heh, _m.face(heh) );
		_m.set_halfedge( vh, new_heh);
		_m.set_halfedge( _m.face(heh), heh );
		_m.set_halfedge( vh1, opp_new_heh );

		// Never forget this, when playing with the topology
		_m.adjust_outgoing_halfedge( vh );
		_m.adjust_outgoing_halfedge( vh1 );
	}

private: // geometry helper

	void compute_midpoint(mesh_t& _m, const mesh_t::Edge& _eh)
	{
		mesh_t::Halfedge heh, opp_heh;

		heh      = _m.halfedge( _eh, 0);
		opp_heh  = _m.halfedge( _eh, 1);

		Point pos(0,0,0);

		mesh_t::Vertex a_0(_m.to_vertex(heh));
		mesh_t::Vertex a_1(_m.to_vertex(opp_heh));

		// boundary edge: 4-point scheme
		if (_m.is_boundary(_eh) )
		{
			pos = get_point(a_0);
			pos += get_point(a_1);
			pos *= 9.0/16;
			Point tpos;
			if(_m.is_boundary(heh))
			{
				tpos = get_point(_m.to_vertex(_m.next_halfedge(heh)));
				tpos += get_point(_m.to_vertex(_m.opposite_halfedge(_m.prev_halfedge(heh))));
			}
			else
			{
				assert(_m.is_boundary(opp_heh));
				tpos = get_point(_m.to_vertex(_m.next_halfedge(opp_heh)));
				tpos += get_point(_m.to_vertex(_m.opposite_halfedge(_m.prev_halfedge(opp_heh))));
			}
			tpos *= -1.0/16;
			pos += tpos;
		}
		else
		{
			int valence_a_0 = _m.valence(a_0);
			int valence_a_1 = _m.valence(a_1);
			assert(valence_a_0>2);
			assert(valence_a_1>2);

			if( (valence_a_0==6 && valence_a_1==6) || (_m.is_boundary(a_0) && valence_a_1==6) || (_m.is_boundary(a_1) && valence_a_0==6) || (_m.is_boundary(a_0) && _m.is_boundary(a_1)) )// use 8-point scheme
			{
				real_t alpha    = real_t(1.0/2);
				real_t beta     = real_t(1.0/8);
				real_t gamma    = real_t(-1.0/16);

				//get points
				mesh_t::Vertex b_0, b_1, c_0, c_1, c_2, c_3;
				mesh_t::Halfedge t_he;

				t_he = _m.next_halfedge(_m.opposite_halfedge(heh));
				b_0 = _m.to_vertex(t_he);
				if(!_m.is_boundary(_m.opposite_halfedge(t_he)))
				{
					t_he = _m.next_halfedge(_m.opposite_halfedge(t_he));
					c_0 = _m.to_vertex(t_he);
				}

				t_he = _m.opposite_halfedge(_m.prev_halfedge(heh));
				b_1 = _m.to_vertex(t_he);
				if(!_m.is_boundary(t_he))
				{
					t_he = _m.opposite_halfedge(_m.prev_halfedge(t_he));
					c_1 = _m.to_vertex(t_he);
				}

				t_he = _m.next_halfedge(_m.opposite_halfedge(opp_heh));
				assert(b_1.idx()==_m.to_vertex(t_he).idx());
				if(!_m.is_boundary(_m.opposite_halfedge(t_he)))
				{
					t_he = _m.next_halfedge(_m.opposite_halfedge(t_he));
					c_2 = _m.to_vertex(t_he);
				}

				t_he = _m.opposite_halfedge(_m.prev_halfedge(opp_heh));
				assert(b_0==_m.to_vertex(t_he));
				if(!_m.is_boundary(t_he))
				{
					t_he = _m.opposite_halfedge(_m.prev_halfedge(t_he));
					c_3 = _m.to_vertex(t_he);
				}

				//compute position.
				//a0,a1,b0,b1 must exist.
				assert(a_0.is_valid());
				assert(a_1.is_valid());
				assert(b_0.is_valid());
				assert(b_1.is_valid());
				//The other vertices may be created from symmetry is they are on the other side of the boundary.

				pos = get_point(a_0);
				pos += get_point(a_1);
				pos *= alpha;

				Point tpos ( get_point(b_0) );
				tpos += get_point(b_1);
				tpos *= beta;
				pos += tpos;

				Point pc_0, pc_1, pc_2, pc_3;
				if(c_0.is_valid())
					pc_0 = get_point(c_0);
				else //create the point by symmetry
				{
					pc_0 = get_point(a_1) + get_point(b_0) - get_point(a_0);
				}
				if(c_1.is_valid())
					pc_1 = get_point(c_1);
				else //create the point by symmetry
				{
					pc_1 = get_point(a_1) + get_point(b_1) - get_point(a_0);
				}
				if(c_2.is_valid())
					pc_2 = get_point(c_2);
				else //create the point by symmetry
				{
					pc_2 = get_point(a_0) + get_point(b_1) - get_point(a_1);
				}
				if(c_3.is_valid())
					pc_3 = get_point(c_3);
				else //create the point by symmetry
				{
					pc_3 = get_point(a_0) + get_point(b_0) - get_point(a_1);
				}
				tpos = pc_0;
				tpos += pc_1;
				tpos += pc_2;
				tpos += pc_3;
				tpos *= gamma;
				pos += tpos;
			}
			else //at least one endpoint is [irregular and not in boundary]
			{
				double normFactor = 0.0;

				if(valence_a_0!=6 && !_m.is_boundary(a_0))
				{
					assert((int)weights[valence_a_0].size()==valence_a_0+1);
					mesh_t::Halfedge t_he = opp_heh;
					for(int i = 0; i < valence_a_0 ; t_he=_m.next_halfedge(_m.opposite_halfedge(t_he)), ++i)
					{
						pos += weights[valence_a_0][i] * get_point(_m.to_vertex(t_he));
					}
					assert(t_he==opp_heh);

					//add irregular vertex:
					pos += weights[valence_a_0][valence_a_0] * get_point(a_0);
					++normFactor;
				}

				if(valence_a_1!=6  && !_m.is_boundary(a_1))
				{
					assert((int)weights[valence_a_1].size()==valence_a_1+1);
					mesh_t::Halfedge t_he = heh;
					for(int i = 0; i < valence_a_1 ; t_he=_m.next_halfedge(_m.opposite_halfedge(t_he)), ++i)
					{
						pos += weights[valence_a_1][i] * get_point(_m.to_vertex(t_he));
					}
					assert(t_he==heh);
					//add irregular vertex:
					pos += weights[valence_a_1][valence_a_1] * get_point(a_1);
					++normFactor;
				}

				assert(normFactor>0.1); //normFactor should be 1 or 2

				//if both vertices are irregular, average positions:
				pos /= normFactor;
			}
		}

		ep_pos_ [_eh] = pos;
	}

private: // data

	Surface_mesh::Vertex_property<Point> vp_pos_;
	Surface_mesh::Edge_property<Point> ep_pos_;
	
	Surface_mesh::Vertex_property<Point> points;

	weights_t     weights;
};
