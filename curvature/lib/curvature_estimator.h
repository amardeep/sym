/***************************************************************************
curvature_estimator.h  - curvature tensor estimator for CGAL polyhedron
----------------------------------------------------------------------------
D.Cohen-Steiner and J.-M. Morvan.
Restricted delaunay triangulations and normal cycle. 
In Proc. 19th Annu. ACM Sympos. Comput. Geom. 2003.
----------------------------------------------------------------------------
begin                : nov 2003
copyright            : (C) 2003 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
***************************************************************************/

#ifndef CURVATURE_ESTIMATOR_H
#define CURVATURE_ESTIMATOR_H

#include "Enriched_polyhedron.h"
#include "Curvature.h"

template <class Polyhedron,class kernel>
class CCurvature_estimator
{
	typedef typename kernel::FT FT;
	typedef typename kernel::RT RT;
	typedef typename kernel::Point_3 Point;
	typedef typename kernel::Vector_3 Vector;
	typedef typename kernel::Direction_3 Direction;
	typedef typename Polyhedron::Facet_handle                       Facet_handle;
	typedef typename Polyhedron::Vertex_handle                      Vertex_handle;
	typedef typename Polyhedron::Halfedge_handle                    Halfedge_handle;
	typedef typename Polyhedron::Halfedge_iterator                  Halfedge_iterator;
	typedef typename Polyhedron::Facet_iterator                     Facet_iterator;
	typedef typename Polyhedron::Vertex_iterator                    Vertex_iterator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator  HV_circulator;
	typedef typename Polyhedron::Halfedge_around_facet_circulator   HF_circulator;
	typedef CCurvature<kernel> Curvature;


public:
	// constants
	enum {ONE_RING,SHORTEST_EDGE,SPHERE,MULTIPASS};

private:
	#define HALF_PI 1.57079632679489661;

	// data
	Polyhedron *m_pMesh;

public:

	// life cycle
	CCurvature_estimator(Polyhedron *pMesh)
	{
		m_pMesh = pMesh;
	}
	~CCurvature_estimator() {}

	// estimate curvature tensors for all vertices
	void run(int averaging_region = ONE_RING,
					 FT radius_averaging_region = 0.01)
	{
		Vertex_iterator hVertex = NULL;
		for(hVertex = m_pMesh->vertices_begin();
			  hVertex != m_pMesh->vertices_end();
				hVertex++)
			run(hVertex,averaging_region,radius_averaging_region);
	}

	// estimate curvature tensors for one vertex
	void run(Vertex_handle hVertex,
		       int averaging_region,
					 FT radius_averaging_region)
	{
		switch(averaging_region)
		{
			case ONE_RING:
				run_one_ring(hVertex,false);
				break;
			case SHORTEST_EDGE:
				run_one_ring(hVertex,true);
				break;
			case SPHERE:
				break;
		}
	}

	void run_one_ring(Vertex_handle hVertex,
		                bool shortest_edge)
	{
		Curvature& curvature = hVertex->curvature();
		curvature.init();
		HV_circulator hHalfedge = hVertex->vertex_begin();
		HV_circulator begin = hHalfedge;
		FT min_length = m_pMesh->min_edge_length_around(hVertex);
		CGAL_For_all(hHalfedge,begin)
		{
			// exclude border case
			Facet_handle hFacet1 = hHalfedge->facet();
			Facet_handle hFacet2 = hHalfedge->opposite()->facet();
			if(hFacet1 == NULL || hFacet2 == NULL)
				continue; // border edge

			// build edge vector and normalize it
			const Point& p1 = hHalfedge->vertex()->point();
			const Point& p2 = hHalfedge->opposite()->vertex()->point();
			Vector edge = p1 - p2;
			FT len_edge = (FT)std::sqrt(edge*edge);
			if(len_edge == 0)
				continue;
			Vector normalized_edge = edge / len_edge;

			// compute signed angle between two incident facets
			// and sum up to current 3x3 tensor
			Vector normal1 = hFacet1->normal();
			Vector normal2 = hFacet2->normal();
			FT sinus = CGAL::cross_product(normal1,normal2)*normalized_edge;
			FT beta = (FT)arc_sinus((double)sinus);
			if(shortest_edge)
				add(curvature.tensor3(),beta*min_length,normalized_edge);
			else
				add(curvature.tensor3(),beta*len_edge,normalized_edge);
		}
		curvature.diagonalize();
		curvature.scale(0.75*min_length);
	}

	// beta * len(edge) * edge * edge_transpose
	//**********************************************
	void add(FT* tensor3,
		       const FT& coeff,
					 const Vector& e)
	{
		tensor3[0] += coeff * e.x()*e.x();
		tensor3[1] += coeff * e.y()*e.x();
		tensor3[2] += coeff * e.y()*e.y();
		tensor3[3] += coeff * e.z()*e.x();
		tensor3[4] += coeff * e.z()*e.y();
		tensor3[5] += coeff * e.z()*e.z();
	}

	// arc sinus
	//**********************************************
	double arc_sinus(double sinus)
	{
		if(sinus >= 1.0)
			return HALF_PI;
		if(sinus <= -1.0)
			return -HALF_PI;
		return std::asin(sinus);
	}
};

#endif // CURVATURE_ESTIMATOR_H
