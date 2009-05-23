/***************************************************************************
curvature.h  -  curvature information for CGAL polyhedron
----------------------------------------------------------------------------
begin                : nov 2003
copyright            : (C) 2003 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
***************************************************************************/
#ifndef CURVATURE_H
#define CURVATURE_H

#include "eigen.h"

template <class kernel>
class CCurvature
{
	typedef typename kernel::FT FT;
	typedef typename kernel::Vector_3 Vector;

	enum shape_type {PLANAR,
	                 SPHERICAL,
	                 PARABOLIC,
	                 ELLIPTIC,
	                 HYPERBOLIC};

private:

	// compact storage for 3x3 symmetric tensor
	FT m_pTensor3[6];
	Vector m_direction_kmin;
	Vector m_direction_kmax;
	FT m_kmin;
	FT m_kmax;
	FT m_scale;

// life cycle
public:
	CCurvature()
	{
		init();
	}
	~CCurvature() {}

// data access
public:

	FT* tensor3() { return m_pTensor3; }
	FT tensor3(const int& i) { return m_pTensor3[i]; }
	void tensor3(const int& i,FT value) { m_pTensor3[i] = value; }

	void kmin(const FT& kmin) { m_kmin = kmin; }
	const FT kmin() const { return m_kmin; }
	FT& kmin() { return m_kmin; }

	void kmax(const FT& kmax) { m_kmax = kmax; }
	const FT kmax() const { return m_kmax; }
	FT& kmax() { return m_kmax; }

	void scale(const FT& scale) { m_scale = scale; }
	const FT scale() const { return m_scale; }
	FT& scale() { return m_scale; }

public:
	void init()
	{
		m_pTensor3[0] = m_pTensor3[1] = 
		m_pTensor3[2] = m_pTensor3[3] = 
		m_pTensor3[4] = m_pTensor3[5] = (FT)0.0;
	}

	shape_type shape()
	{
		if(m_kmax * m_kmin >= 0)
			return ELLIPTIC;
		else
			return HYPERBOLIC;
	}

	void trace_tensor()
	{
		TRACE("%g\t%g\t%g\t%g\t%g\t%g\n",m_pTensor3[0],m_pTensor3[1],
			                               m_pTensor3[2],m_pTensor3[3],
																		 m_pTensor3[4],m_pTensor3[5]);
	}

	void diagonalize()
	{
		// extract eigenvalues and eigenvectors
		double eigen_values[3];  // will be given in decreasing order
		double eigen_vectors[9]; // the same
		eigen::semi_definite_symmetric((double *)m_pTensor3,3,
			                             eigen_vectors,eigen_values);

		// reorder them
		int indices[3] = {0,1,2};
		reorder(eigen_values,indices);

		// set member data
		m_kmin = (FT)eigen_values[indices[2]];
		m_kmax = (FT)eigen_values[indices[1]];
		m_direction_kmin = Vector((FT)eigen_vectors[3*indices[2]],
			                        (FT)eigen_vectors[3*indices[2]+1],
															(FT)eigen_vectors[3*indices[2]+2]);
		m_direction_kmax = Vector((FT)eigen_vectors[3*indices[1]],
			                        (FT)eigen_vectors[3*indices[1]+1],
															(FT)eigen_vectors[3*indices[1]+2]);
	}

	// draw kmin
	void gl_draw_kmin(bool exclude_asymptotic)
	{
		if(exclude_asymptotic && shape() == HYPERBOLIC)
			return;
		double x1 = (double) ( m_scale * m_direction_kmin.x());
		double y1 = (double) ( m_scale * m_direction_kmin.y());
		double z1 = (double) ( m_scale * m_direction_kmin.z());
		double x2 = (double) (-m_scale * m_direction_kmin.x());
		double y2 = (double) (-m_scale * m_direction_kmin.y());
		double z2 = (double) (-m_scale * m_direction_kmin.z());
		::glBegin(GL_LINES);
			::glVertex3d(x1,y1,z1);
			::glVertex3d(x2,y2,z2);
		::glEnd();
	}

	// draw kmax
	void gl_draw_kmax(bool exclude_asymptotic)
	{
		if(exclude_asymptotic && shape() == HYPERBOLIC)
			return;
		double x1 = (double) ( m_scale * m_direction_kmax.x());
		double y1 = (double) ( m_scale * m_direction_kmax.y());
		double z1 = (double) ( m_scale * m_direction_kmax.z());
		double x2 = (double) (-m_scale * m_direction_kmax.x());
		double y2 = (double) (-m_scale * m_direction_kmax.y());
		double z2 = (double) (-m_scale * m_direction_kmax.z());
		::glBegin(GL_LINES);
			::glVertex3d(x1,y1,z1);
			::glVertex3d(x2,y2,z2);
		::glEnd();
	}

	// draw asymptotic directions
	void gl_draw_asymptotic()
	{
		if(shape() != HYPERBOLIC)
			return;

		// two asymptotic directions
		Vector v1 = m_kmin * m_direction_kmax + 
			          m_kmax * m_direction_kmin;
		Vector v2 = m_kmin * m_direction_kmax -
			          m_kmax * m_direction_kmin;

		// normalize
		v1 = v1 / std::sqrt(v1*v1);
		v2 = v2 / std::sqrt(v2*v2);

		double x11 = (double) ( m_scale * v1.x());
		double y11 = (double) ( m_scale * v1.y());
		double z11 = (double) ( m_scale * v1.z());
		double x12 = (double) (-m_scale * v1.x());
		double y12 = (double) (-m_scale * v1.y());
		double z12 = (double) (-m_scale * v1.z());

		double x21 = (double) ( m_scale * v2.x());
		double y21 = (double) ( m_scale * v2.y());
		double z21 = (double) ( m_scale * v2.z());
		double x22 = (double) (-m_scale * v2.x());
		double y22 = (double) (-m_scale * v2.y());
		double z22 = (double) (-m_scale * v2.z());
		::glBegin(GL_LINES);
			::glVertex3d(x11,y11,z11);
			::glVertex3d(x12,y12,z12);
			::glVertex3d(x21,y21,z21);
			::glVertex3d(x22,y22,z22);
		::glEnd();
	}

	// reorder: normal, kmin, kmax
	void reorder(double eigen_values[3],
		           int indices[3])
	{
		// extract min abs value
		double abs_eigen[3];
		abs_eigen[0] = std::fabs(eigen_values[0]);
		abs_eigen[1] = std::fabs(eigen_values[1]);
		abs_eigen[2] = std::fabs(eigen_values[2]);

		// min_abs = 0
		if(abs_eigen[0] < abs_eigen[1] && abs_eigen[0] < abs_eigen[2])
		{
			indices[0] = 0;
			if(eigen_values[1]	< eigen_values[2])
			{
				indices[1] = 1;
				indices[2] = 2;
			}
			else
			{
				indices[1] = 2;
				indices[2] = 1;
			}
		}
		else // min_abs = 1
			if(abs_eigen[1] < abs_eigen[0] && abs_eigen[1] < abs_eigen[2])
			{
				indices[0] = 1;
				if(eigen_values[0]	< eigen_values[2])
				{
					indices[1] = 0;
					indices[2] = 2;
				}
				else
				{
					indices[1] = 2;
					indices[2] = 0;
				}
			}
			else // min_abs = 2
			{
				indices[0] = 2;
				if(eigen_values[0]	< eigen_values[1])
				{
					indices[1] = 0;
					indices[2] = 1;
				}
				else
				{
					indices[1] = 1;
					indices[2] = 0;
				}
			}
	}
};

#endif // CURVATURE_H


