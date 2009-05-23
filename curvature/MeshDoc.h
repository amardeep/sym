// MeshDoc.h : interface of the CMeshDoc class
//

#pragma once

// CGAL
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include "lib/enriched_polyhedron.h"

typedef double number_type;
typedef CGAL::Simple_cartesian<number_type> Enriched_kernel;
typedef Enriched_polyhedron<Enriched_kernel,Enriched_items> Mesh;

class CMeshDoc : public CDocument
{
protected: // create from serialization only
	CMeshDoc();
	DECLARE_DYNCREATE(CMeshDoc)

// Attributes
public:

	Mesh *m_pMesh;
	int m_method_curvature;

// Operations
public:

	// status message
	void StatusMessage( char* fmt, ... );
	void ResetMeshProperties();
	void UpdateMeshProperties(bool update_component = false,
		                        bool update_boundary = false);

// Overrides
	public:
	virtual void Serialize(CArchive& ar);

// Implementation
public:
	virtual ~CMeshDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
	virtual BOOL OnSaveDocument(LPCTSTR lpszPathName);
	afx_msg void OnSubdivisionQuad();
	afx_msg void OnUpdateSubdivisionQuad(CCmdUI *pCmdUI);
	afx_msg void OnEstimationOnering();
	afx_msg void OnUpdateEstimationOnering(CCmdUI *pCmdUI);
};


