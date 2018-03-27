// Application4.h: interface for the Application4 class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_APPLICATION5_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_)
#define AFX_APPLICATION5_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Application.h"


// AA kernel
#define AAKERNEL_SIZE 6

static const float AAFilter[AAKERNEL_SIZE][3] =
{
	{ -0.52, 0.38, 0.128 },
	{ 0.41, 0.56, 0.119 },
	{ 0.27, 0.08, 0.294 },
	{ -0.17, -0.29, 0.249 },
	{ 0.58, -0.55, 0.104 },
	{ -0.31, -0.71, 0.106 },
};

inline int ARRAY(int x, int y, int xres) { return (x + y*xres); }	/* simplify fbuf indexing */

class Application5 : public Application  
{
public:
	Application5();
	virtual ~Application5();
	GzRender* aaRenders[AAKERNEL_SIZE];	// I'll use these to compute the offset aa subimage pixel colors and then use the mp_render data memeber 
										// to generate the actual final image.

	int	Initialize();
	virtual int Render(); 
	int Clean();
};

#endif // !defined(AFX_APPLICATION5_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_)
