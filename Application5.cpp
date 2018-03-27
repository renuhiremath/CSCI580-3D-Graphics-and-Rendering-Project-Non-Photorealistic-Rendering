// Application5.cpp: implementation of the Application5 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment #5
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "Application5.h"
#include "Gz.h"
#include "rend.h"
#include "GeneralHelpers.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "ppot.asc"
#define OUTFILE "output.ppm"


extern int tex_fun(float u, float v, GzColor color);		/* image texture function */
extern int ptex_fun(float u, float v, GzColor color);		/* procedural texture function */
extern int toonTex_fun(float u, float v, GzColor color);	/* lookup toon texture */
extern int GzFreeTexture();

void shade(GzCoord norm, GzCoord color);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


Application5::Application5()
{
	
}

Application5::~Application5()
{
	Clean();
}

int Application5::Initialize()
{
	GzCamera	camera;  
	int		    xRes, yRes;	/* display parameters */ 

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	GzToken		nameListAA[3];
	GzPointer	valueListAA[3];

	int			shaderType, interpStyle;
	float		specpower;
	int		status; 
 
	status = 0; 

	/* 
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/* 
	 * initialize the display and the renderer 
	 */ 
 	m_nWidth = 512;		// frame buffer and display width
	m_nHeight = 512;    // frame buffer and display heightZE];

	// Allocate mem and set default color for the final render
	m_pRender = new GzRender(m_nWidth, m_nHeight);
	m_pRender->GzDefault();
	m_pFrameBuffer = m_pRender->framebuffer;

	for (int rendIt = 0; rendIt < AAKERNEL_SIZE; ++rendIt)
	{
		aaRenders[rendIt] = new GzRender(m_nWidth, m_nHeight);
		aaRenders[rendIt]->GzDefault();

		/* Translation matrix */
		GzMatrix	scale =
		{
			3.25,	0.0,	0.0,	0.0,
			0.0,	3.25,	0.0,	-3.25,
			0.0,	0.0,	3.25,	3.5,
			0.0,	0.0,	0.0,	1.0
		};

		GzMatrix	rotateX =
		{
			1.0,	0.0,	0.0,	0.0,
			0.0,	.7071,	.7071,	0.0,
			0.0,	-.7071,	.7071,	0.0,
			0.0,	0.0,	0.0,	1.0
		};

		GzMatrix	rotateY =
		{
			.866,	0.0,	-0.5,	0.0,
			0.0,	1.0,	0.0,	0.0,
			0.5,	0.0,	.866,	0.0,
			0.0,	0.0,	0.0,	1.0
		};

#if 1 	/* set up app-defined camera if desired, else use camera defaults */
		camera.position[X] = -3;
		camera.position[Y] = -25;
		camera.position[Z] = -4;

		camera.lookat[X] = 7.8;
		camera.lookat[Y] = 0.7;
		camera.lookat[Z] = 6.5;

		camera.worldup[X] = -0.2;
		camera.worldup[Y] = 1.0;
		camera.worldup[Z] = 0.0;

		camera.FOV = 63.7;              /* degrees *              /* degrees */

		status |= aaRenders[rendIt]->GzPutCamera(camera);
#endif 

		/* Start Renderer */
		status |= aaRenders[rendIt]->GzBeginRender();

		/* Light */
		GzLight	light1 = { {-0.7071, 0.7071, 0}, {0.5, 0.5, 0.9} };
		GzLight	light2 = { {0, -0.7071, -0.7071}, {0.9, 0.2, 0.3} };
		GzLight	light3 = { {0.7071, 0.0, -0.7071}, {0.2, 0.7, 0.3} };
		GzLight	ambientlight = { {0, 0, 0}, {0.3, 0.3, 0.3} };

		/* Material property */
		GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
		GzColor ambientCoefficient = { 0.1, 0.1, 0.1 };
		GzColor diffuseCoefficient = { 0.7, 0.7, 0.7 };

		/*
		  renderer is ready for frame --- define lights and shader at start of frame
		*/

		/*
		 * Tokens associated with light parameters
		 */
		//nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
		//valueListLights[0] = (GzPointer)&light1;
		//nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
		//valueListLights[1] = (GzPointer)&light2;
		//nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
		//valueListLights[2] = (GzPointer)&light3;

		//just do one light for toon shading testing
		nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[0] = (GzPointer)&light1;

		status |= aaRenders[rendIt]->GzPutAttribute(1, nameListLights, valueListLights);

		nameListLights[0] = GZ_AMBIENT_LIGHT;
		valueListLights[0] = (GzPointer)&ambientlight;
		status |= aaRenders[rendIt]->GzPutAttribute(1, nameListLights, valueListLights);

		/*
		 * Tokens associated with shading
		 */
		nameListShader[0] = GZ_DIFFUSE_COEFFICIENT;
		valueListShader[0] = (GzPointer)diffuseCoefficient;

		/*
		* Select either GZ_COLOR or GZ_NORMALS as interpolation mode
		*/
		nameListShader[1] = GZ_INTERPOLATE;
		// interpStyle = GZ_COLOR;         /* Gouraud shading */
		interpStyle = GZ_NORMALS;         /* Phong shading */
		valueListShader[1] = (GzPointer)&interpStyle;

		nameListShader[2] = GZ_AMBIENT_COEFFICIENT;
		valueListShader[2] = (GzPointer)ambientCoefficient;
		nameListShader[3] = GZ_SPECULAR_COEFFICIENT;
		valueListShader[3] = (GzPointer)specularCoefficient;
		nameListShader[4] = GZ_DISTRIBUTION_COEFFICIENT;
		specpower = 32;
		valueListShader[4] = (GzPointer)&specpower;

		nameListShader[5] = GZ_TEXTURE_MAP;
#if 1  /* set up null texture function or valid pointer */
		valueListShader[5] = (GzPointer)0;
#else
		valueListShader[5] = (GzPointer)(tex_fun);	/* or use ptex_fun */
#endif
		nameListShader[6] = GZ_TOON_TEX_MAP;
		valueListShader[6] = (GzPointer)(toonTex_fun);
		status |= aaRenders[rendIt]->GzPutAttribute(7, nameListShader, valueListShader);

		// Put the AA shift attributes for each of the renderers
		nameListAA[0] = GZ_AASHIFTX;
		float dx = AAFilter[rendIt][X];
		valueListAA[0] = &dx;

		nameListAA[1] = GZ_AASHIFTY;
		float dy = AAFilter[rendIt][Y];
		valueListAA[1] = &dy;

		status |= aaRenders[rendIt]->GzPutAttribute(2, nameListAA, valueListAA);
		
		status |= aaRenders[rendIt]->GzPushMatrix(scale);
		status |= aaRenders[rendIt]->GzPushMatrix(rotateY);
		status |= aaRenders[rendIt]->GzPushMatrix(rotateX);

		if (status) exit(GZ_FAILURE);
	}
	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}

int Application5::Render() 
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */ 
	GzCoord		normalList[3];	/* vertex normals */ 
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */ 
	char		dummy[256]; 
	int			status; 


	/* Initialize Display */
	status |= m_pRender->GzDefault();  /* init for new frame */
	for (int rendIt = 0; rendIt < AAKERNEL_SIZE; ++rendIt)
	{
		status |= aaRenders[rendIt]->GzDefault();
	}
	
	/* 
	* Tokens associated with triangle vertex values 
	*/ 
	nameListTriangle[0] = GZ_POSITION; 
	nameListTriangle[1] = GZ_NORMAL; 
	nameListTriangle[2] = GZ_TEXTURE_INDEX;  

	// I/O File open
	FILE *infile;

	FILE *outfile;
	if( (outfile  = fopen( OUTFILE , "wb" )) == NULL )
	{
         AfxMessageBox( "The output file was not opened\n" );
		 return GZ_FAILURE;
	}

	/* 
	* Walk through the list of triangles, set color 
	* and render each triangle 
	*/ 
	int ii = 0;
	for (ii = 0; ii < AAKERNEL_SIZE; ++ii)
	{
		if ((infile = fopen(INFILE, "r")) == NULL)
		{
			AfxMessageBox("The input file was not opened\n");
			return GZ_FAILURE;
		}
		while( fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
			&(vertexList[0][0]), &(vertexList[0][1]),  
			&(vertexList[0][2]), 
			&(normalList[0][0]), &(normalList[0][1]), 	
			&(normalList[0][2]), 
			&(uvList[0][0]), &(uvList[0][1]) ); 
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
			&(vertexList[1][0]), &(vertexList[1][1]), 	
			&(vertexList[1][2]), 
			&(normalList[1][0]), &(normalList[1][1]), 	
			&(normalList[1][2]), 
			&(uvList[1][0]), &(uvList[1][1]) ); 
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
			&(vertexList[2][0]), &(vertexList[2][1]), 	
			&(vertexList[2][2]), 
			&(normalList[2][0]), &(normalList[2][1]), 	
			&(normalList[2][2]), 
			&(uvList[2][0]), &(uvList[2][1]) ); 

			/* 
			 * Set the value pointers to the first vertex of the 	
			 * triangle, then feed it to the renderer 
			 * NOTE: this sequence matches the nameList token sequence
			 */ 
			 valueListTriangle[0] = (GzPointer)vertexList; 
			 valueListTriangle[1] = (GzPointer)normalList; 
			 valueListTriangle[2] = (GzPointer)uvList; 

			 // m_pRender->GzPutTriangle(3, nameListTriangle, valueListTriangle); 
			 aaRenders[ii]->GzPutTriangle(3, nameListTriangle, valueListTriangle);
		}

		if (fclose(infile))
			AfxMessageBox(_T("The input file was not closed\n"));
	} 

	// Generate accumulated values for the final render
	// every render has a pixel buffer that can be that has the final color values for the pixels

	// Traverse through all the pixels
	for (int xit = 0; xit < m_nWidth; ++xit)
	{
		for (int yit = 0; yit < m_nHeight; ++yit)
		{
			m_pRender->pixelbuffer[ARRAY(xit, yit, m_nWidth)].red = 0;
			m_pRender->pixelbuffer[ARRAY(xit, yit, m_nWidth)].green = 0;
			m_pRender->pixelbuffer[ARRAY(xit, yit, m_nWidth)].blue = 0;

			// Traverse through all the renderers pixel buffers to compute the final weighted r,g,b components
			for (int renIt = 0; renIt < AAKERNEL_SIZE; ++renIt)
			{
				m_pRender->pixelbuffer[ARRAY(xit, yit, m_nWidth)].red += (AAFilter[renIt][2] * aaRenders[renIt]->pixelbuffer[ARRAY(xit, yit, m_nWidth)].red);		// weight * red component value in the renderer
				m_pRender->pixelbuffer[ARRAY(xit, yit, m_nWidth)].green += (AAFilter[renIt][2] * aaRenders[renIt]->pixelbuffer[ARRAY(xit, yit, m_nWidth)].green);	// weight * green component value in the renderer
				m_pRender->pixelbuffer[ARRAY(xit, yit, m_nWidth)].blue += (AAFilter[renIt][2] * aaRenders[renIt]->pixelbuffer[ARRAY(xit, yit, m_nWidth)].blue);
			}
		}
	} 
	
	if (m_pRender->currentToonType == ToonShadingType::Vanilla)
		m_pRender->OutlineForToonShading();

	// Set matlevel and matrix buffers for the final buffer
	m_pRender->matlevel = aaRenders[0]->matlevel;
	for (int it = 0; it <= m_pRender->matlevel; ++it)
	{
		CopyMatrix(aaRenders[0]->Ximage[it], m_pRender->Ximage[it]);
		CopyMatrix(aaRenders[0]->Xnorm[it], m_pRender->Xnorm[it]);
	}

	// Flush final renderer buffer to frame buffer or the screen
	m_pRender->GzFlushDisplay2File(outfile); 	
	m_pRender->GzFlushDisplay2FrameBuffer();	

	/* 
	 * Close file
	 */ 

	if( fclose( outfile ) )
      AfxMessageBox(_T( "The output file was not closed\n" ));
 
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Clean()
{
	/* 
	 * Clean up and exit 
	 */ 
	int	status = 0; 

	free(m_pRender);
	status |= GzFreeTexture();
	
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS);

	// Free all the renderers
	for (int i = 0; i < AAKERNEL_SIZE; ++i)
	{
		free(aaRenders[i]);
	}
}
