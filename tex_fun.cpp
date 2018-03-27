/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

#define arrIndex(x, y) (x+y*xs)

template< typename type >
inline type clampvalue(type val, type min, type max)
{
	return (val > max) ? max : (val < min ? min : val);
}

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */

  u = clampvalue<float>(u, 0, 1);
  v = clampvalue<float>(v, 0, 1);
  u *= (xs - 1);
  v *= (ys - 1);
  int u1 = floor(u), u2 = ceil(u), v1 = floor(v), v2 = ceil(v);
  float s = u - u1,
  t = v - v1;

  // A = {u1, v1},
  // B = {u2, v1},
  // C = {u1, v2},
  // D = {u2, v2};

  // Read in pixel color values
  GzColor Acol = { image[arrIndex(u1, v1)][RED], image[arrIndex(u1, v1)][GREEN] , image[arrIndex(u1, v1)][BLUE] },
	  Bcol = { image[arrIndex(u2, v1)][RED], image[arrIndex(u2, v1)][GREEN] , image[arrIndex(u2, v1)][BLUE] },
	  Dcol = { image[arrIndex(u1, v2)][RED], image[arrIndex(u1, v2)][GREEN] , image[arrIndex(u1, v2)][BLUE] },
	  Ccol = { image[arrIndex(u2, v2)][RED], image[arrIndex(u2, v2)][GREEN] , image[arrIndex(u2, v2)][BLUE] };

  // Bilinear interpolation
  color[RED] = s * t * Ccol[RED] + (1 - s) * t * Dcol[RED] + s * (1 - t) * Bcol[RED] + (1 - s) * (1 - t) * Acol[RED];
  color[GREEN] = s * t * Ccol[GREEN] + (1 - s) * t * Dcol[GREEN] + s * (1 - t) * Bcol[GREEN] + (1 - s) * (1 - t) * Acol[GREEN];
  color[BLUE] = s * t * Ccol[BLUE] + (1 - s) * t * Dcol[BLUE] + s * (1 - t) * Bcol[BLUE] + (1 - s) * (1 - t) * Acol[BLUE];

	return GZ_SUCCESS;
}

void compSquare(float complexX[2])
{
	complexX[0] = (complexX[0] * complexX[0]) - (complexX[1] * complexX[1]);
	complexX[1] = (-2) * complexX[0] * complexX[1];
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{

	u = clampvalue<float>(u, 0, 1);
	v = clampvalue<float>(v, 0, 1);

	// Same color at corner square elements 
	int vf = floor(v * 10);
	int uf = floor(u * 10);
	if (uf > 5)
	{
		uf = 10 - uf;
	}

	if (vf > 5)
	{
		vf = 10 - vf;
	}

	color[RED] = uf;
	color[RED] /= 10;
	color[GREEN] = vf * uf;
	color[GREEN] /= 10;
	color[GREEN] = clampvalue<float>(color[GREEN], 0, 1) / 2;
	color[BLUE] = vf;
	color[BLUE] /= 10;

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}



/* ------------------------------------------------------------------------------------------------------------------------------------- */
/* Image texture function */
int toonTex_fun(float u, float v, GzColor color)
{
	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE			*fd;

	if (reset) {          /* open and load texture file */
		// fd = fopen("fadingTexture.ppm", "rb");
		if (color[RED] == 3)
		{
			fd = fopen("silhouetteAbstractionTexture.ppm", "rb");
		}
		else if (color[RED] == 4)
		{
			fd = fopen("backlightingTexture.ppm", "rb");
		}
		else if (color[RED] == 5)
		{
			fd = fopen("opacityTexture.ppm", "rb");
		}
		else
		{
			fd = fopen("anotherTexture.ppm", "rb");
		}
		if (fd == NULL) {
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}
		fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
		image = (GzColor*)malloc(sizeof(GzColor)*(xs + 1)*(ys + 1));
		if (image == NULL) {
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
			image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
			image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		}

		reset = 0;          /* init is done */
		fclose(fd);
	}

	/* bounds-test u,v to make sure nothing will overflow image array bounds */
	/* determine texture cell corner values and perform bilinear interpolation */
	/* set color to interpolated GzColor value and return */

	u = clampvalue<float>(u, 0, 1);
	v = clampvalue<float>(v, 0, 1);
	u *= (xs - 1);
	v *= (ys - 1);
	int u1 = floor(u), u2 = ceil(u), v1 = floor(v), v2 = ceil(v);
	float s = u - u1,
		  t = v - v1;

	// Read in pixel color values
	GzColor Acol = { image[arrIndex(u1, v1)][RED], image[arrIndex(u1, v1)][GREEN] , image[arrIndex(u1, v1)][BLUE] },
		Bcol = { image[arrIndex(u2, v1)][RED], image[arrIndex(u2, v1)][GREEN] , image[arrIndex(u2, v1)][BLUE] },
		Dcol = { image[arrIndex(u1, v2)][RED], image[arrIndex(u1, v2)][GREEN] , image[arrIndex(u1, v2)][BLUE] },
		Ccol = { image[arrIndex(u2, v2)][RED], image[arrIndex(u2, v2)][GREEN] , image[arrIndex(u2, v2)][BLUE] };

	// Bilinear interpolation
	color[RED] = s * t * Ccol[RED] + (1 - s) * t * Dcol[RED] + s * (1 - t) * Bcol[RED] + (1 - s) * (1 - t) * Acol[RED];
	color[GREEN] = s * t * Ccol[GREEN] + (1 - s) * t * Dcol[GREEN] + s * (1 - t) * Bcol[GREEN] + (1 - s) * (1 - t) * Acol[GREEN];
	color[BLUE] = s * t * Ccol[BLUE] + (1 - s) * t * Dcol[BLUE] + s * (1 - t) * Bcol[BLUE] + (1 - s) * (1 - t) * Acol[BLUE];

	//color[RED] = image[arrIndex(u1, v1)][RED];
	//color[GREEN] = image[arrIndex(u1, v1)][GREEN];
	//color[BLUE] = image[arrIndex(u1, v1)][BLUE];

	return GZ_SUCCESS;
}
