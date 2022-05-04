# include "tsunami.h"

#define CACHE_SIZE	240
#define SIN sin
#define COS cos

struct GLUquadric {
    GLint	normals;
    GLboolean	textureCoords;
    GLint	orientation;
    GLint	drawStyle;
    void	(*errorCallback)( GLint );
};

GLUquadric *gluNewQuadric(void)
{
    GLUquadric *newstate;

    newstate = (GLUquadric *) malloc(sizeof(GLUquadric));
    if (newstate == NULL) {
	/* Can't report an error at this point... */
	return NULL;
    }
    newstate->normals = GLU_SMOOTH;
    newstate->textureCoords = GL_FALSE;
    newstate->orientation = GLU_OUTSIDE;
    newstate->drawStyle = GLU_FILL;
    newstate->errorCallback = NULL;
    return newstate;
}

void gluQuadricNormals(GLUquadric *qobj, GLenum normals)
{
    switch (normals) {
      case GLU_SMOOTH:
      case GLU_FLAT:
      case GLU_NONE:
	break;
      default:
//	gluQuadricError(qobj, GLU_INVALID_ENUM);
	return;
    }
    qobj->normals = normals;
}

static void __gluMakeIdentityd(GLdouble m[16])
{
    m[0+4*0] = 1; m[0+4*1] = 0; m[0+4*2] = 0; m[0+4*3] = 0;
    m[1+4*0] = 0; m[1+4*1] = 1; m[1+4*2] = 0; m[1+4*3] = 0;
    m[2+4*0] = 0; m[2+4*1] = 0; m[2+4*2] = 1; m[2+4*3] = 0;
    m[3+4*0] = 0; m[3+4*1] = 0; m[3+4*2] = 0; m[3+4*3] = 1;
}

static void __gluMakeIdentityf(GLfloat m[16])
{
    m[0+4*0] = 1; m[0+4*1] = 0; m[0+4*2] = 0; m[0+4*3] = 0;
    m[1+4*0] = 0; m[1+4*1] = 1; m[1+4*2] = 0; m[1+4*3] = 0;
    m[2+4*0] = 0; m[2+4*1] = 0; m[2+4*2] = 1; m[2+4*3] = 0;
    m[3+4*0] = 0; m[3+4*1] = 0; m[3+4*2] = 0; m[3+4*3] = 1;
}

static void normalize(float v[3])
{
    float r;

    r = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
    if (r == 0.0) return;

    v[0] /= r;
    v[1] /= r;
    v[2] /= r;
}

static void cross(float v1[3], float v2[3], float result[3])
{
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
}


void gluPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
    GLdouble m[4][4];
    double sine, cotangent, deltaZ;
    double radians = fovy / 2 * PI / 180;

    deltaZ = zFar - zNear;
    sine = sin(radians);
    if ((deltaZ == 0) || (sine == 0) || (aspect == 0)) {
	return;
    }
    cotangent = cos(radians) / sine;

    __gluMakeIdentityd(&m[0][0]);
    m[0][0] = cotangent / aspect;
    m[1][1] = cotangent;
    m[2][2] = -(zFar + zNear) / deltaZ;
    m[2][3] = -1;
    m[3][2] = -2 * zNear * zFar / deltaZ;
    m[3][3] = 0;
    glMultMatrixd(&m[0][0]);
}

void gluLookAt(GLdouble eyex, GLdouble eyey, GLdouble eyez, GLdouble centerx,
	  GLdouble centery, GLdouble centerz, GLdouble upx, GLdouble upy,
	  GLdouble upz)
{
    float forward[3], side[3], up[3];
    GLfloat m[4][4];

    forward[0] = centerx - eyex;
    forward[1] = centery - eyey;
    forward[2] = centerz - eyez;

    up[0] = upx;
    up[1] = upy;
    up[2] = upz;

    normalize(forward);

    /* Side = forward x up */
    cross(forward, up, side);
    normalize(side);

    /* Recompute up as: up = side x forward */
    cross(side, forward, up);

    __gluMakeIdentityf(&m[0][0]);
    m[0][0] = side[0];
    m[1][0] = side[1];
    m[2][0] = side[2];

    m[0][1] = up[0];
    m[1][1] = up[1];
    m[2][1] = up[2];

    m[0][2] = -forward[0];
    m[1][2] = -forward[1];
    m[2][2] = -forward[2];

    glMultMatrixf(&m[0][0]);
    glTranslated(-eyex, -eyey, -eyez);
}


void gluSphere(GLUquadric *qobj, GLdouble radius, GLint slices, GLint stacks)
{
    GLint i,j;
    GLfloat sinCache1a[CACHE_SIZE];
    GLfloat cosCache1a[CACHE_SIZE];
    GLfloat sinCache2a[CACHE_SIZE];
    GLfloat cosCache2a[CACHE_SIZE];
    GLfloat sinCache3a[CACHE_SIZE];
    GLfloat cosCache3a[CACHE_SIZE];
    GLfloat sinCache1b[CACHE_SIZE];
    GLfloat cosCache1b[CACHE_SIZE];
    GLfloat sinCache2b[CACHE_SIZE];
    GLfloat cosCache2b[CACHE_SIZE];
    GLfloat sinCache3b[CACHE_SIZE];
    GLfloat cosCache3b[CACHE_SIZE];
    GLfloat angle;
    GLfloat zLow, zHigh;
    GLfloat sintemp1 = 0.0, sintemp2 = 0.0, sintemp3 = 0.0, sintemp4 = 0.0;
    GLfloat costemp1 = 0.0, costemp2 = 0.0, costemp3 = 0.0, costemp4 = 0.0;
    GLboolean needCache2, needCache3;
    GLint start, finish;

    if (slices >= CACHE_SIZE) slices = CACHE_SIZE-1;
    if (stacks >= CACHE_SIZE) stacks = CACHE_SIZE-1;
    if (slices < 2 || stacks < 1 || radius < 0.0) {
//	gluQuadricError(qobj, GLU_INVALID_VALUE);
	return;
    }

    /* Cache is the vertex locations cache */
    /* Cache2 is the various normals at the vertices themselves */
    /* Cache3 is the various normals for the faces */
    needCache2 = needCache3 = GL_FALSE;

    if (qobj->normals == GLU_SMOOTH) {
	needCache2 = GL_TRUE;
    }

    if (qobj->normals == GLU_FLAT) {
	if (qobj->drawStyle != GLU_POINT) {
	    needCache3 = GL_TRUE;
	}
	if (qobj->drawStyle == GLU_LINE) {
	    needCache2 = GL_TRUE;
	}
    }

    for (i = 0; i < slices; i++) {
	angle = 2 * PI * i / slices;
	sinCache1a[i] = SIN(angle);
	cosCache1a[i] = COS(angle);
	if (needCache2) {
	    sinCache2a[i] = sinCache1a[i];
	    cosCache2a[i] = cosCache1a[i];
	}
    }

    for (j = 0; j <= stacks; j++) {
	angle = PI * j / stacks;
	if (needCache2) {
	    if (qobj->orientation == GLU_OUTSIDE) {
		sinCache2b[j] = SIN(angle);
		cosCache2b[j] = COS(angle);
	    } else {
		sinCache2b[j] = -SIN(angle);
		cosCache2b[j] = -COS(angle);
	    }
	}
	sinCache1b[j] = radius * SIN(angle);
	cosCache1b[j] = radius * COS(angle);
    }
    /* Make sure it comes to a point */
    sinCache1b[0] = 0;
    sinCache1b[stacks] = 0;

    if (needCache3) {
	for (i = 0; i < slices; i++) {
	    angle = 2 * PI * (i-0.5) / slices;
	    sinCache3a[i] = SIN(angle);
	    cosCache3a[i] = COS(angle);
	}
	for (j = 0; j <= stacks; j++) {
	    angle = PI * (j - 0.5) / stacks;
	    if (qobj->orientation == GLU_OUTSIDE) {
		sinCache3b[j] = SIN(angle);
		cosCache3b[j] = COS(angle);
	    } else {
		sinCache3b[j] = -SIN(angle);
		cosCache3b[j] = -COS(angle);
	    }
	}
    }

    sinCache1a[slices] = sinCache1a[0];
    cosCache1a[slices] = cosCache1a[0];
    if (needCache2) {
	sinCache2a[slices] = sinCache2a[0];
	cosCache2a[slices] = cosCache2a[0];
    }
    if (needCache3) {
	sinCache3a[slices] = sinCache3a[0];
	cosCache3a[slices] = cosCache3a[0];
    }

    switch (qobj->drawStyle) {
      case GLU_FILL:
	/* Do ends of sphere as TRIANGLE_FAN's (if not texturing)
	** We don't do it when texturing because we need to respecify the
	** texture coordinates of the apex for every adjacent vertex (because
	** it isn't a constant for that point)
	*/
	if (!(qobj->textureCoords)) {
	    start = 1;
	    finish = stacks - 1;

	    /* Low end first (j == 0 iteration) */
	    sintemp2 = sinCache1b[1];
	    zHigh = cosCache1b[1];
	    switch(qobj->normals) {
	      case GLU_FLAT:
		sintemp3 = sinCache3b[1];
		costemp3 = cosCache3b[1];
		break;
	      case GLU_SMOOTH:
		sintemp3 = sinCache2b[1];
		costemp3 = cosCache2b[1];
		glNormal3f(sinCache2a[0] * sinCache2b[0],
			cosCache2a[0] * sinCache2b[0],
			cosCache2b[0]);
		break;
	      default:
		break;
	    }
	    glBegin(GL_TRIANGLE_FAN);
	    glVertex3f(0.0, 0.0, radius);
	    if (qobj->orientation == GLU_OUTSIDE) {
		for (i = slices; i >= 0; i--) {
		    switch(qobj->normals) {
		      case GLU_SMOOTH:
			glNormal3f(sinCache2a[i] * sintemp3,
				cosCache2a[i] * sintemp3,
				costemp3);
			break;
		      case GLU_FLAT:
			if (i != slices) {
			    glNormal3f(sinCache3a[i+1] * sintemp3,
				    cosCache3a[i+1] * sintemp3,
				    costemp3);
			}
			break;
		      case GLU_NONE:
		      default:
			break;
		    }
		    glVertex3f(sintemp2 * sinCache1a[i],
			    sintemp2 * cosCache1a[i], zHigh);
		}
	    } else {
		for (i = 0; i <= slices; i++) {
		    switch(qobj->normals) {
		      case GLU_SMOOTH:
			glNormal3f(sinCache2a[i] * sintemp3,
				cosCache2a[i] * sintemp3,
				costemp3);
			break;
		      case GLU_FLAT:
			glNormal3f(sinCache3a[i] * sintemp3,
				cosCache3a[i] * sintemp3,
				costemp3);
			break;
		      case GLU_NONE:
		      default:
			break;
		    }
		    glVertex3f(sintemp2 * sinCache1a[i],
			    sintemp2 * cosCache1a[i], zHigh);
		}
	    }
	    glEnd();

	    /* High end next (j == stacks-1 iteration) */
	    sintemp2 = sinCache1b[stacks-1];
	    zHigh = cosCache1b[stacks-1];
	    switch(qobj->normals) {
	      case GLU_FLAT:
		sintemp3 = sinCache3b[stacks];
		costemp3 = cosCache3b[stacks];
		break;
	      case GLU_SMOOTH:
		sintemp3 = sinCache2b[stacks-1];
		costemp3 = cosCache2b[stacks-1];
		glNormal3f(sinCache2a[stacks] * sinCache2b[stacks],
			cosCache2a[stacks] * sinCache2b[stacks],
			cosCache2b[stacks]);
		break;
	      default:
		break;
	    }
	    glBegin(GL_TRIANGLE_FAN);
	    glVertex3f(0.0, 0.0, -radius);
	    if (qobj->orientation == GLU_OUTSIDE) {
		for (i = 0; i <= slices; i++) {
		    switch(qobj->normals) {
		      case GLU_SMOOTH:
			glNormal3f(sinCache2a[i] * sintemp3,
				cosCache2a[i] * sintemp3,
				costemp3);
			break;
		      case GLU_FLAT:
			glNormal3f(sinCache3a[i] * sintemp3,
				cosCache3a[i] * sintemp3,
				costemp3);
			break;
		      case GLU_NONE:
		      default:
			break;
		    }
		    glVertex3f(sintemp2 * sinCache1a[i],
			    sintemp2 * cosCache1a[i], zHigh);
		}
	    } else {
		for (i = slices; i >= 0; i--) {
		    switch(qobj->normals) {
		      case GLU_SMOOTH:
			glNormal3f(sinCache2a[i] * sintemp3,
				cosCache2a[i] * sintemp3,
				costemp3);
			break;
		      case GLU_FLAT:
			if (i != slices) {
			    glNormal3f(sinCache3a[i+1] * sintemp3,
				    cosCache3a[i+1] * sintemp3,
				    costemp3);
			}
			break;
		      case GLU_NONE:
		      default:
			break;
		    }
		    glVertex3f(sintemp2 * sinCache1a[i],
			    sintemp2 * cosCache1a[i], zHigh);
		}
	    }
	    glEnd();
	} else {
	    start = 0;
	    finish = stacks;
	}
	for (j = start; j < finish; j++) {
	    zLow = cosCache1b[j];
	    zHigh = cosCache1b[j+1];
	    sintemp1 = sinCache1b[j];
	    sintemp2 = sinCache1b[j+1];
	    switch(qobj->normals) {
	      case GLU_FLAT:
		sintemp4 = sinCache3b[j+1];
		costemp4 = cosCache3b[j+1];
		break;
	      case GLU_SMOOTH:
		if (qobj->orientation == GLU_OUTSIDE) {
		    sintemp3 = sinCache2b[j+1];
		    costemp3 = cosCache2b[j+1];
		    sintemp4 = sinCache2b[j];
		    costemp4 = cosCache2b[j];
		} else {
		    sintemp3 = sinCache2b[j];
		    costemp3 = cosCache2b[j];
		    sintemp4 = sinCache2b[j+1];
		    costemp4 = cosCache2b[j+1];
		}
		break;
	      default:
		break;
	    }

	    glBegin(GL_QUAD_STRIP);
	    for (i = 0; i <= slices; i++) {
		switch(qobj->normals) {
		  case GLU_SMOOTH:
		    glNormal3f(sinCache2a[i] * sintemp3,
			    cosCache2a[i] * sintemp3,
			    costemp3);
		    break;
		  case GLU_FLAT:
		  case GLU_NONE:
		  default:
		    break;
		}
		if (qobj->orientation == GLU_OUTSIDE) {
		    if (qobj->textureCoords) {
			glTexCoord2f(1 - (float) i / slices,
				1 - (float) (j+1) / stacks);
		    }
		    glVertex3f(sintemp2 * sinCache1a[i],
			    sintemp2 * cosCache1a[i], zHigh);
		} else {
		    if (qobj->textureCoords) {
			glTexCoord2f(1 - (float) i / slices,
				1 - (float) j / stacks);
		    }
		    glVertex3f(sintemp1 * sinCache1a[i],
			    sintemp1 * cosCache1a[i], zLow);
		}
		switch(qobj->normals) {
		  case GLU_SMOOTH:
		    glNormal3f(sinCache2a[i] * sintemp4,
			    cosCache2a[i] * sintemp4,
			    costemp4);
		    break;
		  case GLU_FLAT:
		    glNormal3f(sinCache3a[i] * sintemp4,
			    cosCache3a[i] * sintemp4,
			    costemp4);
		    break;
		  case GLU_NONE:
		  default:
		    break;
		}
		if (qobj->orientation == GLU_OUTSIDE) {
		    if (qobj->textureCoords) {
			glTexCoord2f(1 - (float) i / slices,
				1 - (float) j / stacks);
		    }
		    glVertex3f(sintemp1 * sinCache1a[i],
			    sintemp1 * cosCache1a[i], zLow);
		} else {
		    if (qobj->textureCoords) {
			glTexCoord2f(1 - (float) i / slices,
				1 - (float) (j+1) / stacks);
		    }
		    glVertex3f(sintemp2 * sinCache1a[i],
			    sintemp2 * cosCache1a[i], zHigh);
		}
	    }
	    glEnd();
	}
	break;
      case GLU_POINT:
	glBegin(GL_POINTS);
	for (j = 0; j <= stacks; j++) {
	    sintemp1 = sinCache1b[j];
	    costemp1 = cosCache1b[j];
	    switch(qobj->normals) {
	      case GLU_FLAT:
	      case GLU_SMOOTH:
		sintemp2 = sinCache2b[j];
		costemp2 = cosCache2b[j];
		break;
	      default:
		break;
	    }
	    for (i = 0; i < slices; i++) {
		switch(qobj->normals) {
		  case GLU_FLAT:
		  case GLU_SMOOTH:
		    glNormal3f(sinCache2a[i] * sintemp2,
			    cosCache2a[i] * sintemp2,
			    costemp2);
		    break;
		  case GLU_NONE:
		  default:
		    break;
		}

		zLow = j * radius / stacks;

		if (qobj->textureCoords) {
		    glTexCoord2f(1 - (float) i / slices,
			    1 - (float) j / stacks);
		}
		glVertex3f(sintemp1 * sinCache1a[i],
			sintemp1 * cosCache1a[i], costemp1);
	    }
	}
	glEnd();
	break;
      case GLU_LINE:
      case GLU_SILHOUETTE:
	for (j = 1; j < stacks; j++) {
	    sintemp1 = sinCache1b[j];
	    costemp1 = cosCache1b[j];
	    switch(qobj->normals) {
	      case GLU_FLAT:
	      case GLU_SMOOTH:
		sintemp2 = sinCache2b[j];
		costemp2 = cosCache2b[j];
		break;
	      default:
		break;
	    }

	    glBegin(GL_LINE_STRIP);
	    for (i = 0; i <= slices; i++) {
		switch(qobj->normals) {
		  case GLU_FLAT:
		    glNormal3f(sinCache3a[i] * sintemp2,
			    cosCache3a[i] * sintemp2,
			    costemp2);
		    break;
		  case GLU_SMOOTH:
		    glNormal3f(sinCache2a[i] * sintemp2,
			    cosCache2a[i] * sintemp2,
			    costemp2);
		    break;
		  case GLU_NONE:
		  default:
		    break;
		}
		if (qobj->textureCoords) {
		    glTexCoord2f(1 - (float) i / slices,
			    1 - (float) j / stacks);
		}
		glVertex3f(sintemp1 * sinCache1a[i],
			sintemp1 * cosCache1a[i], costemp1);
	    }
	    glEnd();
	}
	for (i = 0; i < slices; i++) {
	    sintemp1 = sinCache1a[i];
	    costemp1 = cosCache1a[i];
	    switch(qobj->normals) {
	      case GLU_FLAT:
	      case GLU_SMOOTH:
		sintemp2 = sinCache2a[i];
		costemp2 = cosCache2a[i];
		break;
	      default:
		break;
	    }

	    glBegin(GL_LINE_STRIP);
	    for (j = 0; j <= stacks; j++) {
		switch(qobj->normals) {
		  case GLU_FLAT:
		    glNormal3f(sintemp2 * sinCache3b[j],
			    costemp2 * sinCache3b[j],
			    cosCache3b[j]);
		    break;
		  case GLU_SMOOTH:
		    glNormal3f(sintemp2 * sinCache2b[j],
			    costemp2 * sinCache2b[j],
			    cosCache2b[j]);
		    break;
		  case GLU_NONE:
		  default:
		    break;
		}

		if (qobj->textureCoords) {
		    glTexCoord2f(1 - (float) i / slices,
			    1 - (float) j / stacks);
		}
		glVertex3f(sintemp1 * sinCache1b[j],
			costemp1 * sinCache1b[j], cosCache1b[j]);
	    }
	    glEnd();
	}
	break;
      default:
	break;
    }
}