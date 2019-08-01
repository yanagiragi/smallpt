#ifndef _DISPLAY_HPP
#define _DISPLAY_HPP


#ifdef _MSC_VER
	#include <GL\glut.h>
#else
	#include <GL/freeglut.h>
#endif // _MSC_VER

#include "Ray.hpp"
#include "Vec.hpp"
#include "Utils.hpp"

#include "Config.hpp"

namespace smallPT {
	
	void idleFunc(void) {
		if (currentSpp < limitSpp) {
			UpdateRendering();
		}
		glutPostRedisplay();
	}

	void displayFunc(void) {
		glClear(GL_COLOR_BUFFER_BIT);
		glRasterPos2i(-1, -1);
		glDrawPixels(width, height, GL_RGB, GL_FLOAT, pixels);
		glutSwapBuffers();
	}

	void reshapeFunc(int newWidth, int newHeight) {
		// no reshape!
		/*
			globalConfig::width = newWidth;
			globalConfig::height = newHeight;
			glViewport(0, 0, globalConfig::width, globalConfig::height);
			glLoadIdentity();
			glutPostRedisplay();
		*/
	}

	void SetupGlutDisplay(int argc, char **argv)
	{
		glutInit(&argc, argv);
		glutInitWindowSize(width, height);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);

		glutCreateWindow("Small PT");

		glutReshapeFunc(reshapeFunc);
		glutDisplayFunc(displayFunc);
		glutIdleFunc(idleFunc);
	}

}

#endif