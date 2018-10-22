#include <GL/freeglut.h>

int width = 800;
int height = 600;

void idleFunc(void) {
	
	glutPostRedisplay();
}

void displayFunc(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	
	glutSwapBuffers();
}

void reshapeFunc(int newWidth, int newHeight) {
	width = newWidth;
	height = newHeight;

	glViewport(0, 0, width, height);
	glLoadIdentity();
	
	glutPostRedisplay();
}

int main(int argc, char **argv)
{
    glutInitWindowSize(width, height);
    glutInitWindowPosition(0,0);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInit(&argc, argv);

	glutCreateWindow("Small PT");

    glutReshapeFunc(reshapeFunc);
    glutDisplayFunc(displayFunc);
	glutIdleFunc(idleFunc);

    //glutKeyboardFunc(keyFunc);
    //glutSpecialFunc(specialFunc);
    

    glutMainLoop();
	
    return 0;
}