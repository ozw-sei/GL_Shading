#include <iostream>
#include <stdio.h>
#include <GLUT/glut.h>

//画面サイズ
const int g_SCREEN_WIDTH = 640;
const int g_SCREEN_HEIGHT = 480;

//関数プロトタイプ宣言
void init();
void resize(int w, int h);
void display();
void keyboard(unsigned char key,int x,int y);

void drawRectangle(
        float x0, float y0,
        float x1, float y1,
        float x2 ,float y2,
        float x3 ,float y3 );

void init(){
    //背景色
    glClearColor(0.2, 0.2, 0.2, 1.0);
}
void resize( int w, int h ){

}
void drawRectangle(
        const float x0, const float y0,
        const float x1, const float y1,
        const float x2 ,const float y2,
        const float x3 ,const float y3 ){
    //緑色
    glColor3f(0.0, 1.0, 1.0);
    glLineWidth(0.5);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_QUADS);

        glVertex2f(x0, y0);
        glVertex2f(x1, y1);
        glVertex2f(x2, y2);
        glVertex2f(x3, y3);

    glEnd();

}


void display()
{
    //カラーバッファのクリア
    glClear(GL_COLOR_BUFFER_BIT);


    drawRectangle(
            -0.5, -0.5,
            0.5, -0.5,
            0.5, 0.5,
            -0.5, 0.5);

    glFlush();


}

void keyboard(unsigned char key,int x,int y){
    switch ((unsigned char)key) {
        case 27://ESC
            exit(0);
            break;
        default:
            break;
    }
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGBA);
    
    glutInitWindowSize(g_SCREEN_WIDTH,g_SCREEN_HEIGHT);

    glutInitWindowPosition(100, 100);

    glutCreateWindow("GLUT Program");

    glutKeyboardFunc(keyboard);

    glutDisplayFunc(display);

    init();

    glutMainLoop();
    return 0;
}