//
//  main.cpp
//  gl
//
//  Created by Ozawa Seijiro on 2013/05/26.
//  Copyright (c) 2013年 Ozawa Seijiro. All rights reserved.
//

#include <iostream>
#include <GLUT/GLUT.h>

// コールバック関数
void display(void)
{
    return;
}

// メイン関数
int main(int argc, const char* argv[])
{
    // 変数を定義＆初期化する
    char title[80] = "sample";
    // GLUTライブラリを初期化する
    glutInit(&argc, (char**)argv);
    // ウィンドウを生成する
    glutCreateWindow(title);
    // カレントウィンドウに対するコールバック関数を設定する
    glutDisplayFunc(display);
    // メインループを開始する
    glutMainLoop();
    // 終了
    return 0;
}