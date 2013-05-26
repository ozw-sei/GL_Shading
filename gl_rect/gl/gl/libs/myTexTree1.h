//myTexTree1.h
//樹木
//y軸が上下方向

//最大個数
//#define MAX_REPEAT   5//幹および各枝の世代数
#define MAX_TRUNK    6//1つの幹は3個のlinkからなる
#define MAX_BRANCH0  4//各幹の親枝の個数
#define MAX_BRANCH 250//各親枝から分岐する子枝，孫枝の総数

double random(double minV, double maxV)
{//minV～maxVの一様乱数乱数
	return( (maxV - minV) * (double)rand() / (double)RAND_MAX ) + minV;
}
class CmyTexTree1
{
public:
	CmyTexTree1(void);
	~CmyTexTree1(void) { };
	void create();
	void draw();
	//メンバ変数
	CVector vPos0;//　木（幹の下端位置）
	//幹
	int numRepeatTrunk ;//世代数
	float lenTrunk0;//最下端の長さ(3個分の長さ，link１個分は1/3)
	float radTrunk0;//最下端の半径
	float rateTrunk;//成長率(link１個分）
	//枝
	float alpha0;  //分岐角
	float beta0 ;  //水平面からの偏角
	float length0; //親枝の長さ(イニシエータ）
	float radius0; //親枝の半径
	float rate   ; //成長率
	float time;//経過時間

private:
	//メンバ関数
  void drawBranch(int noK, int noI, int noBranch);
	void drawLeaf(int noK, int noJ, int noBranch);
	//メンバ変数
	float* diffuse ;//拡散光反射率
	float* specular;//鏡面光反射率
	float* ambient; //環境光反射率
	float* shadowDiffuse; //影の拡散光
	float* shadowSpecular; //影の鏡面光
	float highlight;//光沢

	//幹
	CVector vAxisTrunk[MAX_TRUNK * 3];
	float angTrunk[MAX_TRUNK * 3];//回転角
	float lenTrunk[MAX_TRUNK * 3];//長さ
	float radTrunk[MAX_TRUNK * 3];//半径

	//枝
	float alpha[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//分岐角（y軸回転）
	float beta[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH]; //親枝水平面からの偏角（x軸回転）
	float gamma[MAX_TRUNK][MAX_BRANCH0]; //幹から枝分かれする親枝の回転角（y軸回転）
	float lengthBranch[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];
	float radiusBranch[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];
	int noGeneration[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];
	int numBranch0;//幹の各世代における親枝の個数
	int numBranchAll;//各親枝に含まれる総枝数
	int numRepeat   ;//枝の分岐回数

	//葉
	float lengthLeaf[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//z軸方向サイズ
	float widthLeaf[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//x軸方向サイズ
	float betaLeaf[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//親枝水平面からの偏角（x軸回転）
	float amp[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//そよ風振動の振幅（角度）
	float period[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//そよ風振動の周期
	float phase[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//そよ風振動の位相

};
//----------------------------------------------------------------
CmyTexTree1::CmyTexTree1()
{
	static float diffuse0[] = {0.32f, 0.3f, 0.3f, 1.0f};
  static float ambient0[] = {0.2f, 0.2f, 0.2f, 1.0f};
	static float specular0[] = {0.1f, 0.1f, 0.1f, 1.0f};
//	static float shadowD[] = {0.2f, 0.2f, 0.2f, 0.5f};
//	static float shadowS[] = {0.0f, 0.0f, 0.0f, 1.0f};

	//マテリアル特性
	diffuse = diffuse0;
  ambient = ambient0;
  specular = specular0;
  highlight = 10.0f;

	//幹
	numRepeatTrunk = 4;
	rateTrunk = 0.95 + random(-0.02, 0.02);
	lenTrunk0 = 1.0 + random(-0.2, 0.2);
	radTrunk0 = 0.1 + random(-0.01, 0.01);
	for(int k = 0; k < numRepeatTrunk; k++)
	{
		vAxisTrunk[3*k]   = CVector(0.0, 0.0, -1.0);
		vAxisTrunk[3*k+1] = CVector(0.0, 0.0, -1.0);
		vAxisTrunk[3*k+2] = CVector(0.0, 0.0, -1.0);
		angTrunk[3*k] = 0.0;
		angTrunk[3*k+1] = random(-3.0, 3.0);
		angTrunk[3*k+2] = 2.0 + random(-3.0, 3.0);
	}
	
	//枝
	numRepeat = 4;
	rate = 0.8 + random(-0.05, 0.05);
	alpha0 = 40.0 + random(-10.0, 10.0) ;
	beta0  =  -10.0 + random(-5.0, 5.0);
	length0 = 0.6 + random(-0.1, 0.1);
	radius0 = 0.04 + random(-0.005, 0.005);
	numBranch0 = 3;//幹から分かれる親枝個数
	//葉

}
//-----------------------------------------------------------------------------
void CmyTexTree1::create()
{
	int i, j, k, m;

	//幹
	lenTrunk[0] = lenTrunk0 / 3.0;
	radTrunk[0] = radTrunk0;
	for(k = 1; k < numRepeatTrunk * 3; k++)
	{
		lenTrunk[k] = lenTrunk[k-1] * rateTrunk ;
		radTrunk[k] = radTrunk[k-1] * rateTrunk ;
	}

	//枝
	int numBranch = 0;//各世代の枝数
	float len0 = length0;//親枝の長さ
	float rad0 = radius0;//親枝の半径
	float len, rad;

	//各世代の枝のデータ
	for(k = 0; k < numRepeatTrunk; k++)
	{
		len0 *= rateTrunk * (0.8 + random(-0.05, 0.05));
		rad0 *= rateTrunk * (0.8 + random(-0.05, 0.05));
		for(j = 0; j < numBranch0; j++)
		{
			//親枝
			len = len0;
			rad = rad0;
			lengthBranch[k][j][0] = len;
			radiusBranch[k][j][0] = rad;
			alpha[k][j][0] = random(-10.0, 10.0);
			beta[k][j][0] = 30.0 + random(-10.0, 10.0);
			gamma[k][j] = (360.0 / (float)numBranch0) * float(j) + random(-20.0, 20.0) ;
			noGeneration[k][j][0] = 0;
			//子枝，孫枝
			numBranch = 1;
			numBranchAll = 1;
			for(i = 1; i < numRepeat; i++)
			{
				len *= rate;
				rad *= rate;
				for(m = 0; m < numBranch; m++)
				{
					//中央
					lengthBranch[k][j][numBranchAll] = len * random(0.8, 1.2);
					radiusBranch[k][j][numBranchAll] = rad * random(0.8, 1.2);
					alpha[k][j][numBranchAll] = random(-10, 10.0);
					beta[k][j][numBranchAll] = beta0 * random(0.8, 1.2);
					noGeneration[k][j][numBranchAll] = i;
					numBranchAll ++;
					if(numBranchAll > MAX_BRANCH) { printf("枝総数が多すぎます \n"); return;}
					//左側
					lengthBranch[k][j][numBranchAll] = len * random(0.8, 1.2);
					radiusBranch[k][j][numBranchAll] = rad * random(0.8, 1.2);
					alpha[k][j][numBranchAll] = -alpha0 * random(0.8, 1.2);
					beta[k][j][numBranchAll] = beta0 * random(0.8, 1.2);
					noGeneration[k][j][numBranchAll] = i;
					numBranchAll ++;
					if(numBranchAll > MAX_BRANCH) { printf("枝総数が多すぎます \n"); return;}
					//右側
					lengthBranch[k][j][numBranchAll] = len * random(0.8, 1.2);
					radiusBranch[k][j][numBranchAll] = rad * random(0.8, 1.2);
					alpha[k][j][numBranchAll] = alpha0 * random(0.8, 1.2);
					beta[k][j][numBranchAll] = beta0 * random(0.8, 1.2);
					noGeneration[k][j][numBranchAll] = i;
					numBranchAll ++;
					if(numBranchAll > MAX_BRANCH) { printf("枝総数が多すぎます \n"); return;}
				}
				numBranch *= 3;//次世代の枝数
			}//i
//printf("k=%d, j=%d, nnumBranch=%d, numBranchAll=%d \n",k, j, numBranch, numBranchAll);
		}//j
	}//k

	//葉のデータ
	for(k = 0; k < numRepeatTrunk; k++)
	{
		for(j = 0; j < numBranch0; j++)
		{
			for(i = 0; i < numBranchAll; i++)
			{
				lengthLeaf[k][j][i] = 0.7 + random(-0.2, 0.2);
				widthLeaf[k][j][i]  = 0.5 + random(-0.1, 0.1);
				betaLeaf[k][j][i]   = 30.0 + random(-10.0, 10.0);
				amp[k][j][i]        = 3.0 + random(-2.0, 2.0);
				period[k][j][i]     = 2.0 + random(-0.1, 0.1);
				phase[k][j][i]      = random(-10.0, 10.0);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void CmyTexTree1::draw()
{
	int j, k;

	//木全体
	glPushMatrix();
	glTranslatef(vPos0.x, vPos0.y, vPos0.z);
	//幹
	for(k = 0; k < numRepeatTrunk; k++)
	{
		glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
		glMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
		glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
		glMaterialf(GL_FRONT,GL_SHININESS,highlight);

		glRotatef(angTrunk[3*k], vAxisTrunk[3*k].x, vAxisTrunk[3*k].y, vAxisTrunk[3*k].z);
		glPushMatrix();
			glScalef(radTrunk[3*k], lenTrunk[3*k], radTrunk[3*k]);
			glTranslatef(0.0, 0.5, 0.0);
			drawCylinderY(1.0, rateTrunk, 1.1, 6, 2);
		glPopMatrix();
		glTranslatef(0.0, lenTrunk[3*k], 0.0);

		glRotatef(angTrunk[3*k+1], vAxisTrunk[3*k+1].x, vAxisTrunk[3*k+1].y, vAxisTrunk[3*k+1].z);
		glPushMatrix();
			glScalef(radTrunk[3*k+1], lenTrunk[3*k+1], radTrunk[3*k+1]);
			glTranslatef(0.0, 0.5, 0.0);
			drawCylinderY(1.0, rateTrunk, 1.1, 6, 2);
		glPopMatrix();
		glTranslatef(0.0, lenTrunk[3*k+1], 0.0);

		glRotatef(angTrunk[3*k+2], vAxisTrunk[3*k+2].x, vAxisTrunk[3*k+2].y, vAxisTrunk[3*k+2].z);
		glPushMatrix();
			glScalef(radTrunk[3*k+2], lenTrunk[3*k+2], radTrunk[3*k+2]);
			glTranslatef(0.0, 0.5, 0.0);
			drawCylinderY(1.0, rateTrunk, 1.1, 6, 2);
		glPopMatrix();
		glTranslatef(0.0, lenTrunk[3*k+2], 0.0);

		//枝
		for(j = 0; j < numBranch0; j++) 
		{
			
//printf("k =%d, j=%d, gamma=%f \n", k, j, gamma[k][j]);

			glPushMatrix();
				glRotatef(gamma[k][j], 0.0, 1.0, 0.0);
				drawBranch(k, j, 0) ;
			glPopMatrix();
		}
	}
	glPopMatrix();
}

//-----------------------------------------------------------------------------
void CmyTexTree1::drawBranch(int noK, int noJ, int noBranch)
{
	//noK;幹の世代番号
	//noJ:幹から分かれる親枝番号
	//noBranch:1つの親枝集団に含まれる枝番号

	if(noBranch >= numBranchAll) return;
	glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
	glMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
	glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
	glMaterialf(GL_FRONT,GL_SHININESS,highlight);

/*
if(noK == 0 && noJ == 0)
{
printf("noBranch=%d, noGeneration=%d \n",noBranch, noGeneration[noK][noJ][noBranch]);
if(noGeneration[noK][noJ][noBranch] >= numRepeat-1)
{
	float diffuseA[] = {0.99f, 0.0f, 0.0f, 1.0f};
	glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuseA);
}
else
{
	glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
}
}*/
	glPushMatrix();
	//親
	//glRotatef(ang[3*k], vAxis[3*k].x, vAxis[3*k].y, vAxis[3*k].z);
	glRotatef(beta[noK][noJ][noBranch], -1.0, 0.0, 0.0);
	glRotatef(alpha[noK][noJ][noBranch], 0.0, 1.0, 0.0);
	glPushMatrix();
		glScalef(radiusBranch[noK][noJ][noBranch], radiusBranch[noK][noJ][noBranch], lengthBranch[noK][noJ][noBranch]);
		glTranslatef(0.0, 0.0, 0.5);
		drawCylinder(1.0, rate, 1.1, 4, 1);
	glPopMatrix();
	glTranslatef(0.0, 0.0, lengthBranch[noK][noJ][noBranch]);
//if(noGeneration[noK][noJ][noBranch] >= numRepeat-2) 
	drawLeaf(noK, noJ, noBranch);

	//子
	drawBranch(noK, noJ, 3*noBranch+1);
	drawBranch(noK, noJ, 3*noBranch+2);
	drawBranch(noK, noJ, 3*noBranch+3);
	glPopMatrix();
}

void CmyTexTree1::drawLeaf(int noK, int noJ, int noBranch)
{
//	static float diffuseL[] = {0.6f, 0.6f, 0.6f, 1.0f};
//  static float ambientL[] = {0.6f, 0.6f, 0.6f, 1.0f};
	static float diffuseL[] = {0.3f, 0.5f, 0.2f, 1.0f};
  static float ambientL[] = {0.2f, 0.4f, 0.1f, 1.0f};
	static float specularL[] = {0.1f, 0.1f, 0.1f, 1.0f};
	static float highlightL = 5.0;
	glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuseL);
	glMaterialfv(GL_FRONT,GL_AMBIENT,ambientL);
	glMaterialfv(GL_FRONT,GL_SPECULAR,specularL);
	glMaterialf(GL_FRONT,GL_SHININESS,highlightL);

	float amp0 = amp[noK][noJ][noBranch];
	float period0 = period[noK][noJ][noBranch];
	float phase0 = phase[noK][noJ][noBranch];
	float ang = amp0 * sin(2.0 * M_PI * time / period0 + phase0);
	glPushMatrix();
		glRotatef(betaLeaf[noK][noJ][noBranch] + ang, 1.0, 0.0, 0.0);
		glScalef(widthLeaf[noK][noJ][noBranch], 1.0, lengthLeaf[noK][noJ][noBranch]);
		glTranslatef(0.0, 0.0, 0.5);
		//drawPlateY(1.0);
		drawTexPlateY(1.0, 1, 1);
	glPopMatrix();


}


