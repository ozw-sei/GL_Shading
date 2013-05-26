//#include <math.h>
#include <stdio.h>
//#include "../../spline.h"
//#include "../../myMath3.h"

#define M_PI 3.14159265358979323846

//法線方向計算ルーチン
void calcNormal(float* p1,float* p2,float* p3,float* nn)
{
	double mag;
	float a[3], b[3];
	int k, k1, k2;


	for(k = 0; k < 3; k++)
	{
		a[k] = p2[k] - p1[k];
		b[k] = p3[k] - p2[k];
	}

	for(k = 0; k < 3; k++)
	{
		k1 = k + 1; k2 = k + 2;
		if(k1 >= 3) k1 = k1 -3;
		if(k2 >= 3) k2 = k2 -3;
		nn[k] = a[k1]*b[k2] - a[k2]*b[k1];
	}
	//normalize
	mag = nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2];
	for(k = 0; k < 3; k++) nn[k] = nn[k] / mag; 
}

//------------------------------------------------------------------------------
void drawPlateZ(float s)//x-y平面
{//面の法線方向がｚ軸
	float p[4][3] = 
	{ { s/2.0,-s/2.0, 0.0}, { s/2.0, s/2.0, 0.0}, {-s/2.0, s/2.0, 0.0}, 
	  {-s/2.0,-s/2.0, 0.0}};

	glBegin(GL_QUADS);
		glNormal3f(0.0, 0.0, 1.0); //z方向の法線
		glVertex3fv(p[0]);
		glVertex3fv(p[1]);
		glVertex3fv(p[2]);
		glVertex3fv(p[3]);
	glEnd();
}
//------------------------------------------------------------------------------
void drawPlateY(float s)//x-z平面
{//面の法線方向がy軸
	float p[4][3] = 
	{ { s/2.0, 0.0, s/2.0}, { s/2.0, 0.0,-s/2.0}, {-s/2.0, 0.0,-s/2.0}, {-s/2.0, 0.0, s/2.0}};

	glBegin(GL_QUADS);
		glNormal3f(0.0, 1.0, 0.0); //y方向の法線
		glVertex3fv(p[0]);
		glVertex3fv(p[1]);
		glVertex3fv(p[2]);
		glVertex3fv(p[3]);
	glEnd();
}

//-----------------------------------------------------------------------

void drawCube(float s)
{
	float p[8][3] = 
	{ {0.5*s,0.5*s,0.5*s}, {-0.5*s,0.5*s,0.5*s}, {-0.5*s,-0.5*s,0.5*s}, 
	  {0.5*s,-0.5*s,0.5*s},{0.5*s,0.5*s,-0.5*s}, {-0.5*s,0.5*s,-0.5*s},
	  {-0.5*s,-0.5*s,0-0.5*s}, {0.5*s,-0.5*s,-0.5*s}
	};

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,1.0f); //z方向
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);

		glNormal3f(1.0f,0.0f,0.0f); //x方向(正面）
		glVertex3fv(p[0]); glVertex3fv(p[3]);
		glVertex3fv(p[7]); glVertex3fv(p[4]);

		glNormal3f(0.0f,1.0f,0.0f); //y方向
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[5]); glVertex3fv(p[1]);

	 	glNormal3f(-1.0f,0.0f,0.0f); //-x方向
		glVertex3fv(p[1]); glVertex3fv(p[5]);
		glVertex3fv(p[6]); glVertex3fv(p[2]);

		glNormal3f(0.0f,-1.0f,0.0f); //-y方向
		glVertex3fv(p[2]); glVertex3fv(p[6]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);

		glNormal3f(0.0f,0.0f,-1.0f); //-z方向
		glVertex3fv(p[4]); glVertex3fv(p[7]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
}
//-----------------------------------------------------------------

void drawGridPlate(float sizeX, float sizeY, int nSliceX, int nSliceY)
{//xy平面，中心は原点
	int i, j;
	//double ;
	float p[2][3];
	float pitchX = sizeX / (float)nSliceX;
	float pitchY = sizeY / (float)nSliceY;

	for(j = 0; j < nSliceY; j++)
	{
		//頂点z座標
		p[0][2] = 0.0;
		p[1][2] = 0.0;

		glBegin(GL_QUAD_STRIP);
		for(i = 0; i <= nSliceX; i++)
		{
			p[0][0] = (float)i * pitchX - sizeX / 2.0;//x座標
			p[0][1] = (float)j * pitchY - sizeY / 2.0;//y座標
			p[1][0] = (float)i * pitchX - sizeX / 2.0;//x座標
			p[1][1] = (float)(j+1) * pitchY - sizeY / 2.0;//y座標

			glNormal3f(0.0, 0.0, 1.0);//法線ベクトル
			glVertex3fv(p[0]);//頂点座標
			glVertex3fv(p[1]);//頂点座標			
		}
		glEnd();
	}
}

//-----------------------------------------------------------------------
void drawSphere(float radius, int nSlice, int nStack)
{
	int i, j;
	double r0, r1, th0, th1, phi;
	double p[2][3];

	for(j = 0; j < nStack; j++)
	{
		//j=0は北極点(x=0, y=0, z=radius)
		//天頂角
		th0 = M_PI * (double)j / (double)nStack;
		th1 = M_PI * (double)(j+1) / (double)nStack;
		//x-y平面に投影した半径
		r0 = radius * sin(th0);
		r1 = radius * sin(th1);
		//頂点z座標
		p[0][2] = radius * cos(th0);
		p[1][2] = radius * cos(th1);

		glBegin(GL_QUAD_STRIP);
		for(i = 0; i <= nSlice; i++)
		{
			phi = 2.0 * M_PI * (double)i / (double)nSlice;
			//頂点のxy座標(i=0をobjectからみて右端) 
			p[0][0] =   r0 * sin(phi);//x座標
			p[0][1] = - r0 * cos(phi);//y座標
			p[1][0] =   r1 * sin(phi);//x座標
			p[1][1] = - r1 * cos(phi);//y座標

			glNormal3dv(p[0]);//法線ベクトル
			glVertex3dv(p[0]);//頂点座標

			glNormal3dv(p[1]);//法線ベクトル
			glVertex3dv(p[1]);//頂点座標			
		}
		glEnd();
	}
}
//----------------------------------------------------------------------
void drawCylinder(float rBottom, float rTop, float height, int nSlice, int nStack)
{ //上方向はz軸
	//円柱(rBottom=rTop))、円錐台、円錐(rTop = 0.0)
	//rBottom:下底半径, rTop:上底半径
	//nSlice--xy断面分割数
	//nStack---ｚ方向分割数
	//物体の中心は下底と上底の中間

	int i, j;
	double x, y, z, z0, z1;
	double theta;
	double theta0 = 2.0*M_PI/(double)nSlice;
	//上底(Top)
	glBegin(GL_TRIANGLE_FAN);
	glNormal3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, height/2.0);
	for(i = 0; i <= nSlice; i++) 
	{ 
		theta = theta0 * (double)i;
		x = (float)(rTop * cos(theta)); //x成分
		y = (float)(rTop * sin(theta)); //ｙ成分
		z = height/2.0;
		glVertex3f(x, y, z);
	}
	glEnd();

	//下底(Bottom)
	glBegin(GL_TRIANGLE_FAN);
	glNormal3f(0.0, 0.0, -1.0);
	glVertex3f(0.0, 0.0, -height/2.0);
	for(i = 0; i <= nSlice; i++) 
	{ 
		theta = theta0 * (double)(nSlice - i);
		x = (float)(rBottom * cos(theta)); //x成分
		y = (float)(rBottom * sin(theta)); //ｙ成分
		z = -height/2.0;
		glVertex3f(x, y, z);
	}
	glEnd();

	double s, t0, t1, r0, r1, phi;
	double p[2][3], n[3];

	double rr = rBottom - rTop;
	double nz = rr / sqrt(rr*rr + height*height);
  double nxy = sqrt(1.0 - nz * nz);
	//側面
	for(j = 0; j < nStack; j++)
	{
		//j=0は上底(x=0, y=0, z=height/2)
		//2つのt座標
		t0 = (double)j / (double)nStack;
		t1 = (double)(j+1) / (double)nStack;
		//底面からの高さ
		z0 = height * (1.0 - t0);
		z1 = height * (1.0 - t1);
		//半径
		r0 = rBottom + (rTop - rBottom) * z0 / height;
		r1 = rBottom + (rTop - rBottom) * z1 / height;

		//頂点z座標
		p[0][2] = z0 - height / 2.0;
		p[1][2] = z1 - height / 2.0;

		glBegin(GL_QUAD_STRIP);
		for(i = 0; i <= nSlice; i++)
		{
			//s座標
			s = (double)i / (double)nSlice;
			phi = 2.0 * M_PI * s + M_PI / 6.0;
			//頂点のxy座標
			p[0][0] = r0 * cos(phi);//x座標
			p[0][1] = r0 * sin(phi);//y座標
			p[1][0] = r1 * cos(phi);//x座標
			p[1][1] = r1 * sin(phi);//y座標
			//法線ベクトル
			n[0] = nxy * cos(phi);
			n[1] = nxy * sin(phi);
			n[2] = nz;

			glNormal3dv(n);//法線ベクトル
			glVertex3dv(p[0]);//頂点座標
			glVertex3dv(p[1]);//頂点座標
		}
		glEnd();
	}
}
//----------------------------------------------------------------------
void drawCylinderY(float rBottom, float rTop, float height, int nSlice, int nStack)
{ //上方向はY軸
	//円柱(rBottom=rTop))、円錐台、円錐(rTop = 0.0)
	//rBottom:下底半径, rTop:上底半径
	//nSlice--ZX断面分割数
	//nStack---Y方向分割数
	//物体の中心は下底と上底の中間

	int i, j;
	double Z, X, Y, Y0, Y1;
	double theta;
	double theta0 = 2.0*M_PI/(double)nSlice;
	//上底(Top)
	glBegin(GL_TRIANGLE_FAN);
	glNormal3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, height/2.0, 0.0);
	for(i = 0; i <= nSlice; i++) 
	{ 
		theta = theta0 * (double)i;
		Z = (float)(rTop * cos(theta)); //Z成分
		X = (float)(rTop * sin(theta)); //X成分
		Y = height/2.0;
		glVertex3f(X, Y, Z);
	}
	glEnd();

	//下底(Bottom)
	glBegin(GL_TRIANGLE_FAN);
	glNormal3f(0.0, -1.0, 0.0);
	glVertex3f(0.0, -height/2.0, 0.0);
	for(i = 0; i <= nSlice; i++) 
	{ 
		theta = theta0 * (double)(nSlice - i);
		Z = (float)(rBottom * cos(theta)); //Z成分
		X = (float)(rBottom * sin(theta)); //ｙ成分
		Y = -height/2.0;
		glVertex3f(X, Y, Z);
	}
	glEnd();

	double s, t0, t1, r0, r1, phi;
	double p[2][3], n[3];

	double rr = rBottom - rTop;
	double nY = rr / sqrt(rr*rr + height*height);
  double nZX = sqrt(1.0 - nY * nY);
	//側面
	for(j = 0; j < nStack; j++)
	{
		//j=0は上底(Z=0, X=0, Y=height/2)
		//2つのt座標
		t0 = (double)j / (double)nStack;
		t1 = (double)(j+1) / (double)nStack;
		//底面からの高さ
		Y0 = height * (1.0 - t0);
		Y1 = height * (1.0 - t1);
		//半径
		r0 = rBottom + (rTop - rBottom) * Y0 / height;
		r1 = rBottom + (rTop - rBottom) * Y1 / height;

		//頂点Y座標
		p[0][1] = Y0 - height / 2.0;
		p[1][1] = Y1 - height / 2.0;

		glBegin(GL_QUAD_STRIP);
		for(i = 0; i <= nSlice; i++)
		{
			//s座標
			s = (double)i / (double)nSlice;
			phi = 2.0 * M_PI * s + M_PI / 6.0;
			//頂点のZX座標
			p[0][2] = r0 * cos(phi);//Z座標
			p[0][0] = r0 * sin(phi);//X座標
			p[1][2] = r1 * cos(phi);//Z座標
			p[1][0] = r1 * sin(phi);//X座標
			//法線ベクトル
			n[2] = nZX * cos(phi);
			n[0] = nZX * sin(phi);
			n[1] = nY;

			glNormal3dv(n);//法線ベクトル
			glVertex3dv(p[0]);//頂点座標
			glVertex3dv(p[1]);//頂点座標
		}
		glEnd();
	}
}

//-----------------------------------------------------------------------------------------

void drawTorus(float radius1, float radius2, int nSide, int nRing)
{	
	//radius1:円環断面半径
	//radius2:円環の中心軸半径
	//nSide:円環断面における表面分割点数
	//nRing:円環の分割数
	if(radius1 > radius2) { printf("radius1 < radius2としてください \n "); return;}

	int i, j;
	double rr, zz;
	double phi0, phi1, theta;
	double p[2][3];

	for(i = 0; i < nRing; i++)
	{
		//i=0は基本断面(x=radius2を中心とする円, y=0）
		phi0 = 2.0 * M_PI * (double)i / (double)nRing;
		phi1 = 2.0 * M_PI * (double)(i+1) / (double)nRing;

		glBegin(GL_QUAD_STRIP);
		for(j = 0; j <= nSide; j++)
		{
			theta = M_PI - 2.0 * M_PI * (double)j / (double)nSide;
			rr = radius2 + radius1 * cos(theta);//z軸からの距離
			zz = radius1 * sin(theta);
			//頂点のxyz座標(j=0を内側xy平面)
			p[0][0] = rr * cos(phi0);//x座標
			p[0][1] = rr * sin(phi0);//y
			p[0][2] = zz;            //z
			p[1][0] = rr * cos(phi1);//x座標
			p[1][1] = rr * sin(phi1);//y
			p[1][2] = zz;            //z      

			glNormal3d(cos(theta)*cos(phi0),cos(theta)*sin(phi0),sin(theta));
			glVertex3dv(p[0]);//頂点座標

			glNormal3d(cos(theta)*cos(phi1),cos(theta)*sin(phi1),sin(theta));
			glVertex3dv(p[1]);//頂点座標
		}
		glEnd();
	}
}
//-----------------------------------------------------------------
void drawSuper(float r, int nSlice, int nStack, double eps1, double eps2)
{
	//上下の中心が原点
	int i,j,ip,im,np,npL,npR,npU,npD,k1,k2;
	double ct,theta,phi,z,cc;
	float a[31][31], b[31][31], c[31][31];
	float n1[3], n2[3], n3[3], n4[3];
	float pd[31*31][3];

	if(nSlice > 30) nSlice = 30;
	if(nStack > 30) nStack = 30;

	for(j = 0 ;j <= nStack;j++)
	{
		theta = (M_PI/(double)nStack) * ((double)nStack / 2.0 - (double)j);
		                //thetaはx-y平面からの偏角となっている

		if(theta >= 0.0) //上半分
		{
			if(theta == 0.0) z = 0.0;//pow(a,b)のaがa<=0.0でエラー
			else z = pow(sin(fabs(theta)),eps1);//z
		}
		else  //下半分        
		{
			z = - pow(sin(fabs(theta)), eps1);
		}
		for (i = 0 ;i <= nSlice / 2;i++)
		{
			k1 = nSlice * j + i;//objectから見て左側
			k2 = nSlice * j + nSlice - i;//右側
			phi = 2.0 * M_PI * (double)i/(double)nSlice;
			ct = cos(phi);
			if( ct == 0.0 ) cc = 0.0;
			else if (ct > 0) { cc = pow(ct, eps2);}
			else         { cc = -pow(fabs(ct),eps2); }
			if(j == 0 || j == nStack) 
			{
				pd[k1][0] = 0.0f;
				pd[k1][1] = 0.0f;
			}

			else 
			{
				pd[k1][0] = r * (float)(pow(cos(theta),eps1)*cc);
				if(sin(phi) == 0.0) pd[k1][1] = 0.0f;
				else pd[k1][1] = r * (float)(pow(cos(theta),eps1)*pow(fabs(sin(phi)),eps2));
			}
			if(i == 0) k2 = k1;
			pd[k2][0] = pd[k1][0];
			pd[k2][1] = -pd[k1][1];
			pd[k1][2] = r * (float)z;
			pd[k2][2] = r * (float)z;
		}
	}

	//側面の法線成分
	for(i = 0;i < nSlice;i++)
	{
		ip = i+1;
		if(ip == nSlice) ip = 0;
		im = i-1;
		if(i == 0) im = nSlice-1;

		//真上(Top)
		a[i][0] = 0.0f; b[i][0] = 0.0f; c[i][0] = 1.0f;
		//真下（Bottom)
		a[i][nStack] = 0.0f; b[i][nStack] = 0.0f; c[i][nStack] = -1.0f;

		for(j=1;j<nStack;j++)//隣り合う4個の三角形の法線ベクトルを平均化
		{
			np = j*nSlice+i;//注目点
			npL = j*nSlice+im;//左側
			npR = j*nSlice+ip;//右側
			npU = np-nSlice;//上
			npD = np+nSlice;//下
			if(j == 1) 
			{
				n1[0]=0.0f; n1[1]=0.0f; n1[2]=1.0f;//Top
				n2[0]=0.0f; n2[1]=0.0f; n2[2]=1.0f;//Top
				calcNormal(pd[np],pd[npL],pd[npD],n3);//外から見て左下
				calcNormal(pd[np],pd[npD],pd[npR],n4);//右下
			}
			if(j == nStack-1)
			{
				calcNormal(pd[np],pd[npU],pd[npL],n1);//外から見て左上
				calcNormal(pd[np],pd[npR],pd[npU],n2);//右上
				n3[0]=0.0f; n3[1]=0.0f; n3[2]=-1.0f;//Bottom
				n4[0]=0.0f; n4[1]=0.0f; n4[2]=-1.0f;//Bottom
			}
			else 
			{
				calcNormal(pd[np],pd[npU],pd[npL],n1);//外から見て左上
				calcNormal(pd[np],pd[npR],pd[npU],n2);//右上
				calcNormal(pd[np],pd[npL],pd[npD],n3);//外から見て左下
				calcNormal(pd[np],pd[npD],pd[npR],n4);//右下
			}
			a[i][j] = (float)((n1[0]+n2[0]+n3[0]+n4[0])/4.0f);//ｘ方向
			b[i][j] = (float)((n1[1]+n2[1]+n3[1]+n4[1])/4.0f);//ｙ
			c[i][j] = (float)((n1[2]+n2[2]+n3[2]+n4[2])/4.0f);//ｚ
		}
	}

	glBegin(GL_QUADS);
	for(i = 0;i < nSlice;i++)
	{
		ip = i + 1;
		if(ip == nSlice) ip = 0;
		for(j = 0;j < nStack; j++)
		{
			glNormal3f(a[i][j],b[i][j],c[i][j]);
			glVertex3fv(pd[i+j*nSlice]);
			glNormal3f(a[i][j+1],b[i][j+1],c[i][j+1]);
			glVertex3fv(pd[i+(j+1)*nSlice]);
			glNormal3f(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
			glVertex3fv(pd[ip+(j+1)*nSlice]);
			glNormal3f(a[ip][j],b[ip][j],c[ip][j]);
			glVertex3fv(pd[ip +j*nSlice]);
		}
	}
	glEnd();
}

//-----------------------------------------------------------------------------
//四角形のメッシュ（x-y平面,中心が原点）
//ｘ軸方向，ｙ軸方向の幅を固定
void drawElevation(int Nx, int Ny, float sizeX, float sizeY, int sideFlag, float* data)
{
	//全体の幅,長さsizeX, sizeY
	//sideFlag = 0:側面表示せず
	//sideFlag = 1:側面表示する

	const int NMAX = 130;
	int i, j;
	float p[NMAX][NMAX][3]; //頂点座標
	float a[NMAX][NMAX], b[NMAX][NMAX], c[NMAX][NMAX];//頂点の法線
	float pitchX, pitchY;
	float n1[3], n2[3], n3[3], n4[3];

	if(Nx > NMAX) printf("NxがNMAXを超えています(drawElevation1) \n");
	if(Ny > NMAX) printf("NyがNMAXを超えています(drawElevation1) \n");

	//セルのサイズ
	pitchX = sizeX / (float)Nx;
	pitchY = sizeY / (float)Ny;

	//各頂点の座標
	for(j = 0; j <= Ny; j++){
		for(i = 0; i <= Nx; i++){
			p[i][j][0] = (float)(i - Nx / 2) * pitchX;
			p[i][j][1] = (float)(j - Ny / 2) * pitchY;
			p[i][j][2] = data[j * (Nx+1) + i];
		}
	}
	
	//法線成分
	for(i = 0;i <= Nx;i++)
		for(j = 0;j <= Ny;j++)
		{			
			if(j == 0 )
			{
				if(i == 0) 
				{
					calcNormal(p[0][0],p[1][0],p[0][1],n1);
					a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
				else if(i == Nx) 
				{
					calcNormal(p[Nx-1][0],p[Nx][0],p[Nx][1],n1);
					a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
				else 
				{
					calcNormal(p[i][0],p[i][1],p[i-1][0],n1);//左側
					calcNormal(p[i][0],p[i+1][0],p[i][1],n2);//右側
					a[i][j] = (n1[0]+n2[0])/2.0f;
					b[i][j] = (n1[1]+n2[1])/2.0f;
					c[i][j] = (n1[2]+n2[2])/2.0f; }
			}
			else if(j == Ny)
			{
				if(i == 0) 
				{
					calcNormal(p[0][Ny],p[0][Ny-1],p[1][Ny],n1);
					a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
				else if(i == Nx) 
				{
					calcNormal(p[Nx][Ny],p[Nx-1][Ny],p[Nx][Ny-1],n1);
					a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
				else 
				{
					calcNormal(p[i][Ny],p[i-1][Ny],p[i][Ny-1],n1);//左側
					calcNormal(p[i][Ny],p[i][Ny-1],p[i+1][Ny],n2);//右側
					a[i][j] = (n1[0]+n2[0])/2.0f;
					b[i][j] = (n1[1]+n2[1])/2.0f;
					c[i][j] = (n1[2]+n2[2])/2.0f; }
			}
			else
			{
				if(i == 0) 
				{
					calcNormal(p[0][j],p[1][j],p[0][j+1],n1);//上
					calcNormal(p[0][j],p[0][j-1],p[0][1],n2);//下
					a[i][j] = (n1[0]+n2[0])/2.0f;
					b[i][j] = (n1[1]+n2[1])/2.0f;
					c[i][j] = (n1[2]+n2[2])/2.0f; }
				else if(i == Nx) 
				{
					calcNormal(p[Nx][j],p[Nx][j+1],p[Nx-1][j],n1);//上
					calcNormal(p[Nx][j],p[Nx-1][j],p[Nx][j-1],n2);//下
					a[i][j] = (n1[0]+n2[0])/2.0f;
					b[i][j] = (n1[1]+n2[1])/2.0f;
					c[i][j] = (n1[2]+n2[2])/2.0f; }
				else 
				{//上下左右４個の三角形の平均
					calcNormal(p[i][j],p[i][j+1],p[i-1][j],n1);//左上
					calcNormal(p[i][j],p[i+1][j],p[i][j+1],n2);//右上
					calcNormal(p[i][j],p[i-1][j],p[i][j-1],n3);//左下
					calcNormal(p[i][j],p[i][j-1],p[i+1][j],n4);//右下
					a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
					b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
					c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
			}
		}

//	int nC;
	//三角形で面を定義
	glBegin(GL_TRIANGLES);
	for(j = 0;j < Ny;j++)
		for(i = 0;i < Nx;i++)
		{
			//左下の三角形
			//各頂点の法線方向,ﾃｸｽﾁｬｰ座標,頂点座標を与える。
			glNormal3f(a[i][j],b[i][j],c[i][j]);//法線方向
			glVertex3fv(p[i][j]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
			glNormal3f(a[i+1][j],b[i+1][j],c[i+1][j]);
			glVertex3fv(p[i+1][j]);
			glNormal3f(a[i][j+1],b[i][j+1],c[i][j+1]);
			glVertex3fv(p[i][j+1]);
			//右上の三角形
			glNormal3f(a[i+1][j],b[i+1][j],c[i+1][j]);
			glVertex3fv(p[i+1][j]);
			glNormal3f(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glVertex3fv(p[i+1][j+1]);
			glNormal3f(a[i][j+1],b[i][j+1],c[i][j+1]);
			glVertex3fv(p[i][j+1]);
		}
	glEnd();

	if(sideFlag == 1)//側面描画
	{
		glBegin(GL_QUADS);
		//+x方向（i=Nx)
		glNormal3f(1.0, 0.0, 0.0);
		for(j = 0; j < Ny; j++)
		{
			glVertex3f(p[Nx][j][0], p[Nx][j][1], 0.0f);
			glVertex3f(p[Nx][j+1][0], p[Nx][j+1][1], 0.0f);
			glVertex3f(p[Nx][j+1][0], p[Nx][j+1][1], p[Nx][j+1][2]);
			glVertex3f(p[Nx][j][0], p[Nx][j][1], p[Nx][j][2]);
		}
		//-x方向（i=0)
		glNormal3f(-1.0, 0.0, 0.0);
		for(j = 0; j < Ny; j++)
		{
			glVertex3f(p[0][j][0], p[0][j][1], 0.0f);
			glVertex3f(p[0][j][0], p[0][j][1], p[0][j][2]);
			glVertex3f(p[0][j+1][0], p[0][j+1][1], p[0][j+1][2]);
			glVertex3f(p[0][j+1][0], p[0][j+1][1], 0.0f);
		}
		//+y方向（j=Ny)
		glNormal3f(0.0, 1.0, 0.0);
		for(i = 0; i < Nx; i++)
		{
			glVertex3f(p[i][Ny][0], p[i][Ny][1], 0.0f);
			glVertex3f(p[i][Ny][0], p[i][Ny][1], p[i][Ny][2]);
			glVertex3f(p[i+1][Ny][0], p[i+1][Ny][1], p[i+1][Ny][2]);
			glVertex3f(p[i+1][Ny][0], p[i+1][Ny][1], 0.0f);
		}
		//-y方向（j=0)
		glNormal3f(0.0, -1.0, 0.0);
		for(i = 0; i < Nx; i++)
		{
			glVertex3f(p[i][0][0], p[i][0][1], 0.0f);
			glVertex3f(p[i+1][0][0], p[i+1][0][1], 0.0f);
			glVertex3f(p[i+1][0][0], p[i+1][0][1], p[i+1][0][2]);
			glVertex3f(p[i][0][0], p[i][0][1], p[i][0][2]);
		}
		glEnd();
	}
}

//-------------------------------------------------------------------
void drawFloor(float widthX, float widthZ, int nx, int nz)
{
  int i, j;
  //Floor１枚当たりの幅
  float wX = widthX / (float)nx;
  float wZ = widthZ / (float)nz;

  float diffuse[][4] = {
	{ 0.7, 0.7, 0.6, 1.0}, { 0.3f, 0.4, 0.4, 1.0} };
  float ambient[] = { 0.2, 0.2, 0.2, 1.0};
  float specular[]= { 0.5, 0.5, 0.5, 1.0};
  glMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
  glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
  glMaterialf(GL_FRONT,GL_SHININESS,100);

  glNormal3f(0.0, 1.0, 0.0);
  glPushMatrix();
  for (j = 0; j < nz; j++) 
  {
    float z1 = -widthZ / 2.0 + wZ * j; float z2 = z1 + wZ;
    for (i = 0; i < nx; i++) {
      float x1 = -widthX / 2.0 + wX * i; float x2 = x1 + wX;

      glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse[(i + j) & 1]);
			glBegin(GL_QUADS);
      glVertex3f(x1, 0.0, z1);
      glVertex3f(x1, 0.0, z2);
      glVertex3f(x2, 0.0, z2);
      glVertex3f(x2, 0.0, z1);
			glEnd();
    }
  }
  glPopMatrix();
}

//-----------------------------------------------------------------------------
/*
//四角形のメッシュ（x-y平面,中心が原点）
//ｘ軸方向，ｙ軸方向の幅を固定
void drawGrid(int Nx, int Ny, float sizeX, float sizeY, int sideFlag)
{
	//全体の幅,長さsizeX, sizeY
	//sideFlag = 0:側面表示せず
	//sideFlag = 1:側面表示する
	const int NMAX = 130;
	int i, j;
	float p[NMAX][NMAX][3]; //頂点座標
	float pitchX, pitchY;
//	float n1[3], n2[3], n3[3], n4[3];

	if(Nx > NMAX) printf("NxがNMAXを超えています(drawGrid) \n");
	if(Ny > NMAX) printf("NyがNMAXを超えています(drawGrid) \n");

	//セルのサイズ
	pitchX = sizeX / (float)Nx;
	pitchY = sizeY / (float)Ny;

	//各頂点の座標
	for(j = 0; j <= Ny; j++){
		for(i = 0; i <= Nx; i++){
			p[i][j][0] = (float)(i - Nx / 2) * pitchX;
			p[i][j][1] = (float)(j - Ny / 2) * pitchY;
			p[i][j][2] = 0.0;
		}
	}
		
	//三角形で面を定義
	glNormal3f(0.0, 0.0, 1.0);//すべてz方向
	glBegin(GL_TRIANGLES);
	for(j = 0;j < Ny;j++)
		for(i = 0;i < Nx;i++)
		{
			//左下の三角形
			//各頂点の法線方向,ﾃｸｽﾁｬｰ座標,頂点座標を与える。
			glVertex3fv(p[i][j]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
			glVertex3fv(p[i+1][j]);
			glVertex3fv(p[i][j+1]);
			//右上の三角形
			glVertex3fv(p[i+1][j]);
			glVertex3fv(p[i+1][j+1]);
			glVertex3fv(p[i][j+1]);
		}
	glEnd();

	if(sideFlag == 1)//側面描画
	{
		glBegin(GL_QUADS);
		//+x方向（i=Nx)
		glNormal3f(1.0, 0.0, 0.0);
		for(j = 0; j < Ny; j++)
		{
			glVertex3f(p[Nx][j][0], p[Nx][j][1], 0.0f);
			glVertex3f(p[Nx][j+1][0], p[Nx][j+1][1], 0.0f);
			glVertex3f(p[Nx][j+1][0], p[Nx][j+1][1], p[Nx][j+1][2]);
			glVertex3f(p[Nx][j][0], p[Nx][j][1], p[Nx][j][2]);
		}
		//-x方向（i=0)
		glNormal3f(-1.0, 0.0, 0.0);
		for(j = 0; j < Ny; j++)
		{
			glVertex3f(p[0][j][0], p[0][j][1], 0.0f);
			glVertex3f(p[0][j][0], p[0][j][1], p[0][j][2]);
			glVertex3f(p[0][j+1][0], p[0][j+1][1], p[0][j+1][2]);
			glVertex3f(p[0][j+1][0], p[0][j+1][1], 0.0f);
		}
		//+y方向（j=Ny)
		glNormal3f(0.0, 1.0, 0.0);
		for(i = 0; i < Nx; i++)
		{
			glVertex3f(p[i][Ny][0], p[i][Ny][1], 0.0f);
			glVertex3f(p[i][Ny][0], p[i][Ny][1], p[i][Ny][2]);
			glVertex3f(p[i+1][Ny][0], p[i+1][Ny][1], p[i+1][Ny][2]);
			glVertex3f(p[i+1][Ny][0], p[i+1][Ny][1], 0.0f);
		}
		//-y方向（j=0)
		glNormal3f(0.0, -1.0, 0.0);
		for(i = 0; i < Nx; i++)
		{
			glVertex3f(p[i][0][0], p[i][0][1], 0.0f);
			glVertex3f(p[i+1][0][0], p[i+1][0][1], 0.0f);
			glVertex3f(p[i+1][0][0], p[i+1][0][1], p[i+1][0][2]);
			glVertex3f(p[i][0][0], p[i][0][1], p[i][0][2]);
		}
		glEnd();
	}
}
*/