//myTexTree1.h
//����
//y�����㉺����

//�ő��
//#define MAX_REPEAT   5//������ъe�}�̐��㐔
#define MAX_TRUNK    6//1�̊���3��link����Ȃ�
#define MAX_BRANCH0  4//�e���̐e�}�̌�
#define MAX_BRANCH 250//�e�e�}���番�򂷂�q�}�C���}�̑���

double random(double minV, double maxV)
{//minV�`maxV�̈�l��������
	return( (maxV - minV) * (double)rand() / (double)RAND_MAX ) + minV;
}
class CmyTexTree1
{
public:
	CmyTexTree1(void);
	~CmyTexTree1(void) { };
	void create();
	void draw();
	//�����o�ϐ�
	CVector vPos0;//�@�؁i���̉��[�ʒu�j
	//��
	int numRepeatTrunk ;//���㐔
	float lenTrunk0;//�ŉ��[�̒���(3���̒����Clink�P����1/3)
	float radTrunk0;//�ŉ��[�̔��a
	float rateTrunk;//������(link�P���j
	//�}
	float alpha0;  //����p
	float beta0 ;  //�����ʂ���̕Ίp
	float length0; //�e�}�̒���(�C�j�V�G�[�^�j
	float radius0; //�e�}�̔��a
	float rate   ; //������
	float time;//�o�ߎ���

private:
	//�����o�֐�
  void drawBranch(int noK, int noI, int noBranch);
	void drawLeaf(int noK, int noJ, int noBranch);
	//�����o�ϐ�
	float* diffuse ;//�g�U�����˗�
	float* specular;//���ʌ����˗�
	float* ambient; //�������˗�
	float* shadowDiffuse; //�e�̊g�U��
	float* shadowSpecular; //�e�̋��ʌ�
	float highlight;//����

	//��
	CVector vAxisTrunk[MAX_TRUNK * 3];
	float angTrunk[MAX_TRUNK * 3];//��]�p
	float lenTrunk[MAX_TRUNK * 3];//����
	float radTrunk[MAX_TRUNK * 3];//���a

	//�}
	float alpha[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//����p�iy����]�j
	float beta[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH]; //�e�}�����ʂ���̕Ίp�ix����]�j
	float gamma[MAX_TRUNK][MAX_BRANCH0]; //������}�����ꂷ��e�}�̉�]�p�iy����]�j
	float lengthBranch[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];
	float radiusBranch[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];
	int noGeneration[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];
	int numBranch0;//���̊e����ɂ�����e�}�̌�
	int numBranchAll;//�e�e�}�Ɋ܂܂�鑍�}��
	int numRepeat   ;//�}�̕����

	//�t
	float lengthLeaf[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//z�������T�C�Y
	float widthLeaf[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//x�������T�C�Y
	float betaLeaf[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//�e�}�����ʂ���̕Ίp�ix����]�j
	float amp[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//���敗�U���̐U���i�p�x�j
	float period[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//���敗�U���̎���
	float phase[MAX_TRUNK][MAX_BRANCH0][MAX_BRANCH];//���敗�U���̈ʑ�

};
//----------------------------------------------------------------
CmyTexTree1::CmyTexTree1()
{
	static float diffuse0[] = {0.32f, 0.3f, 0.3f, 1.0f};
  static float ambient0[] = {0.2f, 0.2f, 0.2f, 1.0f};
	static float specular0[] = {0.1f, 0.1f, 0.1f, 1.0f};
//	static float shadowD[] = {0.2f, 0.2f, 0.2f, 0.5f};
//	static float shadowS[] = {0.0f, 0.0f, 0.0f, 1.0f};

	//�}�e���A������
	diffuse = diffuse0;
  ambient = ambient0;
  specular = specular0;
  highlight = 10.0f;

	//��
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
	
	//�}
	numRepeat = 4;
	rate = 0.8 + random(-0.05, 0.05);
	alpha0 = 40.0 + random(-10.0, 10.0) ;
	beta0  =  -10.0 + random(-5.0, 5.0);
	length0 = 0.6 + random(-0.1, 0.1);
	radius0 = 0.04 + random(-0.005, 0.005);
	numBranch0 = 3;//�����番�����e�}��
	//�t

}
//-----------------------------------------------------------------------------
void CmyTexTree1::create()
{
	int i, j, k, m;

	//��
	lenTrunk[0] = lenTrunk0 / 3.0;
	radTrunk[0] = radTrunk0;
	for(k = 1; k < numRepeatTrunk * 3; k++)
	{
		lenTrunk[k] = lenTrunk[k-1] * rateTrunk ;
		radTrunk[k] = radTrunk[k-1] * rateTrunk ;
	}

	//�}
	int numBranch = 0;//�e����̎}��
	float len0 = length0;//�e�}�̒���
	float rad0 = radius0;//�e�}�̔��a
	float len, rad;

	//�e����̎}�̃f�[�^
	for(k = 0; k < numRepeatTrunk; k++)
	{
		len0 *= rateTrunk * (0.8 + random(-0.05, 0.05));
		rad0 *= rateTrunk * (0.8 + random(-0.05, 0.05));
		for(j = 0; j < numBranch0; j++)
		{
			//�e�}
			len = len0;
			rad = rad0;
			lengthBranch[k][j][0] = len;
			radiusBranch[k][j][0] = rad;
			alpha[k][j][0] = random(-10.0, 10.0);
			beta[k][j][0] = 30.0 + random(-10.0, 10.0);
			gamma[k][j] = (360.0 / (float)numBranch0) * float(j) + random(-20.0, 20.0) ;
			noGeneration[k][j][0] = 0;
			//�q�}�C���}
			numBranch = 1;
			numBranchAll = 1;
			for(i = 1; i < numRepeat; i++)
			{
				len *= rate;
				rad *= rate;
				for(m = 0; m < numBranch; m++)
				{
					//����
					lengthBranch[k][j][numBranchAll] = len * random(0.8, 1.2);
					radiusBranch[k][j][numBranchAll] = rad * random(0.8, 1.2);
					alpha[k][j][numBranchAll] = random(-10, 10.0);
					beta[k][j][numBranchAll] = beta0 * random(0.8, 1.2);
					noGeneration[k][j][numBranchAll] = i;
					numBranchAll ++;
					if(numBranchAll > MAX_BRANCH) { printf("�}�������������܂� \n"); return;}
					//����
					lengthBranch[k][j][numBranchAll] = len * random(0.8, 1.2);
					radiusBranch[k][j][numBranchAll] = rad * random(0.8, 1.2);
					alpha[k][j][numBranchAll] = -alpha0 * random(0.8, 1.2);
					beta[k][j][numBranchAll] = beta0 * random(0.8, 1.2);
					noGeneration[k][j][numBranchAll] = i;
					numBranchAll ++;
					if(numBranchAll > MAX_BRANCH) { printf("�}�������������܂� \n"); return;}
					//�E��
					lengthBranch[k][j][numBranchAll] = len * random(0.8, 1.2);
					radiusBranch[k][j][numBranchAll] = rad * random(0.8, 1.2);
					alpha[k][j][numBranchAll] = alpha0 * random(0.8, 1.2);
					beta[k][j][numBranchAll] = beta0 * random(0.8, 1.2);
					noGeneration[k][j][numBranchAll] = i;
					numBranchAll ++;
					if(numBranchAll > MAX_BRANCH) { printf("�}�������������܂� \n"); return;}
				}
				numBranch *= 3;//������̎}��
			}//i
//printf("k=%d, j=%d, nnumBranch=%d, numBranchAll=%d \n",k, j, numBranch, numBranchAll);
		}//j
	}//k

	//�t�̃f�[�^
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

	//�ؑS��
	glPushMatrix();
	glTranslatef(vPos0.x, vPos0.y, vPos0.z);
	//��
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

		//�}
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
	//noK;���̐���ԍ�
	//noJ:�����番�����e�}�ԍ�
	//noBranch:1�̐e�}�W�c�Ɋ܂܂��}�ԍ�

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
	//�e
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

	//�q
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


