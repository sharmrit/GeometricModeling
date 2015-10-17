/*

Functions for learnply

Eugene Zhang, 2005
*/
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "glut.h"
#include <string.h>
#include <string>
#include <fstream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "learnply_io.h"
#include "trackball.h"
#include "tmatrix.h"
#include "eig3.h"
#include "SOIL.h"
void calculate_eigens(float a[][3], int n, float d[], float v[][3]);

static PlyFile *in_ply;
char * filename[11];
unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width=1024;
const int win_height=1024;

double radius_factor = 0.9;
int FirstTime =0;
int display_mode = 0; 
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;	
int ACSIZE = 1; // for antialiasing
int view_mode=0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;
int eulercharacteristics;
int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;
static GLuint texName;
unsigned char *checkImage;
float eigenvector[3][3];
float eigval[3];
icVector3 center0fGravity;
int countXmax,countXmin,countYmax,countYmin,countZmax,countZmin;
float boundbox1[3][8],boundbox2[3][8],boundbox3[3][8];
float L = 0.0f;
int previousColor = 0;
int isSecondTime=0;

struct jitter_struct{
	double x;
	double y;
} jitter_para;

jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, {0.875, 0.125}, 
{0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375}, 
{0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625}, 
{0.125, 0.875}, {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}, };

Polyhedron *poly = NULL;
Polyhedron *polyRegular = NULL;
Polyhedron *polyIrregular = NULL;
Polyhedron *PolyInitial= NULL;
float upperBound;

int selection=0;
float dt = 1;
int iteration=1;


void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);
void CreateTexture();


/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char *argv[])
{
	char *progname;
	int num = 1;
	FILE *this_file;
	FILE *this_file1;

	int choice=0;
	char *filename[11]={"../tempmodels/bunny.ply","../tempmodels/dodecahedron.ply",
		"../tempmodels/dragon.ply","../tempmodels/feline.ply",
		"../tempmodels/happy.ply","../tempmodels/hexahedron.ply",
		"../tempmodels/icosahedron.ply","../tempmodels/octahedron.ply",
		"../tempmodels/sphere.ply","../tempmodels/tetrahedron.ply","../tempmodels/torus.ply"};
	progname = argv[0];

	printf("\n 1. Subdivision\n 2. Smoothing\n Enter Your Selection:");
	scanf("%d",&selection);
	if(selection==1)
	{
		printf("Press the number as shown against the choices :\n 1. Original Mesh\n 2. Loop Subdivision(Regular Subdivision)\n 3. Irregular Subdivision based on threshold value=%f"
			"\n 4. 3D checker board color scheme\n 5. 3D checker board color scheme with Loop Subdivision\n 6. 3D checker board color scheme with Irregular subdivision"
			"\n 7. Highlight valence deficit on the original mesh\n 8. Highlight valence deficit after Loop Subdivision(Reguular Subdivision)\n 9. Highlight valence deficit after Irregular Subdivision\n",upperBound);

		printf("Enter the threshold value of U for Irregular Subdivision");
		scanf("%f",&upperBound);
	}
	else if(selection==2)
	{
		printf("\n Press the number as shown against the choices :"
			"\n 1. Mesh Smoothing with uniform scheme"
			"\n 2. Mesh Smoothing with Cord Scheme"
			"\n 3. Mesh Smoothing with Mean Curvature Flow Scheme"
			"\n 4. Mesh smoothing with Mean Value Coordinate Scheme"
			"\n 5. Mesh Smoothing with Cord Scheme with initial mesh"
			"\n 6. Mesh Smoothing with Mean Curvature Flow Scheme with initial mesh"
			"\n 7. Mesh Smoothing with Mean Value Coordinate Scheme with initial mesh"
			"\n 8. Morse Design based on one point function using Explicit Method "
			"\n 9. Morse Design based on one point function using Gauss Seidal Method"
			"\n a. Morse Design based on two points using Explicit Method "
			"\n b. Morse Design based on two points using Gauss Seidal Method "
			"\n c. Uniform Mesh Smoothing with 3D checkerboard color scheme");
		printf("\n Enter the value of dt:");
		scanf("%f",&dt);
	}

	printf(" Models:\n 1. bunny\n 2. dodecahedron\n 3. dragon\n 4. feline\n 5. happy\n 6. hexahedron\n 7. icosahedron\n 8. octahedron\n 9. sphere\n 10. tetrahedron\n 11. torus\n");


	while(choice>11 || choice==0)
	{
		printf("\nEnter the Model id:");
		scanf("%d", &choice);
		if(choice>11 || choice==0)
		{
			printf("wrong choice. Select Again!!\n");
		}
	}

	this_file = fopen(filename[choice-1], "r");
	poly = new Polyhedron (this_file);

	fclose(this_file);


	mat_ident( rotmat );	

	poly->initialize(); // initialize everything
	//polynew->initialize();

	//polynew->loopsubdivision(poly);

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();



	//int choiceCornerListConstruction;
	//printf("Do you want to construct the Corner List:(Enter 1 for Yes/ 0 for No)");
	//scanf("%d",&choiceCornerListConstruction);
	//if(choiceCornerListConstruction==1)
	//{
	//	printf("\nWait....");
	poly->ConstructCornerList();
	//	printf("\nCornerList Created");
	//	
	//}



	printf("\nWait.... The model will appear soon !!\n\n");

	polyRegular = new Polyhedron(this_file,1);
	polyRegular->initialize();
	fclose(this_file);
	mat_ident( rotmat );

	polyIrregular = new Polyhedron(this_file,2);
	polyIrregular->initialize();
	fclose(this_file);
	mat_ident( rotmat );


	//const char *csFile1 = sFile.c_str();
	this_file1 = fopen(filename[choice-1], "r");
	PolyInitial = new Polyhedron(this_file1);    
	fclose(this_file1);
	mat_ident( rotmat );
	PolyInitial->initialize();




	int newEulerCharacteristics = -1;
	int numOfVertices = -1;
	int numOfEdges = -1;
	int numOfFaces = -1;

	if(poly)
	{
		numOfVertices = poly->nverts;
		numOfEdges = poly->nedges;
		numOfFaces = poly->ntris;
	}

	newEulerCharacteristics = numOfVertices - numOfEdges + numOfFaces;

	printf("The Euler characteristics of the mesh in the model without subdivision is:  %d - %d + %d = %d",numOfVertices,numOfEdges,numOfFaces,newEulerCharacteristics);

	if(polyRegular)
	{
		numOfVertices = polyRegular->nverts;
		numOfEdges = polyRegular->nedges;
		numOfFaces = polyRegular->ntris;
	}

	newEulerCharacteristics = numOfVertices - numOfEdges + numOfFaces;

	printf("\nThe Euler characteristics of the mesh in the model after loop subdivision(Regular) is:  %d - %d + %d = %d",numOfVertices,numOfEdges,numOfFaces,newEulerCharacteristics);

	if(polyIrregular)
	{
		numOfVertices = polyIrregular->nverts;
		numOfEdges = polyIrregular->nedges;
		numOfFaces = polyIrregular->ntris;
	}

	newEulerCharacteristics = numOfVertices - numOfEdges + numOfFaces;

	printf("\nThe Eulercharacteristics of the mesh in the model after Irregular subdivision is:  %d - %d + %d = %d",numOfVertices,numOfEdges,numOfFaces,newEulerCharacteristics);

	poly->calc_bounding_box();
	poly->calc_bounding_box_withnormals();

	eulercharacteristics=poly->max_verts-poly->nedges+poly->max_tris;
	printf("\nEuler Characteristics=%d-%d+%d=%d",poly->max_verts,poly->nedges,poly->max_tris,eulercharacteristics);
	int handles=(1-eulercharacteristics*0.5);
	printf("\nNumber of handles:%d",handles);







	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition (20, 20);
	glutInitWindowSize (win_width, win_height); 
	glutCreateWindow ("Geometric Modeling");
	init ();
	glutKeyboardFunc (keyboard);
	glutDisplayFunc(display); 
	glutMotionFunc (motion);
	glutMouseFunc (mouse);
	glutMainLoop(); 
	poly->finalize();  // finalize everything

	return 0;    /* ANSI C requires main to return int. */
}





void color_mapping(double percentage, double col[3])
{
	if (percentage == 0.0){
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 1.0;
	}
	else if (percentage <= 1.0/3){
		col[0] = 1.0;
		col[1] = 1.0-percentage*3.0;
		col[2] = 1.0-percentage*3.0;
	}
	else if (percentage <= 2.0/3){
		col[0] = 1.0;
		col[1] = percentage*3.0-1.0;
		col[2] = 0.0;
	}
	else if (percentage <= 3.0/3){
		col[0] = 3.0-percentage*3.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
	else {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
}

/******************************************************************************
Read in a polyhedron from a file.
******************************************************************************/

Polyhedron::Polyhedron(FILE *file,int type)
{
	// Normal Mesh
	if(0 == type)
	{
		int i,j;
		int elem_count;
		char *elem_name;

		/*** Read in the original PLY object ***/
		in_ply = read_ply (file);

		for (i = 0; i < in_ply->num_elem_types; i++) {

			/* prepare to read the i'th list of elements */
			elem_name = setup_element_read_ply (in_ply, i, &elem_count);
			if (equal_strings ("vertex", elem_name)) {

				/* create a vertex list to hold all the vertices */
				nverts = max_verts = elem_count;
				vlist = new Vertex *[nverts];

				/* set up for getting vertex elements */

				setup_property_ply (in_ply, &vert_props[0]);
				setup_property_ply (in_ply, &vert_props[1]);
				setup_property_ply (in_ply, &vert_props[2]);
				vert_other = get_other_properties_ply (in_ply, 
					offsetof(Vertex_io,other_props));

				/* grab all the vertex elements */
				for (j = 0; j < nverts; j++) {
					Vertex_io vert;
					get_element_ply (in_ply, (void *) &vert);

					/* copy info from the "vert" structure */
					vlist[j] = new Vertex (vert.x, vert.y, vert.z);
					vlist[j]->other_props = vert.other_props;
				}
			}
			else if (equal_strings ("face", elem_name)) {

				/* create a list to hold all the face elements */
				ntris = max_tris = elem_count;
				tlist = new Triangle *[ntris];

				/* set up for getting face elements */
				setup_property_ply (in_ply, &face_props[0]);
				face_other = get_other_properties_ply (in_ply, offsetof(Face_io,other_props));

				/* grab all the face elements */
				for (j = 0; j < elem_count; j++) {
					Face_io face;
					get_element_ply (in_ply, (void *) &face);

					if (face.nverts != 3) {
						fprintf (stderr, "Face has %d vertices (should be three).\n",
							face.nverts);
						exit (-1);
					}

					/* copy info from the "face" structure */
					tlist[j] = new Triangle;
					tlist[j]->nverts = 3;
					tlist[j]->verts[0] = (Vertex *) face.verts[0];
					tlist[j]->verts[1] = (Vertex *) face.verts[1];
					tlist[j]->verts[2] = (Vertex *) face.verts[2];
					tlist[j]->other_props = face.other_props;
				}
			}
			else
				get_other_element_ply (in_ply);
		}

		/* close the file */
		close_ply (in_ply);

		/* fix up vertex pointers in triangles */
		for (i = 0; i < ntris; i++) {
			tlist[i]->verts[0] = vlist[(int) tlist[i]->verts[0]];
			tlist[i]->verts[1] = vlist[(int) tlist[i]->verts[1]];
			tlist[i]->verts[2] = vlist[(int) tlist[i]->verts[2]];
		}

		/* get rid of triangles that use the same vertex more than once */

		for (i = ntris-1; i >= 0; i--) {

			Triangle *tri = tlist[i];
			Vertex *v0 = tri->verts[0];
			Vertex *v1 = tri->verts[1];
			Vertex *v2 = tri->verts[2];

			if (v0 == v1 || v1 == v2 || v2 == v0) {
				free (tlist[i]);
				ntris--;
				tlist[i] = tlist[ntris];
			}
		}

	}

	// Regular subdivision
	if(type==1)
	{
		int i,j;
		int elem_count;
		char *elem_name;


		this->tlist = new Triangle *[4*poly->ntris];
		ntris = max_tris = 4*poly->ntris;

		Vertex** newvList = new Vertex *[4*poly->nverts];

		nverts = 0;
		int VertexCount = 0;
		for(int i = 0; i<poly->ntris; i++)
		{
			Triangle *tri = poly->tlist[i];
			Vertex *newV0,*newV1,*newV2,*newV3,*newV4,*newV5;

			newV0 = tri->verts[0];
			newV1 = tri->verts[1];
			newV2 = tri->verts[2];
			newV3 = new Vertex((newV0->x + newV1->x)/2 , (newV0->y + newV1->y)/2 ,(newV0->z + newV1->z)/2);
			newV4 = new Vertex((newV1->x + newV2->x)/2 , (newV1->y + newV2->y)/2 ,(newV1->z + newV2->z)/2);
			newV5 = new Vertex((newV2->x + newV0->x)/2 , (newV2->y + newV0->y)/2 ,(newV2->z + newV0->z)/2);


			newV3->ntris = 0;
			newV4->ntris = 0;
			newV5->ntris = 0;

			Vertex *v0n = new Vertex(newV0->x , newV0->y , newV0->z);
			Vertex *v1n = new Vertex(newV1->x , newV1->y , newV1->z);
			Vertex *v2n = new Vertex(newV2->x , newV2->y , newV2->z);

			v0n->ntris = 0;
			v1n->ntris = 0;
			v2n->ntris = 0;

			Vertex *v0=NULL,*v1=NULL,*v2=NULL,*v3=NULL,*v4=NULL,*v5=NULL;


			for(int j=0; j<VertexCount ; j++)
			{
				Vertex *currentV = newvList[j];
				if(currentV->x == newV3->x && currentV->y == newV3->y && currentV->z == newV3->z)
				{
					v3 = currentV;
				}
				if(currentV->x == newV4->x && currentV->y == newV4->y && currentV->z == newV4->z)
				{
					v4 = currentV;
				}
				if(currentV->x == newV5->x && currentV->y == newV5->y && currentV->z == newV5->z)
				{
					v5 = currentV;
				}

				if(currentV->x == v0n->x && currentV->y == v0n->y && currentV->z == v0n->z)
				{
					v0 = currentV;
				}
				if(currentV->x == v1n->x && currentV->y == v1n->y && currentV->z == v1n->z)
				{
					v1 = currentV;
				}
				if(currentV->x == v2n->x && currentV->y == v2n->y && currentV->z == v2n->z)
				{
					v2 = currentV;
				}
			}

			if(v3==NULL)
			{
				newvList[VertexCount++] = newV3;		
				nverts++;	   
			}

			if(v4==NULL)
			{
				newvList[VertexCount++] = newV4;
				nverts++;
			}

			if(v5==NULL)
			{
				newvList[VertexCount++] = newV5;
				nverts++;
			}

			if(v0==NULL)
			{
				newvList[VertexCount++] = v0n;
				nverts++;
			}

			if(v1==NULL)
			{
				newvList[VertexCount++] = v1n;
				nverts++;
			}

			if(v2==NULL)
			{
				newvList[VertexCount++] = v2n;
				nverts++;
			}
		}

		this->vlist = new Vertex *[nverts];

		VertexCount = 0;
		for(int i = 0; i<poly->ntris; i++)
		{
			Triangle *tri = poly->tlist[i];
			Vertex *newV0 = tri->verts[0];
			Vertex *newV1 = tri->verts[1];
			Vertex *newV2 = tri->verts[2];
			Vertex *newV3 = new Vertex((newV0->x + newV1->x)/2 , (newV0->y + newV1->y)/2 ,(newV0->z + newV1->z)/2);
			Vertex *newV4 = new Vertex((newV1->x + newV2->x)/2 , (newV1->y + newV2->y)/2 ,(newV1->z + newV2->z)/2);
			Vertex *newV5 = new Vertex((newV2->x + newV0->x)/2 , (newV2->y + newV0->y)/2 ,(newV2->z + newV0->z)/2);
			newV3->ntris = 0;
			newV4->ntris = 0;
			newV5->ntris = 0;

			Vertex *v0n,*v1n,*v2n;
			v0n = new Vertex(newV0->x , newV0->y , newV0->z);
			v1n = new Vertex(newV1->x , newV1->y , newV1->z);
			v2n = new Vertex(newV2->x , newV2->y , newV2->z);

			v0n->ntris = 0;
			v1n->ntris = 0;
			v2n->ntris = 0;



			Vertex *v3 = NULL,*v4 = NULL,*v5 = NULL,*v0 = NULL,*v1 = NULL,*v2 = NULL;

			for(int j=0; j<VertexCount ; j++)
			{
				Vertex *currentV = vlist[j];
				if(currentV->x == newV3->x && currentV->y == newV3->y && currentV->z == newV3->z)
				{
					v3 = currentV;
				}
				if(currentV->x == newV4->x && currentV->y == newV4->y && currentV->z == newV4->z)
				{
					v4 = currentV;
				}
				if(currentV->x == newV5->x && currentV->y == newV5->y && currentV->z == newV5->z)
				{
					v5 = currentV;
				}

				if(currentV->x == v0n->x && currentV->y == v0n->y && currentV->z == v0n->z)
				{
					v0 = currentV;
				}
				if(currentV->x == v1n->x && currentV->y == v1n->y && currentV->z == v1n->z)
				{
					v1 = currentV;
				}
				if(currentV->x == v2n->x && currentV->y == v2n->y && currentV->z == v2n->z)
				{
					v2 = currentV;
				}
			}

			if(v3==NULL)
			{
				vlist[VertexCount++] = newV3;
				v3 = newV3;
			}
			else
			{
				delete newV3;
			}

			if(v4==NULL)
			{
				vlist[VertexCount++] = newV4;
				v4 = newV4;
			}	  
			else
			{
				delete newV4;
			}

			if(v5==NULL)
			{
				vlist[VertexCount++] = newV5;

				v5 = newV5;
			}
			else
			{
				delete newV5;
			}
			if(v0==NULL)
			{
				vlist[VertexCount++] = v0n;
				v0 = v0n;
			}
			else
			{
				delete v0n;
			}

			if(v1==NULL)
			{
				vlist[VertexCount++] = v1n;
				v1 = v1n;
			}
			else
			{
				delete v1n;
			}
			if(v2 ==NULL)
			{
				vlist[VertexCount++] = v2n;
				v2 = v2n;
			}
			else
			{
				delete v2n;
			}
			this->tlist[4*i] = new Triangle;
			this->tlist[4*i]->nverts = 3;
			this->tlist[4*i]->verts[0] = v0;
			this->tlist[4*i]->verts[1] = v3;
			this->tlist[4*i]->verts[2] = v5;

			this->tlist[4*i +1] = new Triangle;
			this->tlist[4*i+1]->nverts = 3;
			this->tlist[4*i+1]->verts[0] = v3;
			this->tlist[4*i+1]->verts[1] = v1;
			this->tlist[4*i+1]->verts[2] = v4;

			this->tlist[4*i+2] = new Triangle;
			this->tlist[4*i+2]->nverts = 3;
			this->tlist[4*i+2]->verts[0] = v4;
			this->tlist[4*i+2]->verts[1] = v2;
			this->tlist[4*i+2]->verts[2] = v5;

			this->tlist[4*i+3] = new Triangle;
			this->tlist[4*i+3]->nverts = 3;
			this->tlist[4*i+3]->verts[0] = v3;
			this->tlist[4*i+3]->verts[1] = v4;
			this->tlist[4*i+3]->verts[2] = v5;
		}
	}

	// Irregular Subdivision
	if(type==2)
	{

		int elem_count;
		char *elem_name;

		// Tag all the edge whose length < upperbound

		ntris = poly->ntris;
		nverts = 0;

		Vertex** newvList = new Vertex *[4*poly->nverts];
		for(int i = 0; i<poly->nedges; i++)
		{
			Edge *currentEdge = poly->elist[i];

			if(currentEdge->length > upperBound)
			{
				currentEdge->tag = 1;
			}
		}

		for(int i = 0; i<poly->ntris; i++)
		{
			int ShorterEdgeCount=0;
			Triangle *CurrentTriangle =  poly->tlist[i];
			for(int j=0;j<3;j++)
			{
				if( 1 == CurrentTriangle->edges[j]->tag)
				{
					ShorterEdgeCount++;
				}
			}   

			if( ShorterEdgeCount==1)
			{
				ntris += 1;
			}

			if(ShorterEdgeCount==2)
			{
				ntris += 2;
			}

			if(ShorterEdgeCount==3)
			{
				ntris += 3;
			}
		}

		tlist = new Triangle *[ntris];
		max_tris = ntris;
		int VertexCount = 0;
		int triangleCount = 0;
		for(int i = 0; i<poly->ntris; i++)
		{
			int ShorterEdgeCount=0;
			Triangle *CurrentTriangle =  poly->tlist[i];

			Edge* tEdge0 = NULL,*tEdge1 = NULL,*tEdge2 = NULL;

			Vertex *v0=NULL,*v1=NULL,*v2=NULL,*v3=NULL,*v4=NULL,*v5=NULL;

			Vertex *newV0 = CurrentTriangle->verts[0];
			Vertex *newV1 = CurrentTriangle->verts[1];
			Vertex *newV2 = CurrentTriangle->verts[2];

			Vertex *newnV0 = newV0;
			Vertex *newnV1 = newV1;
			Vertex *newnV2 = newV2;


			for(int j=0;j<3;j++)
			{
				if( 1 == CurrentTriangle->edges[j]->tag)
				{
					ShorterEdgeCount++;
					if(ShorterEdgeCount==1)
					{
						tEdge0 = CurrentTriangle->edges[j];
						newnV1 = tEdge0->verts[0];
						newnV2 = tEdge0->verts[1];
						for(int dd=0;dd<3;dd++)
						{
							int vin1 = 0;
							int vin2 = 0;
							if(CurrentTriangle->verts[dd]->x == newnV1->x && CurrentTriangle->verts[dd]->y == newnV1->y && CurrentTriangle->verts[dd]->z == newnV1->z)
							{
								vin1= 1;
							}
							else if(CurrentTriangle->verts[dd]->x == newnV2->x && CurrentTriangle->verts[dd]->y == newnV2->y && CurrentTriangle->verts[dd]->z == newnV2->z)
							{
								vin2 = 1;
							}
							if(0 == vin1 && 0 == vin2)
							{
								newnV0 = CurrentTriangle->verts[dd];
							}		  
						}
					}
					if(ShorterEdgeCount==2)
					{
						tEdge1 = CurrentTriangle->edges[j];
						if(tEdge0->verts[0] == tEdge1->verts[0])
						{
							newnV1 = tEdge1->verts[0];
							newnV0 = tEdge0->verts[1];
							newnV2 = tEdge1->verts[1];
						}
						if(tEdge0->verts[1] == tEdge1->verts[1])
						{
							newnV1 = tEdge1->verts[1];
							newnV0 = tEdge0->verts[0];
							newnV2 = tEdge1->verts[0];
						}
						if(tEdge0->verts[1] == tEdge1->verts[0])
						{
							newnV1 = tEdge1->verts[0];
							newnV0 = tEdge0->verts[0];
							newnV2 = tEdge1->verts[1];
						}
						if(tEdge0->verts[0] == tEdge1->verts[1])
						{
							newnV1 = tEdge1->verts[1];
							newnV0 = tEdge0->verts[1];
							newnV2 = tEdge1->verts[0];
						}

					}
					if(ShorterEdgeCount ==3)
					{
						tEdge2 = CurrentTriangle->edges[j];
						newnV0 = newV0;
						newnV1 = newV1;
						newnV2 = newV2;
					}
				}
			}   

			if(1 == ShorterEdgeCount) 	// if shorterEdgeCount equals 1 then divide the current Triangle into two triangles
			{
				Vertex *v3n = new Vertex(newnV0->x , newnV0->y , newnV0->z);
				Vertex *v4n = new Vertex(newnV1->x , newnV1->y , newnV1->z);
				Vertex *v5n = new Vertex(newnV2->x , newnV2->y , newnV2->z);

				if(tEdge0 != NULL)
				{
					Vertex *v0n = new Vertex((tEdge0->verts[0]->x + tEdge0->verts[1]->x)/2 , (tEdge0->verts[0]->y + tEdge0->verts[1]->y)/2 ,(tEdge0->verts[0]->z + tEdge0->verts[1]->z)/2);

					v0n->ntris = 0;
					for(int k=0; k<VertexCount ; k++)
					{
						Vertex *currentV = newvList[k];
						if(currentV->x == v0n->x && currentV->y == v0n->y && currentV->z == v0n->z)
						{
							v0 = currentV;
						}
						if(currentV->x == v3n->x && currentV->y == v3n->y && currentV->z == v3n->z)
						{
							v3 = currentV;
						}
						if(currentV->x == v4n->x && currentV->y == v4n->y && currentV->z == v4n->z)
						{
							v4 = currentV;
						}
						if(currentV->x == v5n->x && currentV->y == v5n->y && currentV->z == v5n->z)
						{
							v5 = currentV;
						}
					}
					if(v0 ==NULL)
					{
						newvList[VertexCount++] = v0n;
						v0 = v0n;
						nverts++;
					}
					if(v3 == NULL)
					{
						newvList[VertexCount++] = v3n;
						v3 = v3n;
						nverts++;
					}
					if(v4 == NULL)
					{
						newvList[VertexCount++] = v4n;
						v4 = v4n;
						nverts++;
					}
					if(v5 ==NULL)
					{
						newvList[VertexCount++] = v5n;
						v5 = v5n;
						nverts++;
					}
				}
				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v3;
				this->tlist[triangleCount]->verts[1] = v4;
				this->tlist[triangleCount]->verts[2] = v0;
				triangleCount++;

				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v3;
				this->tlist[triangleCount]->verts[1] = v0;
				this->tlist[triangleCount]->verts[2] = v5;
				triangleCount++;
			}

			if(2 == ShorterEdgeCount) // if shorter edge count equals 2 then divide the triangle into 3 triangle
			{

				Vertex *v3n = new Vertex(newnV0->x , newnV0->y , newnV0->z);
				Vertex *v4n = new Vertex(newnV1->x , newnV1->y , newnV1->z);
				Vertex *v5n = new Vertex(newnV2->x , newnV2->y , newnV2->z);

				if(tEdge0 != NULL && tEdge1 != NULL)	// Divide the Triangle into two triangle
				{
					Vertex *v0n = new Vertex((tEdge0->verts[0]->x + tEdge0->verts[1]->x)/2 , (tEdge0->verts[0]->y + tEdge0->verts[1]->y)/2 ,(tEdge0->verts[0]->z + tEdge0->verts[1]->z)/2);

					Vertex *v1n = new Vertex((tEdge1->verts[0]->x + tEdge1->verts[1]->x)/2 , (tEdge1->verts[0]->y + tEdge1->verts[1]->y)/2 ,(tEdge1->verts[0]->z + tEdge1->verts[1]->z)/2);

					v0n->ntris = 0;
					v1n->ntris = 0;

					for(int k=0; k<VertexCount ; k++)
					{
						Vertex *currentV = newvList[k];
						if(currentV->x == v0n->x && currentV->y == v0n->y && currentV->z == v0n->z)
						{
							v0 = currentV;
						}
						if(currentV->x == v1n->x && currentV->y == v1n->y && currentV->z == v1n->z)
						{
							v1 = currentV;
						}

						if(currentV->x == v3n->x && currentV->y == v3n->y && currentV->z == v3n->z)
						{
							v3 = currentV;
						}
						if(currentV->x == v4n->x && currentV->y == v4n->y && currentV->z == v4n->z)
						{
							v4 = currentV;
						}
						if(currentV->x == v5n->x && currentV->y == v5n->y && currentV->z == v5n->z)
						{
							v5 = currentV;
						}
					}
					if(v0  == NULL)
					{
						newvList[VertexCount++] = v0n;
						v0 = v0n;
						nverts++;
					}
					if(v1  ==NULL)
					{
						newvList[VertexCount++] = v1n;
						v1 = v1n;
						nverts++;
					}
					if( v3  ==NULL)
					{
						newvList[VertexCount++] = v3n;
						v3 = v3n;
						nverts++;
					}
					if( v4 == NULL)
					{
						newvList[VertexCount++] = v4n;
						v4 = v4n;
						nverts++;
					}
					if(v5  == NULL)
					{
						newvList[VertexCount++] = v5n;
						v5 = v5n;
						nverts++;
					}
				}


				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v0;
				this->tlist[triangleCount]->verts[1] = v4;
				this->tlist[triangleCount]->verts[2] = v1;
				triangleCount++;

				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v0;
				this->tlist[triangleCount]->verts[1] = v1;
				this->tlist[triangleCount]->verts[2] = v3;
				triangleCount++;

				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v3;
				this->tlist[triangleCount]->verts[1] = v1;
				this->tlist[triangleCount]->verts[2] = v5;
				triangleCount++;
			}

			if(3 == ShorterEdgeCount) // if shorterEdgeCount equals to 3 then divide the Triangle into four triangles
			{

				Vertex *v3n = new Vertex(newnV0->x , newnV0->y , newnV0->z);
				Vertex *v4n = new Vertex(newnV1->x , newnV1->y , newnV1->z);
				Vertex *v5n = new Vertex(newnV2->x , newnV2->y , newnV2->z);

				if(tEdge0 != NULL && tEdge1 != NULL && tEdge2 != NULL) // divide the current Triangle into two triangle
				{
					Vertex *v0n = new Vertex((tEdge0->verts[0]->x + tEdge0->verts[1]->x)/2 , (tEdge0->verts[0]->y + tEdge0->verts[1]->y)/2 ,(tEdge0->verts[0]->z + tEdge0->verts[1]->z)/2);
					Vertex *v1n = new Vertex((tEdge1->verts[0]->x + tEdge1->verts[1]->x)/2 , (tEdge1->verts[0]->y + tEdge1->verts[1]->y)/2 ,(tEdge1->verts[0]->z + tEdge1->verts[1]->z)/2);
					Vertex *v2n = new Vertex((tEdge2->verts[0]->x + tEdge2->verts[1]->x)/2 , (tEdge2->verts[0]->y + tEdge2->verts[1]->y)/2 ,(tEdge2->verts[0]->z + tEdge2->verts[1]->z)/2);

					v0n->ntris = 0;
					v1n->ntris = 0;
					v2n->ntris = 0;

					for(int k=0; k<VertexCount ; k++)
					{
						Vertex *currentV = newvList[k];
						if(currentV->x == v0n->x && currentV->y == v0n->y && currentV->z == v0n->z)
						{
							v0 = currentV;
						}
						if(currentV->x == v1n->x && currentV->y == v1n->y && currentV->z == v1n->z)
						{
							v1 = currentV;
						}
						if(currentV->x == v2n->x && currentV->y == v2n->y && currentV->z == v2n->z)
						{
							v2 = currentV;
						}
						if(currentV->x == v3n->x && currentV->y == v3n->y && currentV->z == v3n->z)
						{
							v3 = currentV;
						}
						if(currentV->x == v4n->x && currentV->y == v4n->y && currentV->z == v4n->z)
						{
							v4 = currentV;
						}
						if(currentV->x == v5n->x && currentV->y == v5n->y && currentV->z == v5n->z)
						{
							v5 = currentV;
						}
					}
					if(v0 ==NULL)
					{
						newvList[VertexCount++] = v0n;
						v0 = v0n;
						nverts++;
					}
					if(v1  == NULL)
					{
						newvList[VertexCount++] = v1n;
						v1 = v1n;
						nverts++;
					}
					if( v2 ==NULL )
					{
						newvList[VertexCount++] = v2n;
						v2 = v2n;
						nverts++;
					}
					if( v3 ==NULL )
					{
						newvList[VertexCount++] = v3n;
						v3 = v3n;
						nverts++;
					}
					if( v4 ==NULL )
					{
						newvList[VertexCount++] = v4n;
						v4 = v4n;
						nverts++;
					}
					if(v5 == NULL )
					{
						newvList[VertexCount++] = v5n;
						v5 = v5n;
						nverts++;
					}
				}


				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v3;
				this->tlist[triangleCount]->verts[1] = v0;
				this->tlist[triangleCount]->verts[2] = v2;
				triangleCount++;

				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v0;
				this->tlist[triangleCount]->verts[1] = v4;
				this->tlist[triangleCount]->verts[2] = v1;
				triangleCount++;

				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v0;
				this->tlist[triangleCount]->verts[1] = v1;
				this->tlist[triangleCount]->verts[2] = v2;
				triangleCount++;

				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v1;
				this->tlist[triangleCount]->verts[1] = v5;
				this->tlist[triangleCount]->verts[2] = v2;
				triangleCount++;
			}

			if(0 == ShorterEdgeCount)
			{
				Vertex *v3n = new Vertex(newnV0->x , newnV0->y , newnV0->z);
				Vertex *v4n = new Vertex(newnV1->x , newnV1->y , newnV1->z);
				Vertex *v5n = new Vertex(newnV2->x , newnV2->y , newnV2->z);

				for(int k=0; k<VertexCount ; k++)
				{
					Vertex *currentV = newvList[k];
					if(currentV->x == v3n->x && currentV->y == v3n->y && currentV->z == v3n->z)
					{
						v3 = currentV;
					}
					if(currentV->x == v4n->x && currentV->y == v4n->y && currentV->z == v4n->z)
					{
						v4 = currentV;
					}
					if(currentV->x == v5n->x && currentV->y == v5n->y && currentV->z == v5n->z)
					{
						v5 = currentV;
					}
				}
				if( v3  == NULL)
				{
					newvList[VertexCount++] = v3n;
					v3 = v3n;
					nverts++;
				}
				if(v4 == NULL)
				{
					newvList[VertexCount++] = v4n;
					v4 = v4n;
					nverts++;
				}
				if(v5  == NULL)
				{
					newvList[VertexCount++] = v5n;
					v5 = v5n;
					nverts++;
				}

				this->tlist[triangleCount] = new Triangle;
				this->tlist[triangleCount]->nverts = 3;
				this->tlist[triangleCount]->verts[0] = v3;
				this->tlist[triangleCount]->verts[1] = v4;
				this->tlist[triangleCount]->verts[2] = v5;
				triangleCount++;

			}
		}
		vlist = new Vertex *[nverts];
		for(int i = 0; i<nverts; i++)
		{
			vlist[i] = newvList[i];
		}
	}
}

void Polyhedron::ConstructCornerList()
{
	int i=0;
	clist = new Corner*[3*ntris];
	nCorners = 3*ntris;
	int cListriangleCount = 0;
	for (i = 0; i<ntris; i++) 
	{
		Triangle *tri = tlist[i];
		Vertex *v0 = tri->verts[0];
		Vertex *v1 = tri->verts[1];
		Vertex *v2 = tri->verts[2];

		Corner *corner0 = NULL;
		Corner *corner1 = NULL;
		Corner *corner2 = NULL;

		Edge *sEdge = NULL;

		Edge *edge0 = NULL;
		Edge *edge1 = NULL;
		Edge *edge2 = NULL;

		// First corner
		clist[cListriangleCount] = new Corner;

		clist[cListriangleCount]->tri = tri;
		clist[cListriangleCount]->vert = v0;
		corner0 = clist[cListriangleCount];

		edge0 = tri->edges[0];
		edge1 = tri->edges[1];
		edge2 = tri->edges[2];

		int pointInEdge0 = 0;
		int pointInEdge1 = 0;
		int pointInEdge2 = 0;

		if( edge0->verts[0]->x == v0->x && edge0->verts[0]->y == v0->y && edge0->verts[0]->z == v0->z)
		{
			pointInEdge0 = 1;
		}

		if( edge0->verts[1]->x == v0->x && edge0->verts[1]->y == v0->y && edge0->verts[1]->z == v0->z)
		{
			pointInEdge0 = 1;
		}



		if( edge1->verts[0]->x == v0->x && edge1->verts[0]->y == v0->y && edge1->verts[0]->z == v0->z)
		{
			pointInEdge1 = 1;
		}

		if( edge1->verts[1]->x == v0->x && edge1->verts[1]->y == v0->y && edge1->verts[1]->z == v0->z)
		{
			pointInEdge1 = 1;
		}


		if( edge2->verts[0]->x == v0->x && edge2->verts[0]->y == v0->y && edge2->verts[0]->z == v0->z)
		{
			pointInEdge2= 1;
		}

		if( edge2->verts[1]->x == v0->x && edge2->verts[1]->y == v0->y && edge2->verts[1]->z == v0->z)
		{
			pointInEdge2 = 1;
		}


		if(pointInEdge0 ==0)
		{
			sEdge = edge0;
		}


		if(pointInEdge1 ==0)
		{
			sEdge = edge1;
		}

		if(pointInEdge2 ==0)
		{
			sEdge = edge2;
		}

		clist[cListriangleCount]->edge = sEdge;

		cListriangleCount++;

		// Second corner

		pointInEdge0 = 0;
		pointInEdge1 = 0;
		pointInEdge2 = 0;
		clist[cListriangleCount] = new Corner;

		clist[cListriangleCount]->tri = tri;
		clist[cListriangleCount]->vert = v1;
		corner1 = clist[cListriangleCount];


		if( edge0->verts[0]->x == v1->x && edge0->verts[0]->y == v1->y && edge0->verts[0]->z == v1->z)
		{
			pointInEdge0 = 1;
		}

		if( edge0->verts[1]->x == v1->x && edge0->verts[1]->y == v1->y && edge0->verts[1]->z == v1->z)
		{
			pointInEdge0 = 1;
		}



		if( edge1->verts[0]->x == v1->x && edge1->verts[0]->y == v1->y && edge1->verts[0]->z == v1->z)
		{
			pointInEdge1 = 1;
		}

		if( edge1->verts[1]->x == v1->x && edge1->verts[1]->y == v1->y && edge1->verts[1]->z == v1->z)
		{
			pointInEdge1 = 1;
		}


		if( edge2->verts[0]->x == v1->x && edge2->verts[0]->y == v1->y && edge2->verts[0]->z == v1->z)
		{
			pointInEdge2= 1;
		}

		if( edge2->verts[1]->x == v1->x && edge2->verts[1]->y == v1->y && edge2->verts[1]->z == v1->z)
		{
			pointInEdge2 = 1;
		}

		if(pointInEdge0  == 0)
		{
			sEdge = edge0;
		}


		if( pointInEdge1 ==0)
		{
			sEdge = edge1;
		}

		if(pointInEdge2==0)
		{
			sEdge = edge2;
		}

		clist[cListriangleCount]->edge = sEdge;
		clist[cListriangleCount]->previous = corner0;
		cListriangleCount++;

		// Third corner
		pointInEdge0 = 0;
		pointInEdge1 = 0;
		pointInEdge2 = 0;
		clist[cListriangleCount] = new Corner;

		clist[cListriangleCount]->tri = tri;
		clist[cListriangleCount]->vert = v2;
		corner2 = clist[cListriangleCount];


		edge0 = tri->edges[0];
		edge1 = tri->edges[1];
		edge2 = tri->edges[2];


		if( edge0->verts[0]->x == v2->x && edge0->verts[0]->y == v2->y && edge0->verts[0]->z == v2->z)
		{
			pointInEdge0 = 1;
		}

		if( edge0->verts[1]->x == v2->x && edge0->verts[1]->y == v2->y && edge0->verts[1]->z == v2->z)
		{
			pointInEdge0 = 1;
		}



		if( edge1->verts[0]->x == v2->x && edge1->verts[0]->y == v2->y && edge1->verts[0]->z == v2->z)
		{
			pointInEdge1 = 1;
		}

		if( edge1->verts[1]->x == v2->x && edge1->verts[1]->y == v2->y && edge1->verts[1]->z == v2->z)
		{
			pointInEdge1 = 1;
		}


		if( edge2->verts[0]->x == v2->x && edge2->verts[0]->y == v2->y && edge2->verts[0]->z == v2->z)
		{
			pointInEdge2= 1;
		}

		if( edge2->verts[1]->x == v2->x && edge2->verts[1]->y == v2->y && edge2->verts[1]->z == v2->z)
		{
			pointInEdge2 = 1;
		}

		if( pointInEdge0 ==0)
		{
			sEdge = edge0;
		}


		if(pointInEdge1 ==0)
		{
			sEdge = edge1;
		}

		if( pointInEdge2  ==0)
		{
			sEdge = edge2;
		}

		clist[cListriangleCount]->edge = sEdge;
		clist[cListriangleCount]->previous = corner1;
		clist[cListriangleCount]->next = corner0;

		corner1->next = corner2;

		corner0->next = corner1;
		corner0->previous = corner2;
		cListriangleCount++;

	}

	for(i = 0; i<3*ntris; i++)
	{
		Corner *CurrentCorner = clist[i];

		for(int j = 0; j<3*ntris; j++)
		{
			Corner *CheckCorner = clist[j];
			if(i != j)
			{
				if(CurrentCorner->edge->verts[0]->x == CheckCorner->edge->verts[0]->x &&
					CurrentCorner->edge->verts[0]->y == CheckCorner->edge->verts[0]->y &&
					CurrentCorner->edge->verts[0]->z == CheckCorner->edge->verts[0]->z &&
					CurrentCorner->edge->verts[1]->x == CheckCorner->edge->verts[1]->x &&
					CurrentCorner->edge->verts[1]->y == CheckCorner->edge->verts[1]->y &&
					CurrentCorner->edge->verts[1]->z == CheckCorner->edge->verts[1]->z)
				{
					CurrentCorner->opposite = CheckCorner;
					//break;
				}
				if(CurrentCorner->edge->verts[0]->x == CheckCorner->edge->verts[1]->x &&
					CurrentCorner->edge->verts[0]->y == CheckCorner->edge->verts[1]->y &&
					CurrentCorner->edge->verts[0]->z == CheckCorner->edge->verts[1]->z &&
					CurrentCorner->edge->verts[1]->x == CheckCorner->edge->verts[0]->x &&
					CurrentCorner->edge->verts[1]->y == CheckCorner->edge->verts[0]->y &&
					CurrentCorner->edge->verts[1]->z == CheckCorner->edge->verts[0]->z)
				{
					CurrentCorner->opposite = CheckCorner; 
					//break;
				}
			}
		}
	}
}

void Polyhedron::MeshSmoothingUniformScheme()
{
	float Xnew = 0.0,Ynew = 0.0,Znew = 0.0;


	for(int i = 0; i<PolyInitial->nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		int SizeOfEEdge =0;
		Vertex* CurrVertex = PolyInitial->vlist[i];
		float Xold = CurrVertex->x;
		float Yold = CurrVertex->y;
		float Zold = CurrVertex->z;

	
		for(int kk = 0 ; kk <nedges; kk++)
		{
			if(CurrVertex->x ==  PolyInitial->elist[kk]->verts[0]->x && CurrVertex->y ==  PolyInitial->elist[kk]->verts[0]->y && CurrVertex->z ==  PolyInitial->elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = PolyInitial->elist[kk];
			}
			if(CurrVertex->x ==  PolyInitial->elist[kk]->verts[1]->x && CurrVertex->y ==  PolyInitial->elist[kk]->verts[1]->y && CurrVertex->z ==  PolyInitial->elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = PolyInitial->elist[kk];
			}
		}


		float SummationX = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x && TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y && TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x && TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y && TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInX1 = selectedVertex->x - CurrVertex->x;
			if(0 != DifferenceInX1)
			{
				SummationX = SummationX + DifferenceInX1;
			}
		}
		
		float averageX = SummationX/SizeOfEEdge;
		Xnew = Xold + dt * averageX;

		for(int yy = 0; yy< poly->nverts;yy++)
		{
			if( poly->vlist[yy]->x == Xold  && poly->vlist[yy]->y == Yold && poly->vlist[yy]->z == Zold
				)
			{
				poly->vlist[yy]->x = Xnew;
			}
		}


		float SummationY = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x && TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y && TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x && TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y && TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInY1 = selectedVertex->y - CurrVertex->y;
			if(0 != DifferenceInY1)
			{
				SummationY = SummationY + DifferenceInY1;
			}

		}

		float averageY = SummationY/SizeOfEEdge;
		Ynew = Yold + dt * averageY;

		for(int yy = 0; yy< poly->nverts;yy++)
		{
			if( poly->vlist[yy]->x == Xold  && poly->vlist[yy]->y == Yold && poly->vlist[yy]->z == Zold)
			{
				poly->vlist[yy]->y = Ynew;
			}
		}


		float SummationZ = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x && TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y && TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x && TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y && TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInZ1 = selectedVertex->z - CurrVertex->z;
			if(0 != DifferenceInZ1)
			{
				SummationZ = SummationZ + DifferenceInZ1;
			}

		}

		float averageZ = SummationZ/SizeOfEEdge;
		Znew = Zold + dt * averageZ;
		for(int yy = 0; yy< poly->nverts;yy++)
		{
			if( poly->vlist[yy]->x == Xold  && poly->vlist[yy]->y == Yold && poly->vlist[yy]->z == Zold )
			{
				poly->vlist[yy]->z = Znew;
			}
		}

	}

	for(int ll = 0 ; ll<nverts; ll++)
	{
		PolyInitial->vlist[ll]->x =  poly->vlist[ll]->x;
		PolyInitial->vlist[ll]->y =  poly->vlist[ll]->y;
		PolyInitial->vlist[ll]->z =  poly->vlist[ll]->z;
	}
	PolyInitial->initialize();

}


void Polyhedron::MeshSmoothingCordScheme()
{
	float Xnew = 0.0;
	float Ynew = 0.0;
	float Znew = 0.0;

	Vertex** vtemplist = new Vertex*[nverts];
	for(int pp=0; pp<nverts ; pp++)
	{
		vtemplist[pp] = new Vertex(0,0,0);
	}


	for(int i = 0; i<nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		int SizeOfEEdge =0;
		Vertex* CurrVertex = vlist[i];
		float Xold = CurrVertex->x;
		float Yold = CurrVertex->y;
		float Zold = CurrVertex->z;


		for(int kk = 0 ; kk <PolyInitial->nedges; kk++)
		{
			if(CurrVertex->x ==  elist[kk]->verts[0]->x && 
				CurrVertex->y ==  elist[kk]->verts[0]->y &&
				CurrVertex->z ==  elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
			if(CurrVertex->x ==  elist[kk]->verts[1]->x && 
				CurrVertex->y == elist[kk]->verts[1]->y &&
				CurrVertex->z == elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
		}
		double sumLength=0;
		for(int g= 0; g<SizeOfEEdge; g++)
		{

			sumLength=sumLength+(1/TemperoryEdgeList[g]->length);

		}


		float SummationX = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInX1 = selectedVertex->x - CurrVertex->x;
			if(0 != DifferenceInX1)
			{
				SummationX = SummationX + DifferenceInX1*(1/TemperoryEdgeList[hh]->length);
			}
		}
	

		float averageX = SummationX/sumLength;
		Xnew = Xold + dt * averageX;


		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				
				vtemplist[yy]->x = Xnew;
			}
		}


	
		float SummationY = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInY1 = selectedVertex->y - CurrVertex->y;
			if(0 != DifferenceInY1)
			{
				SummationY = SummationY + DifferenceInY1*(1/TemperoryEdgeList[hh]->length);
			}
		}
		

		float averageY = SummationY/sumLength;
		Ynew = Yold + dt * averageY;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				
				vtemplist[yy]->y = Ynew;
			}
		}

	
		float SummationZ = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInZ1 = selectedVertex->z - CurrVertex->z;
			if(0 != DifferenceInZ1)
			{
				SummationZ = SummationZ + DifferenceInZ1*(1/TemperoryEdgeList[hh]->length);
			}
		}
		
		float averageZ = SummationZ/sumLength;
		Znew = Zold + dt * averageZ;


		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
			
				vtemplist[yy]->z = Znew;
			}
		}

	}

	for(int ll = 0 ; ll<nverts; ll++)
	{

		poly->vlist[ll]->x = vtemplist[ll]->x; 
		poly->vlist[ll]->y = vtemplist[ll]->y; 
		poly->vlist[ll]->z = vtemplist[ll]->z; 
	}
	
	poly->initialize();

	delete []vtemplist;
}
void Polyhedron::MeshSmoothingCordScheme_WithInitialMesh()
{

	float Xnew = 0.0;
	float Ynew = 0.0;
	float Znew = 0.0;



	Vertex** vtemplist = new Vertex*[nverts];
	for(int pp=0; pp<nverts ; pp++)
	{
		vtemplist[pp] = new Vertex(0,0,0);
	}

	for(int i = 0; i<nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		Edge** TemperoryEdgeListOld = new Edge*[20]; 

		int SizeOfEEdge =0;
		int SizeOfEEdgeOld = 0;
		Vertex* CurrVertex = vlist[i];
		float Xold = CurrVertex->x;
		float Yold = CurrVertex->y;
		float Zold = CurrVertex->z;


		for(int kk = 0 ; kk <nedges; kk++)
		{
			if(CurrVertex->x ==  elist[kk]->verts[0]->x && 
				CurrVertex->y ==  elist[kk]->verts[0]->y &&
				CurrVertex->z ==  elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
				TemperoryEdgeListOld[SizeOfEEdgeOld++] = PolyInitial->elist[kk];
			}
			if(CurrVertex->x == elist[kk]->verts[1]->x && 
				CurrVertex->y == elist[kk]->verts[1]->y &&
				CurrVertex->z == elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
				TemperoryEdgeListOld[SizeOfEEdgeOld++] = PolyInitial->elist[kk];
			}
		}
		double sumLength=0;
		for(int g= 0; g<SizeOfEEdgeOld; g++)
		{
			sumLength=sumLength+(1/TemperoryEdgeListOld[g]->length);
			
		}

	
		float SummationX = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInX1 = selectedVertex->x - CurrVertex->x;
			if(0 != DifferenceInX1)
			{
				
				SummationX = SummationX + DifferenceInX1*(1/TemperoryEdgeListOld[hh]->length);
			}
		}
		
		float averageX = SummationX/sumLength;
		Xnew = Xold + dt * averageX;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				//poly->vlist[yy]->x = Xnew;
				vtemplist[yy]->x = Xnew;
			}
		}


		// For Y
		float SummationY = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInY1 = selectedVertex->y - CurrVertex->y;
			if(0 != DifferenceInY1)
			{
				//SummationY = SummationY + DifferenceInY1*TemperoryEdgeList[hh]->length;
				SummationY = SummationY + DifferenceInY1*(1/TemperoryEdgeListOld[hh]->length);

			}
		}
		// here we update the New Y
		float averageY = SummationY/sumLength;
		Ynew = Yold + dt * averageY;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				//poly->vlist[yy]->y = Ynew;
				vtemplist[yy]->y = Ynew;
			}
		}

		// For Z
		float SummationZ = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float DifferenceInZ1 = selectedVertex->z - CurrVertex->z;
			if(0 != DifferenceInZ1)
			{
				//SummationZ = SummationZ + DifferenceInZ1*TemperoryEdgeList[hh]->length;
				SummationZ = SummationZ + DifferenceInZ1*(1/TemperoryEdgeListOld[hh]->length);
			}
		}
		// here we update the New Z
		float averageZ = SummationZ/sumLength;
		Znew = Zold + dt * averageZ;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				
				vtemplist[yy]->z = Znew;
			}
		}

	}

	for(int yy = 0; yy<poly->nverts;yy++)
	{
		poly->vlist[yy]->x = vtemplist[yy]->x;
		poly->vlist[yy]->y = vtemplist[yy]->y;
		poly->vlist[yy]->z = vtemplist[yy]->z;
	}

	delete []vtemplist;
	poly->initialize();

}

void Polyhedron::MeshSmoothingMeanCurvatureFlowScheme()
{
	float Xnew = 0.0;
	float Ynew = 0.0;
	float Znew = 0.0;
	//float dt = 0.001;
	Vertex** vtemplist = new Vertex*[nverts];
	for(int pp=0; pp<nverts ; pp++)
	{
		vtemplist[pp] = new Vertex(0,0,0);
	}

	for(int i = 0; i<nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		int SizeOfEEdge =0;
		Vertex* CurrVertex = vlist[i];
		float Xold = CurrVertex->x;
		float Yold = CurrVertex->y;
		float Zold = CurrVertex->z;

		// Here , we find the incident edges
		for(int kk = 0 ; kk <nedges; kk++)
		{
			if(CurrVertex->x ==  elist[kk]->verts[0]->x && 
				CurrVertex->y ==  elist[kk]->verts[0]->y &&
				CurrVertex->z ==  elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
			if(CurrVertex->x ==  elist[kk]->verts[1]->x && 
				CurrVertex->y ==  elist[kk]->verts[1]->y &&
				CurrVertex->z ==  elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
		}

		// For X
		float SummationX = 0.0f;
		float SummationY = 0.0f;
		float SummationZ = 0.0f;
		float SumOfWeights = 0.0f;

		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;
			Edge* CurrentEdge = TemperoryEdgeList[hh];

			Edge* Edge1 = NULL;
			Edge* Edge2 = NULL;
			Edge* Edge3 = NULL;
			Edge* Edge4 = NULL;
			int firstfound =0;

			for(int bb = 0; bb<CurrVertex->ntris;bb++)
			{
				Triangle* CurrentTriangle = CurrVertex->tris[bb];
				for(int ll=0;ll<3;ll++)
				{
					if(((CurrentTriangle->edges[ll]->verts[0]->x == CurrentEdge->verts[0]->x &&
						CurrentTriangle->edges[ll]->verts[0]->y == CurrentEdge->verts[0]->y &&
						CurrentTriangle->edges[ll]->verts[0]->z == CurrentEdge->verts[0]->z) &&
						(CurrentTriangle->edges[ll]->verts[1]->x == CurrentEdge->verts[1]->x &&
						CurrentTriangle->edges[ll]->verts[1]->y == CurrentEdge->verts[1]->y &&
						CurrentTriangle->edges[ll]->verts[1]->z == CurrentEdge->verts[1]->z)) ||
						((CurrentTriangle->edges[ll]->verts[0]->x == CurrentEdge->verts[1]->x &&
						CurrentTriangle->edges[ll]->verts[0]->y == CurrentEdge->verts[1]->y &&
						CurrentTriangle->edges[ll]->verts[0]->z == CurrentEdge->verts[1]->z) &&
						(CurrentTriangle->edges[ll]->verts[1]->x == CurrentEdge->verts[0]->x &&
						CurrentTriangle->edges[ll]->verts[1]->y == CurrentEdge->verts[0]->y &&
						CurrentTriangle->edges[ll]->verts[1]->z == CurrentEdge->verts[0]->z)))
					{
						if(0 == firstfound)
						{
							if(0 == ll)
							{
								Edge1 = CurrentTriangle->edges[1];
								Edge2 = CurrentTriangle->edges[2];
							}
							if(1 == ll)
							{
								Edge1 = CurrentTriangle->edges[0];
								Edge2 = CurrentTriangle->edges[2];
							}
							if(2 == ll)
							{
								Edge1 = CurrentTriangle->edges[0];
								Edge2 = CurrentTriangle->edges[1];
							}
							firstfound = 1;
						}
						else if(1 == firstfound)
						{
							if(0 == ll)
							{
								Edge3 = CurrentTriangle->edges[1];
								Edge4 = CurrentTriangle->edges[2];
							}
							if(1 == ll)
							{
								Edge3 = CurrentTriangle->edges[0];
								Edge4 = CurrentTriangle->edges[2];
							}
							if(2 == ll)
							{
								Edge3 = CurrentTriangle->edges[0];
								Edge4 = CurrentTriangle->edges[1];
							}
						}
					}
				}
			}

			// here we got all four edges needed
			// for edge1 and edge2 get to get the proper angle
			float commonx =0.0f;
			float commony =0.0f;
			float commonz =0.0f;

			float edge1x = 0.0f;
			float edge1y = 0.0f;
			float edge1z = 0.0f;

			float edge2x = 0.0f;
			float edge2y = 0.0f;
			float edge2z = 0.0f;

			if(Edge1->verts[0]->x == Edge2->verts[0]->x && 
				Edge1->verts[0]->y == Edge2->verts[0]->y &&
				Edge1->verts[0]->z == Edge2->verts[0]->z)
			{
				commonx = Edge1->verts[0]->x;
				commony = Edge1->verts[0]->y;
				commonz = Edge1->verts[0]->z;

				edge1x = Edge1->verts[1]->x;
				edge1y = Edge1->verts[1]->y;
				edge1z = Edge1->verts[1]->z;

				edge2x = Edge2->verts[1]->x;
				edge2y = Edge2->verts[1]->y;
				edge2z = Edge2->verts[1]->z;
			}
			else if(Edge1->verts[0]->x == Edge2->verts[1]->x && 
				Edge1->verts[0]->y == Edge2->verts[1]->y &&
				Edge1->verts[0]->z == Edge2->verts[1]->z)
			{
				commonx = Edge1->verts[0]->x;
				commony = Edge1->verts[0]->y;
				commonz = Edge1->verts[0]->z;

				edge1x = Edge1->verts[1]->x;
				edge1y = Edge1->verts[1]->y;
				edge1z = Edge1->verts[1]->z;

				edge2x = Edge2->verts[0]->x;
				edge2y = Edge2->verts[0]->y;
				edge2z = Edge2->verts[0]->z;
			}
			else if(Edge1->verts[1]->x == Edge2->verts[0]->x && 
				Edge1->verts[1]->y == Edge2->verts[0]->y &&
				Edge1->verts[1]->z == Edge2->verts[0]->z)
			{
				commonx = Edge1->verts[1]->x;
				commony = Edge1->verts[1]->y;
				commonz = Edge1->verts[1]->z;

				edge1x = Edge1->verts[0]->x;
				edge1y = Edge1->verts[0]->y;
				edge1z = Edge1->verts[0]->z;

				edge2x = Edge2->verts[1]->x;
				edge2y = Edge2->verts[1]->y;
				edge2z = Edge2->verts[1]->z;
			}
			else if(Edge1->verts[1]->x == Edge2->verts[1]->x && 
				Edge1->verts[1]->y == Edge2->verts[1]->y &&
				Edge1->verts[1]->z == Edge2->verts[1]->z)
			{
				commonx = Edge1->verts[1]->x;
				commony = Edge1->verts[1]->y;
				commonz = Edge1->verts[1]->z;

				edge1x = Edge1->verts[0]->x;
				edge1y = Edge1->verts[0]->y;
				edge1z = Edge1->verts[0]->z;

				edge2x = Edge2->verts[0]->x;
				edge2y = Edge2->verts[0]->y;
				edge2z = Edge2->verts[0]->z;
			}

			icVector3 vector1((edge1x-commonx),(edge1y-commony),(edge1z-commonz));
			icVector3 vector2((edge2x-commonx),(edge2y-commony),(edge2z-commonz));
			//icVector3 vector1( (Edge1->verts[0]->x -Edge1->verts[1]->x),(Edge1->verts[0]->y -Edge1->verts[1]->y),(Edge1->verts[0]->z -Edge1->verts[1]->z));
			//icVector3 vector2( (Edge2->verts[0]->x -Edge2->verts[1]->x),(Edge2->verts[0]->y -Edge2->verts[1]->y),(Edge2->verts[0]->z -Edge2->verts[1]->z));

			// for edge3 and edge4 get to get the proper angle

			float edge3x = 0.0f;
			float edge3y = 0.0f;
			float edge3z = 0.0f;

			float edge4x = 0.0f;
			float edge4y = 0.0f;
			float edge4z = 0.0f;

			if(Edge3->verts[0]->x == Edge4->verts[0]->x && 
				Edge3->verts[0]->y == Edge4->verts[0]->y &&
				Edge3->verts[0]->z == Edge4->verts[0]->z)
			{
				commonx = Edge3->verts[0]->x;
				commony = Edge3->verts[0]->y;
				commonz = Edge3->verts[0]->z;

				edge3x = Edge3->verts[1]->x;
				edge3y = Edge3->verts[1]->y;
				edge3z = Edge3->verts[1]->z;

				edge4x = Edge4->verts[1]->x;
				edge4y = Edge4->verts[1]->y;
				edge4z = Edge4->verts[1]->z;
			}
			else if(Edge3->verts[0]->x == Edge4->verts[1]->x && 
				Edge3->verts[0]->y == Edge4->verts[1]->y &&
				Edge3->verts[0]->z == Edge4->verts[1]->z)
			{
				commonx = Edge3->verts[0]->x;
				commony = Edge3->verts[0]->y;
				commonz = Edge3->verts[0]->z;

				edge3x = Edge3->verts[1]->x;
				edge3y = Edge3->verts[1]->y;
				edge3z = Edge3->verts[1]->z;

				edge4x = Edge4->verts[0]->x;
				edge4y = Edge4->verts[0]->y;
				edge4z = Edge4->verts[0]->z;
			}
			else if(Edge3->verts[1]->x == Edge4->verts[0]->x && 
				Edge3->verts[1]->y == Edge4->verts[0]->y &&
				Edge3->verts[1]->z == Edge4->verts[0]->z)
			{
				commonx = Edge3->verts[1]->x;
				commony = Edge3->verts[1]->y;
				commonz = Edge3->verts[1]->z;

				edge3x = Edge3->verts[0]->x;
				edge3y = Edge3->verts[0]->y;
				edge3z = Edge3->verts[0]->z;

				edge4x = Edge4->verts[1]->x;
				edge4y = Edge4->verts[1]->y;
				edge4z = Edge4->verts[1]->z;
			}
			else if(Edge3->verts[1]->x == Edge4->verts[1]->x && 
				Edge3->verts[1]->y == Edge4->verts[1]->y &&
				Edge3->verts[1]->z == Edge4->verts[1]->z)
			{
				commonx = Edge3->verts[1]->x;
				commony = Edge3->verts[1]->y;
				commonz = Edge3->verts[1]->z;

				edge3x = Edge3->verts[0]->x;
				edge3y = Edge3->verts[0]->y;
				edge3z = Edge3->verts[0]->z;

				edge4x = Edge4->verts[0]->x;
				edge4y = Edge4->verts[0]->y;
				edge4z = Edge4->verts[0]->z;
			}

			icVector3 vector3((edge3x-commonx),(edge3y-commony),(edge3z-commonz));
			icVector3 vector4((edge4x-commonx),(edge4y-commony),(edge4z-commonz));

			//icVector3 vector3( (Edge3->verts[0]->x -Edge3->verts[1]->x),(Edge3->verts[0]->y -Edge3->verts[1]->y),(Edge3->verts[0]->z -Edge3->verts[1]->z));
			//icVector3 vector4( (Edge4->verts[0]->x -Edge4->verts[1]->x),(Edge4->verts[0]->y -Edge4->verts[1]->y),(Edge4->verts[0]->z -Edge4->verts[1]->z));

			double dotProduct1 = dot(vector1,vector2);

			double productLength1 = length(vector1)*length(vector2);
			double cosTheta1 = dotProduct1/productLength1;

			double Theta1=acos(cosTheta1);
			double cotTheta1=1/tan(Theta1);


			double dotProduct2 = dot(vector3,vector4);

			double productLength2 = length(vector3)*length(vector4);
			double cosTheta2 = dotProduct2/productLength2;
			double Theta2=acos(cosTheta2);
			double cotTheta2=1/tan(Theta2);


			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float weight = 0.5*(cotTheta1+cotTheta2);

			SumOfWeights = SumOfWeights+weight;


			float DifferenceInX1 = selectedVertex->x - CurrVertex->x;				


			DifferenceInX1 = weight*DifferenceInX1;

			if(0 != DifferenceInX1)
			{
				SummationX = SummationX + DifferenceInX1;
			}

			float DifferenceInY1 = selectedVertex->y - CurrVertex->y;				


			DifferenceInY1 = weight*DifferenceInY1;
			if(0 != DifferenceInY1)
			{
				SummationY = SummationY + DifferenceInY1;
			}

			float DifferenceInZ1 = selectedVertex->z - CurrVertex->z;				

			DifferenceInZ1 = weight*DifferenceInZ1;
			if(0 != DifferenceInZ1)
			{
				SummationZ = SummationZ + DifferenceInZ1;
			}

		}

		
		Xnew = Xold + (dt * SummationX)/SumOfWeights;


		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				vtemplist[yy]->x = Xnew;
			}
		}

		Ynew = Yold + dt * SummationY/SumOfWeights;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				vtemplist[yy]->y = Ynew;
			}
		}

		Znew = Zold + (dt * SummationZ)/SumOfWeights;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				vtemplist[yy]->z = Znew;
			}
		}

	}		
	for(int ll = 0 ; ll<nverts; ll++)
	{
		//PolyInitial->vlist[ll]->x = vtemplist[ll]->x; 
		//PolyInitial->vlist[ll]->y = vtemplist[ll]->y; 
		//PolyInitial->vlist[ll]->z = vtemplist[ll]->z; 

		poly->vlist[ll]->x = vtemplist[ll]->x; 
		poly->vlist[ll]->y = vtemplist[ll]->y; 
		poly->vlist[ll]->z = vtemplist[ll]->z; 
	}
	poly->initialize();

	delete []vtemplist;
}
void Polyhedron::MeshSmoothingMeanCurvatureFlowScheme_WithInitialMesh()
{

	float Xnew = 0.0;
	float Ynew = 0.0;
	float Znew = 0.0;


	Vertex** vtemplist = new Vertex*[nverts];
	for(int pp=0; pp<nverts ; pp++)
	{
		vtemplist[pp] = new Vertex(0,0,0);
	}

	for(int i = 0; i<nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		Edge** TemperoryEdgeListOld = new Edge*[20]; 
		int SizeOfEEdge =0;
		int SizeOfEEdgeOld =0;
		Vertex* CurrVertex = vlist[i];
		Vertex* CurrVertexOld = PolyInitial->vlist[i];

		float Xold = CurrVertex->x;
		float Yold = CurrVertex->y;
		float Zold = CurrVertex->z;

		
		for(int kk = 0 ; kk <nedges; kk++)
		{
			if(CurrVertex->x ==  elist[kk]->verts[0]->x && 
				CurrVertex->y == elist[kk]->verts[0]->y &&
				CurrVertex->z == elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
				TemperoryEdgeListOld[SizeOfEEdgeOld++] = PolyInitial->elist[kk];
			}
			if(CurrVertex->x ==  elist[kk]->verts[1]->x && 
				CurrVertex->y == elist[kk]->verts[1]->y &&
				CurrVertex->z == elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
				TemperoryEdgeListOld[SizeOfEEdgeOld++] = PolyInitial->elist[kk];
			}
		}

	
		float SummationX = 0.0f;
		float SummationY = 0.0f;
		float SummationZ = 0.0f;
		float SumOfWeights = 0.0f;

		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;
			Edge* CurrentEdge = TemperoryEdgeList[hh];
			Edge* CurrentEdgeOld = TemperoryEdgeListOld[hh];

			Edge* Edge1 = NULL;
			Edge* Edge2 = NULL;
			Edge* Edge3 = NULL;
			Edge* Edge4 = NULL;
			int firstfound =0;

			for(int bb = 0; bb<CurrVertex->ntris;bb++)
			{
				Triangle* CurrentTriangle = CurrVertex->tris[bb];
				Triangle* CurrentTriangleOld = CurrVertexOld->tris[bb];

				for(int ll=0;ll<3;ll++)
				{
					if(((CurrentTriangle->edges[ll]->verts[0]->x == CurrentEdge->verts[0]->x &&
						CurrentTriangle->edges[ll]->verts[0]->y == CurrentEdge->verts[0]->y &&
						CurrentTriangle->edges[ll]->verts[0]->z == CurrentEdge->verts[0]->z) &&
						(CurrentTriangle->edges[ll]->verts[1]->x == CurrentEdge->verts[1]->x &&
						CurrentTriangle->edges[ll]->verts[1]->y == CurrentEdge->verts[1]->y &&
						CurrentTriangle->edges[ll]->verts[1]->z == CurrentEdge->verts[1]->z)) ||
						((CurrentTriangle->edges[ll]->verts[0]->x == CurrentEdge->verts[1]->x &&
						CurrentTriangle->edges[ll]->verts[0]->y == CurrentEdge->verts[1]->y &&
						CurrentTriangle->edges[ll]->verts[0]->z == CurrentEdge->verts[1]->z) &&
						(CurrentTriangle->edges[ll]->verts[1]->x == CurrentEdge->verts[0]->x &&
						CurrentTriangle->edges[ll]->verts[1]->y == CurrentEdge->verts[0]->y &&
						CurrentTriangle->edges[ll]->verts[1]->z == CurrentEdge->verts[0]->z)))
					{
						if(0 == firstfound)
						{
							if(0 == ll)
							{
								Edge1 = CurrentTriangleOld->edges[1];
								Edge2 = CurrentTriangleOld->edges[2];
							}
							if(1 == ll)
							{
								Edge1 = CurrentTriangleOld->edges[0];
								Edge2 = CurrentTriangleOld->edges[2];
							}
							if(2 == ll)
							{
								Edge1 = CurrentTriangleOld->edges[0];
								Edge2 = CurrentTriangleOld->edges[1];
							}
							firstfound = 1;
						}
						else if(1 == firstfound)
						{
							if(0 == ll)
							{
								Edge3 = CurrentTriangleOld->edges[1];
								Edge4 = CurrentTriangleOld->edges[2];
							}
							if(1 == ll)
							{
								Edge3 = CurrentTriangleOld->edges[0];
								Edge4 = CurrentTriangleOld->edges[2];
							}
							if(2 == ll)
							{
								Edge3 = CurrentTriangleOld->edges[0];
								Edge4 = CurrentTriangleOld->edges[1];
							}
						}
					}
				}
			}


			float commonx =0.0f;
			float commony =0.0f;
			float commonz =0.0f;

			float edge1x = 0.0f;
			float edge1y = 0.0f;
			float edge1z = 0.0f;

			float edge2x = 0.0f;
			float edge2y = 0.0f;
			float edge2z = 0.0f;

			if(Edge1->verts[0]->x == Edge2->verts[0]->x && 
				Edge1->verts[0]->y == Edge2->verts[0]->y &&
				Edge1->verts[0]->z == Edge2->verts[0]->z)
			{
				commonx = Edge1->verts[0]->x;
				commony = Edge1->verts[0]->y;
				commonz = Edge1->verts[0]->z;

				edge1x = Edge1->verts[1]->x;
				edge1y = Edge1->verts[1]->y;
				edge1z = Edge1->verts[1]->z;

				edge2x = Edge2->verts[1]->x;
				edge2y = Edge2->verts[1]->y;
				edge2z = Edge2->verts[1]->z;
			}
			else if(Edge1->verts[0]->x == Edge2->verts[1]->x && 
				Edge1->verts[0]->y == Edge2->verts[1]->y &&
				Edge1->verts[0]->z == Edge2->verts[1]->z)
			{
				commonx = Edge1->verts[0]->x;
				commony = Edge1->verts[0]->y;
				commonz = Edge1->verts[0]->z;

				edge1x = Edge1->verts[1]->x;
				edge1y = Edge1->verts[1]->y;
				edge1z = Edge1->verts[1]->z;

				edge2x = Edge2->verts[0]->x;
				edge2y = Edge2->verts[0]->y;
				edge2z = Edge2->verts[0]->z;
			}
			else if(Edge1->verts[1]->x == Edge2->verts[0]->x && 
				Edge1->verts[1]->y == Edge2->verts[0]->y &&
				Edge1->verts[1]->z == Edge2->verts[0]->z)
			{
				commonx = Edge1->verts[1]->x;
				commony = Edge1->verts[1]->y;
				commonz = Edge1->verts[1]->z;

				edge1x = Edge1->verts[0]->x;
				edge1y = Edge1->verts[0]->y;
				edge1z = Edge1->verts[0]->z;

				edge2x = Edge2->verts[1]->x;
				edge2y = Edge2->verts[1]->y;
				edge2z = Edge2->verts[1]->z;
			}
			else if(Edge1->verts[1]->x == Edge2->verts[1]->x && 
				Edge1->verts[1]->y == Edge2->verts[1]->y &&
				Edge1->verts[1]->z == Edge2->verts[1]->z)
			{
				commonx = Edge1->verts[1]->x;
				commony = Edge1->verts[1]->y;
				commonz = Edge1->verts[1]->z;

				edge1x = Edge1->verts[0]->x;
				edge1y = Edge1->verts[0]->y;
				edge1z = Edge1->verts[0]->z;

				edge2x = Edge2->verts[0]->x;
				edge2y = Edge2->verts[0]->y;
				edge2z = Edge2->verts[0]->z;
			}

			icVector3 vector1((edge1x-commonx),(edge1y-commony),(edge1z-commonz));
			icVector3 vector2((edge2x-commonx),(edge2y-commony),(edge2z-commonz));


			float edge3x = 0.0f;
			float edge3y = 0.0f;
			float edge3z = 0.0f;

			float edge4x = 0.0f;
			float edge4y = 0.0f;
			float edge4z = 0.0f;

			if(Edge3->verts[0]->x == Edge4->verts[0]->x && 
				Edge3->verts[0]->y == Edge4->verts[0]->y &&
				Edge3->verts[0]->z == Edge4->verts[0]->z)
			{
				commonx = Edge3->verts[0]->x;
				commony = Edge3->verts[0]->y;
				commonz = Edge3->verts[0]->z;

				edge3x = Edge3->verts[1]->x;
				edge3y = Edge3->verts[1]->y;
				edge3z = Edge3->verts[1]->z;

				edge4x = Edge4->verts[1]->x;
				edge4y = Edge4->verts[1]->y;
				edge4z = Edge4->verts[1]->z;
			}
			else if(Edge3->verts[0]->x == Edge4->verts[1]->x && 
				Edge3->verts[0]->y == Edge4->verts[1]->y &&
				Edge3->verts[0]->z == Edge4->verts[1]->z)
			{
				commonx = Edge3->verts[0]->x;
				commony = Edge3->verts[0]->y;
				commonz = Edge3->verts[0]->z;

				edge3x = Edge3->verts[1]->x;
				edge3y = Edge3->verts[1]->y;
				edge3z = Edge3->verts[1]->z;

				edge4x = Edge4->verts[0]->x;
				edge4y = Edge4->verts[0]->y;
				edge4z = Edge4->verts[0]->z;
			}
			else if(Edge3->verts[1]->x == Edge4->verts[0]->x && 
				Edge3->verts[1]->y == Edge4->verts[0]->y &&
				Edge3->verts[1]->z == Edge4->verts[0]->z)
			{
				commonx = Edge3->verts[1]->x;
				commony = Edge3->verts[1]->y;
				commonz = Edge3->verts[1]->z;

				edge3x = Edge3->verts[0]->x;
				edge3y = Edge3->verts[0]->y;
				edge3z = Edge3->verts[0]->z;

				edge4x = Edge4->verts[1]->x;
				edge4y = Edge4->verts[1]->y;
				edge4z = Edge4->verts[1]->z;
			}
			else if(Edge3->verts[1]->x == Edge4->verts[1]->x && 
				Edge3->verts[1]->y == Edge4->verts[1]->y &&
				Edge3->verts[1]->z == Edge4->verts[1]->z)
			{
				commonx = Edge3->verts[1]->x;
				commony = Edge3->verts[1]->y;
				commonz = Edge3->verts[1]->z;

				edge3x = Edge3->verts[0]->x;
				edge3y = Edge3->verts[0]->y;
				edge3z = Edge3->verts[0]->z;

				edge4x = Edge4->verts[0]->x;
				edge4y = Edge4->verts[0]->y;
				edge4z = Edge4->verts[0]->z;
			}

			icVector3 vector3((edge3x-commonx),(edge3y-commony),(edge3z-commonz));
			icVector3 vector4((edge4x-commonx),(edge4y-commony),(edge4z-commonz));


			double dotProduct1 = dot(vector1,vector2);

			double productLength1 = length(vector1)*length(vector2);
			double cosTheta1 = dotProduct1/productLength1;
			double Theta1=acos(cosTheta1);
			double cotTheta1=1/tan(Theta1);


			double dotProduct2 = dot(vector3,vector4);

			double productLength2 = length(vector3)*length(vector4);
			double cosTheta2 = dotProduct2/productLength2;
			double Theta2=acos(cosTheta2);
			double cotTheta2=1/tan(Theta2);


			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}


			float weight = 0.5*(cotTheta1+cotTheta2);

			SumOfWeights = SumOfWeights+weight;


			float DifferenceInX1 = selectedVertex->x - CurrVertex->x;				


			DifferenceInX1 = weight*DifferenceInX1;

			if(0 != DifferenceInX1)
			{
				SummationX = SummationX + DifferenceInX1;
			}

			float DifferenceInY1 = selectedVertex->y - CurrVertex->y;				


			DifferenceInY1 = weight*DifferenceInY1;
			if(0 != DifferenceInY1)
			{
				SummationY = SummationY + DifferenceInY1;
			}

			float DifferenceInZ1 = selectedVertex->z - CurrVertex->z;				

			DifferenceInZ1 = weight*DifferenceInZ1;
			if(0 != DifferenceInZ1)
			{
				SummationZ = SummationZ + DifferenceInZ1;
			}

		}

		
		Xnew = Xold + (dt * SummationX)/SumOfWeights;
		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				vtemplist[yy]->x = Xnew;
			}
		}

		Ynew = Yold + dt * SummationY/SumOfWeights;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				vtemplist[yy]->y = Ynew;
			}
		}

		Znew = Zold + (dt * SummationZ)/SumOfWeights;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				vtemplist[yy]->z = Znew;
			}
		}

	}		
	for(int yy = 0; yy<poly->nverts;yy++)
	{
		poly->vlist[yy]->x = vtemplist[yy]->x;
		poly->vlist[yy]->y = vtemplist[yy]->y;
		poly->vlist[yy]->z = vtemplist[yy]->z;
	}

	delete []vtemplist;
	poly->initialize();

}


void Polyhedron::MeshSmoothingMeanValueScheme()
{

	float Xnew = 0.0;
	float Ynew = 0.0;
	float Znew = 0.0;

	Vertex** vtemplist = new Vertex*[nverts];
	for(int pp=0; pp<nverts ; pp++)
	{
		vtemplist[pp] = new Vertex(0,0,0);
	}


	for(int i = 0; i<nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		int SizeOfEEdge =0;
		Vertex* CurrVertex = vlist[i];
		float Xold = CurrVertex->x;
		float Yold = CurrVertex->y;
		float Zold = CurrVertex->z;


		for(int kk = 0 ; kk <nedges; kk++)
		{
			if(CurrVertex->x ==  elist[kk]->verts[0]->x && 
				CurrVertex->y ==  elist[kk]->verts[0]->y &&
				CurrVertex->z ==  elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
			if(CurrVertex->x ==  elist[kk]->verts[1]->x && 
				CurrVertex->y ==  elist[kk]->verts[1]->y &&
				CurrVertex->z ==  elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
		}




		float SummationX = 0.0f;
		float SummationY = 0.0f;
		float SummationZ = 0.0f;
		float SumOfWeights = 0.0f;

		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;
			Edge* CurrentEdge = TemperoryEdgeList[hh];

			Triangle** tempTriangleList = CurrentEdge->tris;
			Edge* Edge1 = NULL;
			Edge* Edge2 = NULL;

			for(int mm=0;mm<3;mm++)
			{

				if(tempTriangleList[0]->edges[mm]->verts[0]->x == CurrVertex->x &&
					tempTriangleList[0]->edges[mm]->verts[0]->y == CurrVertex->y &&
					tempTriangleList[0]->edges[mm]->verts[0]->z == CurrVertex->z 
					)
				{
					if(tempTriangleList[0]->edges[mm]->verts[1]->x == CurrentEdge->verts[1]->x &&
						tempTriangleList[0]->edges[mm]->verts[1]->y == CurrentEdge->verts[1]->y &&
						tempTriangleList[0]->edges[mm]->verts[1]->z == CurrentEdge->verts[1]->z 
						)
					{
						Edge1 = NULL;
					}
					else
					{
						Edge1 = tempTriangleList[0]->edges[mm];
						break;
					}
				}

				if(tempTriangleList[0]->edges[mm]->verts[1]->x == CurrVertex->x &&
					tempTriangleList[0]->edges[mm]->verts[1]->y == CurrVertex->y &&
					tempTriangleList[0]->edges[mm]->verts[1]->z == CurrVertex->z 
					)
				{
					if(tempTriangleList[0]->edges[mm]->verts[0]->x == CurrentEdge->verts[0]->x &&
						tempTriangleList[0]->edges[mm]->verts[0]->y == CurrentEdge->verts[0]->y &&
						tempTriangleList[0]->edges[mm]->verts[0]->z == CurrentEdge->verts[0]->z 
						)
					{
						Edge1 = NULL;
					}
					else
					{
						Edge1 = tempTriangleList[0]->edges[mm];
						break;
					}
				}
			}

			for(int mm=0;mm<3;mm++)
			{
				if(tempTriangleList[1]->edges[mm]->verts[0]->x == CurrVertex->x &&
					tempTriangleList[1]->edges[mm]->verts[0]->y == CurrVertex->y &&
					tempTriangleList[1]->edges[mm]->verts[0]->z == CurrVertex->z 
					)
				{
					if(tempTriangleList[1]->edges[mm]->verts[1]->x == CurrentEdge->verts[1]->x &&
						tempTriangleList[1]->edges[mm]->verts[1]->y == CurrentEdge->verts[1]->y &&
						tempTriangleList[1]->edges[mm]->verts[1]->z == CurrentEdge->verts[1]->z 
						)
					{
						Edge2 = NULL;
					}
					else
					{
						Edge2 = tempTriangleList[1]->edges[mm];
						break;
					}
				}

				if(tempTriangleList[1]->edges[mm]->verts[1]->x == CurrVertex->x &&
					tempTriangleList[1]->edges[mm]->verts[1]->y == CurrVertex->y &&
					tempTriangleList[1]->edges[mm]->verts[1]->z == CurrVertex->z 
					)
				{
					if(tempTriangleList[1]->edges[mm]->verts[0]->x == CurrentEdge->verts[0]->x &&
						tempTriangleList[1]->edges[mm]->verts[0]->y == CurrentEdge->verts[0]->y &&
						tempTriangleList[1]->edges[mm]->verts[0]->z == CurrentEdge->verts[0]->z 
						)
					{
						Edge2 = NULL;
					}
					else
					{
						Edge2 = tempTriangleList[1]->edges[mm];
						break;
					}
				}
			}


			float commonx =0.0f;
			float commony =0.0f;
			float commonz =0.0f;

			float edge0x = 0.0f;
			float edge0y = 0.0f;
			float edge0z = 0.0f;

			float edge1x = 0.0f;
			float edge1y = 0.0f;
			float edge1z = 0.0f;

			float edge2x = 0.0f;
			float edge2y = 0.0f;
			float edge2z = 0.0f;

			if(Edge1->verts[0]->x == Edge2->verts[0]->x && 
				Edge1->verts[0]->y == Edge2->verts[0]->y &&
				Edge1->verts[0]->z == Edge2->verts[0]->z)
			{
				commonx = Edge1->verts[0]->x;
				commony = Edge1->verts[0]->y;
				commonz = Edge1->verts[0]->z;

				edge1x = Edge1->verts[1]->x;
				edge1y = Edge1->verts[1]->y;
				edge1z = Edge1->verts[1]->z;

				edge2x = Edge2->verts[1]->x;
				edge2y = Edge2->verts[1]->y;
				edge2z = Edge2->verts[1]->z;
				if(CurrentEdge->verts[0]->x == commonx &&
					CurrentEdge->verts[0]->y == commony &&
					CurrentEdge->verts[0]->z == commonz)
				{
					edge0x = CurrentEdge->verts[1]->x;
					edge0y = CurrentEdge->verts[1]->y;
					edge0z = CurrentEdge->verts[1]->z;
				}
				else
				{
					edge0x = CurrentEdge->verts[0]->x;
					edge0y = CurrentEdge->verts[0]->y;
					edge0z = CurrentEdge->verts[0]->z;

				}
			}
			else if(Edge1->verts[0]->x == Edge2->verts[1]->x && 
				Edge1->verts[0]->y == Edge2->verts[1]->y &&
				Edge1->verts[0]->z == Edge2->verts[1]->z)
			{
				commonx = Edge1->verts[0]->x;
				commony = Edge1->verts[0]->y;
				commonz = Edge1->verts[0]->z;

				edge1x = Edge1->verts[1]->x;
				edge1y = Edge1->verts[1]->y;
				edge1z = Edge1->verts[1]->z;

				edge2x = Edge2->verts[0]->x;
				edge2y = Edge2->verts[0]->y;
				edge2z = Edge2->verts[0]->z;
				if(CurrentEdge->verts[0]->x == commonx &&
					CurrentEdge->verts[0]->y == commony &&
					CurrentEdge->verts[0]->z == commonz)
				{
					edge0x = CurrentEdge->verts[1]->x;
					edge0y = CurrentEdge->verts[1]->y;
					edge0z = CurrentEdge->verts[1]->z;
				}
				else
				{
					edge0x = CurrentEdge->verts[0]->x;
					edge0y = CurrentEdge->verts[0]->y;
					edge0z = CurrentEdge->verts[0]->z;

				}
			}
			else if(Edge1->verts[1]->x == Edge2->verts[0]->x && 
				Edge1->verts[1]->y == Edge2->verts[0]->y &&
				Edge1->verts[1]->z == Edge2->verts[0]->z)
			{
				commonx = Edge1->verts[1]->x;
				commony = Edge1->verts[1]->y;
				commonz = Edge1->verts[1]->z;

				edge1x = Edge1->verts[0]->x;
				edge1y = Edge1->verts[0]->y;
				edge1z = Edge1->verts[0]->z;

				edge2x = Edge2->verts[1]->x;
				edge2y = Edge2->verts[1]->y;
				edge2z = Edge2->verts[1]->z;
				if(CurrentEdge->verts[0]->x == commonx &&
					CurrentEdge->verts[0]->y == commony &&
					CurrentEdge->verts[0]->z == commonz)
				{
					edge0x = CurrentEdge->verts[1]->x;
					edge0y = CurrentEdge->verts[1]->y;
					edge0z = CurrentEdge->verts[1]->z;
				}
				else
				{
					edge0x = CurrentEdge->verts[0]->x;
					edge0y = CurrentEdge->verts[0]->y;
					edge0z = CurrentEdge->verts[0]->z;

				}
			}
			else if(Edge1->verts[1]->x == Edge2->verts[1]->x && 
				Edge1->verts[1]->y == Edge2->verts[1]->y &&
				Edge1->verts[1]->z == Edge2->verts[1]->z)
			{
				commonx = Edge1->verts[1]->x;
				commony = Edge1->verts[1]->y;
				commonz = Edge1->verts[1]->z;

				edge1x = Edge1->verts[0]->x;
				edge1y = Edge1->verts[0]->y;
				edge1z = Edge1->verts[0]->z;

				edge2x = Edge2->verts[0]->x;
				edge2y = Edge2->verts[0]->y;
				edge2z = Edge2->verts[0]->z;
				if(CurrentEdge->verts[0]->x == commonx &&
					CurrentEdge->verts[0]->y == commony &&
					CurrentEdge->verts[0]->z == commonz)
				{
					edge0x = CurrentEdge->verts[1]->x;
					edge0y = CurrentEdge->verts[1]->y;
					edge0z = CurrentEdge->verts[1]->z;
				}
				else
				{
					edge0x = CurrentEdge->verts[0]->x;
					edge0y = CurrentEdge->verts[0]->y;
					edge0z = CurrentEdge->verts[0]->z;

				}
			}


			icVector3 vector1((edge1x-commonx),(edge1y-commony),(edge1z-commonz));
			icVector3 vector2((edge2x-commonx),(edge2y-commony),(edge2z-commonz));
			icVector3 vector3((edge0x-commonx),(edge0y-commony),(edge0z-commonz));


			double dotProduct1 = dot(vector1,vector3);

			double productLength1 = length(vector1)*length(vector3);
			double cosTheta1 = dotProduct1/productLength1;
			double Theta1 = acos(cosTheta1);
			double tantheta1 = tan(Theta1/2);


			double dotProduct2 = dot(vector2,vector3);

			double productLength2 = length(vector2)*length(vector3);
			double cosTheta2 = dotProduct2/productLength2;
			double Theta2 = acos(cosTheta2);
			double tantheta2 = tan(Theta2/2);


			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}



			float weight = 0.5*(tantheta1+tantheta2);
			if(0>weight)
			{
				weight = 0;
			}

			SumOfWeights = SumOfWeights+weight;


			float DifferenceInX1 = selectedVertex->x - CurrVertex->x;				


			DifferenceInX1 = weight*DifferenceInX1;

			if(0 != DifferenceInX1)
			{
				SummationX = SummationX + DifferenceInX1;
			}

			float DifferenceInY1 = selectedVertex->y - CurrVertex->y;				

			DifferenceInY1 = weight*DifferenceInY1;
			if(0 != DifferenceInY1)
			{
				SummationY = SummationY + DifferenceInY1;
			}

			float DifferenceInZ1 = selectedVertex->z - CurrVertex->z;				
			DifferenceInZ1 = weight*DifferenceInZ1;
			if(0 != DifferenceInZ1)
			{
				SummationZ = SummationZ + DifferenceInZ1;
			}

		}

		
		Xnew = Xold + (dt * SummationX)/SumOfWeights;
		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				
				vtemplist[yy]->x = Xnew;
			}
		}

	
		Ynew = Yold + (dt * SummationY)/SumOfWeights;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				
				vtemplist[yy]->y = Ynew;
			}
		}

	
		Znew = Zold + (dt * SummationZ)/SumOfWeights;
		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				
				vtemplist[yy]->z = Znew;
			}
		}
	}		

	for(int ll = 0 ; ll<nverts; ll++)
	{

		poly->vlist[ll]->x = vtemplist[ll]->x; 
		poly->vlist[ll]->y = vtemplist[ll]->y; 
		poly->vlist[ll]->z = vtemplist[ll]->z; 
	}

	poly->initialize();

	delete []vtemplist;
}
void Polyhedron::MeshSmoothingMeanValueScheme_WithInitialMesh()
{
	
	float Xnew = 0.0;
	float Ynew = 0.0;
	float Znew = 0.0;

	Vertex** vtemplist = new Vertex*[nverts];
	for(int pp=0; pp<nverts ; pp++)
	{
		vtemplist[pp] = new Vertex(0,0,0);
	}


	for(int i = 0; i<PolyInitial->nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		Edge** TemperoryEdgeListOld = new Edge*[20]; 
		int SizeOfEEdge =0;
		int SizeOfEEdgeOld =0;
		Vertex* CurrVertex = vlist[i];
		Vertex* CurrVertexOld = PolyInitial->vlist[i];

		float Xold = CurrVertex->x;
		float Yold = CurrVertex->y;
		float Zold = CurrVertex->z;

	
		for(int kk = 0 ; kk <nedges; kk++)
		{
			if(CurrVertex->x ==  elist[kk]->verts[0]->x && 
				CurrVertex->y ==  elist[kk]->verts[0]->y &&
				CurrVertex->z ==  elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
				TemperoryEdgeListOld[SizeOfEEdgeOld++] = PolyInitial->elist[kk];
			}
			if(CurrVertex->x ==  elist[kk]->verts[1]->x && 
				CurrVertex->y ==  elist[kk]->verts[1]->y &&
				CurrVertex->z ==  elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
				TemperoryEdgeListOld[SizeOfEEdgeOld++] = PolyInitial->elist[kk];
			}
		}



	
		float SummationX = 0.0f;
		float SummationY = 0.0f;
		float SummationZ = 0.0f;
		float SumOfWeights = 0.0f;

		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;
			Edge* CurrentEdge = TemperoryEdgeList[hh];
			Edge* CurrentEdgeOld = TemperoryEdgeListOld[hh];

			Triangle** tempTriangleList = CurrentEdge->tris;
			Triangle** tempTriangleListOld = CurrentEdgeOld->tris;

			Edge* Edge1 = NULL;
			Edge* Edge2 = NULL;

			for(int mm=0;mm<3;mm++)
			{

				if(tempTriangleList[0]->edges[mm]->verts[0]->x == CurrVertex->x &&
					tempTriangleList[0]->edges[mm]->verts[0]->y == CurrVertex->y &&
					tempTriangleList[0]->edges[mm]->verts[0]->z == CurrVertex->z 
					)
				{
					if(tempTriangleList[0]->edges[mm]->verts[1]->x == CurrentEdge->verts[1]->x &&
						tempTriangleList[0]->edges[mm]->verts[1]->y == CurrentEdge->verts[1]->y &&
						tempTriangleList[0]->edges[mm]->verts[1]->z == CurrentEdge->verts[1]->z 
						)
					{
						Edge1 = NULL;
					}
					else
					{
					
						Edge1 = tempTriangleListOld[0]->edges[mm];
						break;
					}
				}

				if(tempTriangleList[0]->edges[mm]->verts[1]->x == CurrVertex->x &&
					tempTriangleList[0]->edges[mm]->verts[1]->y == CurrVertex->y &&
					tempTriangleList[0]->edges[mm]->verts[1]->z == CurrVertex->z 
					)
				{
					if(tempTriangleList[0]->edges[mm]->verts[0]->x == CurrentEdge->verts[0]->x &&
						tempTriangleList[0]->edges[mm]->verts[0]->y == CurrentEdge->verts[0]->y &&
						tempTriangleList[0]->edges[mm]->verts[0]->z == CurrentEdge->verts[0]->z 
						)
					{
						Edge1 = NULL;
					}
					else
					{
						
						Edge1 = tempTriangleListOld[0]->edges[mm];
						break;
					}
				}
			}

			for(int mm=0;mm<3;mm++)
			{
				if(tempTriangleList[1]->edges[mm]->verts[0]->x == CurrVertex->x &&
					tempTriangleList[1]->edges[mm]->verts[0]->y == CurrVertex->y &&
					tempTriangleList[1]->edges[mm]->verts[0]->z == CurrVertex->z 
					)
				{
					if(tempTriangleList[1]->edges[mm]->verts[1]->x == CurrentEdge->verts[1]->x &&
						tempTriangleList[1]->edges[mm]->verts[1]->y == CurrentEdge->verts[1]->y &&
						tempTriangleList[1]->edges[mm]->verts[1]->z == CurrentEdge->verts[1]->z 
						)
					{
						Edge2 = NULL;
					}
					else
					{
						
						Edge2 = tempTriangleListOld[1]->edges[mm];
						break;
					}
				}

				if(tempTriangleList[1]->edges[mm]->verts[1]->x == CurrVertex->x &&
					tempTriangleList[1]->edges[mm]->verts[1]->y == CurrVertex->y &&
					tempTriangleList[1]->edges[mm]->verts[1]->z == CurrVertex->z 
					)
				{
					if(tempTriangleList[1]->edges[mm]->verts[0]->x == CurrentEdge->verts[0]->x &&
						tempTriangleList[1]->edges[mm]->verts[0]->y == CurrentEdge->verts[0]->y &&
						tempTriangleList[1]->edges[mm]->verts[0]->z == CurrentEdge->verts[0]->z 
						)
					{
						Edge2 = NULL;
					}
					else
					{
						
						Edge2 = tempTriangleListOld[1]->edges[mm];
						break;
					}
				}
			}

		
			float commonx =0.0f;
			float commony =0.0f;
			float commonz =0.0f;

			float edge0x = 0.0f;
			float edge0y = 0.0f;
			float edge0z = 0.0f;

			float edge1x = 0.0f;
			float edge1y = 0.0f;
			float edge1z = 0.0f;

			float edge2x = 0.0f;
			float edge2y = 0.0f;
			float edge2z = 0.0f;

			if(Edge1->verts[0]->x == Edge2->verts[0]->x && 
				Edge1->verts[0]->y == Edge2->verts[0]->y &&
				Edge1->verts[0]->z == Edge2->verts[0]->z)
			{
				commonx = Edge1->verts[0]->x;
				commony = Edge1->verts[0]->y;
				commonz = Edge1->verts[0]->z;

				edge1x = Edge1->verts[1]->x;
				edge1y = Edge1->verts[1]->y;
				edge1z = Edge1->verts[1]->z;

				edge2x = Edge2->verts[1]->x;
				edge2y = Edge2->verts[1]->y;
				edge2z = Edge2->verts[1]->z;
				if(CurrentEdge->verts[0]->x == commonx &&
					CurrentEdge->verts[0]->y == commony &&
					CurrentEdge->verts[0]->z == commonz)
				{
					edge0x = CurrentEdge->verts[1]->x;
					edge0y = CurrentEdge->verts[1]->y;
					edge0z = CurrentEdge->verts[1]->z;
				}
				else
				{
					edge0x = CurrentEdge->verts[0]->x;
					edge0y = CurrentEdge->verts[0]->y;
					edge0z = CurrentEdge->verts[0]->z;

				}
			}
			else if(Edge1->verts[0]->x == Edge2->verts[1]->x && 
				Edge1->verts[0]->y == Edge2->verts[1]->y &&
				Edge1->verts[0]->z == Edge2->verts[1]->z)
			{
				commonx = Edge1->verts[0]->x;
				commony = Edge1->verts[0]->y;
				commonz = Edge1->verts[0]->z;

				edge1x = Edge1->verts[1]->x;
				edge1y = Edge1->verts[1]->y;
				edge1z = Edge1->verts[1]->z;

				edge2x = Edge2->verts[0]->x;
				edge2y = Edge2->verts[0]->y;
				edge2z = Edge2->verts[0]->z;
				if(CurrentEdge->verts[0]->x == commonx &&
					CurrentEdge->verts[0]->y == commony &&
					CurrentEdge->verts[0]->z == commonz)
				{
					edge0x = CurrentEdge->verts[1]->x;
					edge0y = CurrentEdge->verts[1]->y;
					edge0z = CurrentEdge->verts[1]->z;
				}
				else
				{
					edge0x = CurrentEdge->verts[0]->x;
					edge0y = CurrentEdge->verts[0]->y;
					edge0z = CurrentEdge->verts[0]->z;

				}
			}
			else if(Edge1->verts[1]->x == Edge2->verts[0]->x && 
				Edge1->verts[1]->y == Edge2->verts[0]->y &&
				Edge1->verts[1]->z == Edge2->verts[0]->z)
			{
				commonx = Edge1->verts[1]->x;
				commony = Edge1->verts[1]->y;
				commonz = Edge1->verts[1]->z;

				edge1x = Edge1->verts[0]->x;
				edge1y = Edge1->verts[0]->y;
				edge1z = Edge1->verts[0]->z;

				edge2x = Edge2->verts[1]->x;
				edge2y = Edge2->verts[1]->y;
				edge2z = Edge2->verts[1]->z;
				if(CurrentEdge->verts[0]->x == commonx &&
					CurrentEdge->verts[0]->y == commony &&
					CurrentEdge->verts[0]->z == commonz)
				{
					edge0x = CurrentEdge->verts[1]->x;
					edge0y = CurrentEdge->verts[1]->y;
					edge0z = CurrentEdge->verts[1]->z;
				}
				else
				{
					edge0x = CurrentEdge->verts[0]->x;
					edge0y = CurrentEdge->verts[0]->y;
					edge0z = CurrentEdge->verts[0]->z;

				}
			}
			else if(Edge1->verts[1]->x == Edge2->verts[1]->x && 
				Edge1->verts[1]->y == Edge2->verts[1]->y &&
				Edge1->verts[1]->z == Edge2->verts[1]->z)
			{
				commonx = Edge1->verts[1]->x;
				commony = Edge1->verts[1]->y;
				commonz = Edge1->verts[1]->z;

				edge1x = Edge1->verts[0]->x;
				edge1y = Edge1->verts[0]->y;
				edge1z = Edge1->verts[0]->z;

				edge2x = Edge2->verts[0]->x;
				edge2y = Edge2->verts[0]->y;
				edge2z = Edge2->verts[0]->z;
				if(CurrentEdge->verts[0]->x == commonx &&
					CurrentEdge->verts[0]->y == commony &&
					CurrentEdge->verts[0]->z == commonz)
				{
					edge0x = CurrentEdge->verts[1]->x;
					edge0y = CurrentEdge->verts[1]->y;
					edge0z = CurrentEdge->verts[1]->z;
				}
				else
				{
					edge0x = CurrentEdge->verts[0]->x;
					edge0y = CurrentEdge->verts[0]->y;
					edge0z = CurrentEdge->verts[0]->z;

				}
			}


			icVector3 vector1((edge1x-commonx),(edge1y-commony),(edge1z-commonz));
			icVector3 vector2((edge2x-commonx),(edge2y-commony),(edge2z-commonz));
			icVector3 vector3((edge0x-commonx),(edge0y-commony),(edge0z-commonz));


			double dotProduct1 = dot(vector1,vector3);

			double productLength1 = length(vector1)*length(vector3);
			double cosTheta1 = dotProduct1/productLength1;
			double Theta1 = acos(cosTheta1);
			double tantheta1 = tan(Theta1/2);


			double dotProduct2 = dot(vector2,vector3);

			double productLength2 = length(vector2)*length(vector3);
			double cosTheta2 = dotProduct2/productLength2;

			double Theta2 = acos(cosTheta2);
			double tantheta2 = tan(Theta2/2);


			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}



			float weight = 0.5*(tantheta1+tantheta2);
			if(0>weight)
			{
				weight = 0;
			}

			SumOfWeights = SumOfWeights+weight;


			float DifferenceInX1 = selectedVertex->x - CurrVertex->x;				

			DifferenceInX1 = weight*DifferenceInX1;

			if(0 != DifferenceInX1)
			{
				SummationX = SummationX + DifferenceInX1;
			}

			float DifferenceInY1 = selectedVertex->y - CurrVertex->y;				


			DifferenceInY1 = weight*DifferenceInY1;
			if(0 != DifferenceInY1)
			{
				SummationY = SummationY + DifferenceInY1;
			}

			float DifferenceInZ1 = selectedVertex->z - CurrVertex->z;				

			DifferenceInZ1 = weight*DifferenceInZ1;
			if(0 != DifferenceInZ1)
			{
				SummationZ = SummationZ + DifferenceInZ1;
			}

		}

		
		Xnew = Xold + (dt * SummationX)/SumOfWeights;
		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				
				vtemplist[yy]->x = Xnew;
			}
		}

	
		Ynew = Yold + (dt * SummationY)/SumOfWeights;

		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  && poly->vlist[yy]->y == Yold && poly->vlist[yy]->z == Zold
				)
			{
			
				vtemplist[yy]->y = Ynew;
			}
		}

	
		Znew = Zold + (dt * SummationZ)/SumOfWeights;
		for(int yy = 0; yy<poly->nverts;yy++)
		{
			if(poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				
				vtemplist[yy]->z = Znew;
			}
		}


	}		
	for(int yy = 0; yy<poly->nverts;yy++)
	{
		poly->vlist[yy]->x = vtemplist[yy]->x;
		poly->vlist[yy]->y = vtemplist[yy]->y;
		poly->vlist[yy]->z = vtemplist[yy]->z;
	}

	delete []vtemplist;
	poly->initialize();

}

void Polyhedron::MorseDesignUsingExplicitMethod()
{
	Vertex** vtemplist = new Vertex*[nverts];
	for(int pp=0; pp<nverts ; pp++)
	{
		vtemplist[pp] = new Vertex(0,0,0);
	}

	float newMorseVal = 0.0;
	float newMorseVal1 = 0.0;

	for(int i = 0; i<nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		int SizeOfEEdge =0;
		Vertex* CurrVertex = vlist[i];
		float OldMorseVal = CurrVertex->MorseVal;
		float OldMorseVal1 = CurrVertex->MorseVal1;
		float Xold = CurrVertex->x;
		float Yold = CurrVertex->y;
		float Zold = CurrVertex->z;

	
		for(int kk = 0 ; kk <nedges; kk++)
		{
			if(CurrVertex->x ==  elist[kk]->verts[0]->x && CurrVertex->y ==  elist[kk]->verts[0]->y && CurrVertex->z ==  elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
			if(CurrVertex->x ==  elist[kk]->verts[1]->x && CurrVertex->y ==  elist[kk]->verts[1]->y && CurrVertex->z ==  elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
		}

		float SumOf = 0.0f;
		float SumOf1 = 0.0f;
		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x && TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y && TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float diff = selectedVertex->MorseVal - CurrVertex->MorseVal;
			if(0 != diff)
			{
				SumOf = SumOf + diff;
			}
			float diff1 = selectedVertex->MorseVal1 - CurrVertex->MorseVal1;
			if(0 != diff1)
			{
				SumOf1 = SumOf1 + diff1;
			}
		}
	
		double average = SumOf/SizeOfEEdge;
		newMorseVal = OldMorseVal + dt * average;

		double average1 = SumOf1/SizeOfEEdge;
		newMorseVal1 = OldMorseVal1 + dt * average1;

		for(int yy = 0; yy< poly->nverts;yy++)
		{
			if( poly->vlist[yy]->x == Xold  && poly->vlist[yy]->y == Yold && poly->vlist[yy]->z == Zold)
			{
						
				vtemplist[yy]->MorseVal =newMorseVal;
				vtemplist[yy]->MorseVal1 =newMorseVal1;
			}
		}
	}

	
	for(int ll = 0 ; ll<nverts; ll++)
	{
		poly->vlist[ll]->MorseVal =  vtemplist[ll]->MorseVal;
		poly->vlist[ll]->MorseVal1 =  vtemplist[ll]->MorseVal1;
	}
	delete []vtemplist;
}

void Polyhedron::MorseDesignUsingGaussSeidalMethod()
{
	double newMorseVal = 0.0;
	double newMorseVal1 = 0.0;

	for(int i = 0; i<nverts; i++)
	{
		Edge** TemperoryEdgeList = new Edge*[20]; 
		int SizeOfEEdge =0;
		Vertex* CurrVertex = vlist[i];
		float OldMorseVal = CurrVertex->MorseVal;
		float OldMorseVal1 = CurrVertex->MorseVal1;

		float Xold = CurrVertex->x, Yold = CurrVertex->y,Zold = CurrVertex->z;


		for(int kk = 0 ; kk <nedges; kk++)
		{
			if(CurrVertex->x ==  elist[kk]->verts[0]->x && 
				CurrVertex->y ==  elist[kk]->verts[0]->y &&
				CurrVertex->z ==  elist[kk]->verts[0]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
			if(CurrVertex->x ==  elist[kk]->verts[1]->x && 
				CurrVertex->y ==  elist[kk]->verts[1]->y &&
				CurrVertex->z ==  elist[kk]->verts[1]->z)
			{
				TemperoryEdgeList[SizeOfEEdge++] = elist[kk];
			}
		}

		float SumOf = 0.0f, SumOf1 = 0.0f;

		for(int hh= 0; hh<SizeOfEEdge; hh++)
		{
			Vertex* selectedVertex =  NULL;

			if(TemperoryEdgeList[hh]->verts[0]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[0]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[0]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[1];
			}

			if(TemperoryEdgeList[hh]->verts[1]->x == CurrVertex->x &&
				TemperoryEdgeList[hh]->verts[1]->y == CurrVertex->y &&
				TemperoryEdgeList[hh]->verts[1]->z == CurrVertex->z)
			{
				selectedVertex = TemperoryEdgeList[hh]->verts[0];
			}

			float diff = selectedVertex->MorseVal - CurrVertex->MorseVal;
			if(0 != diff)
			{
				SumOf = SumOf + diff;
			}

			float diff1 = selectedVertex->MorseVal1 - CurrVertex->MorseVal1;
			if(0 != diff1)
			{
				SumOf1 = SumOf1 + diff1;
			}
		}
		
		double average = SumOf/SizeOfEEdge;
		newMorseVal = OldMorseVal + dt * average;

		double average1 = SumOf1/SizeOfEEdge;
		newMorseVal1 = OldMorseVal1 + dt * average1;

		for(int yy = 0; yy< poly->nverts;yy++)
		{
			if( poly->vlist[yy]->x == Xold  &&
				poly->vlist[yy]->y == Yold &&
				poly->vlist[yy]->z == Zold
				)
			{
				vlist[yy]->MorseVal = newMorseVal;
				vlist[yy]->MorseVal1 = newMorseVal1;
			}
		}
	}
}




/******************************************************************************
Write out a polyhedron to a file.
******************************************************************************/

void Polyhedron::write_file(FILE *file)
{
	int i;
	PlyFile *ply;
	char **elist;
	int num_elem_types;

	/*** Write out the transformed PLY object ***/

	elist = get_element_list_ply (in_ply, &num_elem_types);
	ply = write_ply (file, num_elem_types, elist, in_ply->file_type);

	/* describe what properties go into the vertex elements */

	describe_element_ply (ply, "vertex", nverts);
	describe_property_ply (ply, &vert_props[0]);
	describe_property_ply (ply, &vert_props[1]);
	describe_property_ply (ply, &vert_props[2]);
	//  describe_other_properties_ply (ply, vert_other, offsetof(Vertex_io,other_props));

	describe_element_ply (ply, "face", ntris);
	describe_property_ply (ply, &face_props[0]);

	//  describe_other_properties_ply (ply, face_other,
	//                                offsetof(Face_io,other_props));

	//  describe_other_elements_ply (ply, in_ply->other_elems);

	copy_comments_ply (ply, in_ply);
	(ply, in_ply);
	char k[1024];
	sprintf(k, "modified by learnply");
	//  append_cokent_ply (ply, "modified by simvizply %f");
	append_comment_ply (ply, k);
	copy_obj_info_ply (ply, in_ply);

	header_complete_ply (ply);

	/* set up and write the vertex elements */
	put_element_setup_ply (ply, "vertex");
	for (i = 0; i < nverts; i++) {
		Vertex_io vert;

		/* copy info to the "vert" structure */
		vert.x = vlist[i]->x;
		vert.y = vlist[i]->y;
		vert.z = vlist[i]->z;
		vert.other_props = vlist[i]->other_props;

		put_element_ply (ply, (void *) &vert);
	}

	/* index all the vertices */
	for (i = 0; i < nverts; i++)
		vlist[i]->index = i;

	/* set up and write the face elements */
	put_element_setup_ply (ply, "face");

	Face_io face;
	face.verts = new int[3];

	for (i = 0; i < ntris; i++) {

		/* copy info to the "face" structure */
		face.nverts = 3;
		face.verts[0] = tlist[i]->verts[0]->index;
		face.verts[1] = tlist[i]->verts[1]->index;
		face.verts[2] = tlist[i]->verts[2]->index;
		face.other_props = tlist[i]->other_props;

		put_element_ply (ply, (void *) &face);
	}
	put_other_elements_ply (ply);

	close_ply (ply);
	free_ply (ply);
}

void Polyhedron::initialize(){
	icVector3 v1, v2;

	create_pointers();
	calc_edge_length();
	seed = -1;
}

void Polyhedron::finalize(){
	int i;

	for (i=0; i<ntris; i++){
		free(tlist[i]->other_props);
		free(tlist[i]);
	}
	for (i=0; i<nedges; i++) {
		free(elist[i]->tris);
		free(elist[i]);
	}
	for (i=0; i<nverts; i++) {
		free(vlist[i]->tris);
		free(vlist[i]->other_props);
		free(vlist[i]);
	}

	free(tlist);
	free(elist);
	free(vlist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
f1    - face that we're looking to share with
v1,v2 - two vertices of f1 that define edge

Exit:
return the matching face, or NULL if there is no such face
******************************************************************************/

Triangle *Polyhedron::find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
	int i,j;
	Triangle *f2;
	Triangle *adjacent = NULL;

	/* look through all faces of the first vertex */

	for (i = 0; i < v1->ntris; i++) {
		f2 = v1->tris[i];
		if (f2 == f1)
			continue;
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < f2->nverts; j++) {

			/* look for a match */
			if (f2->verts[j] == v2) {

#if 0
				/* watch out for triple edges */

				if (adjacent != NULL) {

					fprintf (stderr, "model has triple edges\n");

					fprintf (stderr, "face 1: ");
					for (k = 0; k < f1->nverts; k++)
						fprintf (stderr, "%d ", f1->iverts[k]);
					fprintf (stderr, "\nface 2: ");
					for (k = 0; k < f2->nverts; k++)
						fprintf (stderr, "%d ", f2->iverts[k]);
					fprintf (stderr, "\nface 3: ");
					for (k = 0; k < adjacent->nverts; k++)
						fprintf (stderr, "%d ", adjacent->iverts[k]);
					fprintf (stderr, "\n");

				}

				/* if we've got a match, remember this face */
				adjacent = f2;
#endif

#if 1
				/* if we've got a match, return this face */
				return (f2);
#endif

			}
		}
	}

	return (adjacent);
}


/******************************************************************************
Create an edge.

Entry:
v1,v2 - two vertices of f1 that define edge
******************************************************************************/


void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
	int i,j;
	Triangle *f;

	/* make sure there is enough room for a new edge */

	if (nedges >= max_edges) {

		max_edges += 100;
		Edge **list = new Edge *[max_edges];

		/* copy the old list to the new one */
		for (i = 0; i < nedges; i++)
			list[i] = elist[i];

		/* replace list */
		free (elist);
		elist = list;
	}

	/* create the edge */

	elist[nedges] = new Edge;
	Edge *e = elist[nedges];
	e->index = nedges;
	e->verts[0] = v1;
	e->verts[1] = v2;
	e->tag = 0;
	nedges++;

	/* count all triangles that will share the edge, and do this */
	/* by looking through all faces of the first vertex */

	e->ntris = 0;

	for (i = 0; i < v1->ntris; i++) {
		f = v1->tris[i];
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 3; j++) {
			/* look for a match */
			if (f->verts[j] == v2) {
				e->ntris++;
				break;
			}
		}
	}

	/* make room for the face pointers (at least two) */
	if (e->ntris < 2)
		e->tris = new Triangle *[2];
	else
		e->tris = new Triangle *[e->ntris];

	/* create pointers from edges to faces and vice-versa */

	e->ntris = 0; /* start this out at zero again for creating ptrs to tris */

	for (i = 0; i < v1->ntris; i++) {

		f = v1->tris[i];

		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 3; j++)
			if (f->verts[j] == v2) {

				e->tris[e->ntris] = f;
				e->ntris++;

				if (f->verts[(j+1)%3] == v1)
					f->edges[j] = e;
				else if (f->verts[(j+2)%3] == v1)
					f->edges[(j+2)%3] = e;
				else {
					fprintf (stderr, "Non-recoverable inconsistancy in create_edge()\n");
					exit (-1);
				}

				break;  /* we'll only find one instance of v2 */
			}

	}
}


/******************************************************************************
Create edges.
******************************************************************************/

void Polyhedron::create_edges()
{
	int i,j;
	Triangle *f;
	Vertex *v1,*v2;
	double count = 0;

	/* count up how many edges we may require */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < f->nverts; j++) {
			v1 = f->verts[j];
			v2 = f->verts[(j+1) % f->nverts];
			Triangle *result = find_common_edge (f, v1, v2);
			if (result)
				count += 0.5;
			else
				count += 1;
		}
	}

	/*
	printf ("counted %f edges\n", count);
	*/

	/* create space for edge list */

	max_edges = (int) (count + 10);  /* leave some room for expansion */
	elist = new Edge *[max_edges];
	nedges = 0;

	/* zero out all the pointers from faces to edges */

	for (i = 0; i < ntris; i++)
		for (j = 0; j < 3; j++)
			tlist[i]->edges[j] = NULL;

	/* create all the edges by examining all the triangles */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < 3; j++) {
			/* skip over edges that we've already created */
			if (f->edges[j])
				continue;
			v1 = f->verts[j];
			v2 = f->verts[(j+1) % f->nverts];
			create_edge (v1, v2);
		}
	}
}


/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/

void Polyhedron::vertex_to_tri_ptrs()
{
	int i,j;
	Triangle *f;
	Vertex *v;

	/* zero the count of number of pointers to faces */

	for (i = 0; i < nverts; i++)
		vlist[i]->max_tris = 0;

	/* first just count all the face pointers needed for each vertex */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < f->nverts; j++)
			f->verts[j]->max_tris++;
	}

	/* allocate memory for face pointers of vertices */

	for (i = 0; i < nverts; i++) {
		vlist[i]->tris = (Triangle **)
			malloc (sizeof (Triangle *) * vlist[i]->max_tris);
		vlist[i]->ntris = 0;
	}

	/* now actually create the face pointers */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < f->nverts; j++) {
			v = f->verts[j];
			v->tris[v->ntris] = f;
			v->ntris++;
		}
	}
}


/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
******************************************************************************/

Triangle *Polyhedron::other_triangle(Edge *edge, Triangle *tri)
{
	/* search for any other triangle */

	for (int i = 0; i < edge->ntris; i++)
		if (edge->tris[i] != tri)
			return (edge->tris[i]);

	/* there is no such other triangle if we get here */
	return (NULL);
}


/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
v - vertex whose face list is to be ordered
******************************************************************************/

void Polyhedron::order_vertex_to_tri_ptrs(Vertex *v)
{
	int i,j;
	Triangle *f;
	Triangle *fnext;
	int nf;
	int vindex;
	int boundary;
	int count;

	nf = v->ntris;
	f = v->tris[0];

	/* go backwards (clockwise) around faces that surround a vertex */
	/* to find out if we reach a boundary */

	boundary = 0;

	for (i = 1; i <= nf; i++) {

		/* find reference to v in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[j] == v) {
				vindex = j;
				break;
			}

			/* error check */
			if (vindex == -1) {
				fprintf (stderr, "can't find vertex #1\n");
				exit (-1);
			}

			/* corresponding face is the previous one around v */
			fnext = other_triangle (f->edges[vindex], f);

			/* see if we've reached a boundary, and if so then place the */
			/* current face in the first position of the vertice's face list */

			if (fnext == NULL) {
				/* find reference to f in v */
				for (j = 0; j < v->ntris; j++)
					if (v->tris[j] == f) {
						v->tris[j] = v->tris[0];
						v->tris[0] = f;
						break;
					}
					boundary = 1;
					break;
			}

			f = fnext;
	}

	/* now walk around the faces in the forward direction and place */
	/* them in order */

	f = v->tris[0];
	count = 0;

	for (i = 1; i < nf; i++) {

		/* find reference to vertex in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[(j+1) % f->nverts] == v) {
				vindex = j;
				break;
			}

			/* error check */
			if (vindex == -1) {
				fprintf (stderr, "can't find vertex #2\n");
				exit (-1);
			}

			/* corresponding face is next one around v */
			fnext = other_triangle (f->edges[vindex], f);

			/* break out of loop if we've reached a boundary */
			count = i;
			if (fnext == NULL) {
				break;
			}

			/* swap the next face into its proper place in the face list */
			for (j = 0; j < v->ntris; j++)
				if (v->tris[j] == fnext) {
					v->tris[j] = v->tris[i];
					v->tris[i] = fnext;
					break;
				}

				f = fnext;
	}
}


/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
f - face whose vertex list is to be searched
v - vertex to return reference to

Exit:
returns index in face's list, or -1 if vertex not found
******************************************************************************/

int Polyhedron::face_to_vertex_ref(Triangle *f, Vertex *v)
{
	int j;
	int vindex = -1;

	for (j = 0; j < f->nverts; j++)
		if (f->verts[j] == v) {
			vindex = j;
			break;
		}

		return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/

void Polyhedron::create_pointers()
{
	int i;

	/* index the vertices and triangles */

	for (i = 0; i < nverts; i++)
		vlist[i]->index = i;

	for (i = 0; i < ntris; i++) 
		tlist[i]->index = i;

	/* create pointers from vertices to triangles */
	vertex_to_tri_ptrs();

	/* make edges */
	create_edges();


	/* order the pointers from vertices to faces */
	for (i = 0; i < nverts; i++){
		//		if (i %1000 == 0)
		//			fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
		order_vertex_to_tri_ptrs(vlist[i]);

	}
	/* index the edges */

	for (i = 0; i < nedges; i++){
		//		if (i %1000 == 0)
		//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
		elist[i]->index = i;
	}

}

void Polyhedron::calc_bounding_box_withnormals()
{
	float AA=0.f;
	float BB=0.f;
	float CC=0.f;
	float DD=0.0;
	float EE=0.0;
	float FF=0.0;

	float xSum = 0.f;
	float ySum = 0.f;
	float zSum = 0.f;

	int VertexCount=0;

	for(VertexCount=0; VertexCount<nverts; VertexCount++)
	{
		xSum += vlist[VertexCount]->x;
		ySum += vlist[VertexCount]->y;
		zSum += vlist[VertexCount]->z;
	}

	icVector3 center1(xSum/nverts,ySum/nverts,zSum/nverts);
	for(VertexCount=0; VertexCount<nverts; VertexCount++)
	{
		/*float XX = vlist[VertexCount]->x - center1.x;
		float YY = vlist[VertexCount]->y - center1.y;
		float j = vlist[VertexCount]->z - center1.z;*/

		// these values are from the average normal to the the vertex now
		float XX = vlist[VertexCount]->normal.x;
		float YY = vlist[VertexCount]->normal.y;
		float j = vlist[VertexCount]->normal.z;

		AA = AA + XX*XX;
		BB = BB + XX*YY;
		CC = CC + XX*j;
		DD = DD + YY*YY;
		EE = EE + YY*j;
		FF = FF + j*j;
	}

	// Thus AA,BB, CC, DD, EE and FF are the final entries of the 3 X 3 Matrix
	// Next step is to obtain the Eigenvalues and Eigen vectors from this matrix
	float  inputmatrix[3][3];

	inputmatrix[0][0] = AA/nverts;
	inputmatrix[0][1] = BB/nverts;
	inputmatrix[0][2] = CC/nverts;

	inputmatrix[1][0] = BB/nverts;
	inputmatrix[1][1] = DD/nverts;
	inputmatrix[1][2] = EE/nverts;

	inputmatrix[2][0] = CC/nverts;
	inputmatrix[2][1] = EE/nverts;
	inputmatrix[2][2] = FF/nverts;

	float eigenValues[3];
	float eigenVectors[3][3];
	calculate_eigens(inputmatrix,3,eigenValues,eigenVectors);
	icVector3 median(center1.x,  center1.y,  center1.z);
	icVector3 norm1 (eigenVectors[0][0],eigenVectors[1][0],eigenVectors[2][0]);
	icVector3 norm2 (eigenVectors[0][1],eigenVectors[1][1],eigenVectors[2][1]);
	icVector3 norm3 (eigenVectors[0][2],eigenVectors[1][2],eigenVectors[2][2]);

	float normval1 = eigenValues[0];
	float normval2 = eigenValues[1];
	float normval3 = eigenValues[2];

	// Here we will calculate the extents to which the bounding box has to cover the volume
	// There will be 6 extents (Lengths). This will help us to draw the final bounding box
	// These extents will be the plane perpendicular to the Eigen vectors on either side

	double maxDotXProduct = -10000000;
	double minDotXProduct = 10000000;

	double maxDotYProduct = -10000000;
	double minDotYProduct = 10000000;

	double maxDotZProduct = -10000000;
	double minDotZProduct = 10000000;

	int vertexWithMaxDotX = 0;
	int vertexWithMinDotX = 0;

	int vertexWithMaxDotY = 0;
	int vertexWithMinDotY = 0;

	int vertexWithMaxDotZ = 0;
	int vertexWithMinDotZ = 0;

	for(VertexCount=0; VertexCount<nverts; VertexCount++)
	{
		icVector3 pVector((vlist[VertexCount]->x -median.x),(vlist[VertexCount]->y -median.y),(vlist[VertexCount]->z -median.z));
		double dotProductX = dot(pVector,norm1);
		double dotProductY = dot(pVector,norm2);
		double dotProductZ = dot(pVector,norm3);

		// X Extent
		if(dotProductX>maxDotXProduct)
		{
			maxDotXProduct = dotProductX;
			vertexWithMaxDotX = VertexCount;
		}

		if(dotProductX<minDotXProduct)
		{
			minDotXProduct = dotProductX;
			vertexWithMinDotX = VertexCount;
		}

		// Y Extent

		if(dotProductY>maxDotYProduct)
		{
			maxDotYProduct = dotProductY;
			vertexWithMaxDotY = VertexCount;
		}

		if(dotProductY<minDotYProduct)
		{
			minDotYProduct = dotProductY;
			vertexWithMinDotY = VertexCount;
		}

		// Z Extent
		if(dotProductZ>maxDotZProduct)
		{
			maxDotZProduct = dotProductZ;
			vertexWithMaxDotZ = VertexCount;
		}

		if(dotProductZ<minDotZProduct)
		{
			minDotZProduct = dotProductZ;
			vertexWithMinDotZ = VertexCount;
		}
	}


	double extentX1 = fabs(maxDotXProduct);
	double extentX2 = fabs(minDotXProduct);

	double extentY1 = fabs(maxDotYProduct);
	double extentY2 = fabs(minDotYProduct);

	double extentZ1 = fabs(maxDotZProduct);
	double extentZ2 = fabs(minDotZProduct);

	icVector3 bboxpoint1 = median + norm1*extentX1 + norm2*extentY1 + norm3*extentZ1;
	icVector3 bboxpoint2 = median + norm1*extentX1 + norm2*extentY1 - norm3*extentZ2;
	icVector3 bboxpoint3 = median + norm1*extentX1 - norm2*extentY2 - norm3*extentZ2;
	icVector3 bboxpoint4 = median + norm1*extentX1 - norm2*extentY2 + norm3*extentZ1;
	icVector3 bboxpoint5 = median - norm1*extentX2 + norm2*extentY1 + norm3*extentZ1;
	icVector3 bboxpoint6 = median - norm1*extentX2 + norm2*extentY1 - norm3*extentZ2;
	icVector3 bboxpoint7 = median - norm1*extentX2 - norm2*extentY2 - norm3*extentZ2;
	icVector3 bboxpoint8 = median - norm1*extentX2 - norm2*extentY2 + norm3*extentZ1;

	boundbox3[0][0] = bboxpoint1.x;
	boundbox3[1][0] = bboxpoint1.y;
	boundbox3[2][0] = bboxpoint1.z;

	boundbox3[0][1] = bboxpoint2.x;
	boundbox3[1][1] = bboxpoint2.y;
	boundbox3[2][1] = bboxpoint2.z;

	boundbox3[0][2] = bboxpoint3.x;
	boundbox3[1][2] = bboxpoint3.y;
	boundbox3[2][2] = bboxpoint3.z;

	boundbox3[0][3] = bboxpoint4.x;
	boundbox3[1][3] = bboxpoint4.y;
	boundbox3[2][3] = bboxpoint4.z;

	boundbox3[0][4] = bboxpoint5.x;
	boundbox3[1][4] = bboxpoint5.y;
	boundbox3[2][4] = bboxpoint5.z;

	boundbox3[0][5] = bboxpoint6.x;
	boundbox3[1][5] = bboxpoint6.y;
	boundbox3[2][5] = bboxpoint6.z;

	boundbox3[0][6] = bboxpoint7.x;
	boundbox3[1][6] = bboxpoint7.y;
	boundbox3[2][6] = bboxpoint7.z;

	boundbox3[0][7] = bboxpoint8.x;
	boundbox3[1][7] = bboxpoint8.y;
	boundbox3[2][7] = bboxpoint8.z;
}

void Polyhedron::calc_bounding_box()
{
	unsigned int i;
	float xsum=0,ysum=0,zsum=0;
	double a=0,b=0,c=0,d=0,e=0,f=0;
	icVector3 tempVect;


	for (i=0; i<nverts; i++) 
	{
		xsum=xsum+vlist[i]->x;
		ysum=ysum+vlist[i]->y;
		zsum=zsum+vlist[i]->z;
	}
	center0fGravity.x=xsum/nverts;
	center0fGravity.y=ysum/nverts;
	center0fGravity.z=zsum/nverts;
	for (i=0; i<nverts; i++) 
	{
		a+=(vlist[i]->x-center0fGravity.x)*(vlist[i]->x-center0fGravity.x);
		b+=(vlist[i]->x-center0fGravity.x)*(vlist[i]->y-center0fGravity.y);
		c+=(vlist[i]->x-center0fGravity.x)*(vlist[i]->z-center0fGravity.z);
		d+=(vlist[i]->y-center0fGravity.y)*(vlist[i]->y-center0fGravity.y);
		e+=(vlist[i]->y-center0fGravity.y)*(vlist[i]->z-center0fGravity.z);
		f+=(vlist[i]->z-center0fGravity.z)*(vlist[i]->z-center0fGravity.z);

	}
	float mat[3][3];
	mat[0][0]=a/nverts;
	mat[0][1]=b/nverts;
	mat[0][2]=c/nverts;
	mat[1][0]=b/nverts;
	mat[1][1]=d/nverts;
	mat[1][2]=e/nverts;
	mat[2][0]=c/nverts;
	mat[2][1]=e/nverts;
	mat[2][2]=f/nverts;
	//	float mat[3][3]={a/nverts,b/nverts,c/nverts,b/nverts,d/nverts,e/nverts,c/nverts,e/nverts,f/nverts};
	//	printf("\nfinalmatrix\n%f,%f,%f\n%f,%f,%f,\n%f,%f,%f\n",a,b,c,b,d,e,c,e,f);
	calculate_eigens(mat,3,eigval,eigenvector);

	icVector3 EigVect1(eigenvector[0][0],eigenvector[1][0],eigenvector[2][0]);
	icVector3 EigVect2(eigenvector[0][1],eigenvector[1][1],eigenvector[2][1]);
	icVector3 EigVect3(eigenvector[0][2],eigenvector[1][2],eigenvector[2][2]);

	double dotprodWRTx,dotprodWRTy,dotprodWRTz;
	double maxX=-1000000,maxY=-1000000,maxZ=-1000000,minX=100000,minY=100000,minZ=100000;

	icVector3 vecAlignedX,vecAlignedy,vecAlignedz;
	for (i=0; i<nverts; i++) 
	{
		//tempVect=((vlist[i]->x-center0fGravity.x),(vlist[i]->y-center0fGravity.y),(vlist[i]->z-center0fGravity.z));

		tempVect.x = vlist[i]->x-center0fGravity.x;
		tempVect.y = vlist[i]->y-center0fGravity.y;
		tempVect.z = vlist[i]->z-center0fGravity.z;

		//		dotprodWRTx= eigenvector[0][0]*tempVect.x+eigenvector[1][0]*tempVect.y+eigenvector[2][0]*tempVect.z;
		dotprodWRTx = dot(tempVect,EigVect1);
		dotprodWRTy = dot(tempVect,EigVect2);
		dotprodWRTz = dot(tempVect,EigVect3);


		if(dotprodWRTx>maxX)
		{
			maxX=dotprodWRTx;
			countXmax=i;
		}
		if(dotprodWRTx<minX)
		{
			minX=dotprodWRTx;
			countXmin=i;
		}

		//		dotprodWRTy= eigenvector[0][1]*tempVect.x+eigenvector[1][1]*tempVect.y+eigenvector[2][1]*tempVect.z;

		if(dotprodWRTy>maxY)
		{
			maxY=dotprodWRTy;
			countYmax=i;
		}
		if(dotprodWRTy<minY)
		{
			minY=dotprodWRTy;
			countYmin=i;
		}

		//		dotprodWRTz= eigenvector[0][2]*tempVect.x+eigenvector[1][2]*tempVect.y+eigenvector[2][2]*tempVect.z;

		if(dotprodWRTz>maxZ)
		{
			maxZ=dotprodWRTz;
			countZmax=i;
		}
		if(dotprodWRTz<minZ)
		{
			minZ=dotprodWRTz;
			countZmin=i;
		}


	}
	double extentX1=fabs(maxX);
	double extentX2=fabs(minX);
	double extentY1=fabs(maxY);
	double extentY2=fabs(minY);
	double extentZ1=fabs(maxZ);
	double extentZ2=fabs(minZ);

	icVector3 VectorwrtX(eigenvector[0][0],eigenvector[1][0],eigenvector[2][0]);
	icVector3 VectorwrtY(eigenvector[0][1],eigenvector[1][1],eigenvector[2][1]);
	icVector3 VectorwrtZ(eigenvector[0][2],eigenvector[1][2],eigenvector[2][2]);


	icVector3 bb1= center0fGravity+ VectorwrtX*extentX1 +VectorwrtY*extentY1 +VectorwrtZ*extentZ1;
	icVector3 bb2= center0fGravity+  VectorwrtX*extentX1 +VectorwrtY*extentY1 -VectorwrtZ*extentZ1;
	icVector3 bb3= center0fGravity+  VectorwrtX*extentX1 -VectorwrtY*extentY1 -VectorwrtZ*extentZ1;
	icVector3 bb4= center0fGravity+  VectorwrtX*extentX1 -VectorwrtY*extentY1 +VectorwrtZ*extentZ1;
	icVector3 bb5= center0fGravity-  VectorwrtX*extentX1 +VectorwrtY*extentY1 +VectorwrtZ*extentZ1;
	icVector3 bb6= center0fGravity-  VectorwrtX*extentX1 +VectorwrtY*extentY1 -VectorwrtZ*extentZ1;
	icVector3 bb7= center0fGravity-  VectorwrtX*extentX1 -VectorwrtY*extentY1 -VectorwrtZ*extentZ1;
	icVector3 bb8= center0fGravity-  VectorwrtX*extentX1 -VectorwrtY*extentY1 +VectorwrtZ*extentZ1;

	boundbox1[0][0]=bb1.x;
	boundbox1[1][0]=bb1.y;
	boundbox1[2][0]=bb1.z;

	boundbox1[0][1]=bb2.x;
	boundbox1[1][1]=bb2.y;
	boundbox1[2][1]=bb2.z;

	boundbox1[0][2]=bb3.x;
	boundbox1[1][2]=bb3.y;
	boundbox1[2][2]=bb3.z;


	boundbox1[0][3]=bb4.x;
	boundbox1[1][3]=bb4.y;
	boundbox1[2][3]=bb4.z;


	boundbox1[0][4]=bb5.x;
	boundbox1[1][4]=bb5.y;
	boundbox1[2][4]=bb5.z;

	boundbox1[0][5]=bb6.x;
	boundbox1[1][5]=bb6.y;
	boundbox1[2][5]=bb6.z;

	boundbox1[0][6]=bb7.x;
	boundbox1[1][6]=bb7.y;
	boundbox1[2][6]=bb7.z;

	boundbox1[0][7]=bb8.x;
	boundbox1[1][7]=bb8.y;
	boundbox1[2][7]=bb8.z;




	//	printf("\Eigenvector:\nfirst:%f,%f,%f\nsecond:%f,%f,%f,\nThird:%f,%f,%f\n",eigenvector[0][0],eigenvector[0][1],eigenvector[0][2],eigenvector[1][0],eigenvector[1][1],eigenvector[1][2],eigenvector[2][0],eigenvector[2][1],eigenvector[2][2]);
}


void Polyhedron::calc_bounding_sphere()
{
	unsigned int i;
	icVector3 min, max;

	for (i=0; i<nverts; i++) {
		if (i==0)  {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		}
		else {
			if (vlist[i]->x < min.entry[0])
				min.entry[0] = vlist[i]->x;
			if (vlist[i]->x > max.entry[0])
				max.entry[0] = vlist[i]->x;
			if (vlist[i]->y < min.entry[1])
				min.entry[1] = vlist[i]->y;
			if (vlist[i]->y > max.entry[1])
				max.entry[1] = vlist[i]->y;
			if (vlist[i]->z < min.entry[2])
				min.entry[2] = vlist[i]->z;
			if (vlist[i]->z > max.entry[2])
				max.entry[2] = vlist[i]->z;
		}
	}
	center = (min + max) * 0.5;
	radius = length(center - min);

	icVector3 normalVector1(1,0,0);
	icVector3 normalVector2(0,1,0);
	icVector3 normalVector3(0,0,1);
	icVector3 bb1 ;
	icVector3 bb2 ;
	icVector3 bb3 ;
	icVector3 bb4 ;
	icVector3 bb5 ;
	icVector3 bb6 ;
	icVector3 bb7 ;
	icVector3 bb8 ;

	bb1 = center + normalVector1*fabs(max.entry[0]) + normalVector2*fabs(max.entry[1]) + normalVector3*fabs(max.entry[2]);
	bb2 = center + normalVector1*fabs(max.entry[0]) + normalVector2*fabs(max.entry[1]) - normalVector3*fabs(min.entry[2]);
	bb3 = center + normalVector1*fabs(max.entry[0]) - normalVector2*fabs(min.entry[1]) - normalVector3*fabs(min.entry[2]);
	bb4 = center + normalVector1*fabs(max.entry[0]) - normalVector2*fabs(min.entry[1]) + normalVector3*fabs(max.entry[2]);
	bb5 = center - normalVector1*fabs(min.entry[0]) + normalVector2*fabs(max.entry[1]) + normalVector3*fabs(max.entry[2]);
	bb6 = center - normalVector1*fabs(min.entry[0]) + normalVector2*fabs(max.entry[1]) - normalVector3*fabs(min.entry[2]);
	bb7 = center - normalVector1*fabs(min.entry[0]) - normalVector2*fabs(min.entry[1]) - normalVector3*fabs(min.entry[2]);
	bb8 = center - normalVector1*fabs(min.entry[0]) - normalVector2*fabs(min.entry[1]) + normalVector3*fabs(max.entry[2]);



	boundbox2[0][0] = bb1.x;
	boundbox2[1][0] = bb1.y;
	boundbox2[2][0] = bb1.z;

	boundbox2[0][1] = bb2.x;
	boundbox2[1][1] = bb2.y;
	boundbox2[2][1] = bb2.z;

	boundbox2[0][2] = bb3.x;
	boundbox2[1][2] = bb3.y;
	boundbox2[2][2] = bb3.z;

	boundbox2[0][3] = bb4.x;
	boundbox2[1][3] = bb4.y;
	boundbox2[2][3] = bb4.z;

	boundbox2[0][4] = bb5.x;
	boundbox2[1][4] = bb5.y;
	boundbox2[2][4] = bb5.z;

	boundbox2[0][5] = bb6.x;
	boundbox2[1][5] = bb6.y;
	boundbox2[2][5] = bb6.z;

	boundbox2[0][6] = bb7.x;
	boundbox2[1][6] = bb7.y;
	boundbox2[2][6] = bb7.z;

	boundbox2[0][7] = bb8.x;
	boundbox2[1][7] = bb8.y;
	boundbox2[2][7] = bb8.z;
}

void Polyhedron::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i=0; i<nedges; i++) {
		v1.set(elist[i]->verts[0]->x, elist[i]->verts[0]->y, elist[i]->verts[0]->z);
		v2.set(elist[i]->verts[1]->x, elist[i]->verts[1]->y, elist[i]->verts[1]->z);
		elist[i]->length = length(v1-v2);
	}
}

void Polyhedron::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2;
	Triangle *temp_t;
	double length[3];

	area = 0.0;
	for (i=0; i<ntris; i++){
		for (j=0; j<3; j++)
			length[j] = tlist[i]->edges[j]->length;
		double temp_s = (length[0] + length[1] + length[2])/2.0;
		tlist[i]->area = sqrt(temp_s*(temp_s-length[0])*(temp_s-length[1])*(temp_s-length[2]));

		area += tlist[i]->area;
		temp_t = tlist[i];
		v1.set(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		v2.set(vlist[tlist[i]->verts[1]->index]->x, vlist[tlist[i]->verts[1]->index]->y, vlist[tlist[i]->verts[1]->index]->z);
		v0.set(vlist[tlist[i]->verts[2]->index]->x, vlist[tlist[i]->verts[2]->index]->y, vlist[tlist[i]->verts[2]->index]->z);
		tlist[i]->normal = cross(v0-v1, v2-v1);
		normalize(tlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i=0; i<ntris; i++){
		icVector3 cent(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		signedvolume += dot(test-cent, tlist[i]->normal)*tlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume<0) 
		orientation = 0;
	else {
		orientation = 1;
		for (i=0; i<ntris; i++)
			tlist[i]->normal *= -1.0;
	}
}

void sort(unsigned int *A, unsigned int *B, unsigned int *C, unsigned int sid, unsigned int eid){
	unsigned int i;
	unsigned int *tempA, *tempB, *tempC;
	unsigned int current1, current2, current0;

	if (sid>=eid)
		return;
	sort(A, B, C, sid, (sid+eid)/2);
	sort(A, B, C, (sid+eid)/2+1, eid);
	tempA = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	tempB = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	tempC = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	for (i=0; i<eid-sid+1; i++){
		tempA[i] = A[i+sid];
		tempB[i] = B[i+sid];
		tempC[i] = C[i+sid];
	}
	current1 = sid;
	current2 = (sid+eid)/2+1;
	current0 = sid;
	while ((current1<=(sid+eid)/2) && (current2<=eid)){
		if (tempA[current1-sid] < tempA[current2-sid]) {
			A[current0] = tempA[current1-sid];
			B[current0] = tempB[current1-sid];
			C[current0] = tempC[current1-sid];
			current1++;		
		}
		else if (tempA[current1-sid] > tempA[current2-sid]){
			A[current0] = tempA[current2-sid];
			B[current0] = tempB[current2-sid];
			C[current0] = tempC[current2-sid];
			current2++;		
		}
		else {
			if (tempB[current1-sid] < tempB[current2-sid]) {
				A[current0] = tempA[current1-sid];
				B[current0] = tempB[current1-sid];
				C[current0] = tempC[current1-sid];
				current1++;		
			} else {
				A[current0] = tempA[current2-sid];
				B[current0] = tempB[current2-sid];
				C[current0] = tempC[current2-sid];
				current2++;		
			}
		}
		current0++;
	}
	if (current1<=(sid+eid)/2){
		for (i=current1; i<=(sid+eid)/2; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}
	if (current2<=eid){
		for (i=current2; i<=eid; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}

	free(tempA);
	free(tempB);
	free(tempC);
}

void init(void) {
	/* select clearing color */ 

	glClearColor (0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel (GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	// may need it
	glPixelStorei(GL_PACK_ALIGNMENT,1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0) 
		glFrontFace(GL_CW);
	else 
		glFrontFace(GL_CCW);
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '0':
		display_mode = 0;
		display();
		break;

	case '1':
		display_mode = 1;
		display();
		break;

	case '2':
		display_mode = 2;
		display();
		break;

	case '3':
		display_mode = 3;
		display();
		break;

	case '4':
		display_mode = 4;
		display();
		break;

	case '5':
		display_mode = 5;
		display();
		break;

	case '6':
		display_mode = 6;
		display();
		break;

	case '7':
		display_mode = 7;
		display();
		break;

	case '8':
		display_mode = 8;
		display();
		break;

	case '9':
		display_mode = 9;
		display();
		break;

	case 'a':
		display_mode = 11;
		display();
		break;

	case 'b':
		display_mode = 12;
		display();
		break;

	case 'c':
		display_mode = 13;
		display();
		break;

	case 'x':
		switch(ACSIZE){
		case 1:
			ACSIZE = 16;
			break;

		case 16:
			ACSIZE = 1;
			break;

		default:
			ACSIZE = 1;
			break;
		}
		fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
		display();
		break;
	case '|':
		this_file = fopen("rotmat.txt", "w");
		for (i=0; i<4; i++) 
			fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
		fclose(this_file);
		break;

	case '^':
		this_file = fopen("rotmat.txt", "r");
		for (i=0; i<4; i++) 
			fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
		fclose(this_file);
		display();
		break;

	}
}

Polyhedron::Polyhedron()
{
	nverts = nedges = ntris = 0;
	max_verts = max_tris = 50;

	vlist = new Vertex *[max_verts];
	tlist = new Triangle *[max_tris];		
}


void multmatrix(const Matrix m)
{ 
	int i,j, index = 0;

	GLfloat mat[16];

	for ( i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			mat[index++] = m[i][j];

	glMultMatrixf (mat);
}

void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(0.0, 0.0, -3.0);
	multmatrix( rotmat );

	glScalef(1.0/poly->radius, 1.0/poly->radius, 1.0/poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

void motion(int x, int y) {
	float r[4];
	float xsize, ysize, s, t;

	switch(mouse_mode){
	case -1:

		xsize = (float) win_width;
		ysize = (float) win_height;

		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		mat_to_quat( rotmat, rvec );
		trackball( r, s_old, t_old, s, t );
		add_quats( r, rvec, rvec );
		quat_to_mat( rvec, rotmat );

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth=1.0e+20, current_depth;
	int seed_id=-1; 
	unsigned char need_to_update;

	printf("hits = %d\n", hits);
	ptr = (GLuint *) buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	printf("triangle id = %d\n", seed_id);
	return seed_id;
}
void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		switch(mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float) win_width;
				float ysize = (float) win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_mode = -1;  // down
				mouse_button = button;
				last_x = x;
				last_y = y;
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	} else if (button == GLUT_MIDDLE_BUTTON) {
		if (state == GLUT_DOWN) {  // build up the selection feedback mode

			GLuint selectBuf[win_width];
			GLint hits;
			GLint viewport[4];

			glGetIntegerv(GL_VIEWPORT, viewport);

			glSelectBuffer(win_width, selectBuf);
			(void) glRenderMode(GL_SELECT);

			glInitNames();
			glPushName(0);

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();
			/*  create 5x5 pixel picking region near cursor location */
			gluPickMatrix((GLdouble) x, (GLdouble) (viewport[3] - y),
				1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix ();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
			glPopMatrix();
			glFlush();

			hits = glRenderMode(GL_RENDER);
			poly->seed = processHits(hits, selectBuf);

			// Update the morse value of the selected vertex
			Triangle* ClickedTriangle =  poly->tlist[poly->seed];
			Vertex* SelectedVertex = NULL;
			if(ClickedTriangle)
			{
				SelectedVertex = ClickedTriangle->verts[0];
			}
			if(SelectedVertex)
			{
				if(0 == isSecondTime)
				{
					SelectedVertex->MorseVal = 1.0f;
					isSecondTime = 1;
				}
				else
				{
					SelectedVertex->MorseVal1 = 1.0f;
				}
			}
			// End Update
			display();
		}
	}
}

void CreateTexture(){ 

	/* load an image file directly as a new OpenGL texture */ 
	texName = SOIL_load_OGL_texture 
		( 
		"../sample_textures/images1.jpg", 
		SOIL_LOAD_AUTO, 
		SOIL_CREATE_NEW_ID, 
		SOIL_FLAG_INVERT_Y 
		); 

	if (texName == 0) 
		printf("no picture\n"); 

	// Typical Texture Generation Using Data From The Bitmap 
	glBindTexture(GL_TEXTURE_2D, texName); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 

	//glGenTextures(1, &textureID); 

	// "Bind" the newly created texture : all future texture functions will modify this texture 
	glBindTexture(GL_TEXTURE_2D, texName); 

	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); 
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); 
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); 
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 

	// Give the image to OpenGL 
	//float pixels[] = { 
	// 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 
	// 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f 
	//}; 
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 2, 2, 0, GL_RGB, GL_FLOAT, pixels); 
	// 
	////printf("create the texture\n"); 
	//int width=64; 
	//int height=64; 
	// 
	//unsigned char* image = SOIL_load_image("../sample_textures/107.jpg", &width, &height, 0, SOIL_LOAD_RGB); 
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image); 
	// 
	////raw_texture_load("../sample_textures/107.bmp"); 
	////glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, checkImage); 
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); 

}
void FindSaddlesAndLocalMaxMinWithTexture(Polyhedron *this_poly)
{
	CreateTexture();

	for(int i = 0; i<this_poly->ntris; i++)
	{
		// Draw Texture here
		glEnable(GL_TEXTURE_2D);
		glBegin(GL_POLYGON);
		glTexCoord2d(this_poly->tlist[i]->verts[0]->MorseVal,this_poly->tlist[i]->verts[0]->MorseVal1);
		glVertex3d(this_poly->tlist[i]->verts[0]->x,this_poly->tlist[i]->verts[0]->y,this_poly->tlist[i]->verts[0]->z);

		glTexCoord2d(this_poly->tlist[i]->verts[1]->MorseVal,this_poly->tlist[i]->verts[1]->MorseVal1);
		glVertex3d(this_poly->tlist[i]->verts[1]->x,this_poly->tlist[i]->verts[1]->y,this_poly->tlist[i]->verts[1]->z);

		glTexCoord2d(this_poly->tlist[i]->verts[2]->MorseVal,this_poly->tlist[i]->verts[2]->MorseVal1);
		glVertex3d(this_poly->tlist[i]->verts[2]->x,this_poly->tlist[i]->verts[2]->y,this_poly->tlist[i]->verts[2]->z);
		glEnd();
		glDisable( GL_TEXTURE_2D );
	}
	// End Display Texture
	int nMax     = 0;
	int nMin     = 0;
	int nSaddle = 0;

	//CreateTexture();

	// End Display Texture
	for(int i = 0; i<this_poly->nverts; i++)
	{
		Vertex* currentVertexToUpdate = this_poly->vlist[i];
		int VertexSize = 0;
		Vertex** TempVertexlist = new Vertex*[20];

		// Here , we find the adjacent Vertices
		for(int kk = 0 ; kk <currentVertexToUpdate->ntris; kk++)
		{
			int cCount = 0;
			Triangle* ccTriangle = currentVertexToUpdate->tris[kk];

			for(int jj =0; jj<this_poly->nCorners;jj++)
			{
				if(this_poly->clist[jj]->tri->verts[0]->x == ccTriangle->verts[0]->x &&
					this_poly->clist[jj]->tri->verts[0]->y == ccTriangle->verts[0]->y &&
					this_poly->clist[jj]->tri->verts[0]->z == ccTriangle->verts[0]->z &&
					this_poly->clist[jj]->tri->verts[1]->x == ccTriangle->verts[1]->x &&
					this_poly->clist[jj]->tri->verts[1]->y == ccTriangle->verts[1]->y &&
					this_poly->clist[jj]->tri->verts[1]->z == ccTriangle->verts[1]->z &&
					this_poly->clist[jj]->tri->verts[2]->x == ccTriangle->verts[2]->x &&
					this_poly->clist[jj]->tri->verts[2]->y == ccTriangle->verts[2]->y &&
					this_poly->clist[jj]->tri->verts[2]->z == ccTriangle->verts[2]->z)
				{
					if(this_poly->clist[jj]->vert->x == currentVertexToUpdate->x &&
						this_poly->clist[jj]->vert->y == currentVertexToUpdate->y &&
						this_poly->clist[jj]->vert->z == currentVertexToUpdate->z )
					{
						cCount = jj;
						break;
					}
				}
			}
			TempVertexlist[VertexSize++] = this_poly->clist[cCount]->next->vert;
		}

		int IsMaxMinOrNone = 0; // 0: None, 1: Min, 2: Max
		int SaddleCount = 0;

		int* ListOfMaxMinNone = new int[VertexSize];
		for(int ff=0; ff<VertexSize; ff++)
		{
			ListOfMaxMinNone[ff] = 0;
		}


		for(int hh= 0; hh<VertexSize; hh++)
		{
			Vertex* Correctvertex =  TempVertexlist[hh];

			// Here we have both the current vertex and the correct vertex from the neighbour
			if(currentVertexToUpdate->MorseVal > Correctvertex->MorseVal)
			{
				ListOfMaxMinNone[hh] = 2;
			}
			else if(currentVertexToUpdate->MorseVal < Correctvertex->MorseVal)
			{
				ListOfMaxMinNone[hh] = 1;
			}
			else
			{
				ListOfMaxMinNone[hh] = 0;
			}
		}
		int IsNone = 1;
		int IsMin  = 1;
		int IsMax  = 1;		

		for(int ff=0; ff<VertexSize; ff++)
		{
			if(0 != ListOfMaxMinNone[ff])
			{
				IsNone = 0;
			}	
			if(1 != ListOfMaxMinNone[ff])
			{
				IsMin = 0;
			}
			if(2 != ListOfMaxMinNone[ff])
			{
				IsMax = 0;
			}
		}	


		glPointSize(10);
		glBegin(GL_POINTS);
		if(1 == IsMin)
		{
			nMin = nMin+1;
			glColor3f(0.0, 0.0, 1.0);

		}
		if(1 == IsMax)
		{
			nMax = nMax+1;
			glColor3f(1.0, 0.0, 0.0);

		}
		glEnd();

		// Here we calculate the saddle count 
		// Saddle count 0 and 1 are neglected
		for(int ff=0; ff<VertexSize; ff++)
		{
			if(2 == ListOfMaxMinNone[ff])
			{
				ListOfMaxMinNone[ff] = 0;
			}			
		}	

		int previousSign = ListOfMaxMinNone[VertexSize-1];

		for(int ff=0; ff<VertexSize; ff++)
		{
			int currentCount = VertexSize- 1 -ff;
			int currentSign = ListOfMaxMinNone[ff];
			if((1 == previousSign) && (0 == currentSign))
			{
				SaddleCount++;
				previousSign = currentSign;
			}
			else if((0 == previousSign) && (1 == currentSign))
			{
				previousSign = currentSign;
			}
			// Need to check the first point in case of the last point in the sequence
			// because it is a circular neighbour
			if(ff == VertexSize-1)
			{
				if(1 == currentSign)
				{
					if(0 == ListOfMaxMinNone[VertexSize-1])
					{
						SaddleCount++;
					}
				}
			}
		}

		if(SaddleCount>=2)
		{
			nSaddle = nSaddle + SaddleCount-1;
		}

		glPointSize(5);
		glBegin(GL_POINTS);
		if(2 == SaddleCount)
		{
			glColor3f(1.0, 1.0, 0.0);
			//glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}
		if(3 == SaddleCount)
		{
			glColor3f(1.0, 0.0, 1.0);
			//glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}
		if(4 == SaddleCount)
		{
			glColor3f(0.0, 1.0, 1.0);
			//glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}
		if(5 == SaddleCount)
		{
			glColor3f(1.0, 1.0, 0.5);
			//glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}

		glEnd();
	}
	printf("\n. No of Maxima=%d, Number of Minima=%d,Number of Saddle=%d",nMax,nMin,nSaddle);
}

void FindSaddlesAndLocalMaxMin(Polyhedron *this_poly)
{
	int nMax     = 0;
	int nMin     = 0;
	int nSaddle = 0;


	for(int i = 0; i<this_poly->nverts; i++)
	{
		Vertex* currentVertexToUpdate = this_poly->vlist[i];
		int VertexSize = 0;
		Vertex** TempVertexlist = new Vertex*[20];

		// Here , we find the adjacent Vertices
		for(int kk = 0 ; kk <currentVertexToUpdate->ntris; kk++)
		{
			int cCount = 0;
			Triangle* ccTriangle = currentVertexToUpdate->tris[kk];

			for(int jj =0; jj<this_poly->nCorners;jj++)
			{
				if(this_poly->clist[jj]->tri->verts[0]->x == ccTriangle->verts[0]->x &&
					this_poly->clist[jj]->tri->verts[0]->y == ccTriangle->verts[0]->y &&
					this_poly->clist[jj]->tri->verts[0]->z == ccTriangle->verts[0]->z &&
					this_poly->clist[jj]->tri->verts[1]->x == ccTriangle->verts[1]->x &&
					this_poly->clist[jj]->tri->verts[1]->y == ccTriangle->verts[1]->y &&
					this_poly->clist[jj]->tri->verts[1]->z == ccTriangle->verts[1]->z &&
					this_poly->clist[jj]->tri->verts[2]->x == ccTriangle->verts[2]->x &&
					this_poly->clist[jj]->tri->verts[2]->y == ccTriangle->verts[2]->y &&
					this_poly->clist[jj]->tri->verts[2]->z == ccTriangle->verts[2]->z)
				{
					if(this_poly->clist[jj]->vert->x == currentVertexToUpdate->x &&
						this_poly->clist[jj]->vert->y == currentVertexToUpdate->y &&
						this_poly->clist[jj]->vert->z == currentVertexToUpdate->z )
					{
						cCount = jj;
						break;
					}
				}
			}
			TempVertexlist[VertexSize++] = this_poly->clist[cCount]->next->vert;
		}

		int IsMaxMinOrNone = 0; // 0: None, 1: Min, 2: Max
		int SaddleCount = 0;

		int* ListOfMaxMinNone = new int[VertexSize];
		for(int ff=0; ff<VertexSize; ff++)
		{
			ListOfMaxMinNone[ff] = 0;
		}


		for(int hh= 0; hh<VertexSize; hh++)
		{
			Vertex* Correctvertex =  TempVertexlist[hh];

			// Here we have both the current vertex and the correct vertex from the neighbour
			double DiffVal = fabs(currentVertexToUpdate->MorseVal - Correctvertex->MorseVal);

			if(currentVertexToUpdate->MorseVal > Correctvertex->MorseVal && (DiffVal > 0.0000000000001))
			{
				ListOfMaxMinNone[hh] = 2;
			}
			else if(currentVertexToUpdate->MorseVal < Correctvertex->MorseVal && (DiffVal > 0.0000000000001))
			{
				ListOfMaxMinNone[hh] = 1;
			}
			else
			{
				ListOfMaxMinNone[hh] = 0;
			}
		}
		int IsNone = 1;
		int IsMin  = 1;
		int IsMax  = 1;		

		for(int ff=0; ff<VertexSize; ff++)
		{
			if(0 != ListOfMaxMinNone[ff])
			{
				IsNone = 0;
			}	
			if(1 != ListOfMaxMinNone[ff])
			{
				IsMin = 0;
			}
			if(2 != ListOfMaxMinNone[ff])
			{
				IsMax = 0;
			}
		}	

		glPointSize(10);
		glBegin(GL_POINTS);
		if(1 == IsMin)
		{
			nMin = nMin+1;
			glColor3f(0.0, 0.0, 1.0);
			glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}
		if(1 == IsMax)
		{
			nMax = nMax+1;
			glColor3f(1.0, 0.0, 0.0);
			glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}
		glEnd();

		// Here we calculate the saddle count 
		// Saddle count 0 and 1 are neglected
		for(int ff=0; ff<VertexSize; ff++)
		{
			if(2 == ListOfMaxMinNone[ff])
			{
				ListOfMaxMinNone[ff] = 0;
			}			
		}	

		int previousSign = ListOfMaxMinNone[VertexSize-1];

		for(int ff=0; ff<VertexSize; ff++)
		{
			int currentCount = VertexSize- 1 -ff;
			int currentSign = ListOfMaxMinNone[ff];
			if((1 == previousSign) && (0 == currentSign))
			{
				SaddleCount++;
				previousSign = currentSign;
			}
			else if((0 == previousSign) && (1 == currentSign))
			{
				previousSign = currentSign;
			}
			// Need to check the first point in case of the last point in the sequence
			// because it is a circular neighbour
			if(ff == VertexSize-1)
			{
				if(1 == currentSign)
				{
					if(0 == ListOfMaxMinNone[VertexSize-1])
					{
						SaddleCount++;
					}
				}
			}
		}

		if(SaddleCount>=2)
		{
			nSaddle = nSaddle + SaddleCount-1;
		}

		glPointSize(5);
		glBegin(GL_POINTS);
		if(2 == SaddleCount)
		{
			glColor3f(1.0, 1.0, 0.0);
			glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}
		if(3 == SaddleCount)
		{
			glColor3f(1.0, 0.0, 1.0);
			glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}
		if(4 == SaddleCount)
		{
			glColor3f(0.0, 1.0, 1.0);
			glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}
		if(5 == SaddleCount)
		{
			glColor3f(1.0, 1.0, 0.5);
			glVertex3d(currentVertexToUpdate->x,currentVertexToUpdate->y,currentVertexToUpdate->z);
		}

		glEnd();
	}
	printf("\n No of Maxima=%d, Number of Minima=%d,Number of Saddle=%d",nMax,nMin,nSaddle);
	printf("\n NMax + Nmin - Nsaddle =%d",nMax+nMin-nSaddle);
}

void display_object()
{
	unsigned int i, j;
	Polyhedron *the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (i=0; i<poly->ntris; i++) {
		Triangle *temp_t=poly->tlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		for (j=0; j<3; j++) {
			Vertex *temp_v = temp_t->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];
	int Total_def_Count = 0;
	double Total_Angle_def = 0.0;
	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	if(selection==0)
	{
		L = 0.003f;
		for (i=0; i<this_poly->ntris; i++) {
			if (mode == GL_SELECT)
				glLoadName(i+1);

			Triangle *temp_t=this_poly->tlist[i];
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			} else {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);


			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i==this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
					glColor3f(1.0, 1.0, 0.0);

				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	if(selection==1)
	{
		switch (display_mode) {
		case 0: //Default
			L = 0.003f;
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				if (i == this_poly->seed) {
					mat_diffuse[0] = 0.0;
					mat_diffuse[1] = 0.0;
					mat_diffuse[2] = 1.0;
					mat_diffuse[3] = 1.0;
				} else {
					mat_diffuse[0] = 1.0;
					mat_diffuse[1] = 1.0;
					mat_diffuse[2] = 0.0;
					mat_diffuse[3] = 1.0;
				}

				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);


				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					if (i==this_poly->seed)
						glColor3f(0.0, 0.0, 1.0);
					else
						glColor3f(1.0, 1.0, 0.0);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;


		case 1: // Mesh without subdivision

			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glColor3f(1.0, 1.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

			}
			break;

		case 2: // Mesh with Loop Subdivision(Regular Subdivision)

			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<polyRegular->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=polyRegular->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					if (i==this_poly->seed)
						glColor3f(0.0, 0.0, 1.0);
					else
						glColor3f(1.0, 1.0, 0.0);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);


				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;

		case 3: //  Irregular Subdivision

			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<polyIrregular->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=polyIrregular->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					if (i==this_poly->seed)
						glColor3f(0.0, 0.0, 1.0);
					else
						glColor3f(1.0, 1.0, 0.0);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;

		case 4: //3D ChekerBoard Color scheme
			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
			glDisable(GL_LIGHT1);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			if(4 == previousColor)
			{
				L = L+0.03;
			}
			else
			{
				previousColor = 4;
			}

			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];

					int R = floor(temp_v->x/L);
					if(R%2==0)
					{
						R=1;
					}
					else
					{
						R=0;
					}
					int G =floor(temp_v->y/L);
					if(G%2==0)
					{
						G=1;
					}
					else
					{
						G=0;
					}
					int B =floor(temp_v->z/L);
					if(B%2==0)
					{
						B=1;
					}
					else
					{
						B=0;
					}
					glColor3f(R,G,B);

					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

			}
			break;


		case 5: //  3D ChekerBoard Color scheme for Loop subdivision(regular subdivision)
			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
			glDisable(GL_LIGHT1);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			if(5 == previousColor)
			{
				L = L+0.03;
			}
			else
			{
				previousColor = 5;

			}

			for (i=0; i<polyRegular->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=polyRegular->tlist[i];
				glBegin(GL_POLYGON);

				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];

					int R = floor(temp_v->x/L);
					if(R%2==0)
					{
						R=1;
					}
					else
					{
						R=0;
					}
					int G =floor(temp_v->y/L);
					if(G%2==0)
					{
						G=1;
					}
					else
					{
						G=0;
					}
					int B =floor(temp_v->z/L);
					if(B%2==0)
					{
						B=1;
					}
					else
					{
						B=0;
					}
					glColor3f(R,G,B);

					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

			}
			break;
		case 6: // 3dCheckboard Color scheme for irregular subdivision
			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
			glDisable(GL_LIGHT1);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			if(previousColor==6)
			{
				L = L+0.03;
			}
			else
			{
				previousColor = 6;

			}
			for (i=0; i<polyIrregular->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=polyIrregular->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];

					int R = floor(temp_v->x/L);
					if(R%2==0)
					{
						R=1;
					}
					else
					{
						R=0;
					}
					int G =floor(temp_v->y/L);
					if(G%2==0)
					{
						G=1;
					}
					else
					{
						G=0;
					}
					int B =floor(temp_v->z/L);
					if(B%2==0)
					{
						B=1;
					}
					else
					{
						B=0;
					}
					glColor3f(R,G,B);

					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

			}
			break;

		case 7: // 3D checker board Color Scheme
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			Total_Angle_def = 0.0;
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (int i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glColor3f(1.0, 1.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}


			Total_def_Count = 0;
			glPointSize(8);
			glBegin(GL_POINTS);
			for(int k=0; k<this_poly->nverts;k++)
			{
				int def = 6 - this_poly->vlist[k]->ntris;
				if(def != 0)
				{
					Total_def_Count = Total_def_Count +def;
					if(1 == def)
					{
						glColor3f(0.0, 0.0, 1.0);
						glVertex3d(this_poly->vlist[k]->x, this_poly->vlist[k]->y,this_poly->vlist[k]->z);
					}
					if(2 == def)
					{
						glColor3f(0.0, 1.0, 0.0);
						glVertex3d(this_poly->vlist[k]->x, this_poly->vlist[k]->y,this_poly->vlist[k]->z);
					}
					if(3 == def)
					{
						glColor3f(0.0, 1.0, 1.0);
						glVertex3d(this_poly->vlist[k]->x, this_poly->vlist[k]->y,this_poly->vlist[k]->z);
					}
					if(-1 == def)
					{
						glColor3f(1.0, 0.0, 0.0);
						glVertex3d(this_poly->vlist[k]->x, this_poly->vlist[k]->y,this_poly->vlist[k]->z);
					}
					if(-2 == def)
					{
						glColor3f(1.0, 0.0, 1.0);
						glVertex3d(this_poly->vlist[k]->x, this_poly->vlist[k]->y,this_poly->vlist[k]->z);
					}
					if(-3 == def)
					{
						glColor3f(1.0, 1.0, 0.0);
						glVertex3d(this_poly->vlist[k]->x, this_poly->vlist[k]->y,this_poly->vlist[k]->z);
					}
					if(-3 > def)
					{
						glColor3f(1.0, 1.0, 1.0);
						glVertex3d(this_poly->vlist[k]->x, this_poly->vlist[k]->y,this_poly->vlist[k]->z);
					}
				}

				float totalAngleAtTheVertex = 0.f;
				for(int m = 0; m<this_poly->vlist[k]->ntris;m++)	
				{
					Triangle *currentTriangle = this_poly->vlist[k]->tris[m];
					Vertex *v0 = this_poly->vlist[k];
					Vertex *v1 = NULL;
					Vertex *v2 = NULL;

					if(currentTriangle->verts[0]->x == v0->x &&
						currentTriangle->verts[0]->y == v0->y &&
						currentTriangle->verts[0]->y == v0->y )
					{
						v1 = currentTriangle->verts[1];
						v2 = currentTriangle->verts[2];
					}
					else if(currentTriangle->verts[1]->x == v0->x &&
						currentTriangle->verts[1]->y == v0->y &&
						currentTriangle->verts[1]->y == v0->y )
					{
						v1 = currentTriangle->verts[0];
						v2 = currentTriangle->verts[2];
					}
					if(currentTriangle->verts[2]->x == v0->x &&
						currentTriangle->verts[2]->y == v0->y &&
						currentTriangle->verts[2]->y == v0->y )
					{
						v1 = currentTriangle->verts[0];
						v2 = currentTriangle->verts[1];
					}
					if(NULL != v0 && NULL != v1 && NULL != v2)
					{
						icVector3 vector1( (v1->x - v0->x),(v1->y - v0->y),(v1->z - v0->z));
						icVector3 vector2( (v2->x - v0->x),(v2->y - v0->y),(v2->z - v0->z));
						double dotProduct = dot(vector1,vector2);
						double productLength = length(vector1)*length(vector2);
						double cosTheta = dotProduct/productLength;
						double theta = acos(cosTheta);
						double angle = theta*(180/3.14159265358979323846);
						totalAngleAtTheVertex =totalAngleAtTheVertex + angle;

					}
				}
				double adef = 360 - totalAngleAtTheVertex;
				Total_Angle_def =   Total_Angle_def + adef;
			}

			glEnd();
			printf("\nTotal Valence deficit for the model: %d",Total_def_Count);
			printf("\nTotal Angle deficit for the model is: %f",Total_Angle_def);
			break;

		case 8: // 3D ChekerBoard Color scheme with loop subdivision(Regular subdivision)
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;
			Total_Angle_def = 0.0;
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (int i=0; i<polyRegular->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=polyRegular->tlist[i];
				glBegin(GL_POLYGON);
				for (int j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glColor3f(1.0, 1.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			glPointSize(8);
			glBegin(GL_POINTS);
			for(int k=0; k<polyRegular->nverts;k++)
			{
				int def = 6 - polyRegular->vlist[k]->ntris;
				if(def != 0)
				{
					Total_def_Count = Total_def_Count + def;
					if(1 == def)
					{
						glColor3f(0.0, 0.0, 1.0);
						glVertex3d(polyRegular->vlist[k]->x, polyRegular->vlist[k]->y,polyRegular->vlist[k]->z);
					}
					if(2 == def)
					{
						glColor3f(0.0, 1.0, 0.0);
						glVertex3d(polyRegular->vlist[k]->x, polyRegular->vlist[k]->y,polyRegular->vlist[k]->z);
					}
					if(3 == def)
					{
						glColor3f(0.0, 1.0, 1.0);
						glVertex3d(polyRegular->vlist[k]->x, polyRegular->vlist[k]->y,polyRegular->vlist[k]->z);
					}
					if(-1 == def)
					{
						glColor3f(1.0, 0.0, 0.0);
						glVertex3d(polyRegular->vlist[k]->x, polyRegular->vlist[k]->y,polyRegular->vlist[k]->z);
					}
					if(-2 == def)
					{
						glColor3f(1.0, 0.0, 1.0);
						glVertex3d(polyRegular->vlist[k]->x, polyRegular->vlist[k]->y,polyRegular->vlist[k]->z);
					}
					if(-3 == def)
					{
						glColor3f(1.0, 1.0, 0.0);
						glVertex3d(polyRegular->vlist[k]->x, polyRegular->vlist[k]->y,polyRegular->vlist[k]->z);
					}
					if(-3 > def)
					{
						glColor3f(1.0, 1.0, 1.0);
						glVertex3d(polyRegular->vlist[k]->x, polyRegular->vlist[k]->y,polyRegular->vlist[k]->z);
					}
				}

				float totalAngleAtTheVertex = 0.f;
				for(int m = 0; m<polyRegular->vlist[k]->ntris;m++)	
				{
					Triangle *currentTriangle = polyRegular->vlist[k]->tris[m];
					Vertex *v0 = polyRegular->vlist[k];
					Vertex *v1 = NULL;
					Vertex *v2 = NULL;

					if(currentTriangle->verts[0]->x == v0->x &&
						currentTriangle->verts[0]->y == v0->y &&
						currentTriangle->verts[0]->y == v0->y )
					{
						v1 = currentTriangle->verts[1];
						v2 = currentTriangle->verts[2];
					}
					else if(currentTriangle->verts[1]->x == v0->x &&
						currentTriangle->verts[1]->y == v0->y &&
						currentTriangle->verts[1]->y == v0->y )
					{
						v1 = currentTriangle->verts[0];
						v2 = currentTriangle->verts[2];
					}
					if(currentTriangle->verts[2]->x == v0->x &&
						currentTriangle->verts[2]->y == v0->y &&
						currentTriangle->verts[2]->y == v0->y )
					{
						v1 = currentTriangle->verts[0];
						v2 = currentTriangle->verts[1];
					}
					if(NULL != v0 && NULL != v1 && NULL != v2)
					{
						icVector3 vector1( (v1->x - v0->x),(v1->y - v0->y),(v1->z - v0->z));
						icVector3 vector2( (v2->x - v0->x),(v2->y - v0->y),(v2->z - v0->z));
						double dotProduct = dot(vector1,vector2);
						double productLength = length(vector1)*length(vector2);
						double cosTheta = dotProduct/productLength;
						double theta = acos(cosTheta);
						double angle = theta*(180/3.14159265358979323846);
						totalAngleAtTheVertex =totalAngleAtTheVertex + angle;

					}
				}
				double adef = 360 - totalAngleAtTheVertex;
				Total_Angle_def =   Total_Angle_def + adef;
			}

			glEnd();
			printf("\nTotal Valence deficit for Model with Regular Division is: %d",Total_def_Count);
			printf("\nTotal Angle deficit for Model with Regular Division Mesh is:%f ",Total_Angle_def);
			break;



		case 9: // 3D Checkboard with Irregular subdivision
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;
			Total_Angle_def = 0.0;
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (int i=0; i<polyIrregular->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=polyIrregular->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glColor3f(1.0, 1.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (int j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			glPointSize(8);
			glBegin(GL_POINTS);
			for(int k=0; k<polyIrregular->nverts;k++)
			{
				int def = 6 - polyIrregular->vlist[k]->ntris;
				if(def != 0)
				{
					Total_def_Count = Total_def_Count + def;
					if(1 == def)
					{
						glColor3f(0.0, 0.0, 1.0);
						glVertex3d(polyIrregular->vlist[k]->x, polyIrregular->vlist[k]->y,polyIrregular->vlist[k]->z);
					}
					if(2 == def)
					{
						glColor3f(0.0, 1.0, 0.0);
						glVertex3d(polyIrregular->vlist[k]->x, polyIrregular->vlist[k]->y,polyIrregular->vlist[k]->z);
					}
					if(3 == def)
					{
						glColor3f(0.0, 1.0, 1.0);
						glVertex3d(polyIrregular->vlist[k]->x, polyIrregular->vlist[k]->y,polyIrregular->vlist[k]->z);
					}
					if(-1 == def)
					{
						glColor3f(1.0, 0.0, 0.0);
						glVertex3d(polyIrregular->vlist[k]->x, polyIrregular->vlist[k]->y,polyIrregular->vlist[k]->z);
					}
					if(-2 == def)
					{
						glColor3f(1.0, 0.0, 1.0);
						glVertex3d(polyIrregular->vlist[k]->x, polyIrregular->vlist[k]->y,polyIrregular->vlist[k]->z);
					}
					if(-3 == def)
					{
						glColor3f(1.0, 1.0, 0.0);
						glVertex3d(polyIrregular->vlist[k]->x, polyIrregular->vlist[k]->y,polyIrregular->vlist[k]->z);
					}
					if(-3 > def)
					{
						glColor3f(1.0, 1.0, 1.0);
						glVertex3d(polyIrregular->vlist[k]->x, polyIrregular->vlist[k]->y,polyIrregular->vlist[k]->z);
					}
				}

				float totalAngleAtTheVertex = 0.f;
				for(int m = 0; m<polyIrregular->vlist[k]->ntris;m++)	
				{
					Triangle *currentTriangle = polyIrregular->vlist[k]->tris[m];
					Vertex *v0 = polyIrregular->vlist[k];
					Vertex *v1 = NULL;
					Vertex *v2 = NULL;

					if(currentTriangle->verts[0]->x == v0->x &&
						currentTriangle->verts[0]->y == v0->y &&
						currentTriangle->verts[0]->y == v0->y )
					{
						v1 = currentTriangle->verts[1];
						v2 = currentTriangle->verts[2];
					}
					else if(currentTriangle->verts[1]->x == v0->x &&
						currentTriangle->verts[1]->y == v0->y &&
						currentTriangle->verts[1]->y == v0->y )
					{
						v1 = currentTriangle->verts[0];
						v2 = currentTriangle->verts[2];
					}
					if(currentTriangle->verts[2]->x == v0->x &&
						currentTriangle->verts[2]->y == v0->y &&
						currentTriangle->verts[2]->y == v0->y )
					{
						v1 = currentTriangle->verts[0];
						v2 = currentTriangle->verts[1];
					}
					if(NULL != v0 && NULL != v1 && NULL != v2)
					{
						icVector3 vector1( (v1->x - v0->x),(v1->y - v0->y),(v1->z - v0->z));
						icVector3 vector2( (v2->x - v0->x),(v2->y - v0->y),(v2->z - v0->z));
						double dotProduct = dot(vector1,vector2);
						double productLength = length(vector1)*length(vector2);
						double cosTheta = dotProduct/productLength;
						double theta = acos(cosTheta);
						double angle = theta*(180/3.14159265358979323846);
						totalAngleAtTheVertex =totalAngleAtTheVertex + angle;

					}
				}
				double adef = 360 - totalAngleAtTheVertex;
				Total_Angle_def =   Total_Angle_def + adef;
			}

			glEnd();
			printf("\nTotal Valence deficit for Model with Regular Division is: %d",Total_def_Count);
			printf("\nTotal Angle deficit for Model with Regular Division Mesh is: %f",Total_Angle_def);
			break;


		case 10:
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {
					mat_diffuse[0] = 1.0;
					mat_diffuse[1] = 0.0;
					mat_diffuse[2] = 0.0;
					mat_diffuse[3] = 1.0;

					glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);

					glColor3f(1.0, 0.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;



			//		glDisable(GL_LIGHTING);
			//		glDisable(GL_LIGHT0);
			//	glDisable(GL_LIGHT1);
			//		if (i == this_poly->seed) {
			//			mat_diffuse[0] = 0.0;
			//			mat_diffuse[1] = 0.0;
			//			mat_diffuse[2] = 1.0;
			//			mat_diffuse[3] = 1.0;
			//		} else {
			//			mat_diffuse[0] = 1.0;
			//			mat_diffuse[1] = 1.0;
			//			mat_diffuse[2] = 0.0;
			//			mat_diffuse[3] = 1.0;
			//		}
			//		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			//		glColor3ub((i)%255, ((i)*(i))%255,((i)*(i)*(i))%255);
			//		glBegin(GL_POLYGON);
			//		for (j=0; j<3; j++) {
			//			
			//			Vertex *temp_v = temp_t->verts[j];
			//			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
			//			//glColor3ub((j+i)%255, ((j+i)*(j+i))%255,((j+i)*(j+i)*(j+i))%255);
			//			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			//		}
			//		glEnd();

			//		break;

			//	case 2:
			//		glDisable(GL_LIGHTING);
			//		glDisable(GL_LIGHT0);
			//	glDisable(GL_LIGHT1);
			//		if (i == this_poly->seed) {
			//			mat_diffuse[0] = 0.0;
			//			mat_diffuse[1] = 0.0;
			//			mat_diffuse[2] = 1.0;
			//			mat_diffuse[3] = 1.0;
			//		} else {
			//			mat_diffuse[0] = 1.0;
			//			mat_diffuse[1] = 1.0;
			//			mat_diffuse[2] = 0.0;
			//			mat_diffuse[3] = 1.0;
			//		}
			//		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			//		glBegin(GL_POLYGON);
			//		
			//		
			//		for (j=0; j<3; j++) {
			//			
			//			Vertex *temp_v = temp_t->verts[j];
			//			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
			//
			//			glColor3f(fabs(temp_v->normal.entry[0]),fabs(temp_v->normal.entry[1]),fabs(temp_v->normal.entry[2]));
			//			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			//		}
			//		glEnd();
			//		break;

			//		case 3:
			//		glDisable(GL_LIGHTING);
			//		glDisable(GL_LIGHT0);
			//		glDisable(GL_LIGHT1);
			//		if (i == this_poly->seed) {
			//			mat_diffuse[0] = 0.0;
			//			mat_diffuse[1] = 0.0;
			//			mat_diffuse[2] = 1.0;
			//			mat_diffuse[3] = 1.0;
			//		} else {
			//			mat_diffuse[0] = 1.0;
			//			mat_diffuse[1] = 1.0;
			//			mat_diffuse[2] = 0.0;
			//			mat_diffuse[3] = 1.0;
			//		}
			//		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			//		glBegin(GL_POLYGON);
			//		
			//		
			//		for (j=0; j<3; j++) {
			//			
			//			Vertex *temp_v = temp_t->verts[j];

			//            int R = floor(temp_v->x/0.0001);
			//			if(R%2==0)
			//			{
			//				R=1;
			//			}
			//			else
			//			{
			//				R=0;
			//			}
			//			int G =floor(temp_v->y/0.0001);
			//			if(G%2==0)
			//			{
			//				G=1;
			//			}
			//			else
			//			{
			//				G=0;
			//			}
			//			int B =floor(temp_v->z/0.0001);
			//			if(B%2==0)
			//			{
			//				B=1;
			//			}
			//			else
			//			{
			//				B=0;
			//			}
			//			glColor3f(R,G,B);
			//			
			//			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
			//			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			//		}
			//		glEnd();
			//		break;

			//	case 4:
			//		//glDisable(GL_LIGHTING);
			//		//glDisable(GL_LIGHT0);
			//		//glDisable(GL_LIGHT1);
			//		//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			//		glColor3f(1.0,0,0);
			//		glLineWidth(10);
			//		glBegin(GL_LINES);

			//		glVertex3d(boundbox1[0][0],boundbox1[1][0],boundbox1[2][0]);
			//		glVertex3d(boundbox1[0][1],boundbox1[1][1],boundbox1[2][1]);		

			//		glVertex3d(boundbox1[0][1],boundbox1[1][1],boundbox1[2][1]);
			//		glVertex3d(boundbox1[0][2],boundbox1[1][2],boundbox1[2][2]);		

			//		glVertex3d(boundbox1[0][2],boundbox1[1][2],boundbox1[2][2]);
			//		glVertex3d(boundbox1[0][3],boundbox1[1][3],boundbox1[2][3]);		

			//		glVertex3d(boundbox1[0][3],boundbox1[1][3],boundbox1[2][3]);
			//		glVertex3d(boundbox1[0][0],boundbox1[1][0],boundbox1[2][0]);		

			//		glVertex3d(boundbox1[0][0],boundbox1[1][0],boundbox1[2][0]);
			//		glVertex3d(boundbox1[0][4],boundbox1[1][4],boundbox1[2][4]);		

			//		glVertex3d(boundbox1[0][1],boundbox1[1][1],boundbox1[2][1]);
			//		glVertex3d(boundbox1[0][5],boundbox1[1][5],boundbox1[2][5]);		

			//		glVertex3d(boundbox1[0][2],boundbox1[1][2],boundbox1[2][2]);
			//		glVertex3d(boundbox1[0][6],boundbox1[1][6],boundbox1[2][6]);		

			//		glVertex3d(boundbox1[0][3],boundbox1[1][3],boundbox1[2][3]);
			//		glVertex3d(boundbox1[0][7],boundbox1[1][7],boundbox1[2][7]);	

			//		glVertex3d(boundbox1[0][4],boundbox1[1][4],boundbox1[2][4]);
			//		glVertex3d(boundbox1[0][7],boundbox1[1][7],boundbox1[2][7]);		

			//		glVertex3d(boundbox1[0][7],boundbox1[1][7],boundbox1[2][7]);
			//		glVertex3d(boundbox1[0][6],boundbox1[1][6],boundbox1[2][6]);		

			//		glVertex3d(boundbox1[0][6],boundbox1[1][6],boundbox1[2][6]);
			//		glVertex3d(boundbox1[0][5],boundbox1[1][5],boundbox1[2][5]);		

			//		glVertex3d(boundbox1[0][5],boundbox1[1][5],boundbox1[2][5]);
			//		glVertex3d(boundbox1[0][4],boundbox1[1][4],boundbox1[2][4]);
			//
			//		glEnd();

			//		glBegin(GL_POLYGON);
			//		for (j=0; j<3; j++) {
			//			
			//			Vertex *temp_v = temp_t->verts[j];
			//			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);

			//			glColor3f(1.0, 1.0, 0.0);
			//			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			//		}
			//		glEnd();
			//		break;

			//	case 5:
			///*		glDisable(GL_LIGHTING);
			//		glDisable(GL_LIGHT0);
			//		glDisable(GL_LIGHT1);
			//		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);*/
			//		glColor3f(1.0,0,0);
			//		glLineWidth(10);
			//		glBegin(GL_LINES);

			//		glVertex3d(boundbox2[0][0],boundbox2[1][0],boundbox2[2][0]);
			//		glVertex3d(boundbox2[0][1],boundbox2[1][1],boundbox2[2][1]);		

			//		glVertex3d(boundbox2[0][1],boundbox2[1][1],boundbox2[2][1]);
			//		glVertex3d(boundbox2[0][2],boundbox2[1][2],boundbox2[2][2]);		

			//		glVertex3d(boundbox2[0][2],boundbox2[1][2],boundbox2[2][2]);
			//		glVertex3d(boundbox2[0][3],boundbox2[1][3],boundbox2[2][3]);		

			//		glVertex3d(boundbox2[0][3],boundbox2[1][3],boundbox2[2][3]);
			//		glVertex3d(boundbox2[0][0],boundbox2[1][0],boundbox2[2][0]);		

			//		glVertex3d(boundbox2[0][0],boundbox2[1][0],boundbox2[2][0]);
			//		glVertex3d(boundbox2[0][4],boundbox2[1][4],boundbox2[2][4]);		

			//		glVertex3d(boundbox2[0][1],boundbox2[1][1],boundbox2[2][1]);
			//		glVertex3d(boundbox2[0][5],boundbox2[1][5],boundbox2[2][5]);		

			//		glVertex3d(boundbox2[0][2],boundbox2[1][2],boundbox2[2][2]);
			//		glVertex3d(boundbox2[0][6],boundbox2[1][6],boundbox2[2][6]);		

			//		glVertex3d(boundbox2[0][3],boundbox2[1][3],boundbox2[2][3]);
			//		glVertex3d(boundbox2[0][7],boundbox2[1][7],boundbox2[2][7]);	

			//		glVertex3d(boundbox2[0][4],boundbox2[1][4],boundbox2[2][4]);
			//		glVertex3d(boundbox2[0][7],boundbox2[1][7],boundbox2[2][7]);		

			//		glVertex3d(boundbox2[0][7],boundbox2[1][7],boundbox2[2][7]);
			//		glVertex3d(boundbox2[0][6],boundbox2[1][6],boundbox2[2][6]);		

			//		glVertex3d(boundbox2[0][6],boundbox2[1][6],boundbox2[2][6]);
			//		glVertex3d(boundbox2[0][5],boundbox2[1][5],boundbox2[2][5]);		

			//		glVertex3d(boundbox2[0][5],boundbox2[1][5],boundbox2[2][5]);
			//		glVertex3d(boundbox2[0][4],boundbox2[1][4],boundbox2[2][4]);
			//
			//		glEnd();

			//		glBegin(GL_POLYGON);
			//		for (j=0; j<3; j++) {
			//			
			//			Vertex *temp_v = temp_t->verts[j];
			//			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);

			//			glColor3f(1.0, 1.0, 0.0);
			//			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			//		}
			//		glEnd();

			//		break;

			//		case 6:
			///*		glDisable(GL_LIGHTING);
			//		glDisable(GL_LIGHT0);
			//		glDisable(GL_LIGHT1);
			//		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);*/
			//		glColor3f(1.0,0,0);
			//		glLineWidth(10);
			//		glBegin(GL_LINES);

			//		glVertex3d(boundbox3[0][0],boundbox3[1][0],boundbox3[2][0]);
			//		glVertex3d(boundbox3[0][1],boundbox3[1][1],boundbox3[2][1]);		

			//		glVertex3d(boundbox3[0][1],boundbox3[1][1],boundbox3[2][1]);
			//		glVertex3d(boundbox3[0][2],boundbox3[1][2],boundbox3[2][2]);		

			//		glVertex3d(boundbox3[0][2],boundbox3[1][2],boundbox3[2][2]);
			//		glVertex3d(boundbox3[0][3],boundbox3[1][3],boundbox3[2][3]);		

			//		glVertex3d(boundbox3[0][3],boundbox3[1][3],boundbox3[2][3]);
			//		glVertex3d(boundbox3[0][0],boundbox3[1][0],boundbox3[2][0]);		

			//		glVertex3d(boundbox3[0][0],boundbox3[1][0],boundbox3[2][0]);
			//		glVertex3d(boundbox3[0][4],boundbox3[1][4],boundbox3[2][4]);		

			//		glVertex3d(boundbox3[0][1],boundbox3[1][1],boundbox3[2][1]);
			//		glVertex3d(boundbox3[0][5],boundbox3[1][5],boundbox3[2][5]);		

			//		glVertex3d(boundbox3[0][2],boundbox3[1][2],boundbox3[2][2]);
			//		glVertex3d(boundbox3[0][6],boundbox3[1][6],boundbox3[2][6]);		

			//		glVertex3d(boundbox3[0][3],boundbox3[1][3],boundbox3[2][3]);
			//		glVertex3d(boundbox3[0][7],boundbox3[1][7],boundbox3[2][7]);	

			//		glVertex3d(boundbox3[0][4],boundbox3[1][4],boundbox3[2][4]);
			//		glVertex3d(boundbox3[0][7],boundbox3[1][7],boundbox3[2][7]);		

			//		glVertex3d(boundbox3[0][7],boundbox3[1][7],boundbox3[2][7]);
			//		glVertex3d(boundbox3[0][6],boundbox3[1][6],boundbox3[2][6]);		

			//		glVertex3d(boundbox3[0][6],boundbox3[1][6],boundbox3[2][6]);
			//		glVertex3d(boundbox3[0][5],boundbox3[1][5],boundbox3[2][5]);		

			//		glVertex3d(boundbox3[0][5],boundbox3[1][5],boundbox3[2][5]);
			//		glVertex3d(boundbox3[0][4],boundbox3[1][4],boundbox3[2][4]);
			//
			//		glEnd();

			//		glBegin(GL_POLYGON);
			//		for (j=0; j<3; j++) {
			//			
			//			Vertex *temp_v = temp_t->verts[j];
			//			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);

			//			glColor3f(1.0, 1.0, 0.0);
			//			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			//		}
			//		glEnd();

			//		break;

			//	case 7:
			//		glEnable(GL_LIGHTING);
			//		glEnable(GL_LIGHT0);
			//		glEnable(GL_LIGHT1);
			//		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			//		glBegin(GL_POLYGON);
			//		for (j=0; j<3; j++) {
			//			Vertex *temp_v = temp_t->verts[j];
			//			glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
			//			glColor3f(1.0, 1.0, 1.0);
			//			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			//		}
			//		glEnd();
			//		break;



		}
	}
	else if(selection==2)
	{

		switch (display_mode) {
		case 0: //Default
			L = 0.05f;
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				if (i == this_poly->seed) {
					mat_diffuse[0] = 0.0;
					mat_diffuse[1] = 0.0;
					mat_diffuse[2] = 1.0;
					mat_diffuse[3] = 1.0;
				} else {
					mat_diffuse[0] = 1.0;
					mat_diffuse[1] = 1.0;
					mat_diffuse[2] = 0.0;
					mat_diffuse[3] = 1.0;
				}

				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);


				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					if (i==this_poly->seed)
						glColor3f(0.0, 0.0, 1.0);
					else
						glColor3f(1.0, 1.0, 0.0);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;


		case 1: // Mesh smoothing with uniform scheme

			poly->MeshSmoothingUniformScheme();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					if (i==this_poly->seed)
						glColor3f(0.0, 0.0, 1.0);
					else
						glColor3f(1.0, 1.0, 0.0);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;

		case 2: // Mesh Smoothing with Cord Scheme

			poly->MeshSmoothingCordScheme();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					if (i==this_poly->seed)
						glColor3f(0.0, 0.0, 1.0);
					else
						glColor3f(1.0, 1.0, 0.0);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;


		case 3: // Mesh Smoothing with Mean Curvature Flow Scheme

			poly->MeshSmoothingMeanCurvatureFlowScheme();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					if (i==this_poly->seed)
						glColor3f(0.0, 0.0, 1.0);
					else
						glColor3f(1.0, 1.0, 0.0);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;


		case 4: // Mesh smoothing with Mean Value Coordinate Scheme

			poly->MeshSmoothingMeanValueScheme();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					if (i==this_poly->seed)
						glColor3f(0.0, 0.0, 1.0);
					else
						glColor3f(1.0, 1.0, 0.0);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;


		case 5: // Mesh Smoothing with Cord Scheme	with initial mesh
			poly->MeshSmoothingCordScheme_WithInitialMesh();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glColor3f(1.0, 1.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;

		case 6: // Mesh Smoothing with Mean Curvature Flow Scheme with initial mesh
			poly->MeshSmoothingMeanCurvatureFlowScheme_WithInitialMesh();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glColor3f(1.0, 1.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;

		case 7: // Mesh Smoothing with Mean Value Coordinate Scheme with initial mesh

			poly->MeshSmoothingMeanValueScheme_WithInitialMesh();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glColor3f(1.0, 1.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;

		case 8: // Morse Design based on one point function using Explicit Method 
			this_poly->MorseDesignUsingExplicitMethod();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					//glColor3f(1.0, 1.0, 0.0);
					glColor3f(temp_v->MorseVal,temp_v->MorseVal,temp_v->MorseVal);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0f, 0.0f, 0.0f);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}

			FindSaddlesAndLocalMaxMin(this_poly);

			break;

		case 9: // Morse Design based on one point function using Gauss Seidal Method 
			this_poly->MorseDesignUsingGaussSeidalMethod();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glColor3f(temp_v->MorseVal,temp_v->MorseVal,temp_v->MorseVal);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0f, 0.0f, 0.0f);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}

			FindSaddlesAndLocalMaxMin(this_poly);

			break;



		case 10:
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {
					mat_diffuse[0] = 1.0;
					mat_diffuse[1] = 0.0;
					mat_diffuse[2] = 0.0;
					mat_diffuse[3] = 1.0;

					glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);

					glColor3f(1.0, 0.0, 0.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			break;

		case 11:// Morse Design based on two points using Explicit Method 
			poly->MorseDesignUsingExplicitMethod();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<this_poly->ntris; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];

				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0f, 0.0f, 0.0f);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}
			FindSaddlesAndLocalMaxMinWithTexture(this_poly);

			break;

		case 12:// Morse Design based on two points using Gauss Seidal Method 
			this_poly->MorseDesignUsingGaussSeidalMethod();
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;

			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<this_poly->ntris; i++) 
			{
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=this_poly->tlist[i];

				glLineWidth(1);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
				glDisable(GL_LIGHT1);
				glColor3f(0.0, 0.0, 0.0);

				glBegin(GL_LINE_LOOP);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();
			}

			FindSaddlesAndLocalMaxMinWithTexture(this_poly);
			break;




		case 13: // Mesh smoothing with uniform scheme
			L=0.08;
			poly->MeshSmoothingMeanValueScheme();	
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;
			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
			glDisable(GL_LIGHT1);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			for (i=0; i<poly->ntris; i++) 
			{
				if (mode == GL_SELECT)
					glLoadName(i+1);

				Triangle *temp_t=poly->tlist[i];
				glBegin(GL_POLYGON);
				for (j=0; j<3; j++) {

					Vertex *temp_v = temp_t->verts[j];

					int R = floor(temp_v->x/L);
					if(R%2==0)
					{
						R=1;
					}
					else
					{
						R=0;
					}
					int G =floor(temp_v->y/L);
					if(G%2==0)
					{
						G=1;
					}
					else
					{
						G=0;
					}
					int B =floor(temp_v->z/L);
					if(B%2==0)
					{
						B=1;
					}
					else
					{
						B=0;
					}
					glColor3f(R,G,B);

					glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);

					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glEnd();

			}
			break;

		}

	}


}


void display(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glGetIntegerv (GL_VIEWPORT, viewport);

	glClear(GL_ACCUM_BUFFER_BIT);
	for (jitter = 0; jitter < ACSIZE; jitter++) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);
		glPushMatrix ();
		switch(ACSIZE)
		{
		case 1:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;

		case 16:
			glTranslatef (ji16[jitter].x*2.0/viewport[2], ji16[jitter].y*2.0/viewport[3], 0.0);
			break;

		default:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;
		}
		set_scene(GL_RENDER, poly);
		display_shape(GL_RENDER, poly);
		glPopMatrix ();
		glAccum(GL_ACCUM, 1.0/ACSIZE);
	}
	glAccum (GL_RETURN, 1.0);
	glFlush();
	glutSwapBuffers();
	glFinish();
}

void Polyhedron::average_normals()
{
	int i, j;

	for (i=0; i<nverts; i++) {
		vlist[i]->normal = icVector3(0.0);
		for (j=0; j<vlist[i]->ntris; j++) 
			vlist[i]->normal += vlist[i]->tris[j]->normal;
		normalize(vlist[i]->normal);
	}
}



void calculate_eigens(float a[][3], int n, float d[], float v[][3])
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c;
	float b[3];
	float z[3];
	for (ip=0;ip<n;ip++) { 
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) { 
		b[ip]=d[ip]=a[ip][ip]; 
		z[ip]=0.0;		
	}
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) { 
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) { 			
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]); 
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq);
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq);
					}
					for (j=iq+1;j<n;j++) { 
						ROTATE(a,ip,j,iq,j);
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq);
					}						
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip]; 
			z[ip]=0.0; 
		}
	}

}


