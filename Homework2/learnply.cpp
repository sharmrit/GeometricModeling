
/*

Functions for learnply

Eugene Zhang, 2005
*/
#include <iostream>
using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "glut.h"
#include <string.h>
#include <fstream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "learnply_io.h"
#include "trackball.h"
#include "tmatrix.h"
//#include "NumRecipes.h"

static PlyFile *in_ply;

unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width=1024;
const int win_height=1024;

double radius_factor = 0.9;

int display_mode = 0; 
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;	
int ACSIZE = 1; // for antialiasing
int view_mode=0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;
float L = 0.0f;
int CreateCorner = 0;
int lastcolorscheme = 0;
float boundingBox[3][8];  // Orientation
float boundingBox1[3][8]; // Axis Aligned
float boundingBox2[3][8]; // Normal
float majoreigen[3];
float middleeigen[3];
float minoreigen[3];
int RGB[3]= {0,0,0};

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

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
Polyhedron *polyR = NULL;
Polyhedron *polyIR = NULL;
int isSubdived = 0;
float uLimit;

void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);
void calculate_eigens(float a[][3], int n, float d[], float v[][3]);
void calculate_eigens1(float a[][3], int n, float d[], float v[][3]);
/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char *argv[])
{
  char *progname;
  int num = 1;
	FILE *this_file;

  progname = argv[0];

  string sFile("C:/Users/Prashant/Documents/Visual Studio 2012/Projects/project1/tempmodels/");
  int plyNumber = -1;
  cout<<"1: bunny"<<endl;	  
  cout<<"2: dodecahedron"<<endl;	  
  cout<<"3: dragon"<<endl;	  
  cout<<"4: feline"<<endl;	  
  cout<<"5: happy"<<endl;	  
  cout<<"6: hexahedron"<<endl;	  
  cout<<"7: icosahedron"<<endl;	  
  cout<<"8: octahedron"<<endl;	  
  cout<<"9: sphere"<<endl;	  
  cout<<"10: tetrahedron"<<endl;	  
  cout<<"11: torus"<<endl;
  cout<<"enter the file number from 1 to 11"<<endl;
  cin>>plyNumber;
  //cout<<"enter a real number for checkboard color scheme"<<endl;
  //cin>>L;
  cout<<"enter a value of U for Irregular Subdivision"<<endl;
  cin>>uLimit;
  cout<<"enter 1 to construct CornerList"<<endl;
  cin>>CreateCorner;
  cout<<"Color Scheme: choose 0-6"<<endl;
  cout<<"0: Default"<<endl;
  cout<<"1: Normal Mesh"<<endl;
  cout<<"2: Regular Subdivision"<<endl;
  cout<<"3: Irregular Subdivision"<<endl;
  cout<<"4: 3D CheckerBoard for Normal Mesh"<<endl;
  cout<<"5: 3D CheckerBoard for Regular Division Mesh"<<endl;
  cout<<"6: 3D CheckerBoard for Irregular Division Mesh"<<endl;
  cout<<"7: Valence deficiency in Normal Mesh"<<endl;
  cout<<"8: Valence deficiency in Regular Division Mesh"<<endl;
  cout<<"9: Valence deficiency in Irregular Division Mesh"<<endl;
  
  string sFileAdd;
  switch(plyNumber)
  {
  case 1:
	  sFileAdd = "bunny.ply";	  
	  break;

  case 2:	  
	  sFileAdd = "dodecahedron.ply";
	  break;

  case 3:
	  sFileAdd = "dragon.ply";
	  break;

  case 4:
	  sFileAdd = "feline.ply";
	  break;

  case 5:
	  sFileAdd = "happy.ply";
	  break;

  case 6:
	  sFileAdd = "hexahedron.ply";
	  break;

  case 7:
	  sFileAdd = "icosahedron.ply";
	  break;

  case 8:
	  sFileAdd = "octahedron.ply";
	  break;

  case 9:
	  sFileAdd = "sphere.ply";
	  break;

  case 10:
	  sFileAdd = "tetrahedron.ply";
	  break;

  case 11:
	  sFileAdd = "torus.ply";
	  break;

  default:
	  cout<<"Invalid Entry"<<endl;
	  return 0;
  }  

  const char *csFileAdd = sFileAdd.c_str();
  cout<<"Model Selected: "<<csFileAdd<<endl;
   
    sFile += sFileAdd;
	const char *csFile = sFile.c_str();
	this_file = fopen(csFile, "r");
	poly = new Polyhedron(this_file);
    
	
	fclose(this_file);
	mat_ident( rotmat );	

	poly->initialize(); // initialize everything
	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();
	if(1 == CreateCorner)
	{
	 poly->CreateCornerList();
	}
	polyR = new Polyhedron(this_file,1);
	polyR->initialize();

	polyIR = new Polyhedron(this_file,1.f);
	polyIR->initialize();

	
	// XD8 (1)
	// Here we will get the number of vertices, edges and Faces from the polyhedron 
	// and calculate the eulers number as EU = V -E + F
	int nEU = -1;
	int nV = -1;
	int nE = -1;
	int nF = -1;

	if(poly)
	{
	nV = poly->nverts;
	nE = poly->nedges;
	nF = poly->ntris;
	}

	nEU = nV - nE + nF;

    cout<<"The Euler value of the model is: " <<nV<< "-"<<nE<< "+"<< nF<< "="<< nEU <<endl;

	if(polyR)
	{
	nV = polyR->nverts;
	nE = polyR->nedges;
	nF = polyR->ntris;
	}

	nEU = nV - nE + nF;

	cout<<"The Euler value of the model with Regular Subdivison is: " <<nV<< "-"<<nE<< "+"<< nF<< "="<< nEU <<endl;

   if(polyIR)
	{
	nV = polyIR->nverts;
	nE = polyIR->nedges;
	nF = polyIR->ntris;
	}

	nEU = nV - nE + nF;

	cout<<"The Euler value of the model with Irregular Subdivison is: " <<nV<< "-"<<nE<< "+"<< nF<< "="<< nEU <<endl;
	// To calculate the number of handles of the model we use the genus value g, given by
	// nEU = 2 - 2g
	// or g = 1 - nEU/2;

	// For this count the number of handles by looking at the picture
	// Write the corresponding euler's value and derive the formula based on that


	// End XD8
	poly->calculate_bounding_box_with_moments();
	poly->calculate_bounding_box_with_normals();
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
	polyR->finalize();  
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

Polyhedron::Polyhedron(FILE *file)
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

Polyhedron::Polyhedron(FILE *file, float IRegular)
{
  int i,j;
  int elem_count;
  char *elem_name;

  // First lets tag all the edges with length less than uLimit
  // By Default, this tag is 0
  ntris = poly->ntris;
  nverts = 0;

  Vertex** vlist1 = new Vertex *[4*poly->nverts];
  for(i = 0; i<poly->nedges; i++)
  {
   Edge *currentEdge = poly->elist[i];
   if(currentEdge->length > uLimit)
   {
	   currentEdge->tag = 1;
   }
  }

  for(i = 0; i<poly->ntris; i++)
  {
	  int smallerEdgeCount=0;
	  Triangle *CurrentTriangle =  poly->tlist[i];
	  for(j=0;j<3;j++)
	  {
	  if( 1 == CurrentTriangle->edges[j]->tag)
	  {
	  smallerEdgeCount++;
	  }
	  }   
	  
	  if(1 == smallerEdgeCount)
	  {
	  ntris += 1;
	  }

	  if(2 == smallerEdgeCount)
	  {
	  ntris += 2;
	  }

	  if(3 == smallerEdgeCount)
	  {
	  ntris += 3;
	  }
  }

  tlist = new Triangle *[ntris];
  max_tris = ntris;
  int vCount = 0;
  int tCount = 0;
  for(i = 0; i<poly->ntris; i++)
  {
	  int smallerEdgeCount=0;
	  Triangle *CurrentTriangle =  poly->tlist[i];
	  Edge* tEdge0 = NULL;
	  Edge* tEdge1 = NULL;
	  Edge* tEdge2 = NULL;

	  Vertex *v0 = NULL;
	  Vertex *v1 = NULL;
	  Vertex *v2 = NULL;

	  Vertex *v3 = NULL;
	  Vertex *v4 = NULL;
	  Vertex *v5 = NULL;

      Vertex *v0nn = CurrentTriangle->verts[0];
      Vertex *v1nn = CurrentTriangle->verts[1];
      Vertex *v2nn = CurrentTriangle->verts[2];
	  
	  Vertex *v0nnn = v0nn;
      Vertex *v1nnn = v1nn;
      Vertex *v2nnn = v2nn;
	  

	  for(j=0;j<3;j++)
	  {
	  if( 1 == CurrentTriangle->edges[j]->tag)
	  {
	  smallerEdgeCount++;
	  if(1 == smallerEdgeCount)
	  {
	  tEdge0 = CurrentTriangle->edges[j];
	  v1nnn = tEdge0->verts[0];
	  v2nnn = tEdge0->verts[1];
	  for(int dd=0;dd<3;dd++)
	  {
	  //if(CurrentTriangle->verts[dd]!= v1nnn && CurrentTriangle->verts[dd]!= v2nnn);
	  //{
	 // v0nnn = CurrentTriangle->verts[dd];
	  //}
	  int vin1 = 0;
	  int vin2 = 0;
	  if(CurrentTriangle->verts[dd]->x == v1nnn->x && 
	 CurrentTriangle->verts[dd]->y == v1nnn->y && 
	 CurrentTriangle->verts[dd]->z == v1nnn->z)
	  {
	  vin1= 1;
	  }
	  else if(CurrentTriangle->verts[dd]->x == v2nnn->x && 
	 CurrentTriangle->verts[dd]->y == v2nnn->y && 
	 CurrentTriangle->verts[dd]->z == v2nnn->z)
	  {
	  vin2 = 1;
	  }
	  if(0 == vin1 && 0 == vin2)
	  {
	     v0nnn = CurrentTriangle->verts[dd];
	  }	  
	  }
	  }

	  if(2 == smallerEdgeCount)
	  {
	  tEdge1 = CurrentTriangle->edges[j];
	  // The vertex common in edge0 and edge1 is v1nnn
	  // The other vertex of edge0 is v0nnn
	  // the other vertex of edge1 is v2nnn
	  if(tEdge0->verts[0] == tEdge1->verts[0])
	  {
	  v1nnn = tEdge1->verts[0];
	  v0nnn = tEdge0->verts[1];
	  v2nnn = tEdge1->verts[1];
	  }
	  
	  if(tEdge0->verts[1] == tEdge1->verts[1])
	  {
	  v1nnn = tEdge1->verts[1];
	  v0nnn = tEdge0->verts[0];
	  v2nnn = tEdge1->verts[0];
	  }

	  if(tEdge0->verts[1] == tEdge1->verts[0])
	  {
	  v1nnn = tEdge1->verts[0];
	  v0nnn = tEdge0->verts[0];
	  v2nnn = tEdge1->verts[1];
	  }

	  if(tEdge0->verts[0] == tEdge1->verts[1])
	  {
	  v1nnn = tEdge1->verts[1];
	  v0nnn = tEdge0->verts[1];
	  v2nnn = tEdge1->verts[0];
	  }
	  
	  }
	  if(3 == smallerEdgeCount)
	  {
	  tEdge2 = CurrentTriangle->edges[j];
	  v0nnn = v0nn;
	  v1nnn = v1nn;
	  v2nnn = v2nn;
	  }
	  }
	  }   
	  
	  if(1 == smallerEdgeCount)
	  {
	  Vertex *v3n = new Vertex(v0nnn->x , v0nnn->y , v0nnn->z);
	      Vertex *v4n = new Vertex(v1nnn->x , v1nnn->y , v1nnn->z);
	      Vertex *v5n = new Vertex(v2nnn->x , v2nnn->y , v2nnn->z);
	  // Subdivide the current Triangle into two part
	  if(tEdge0 != NULL)
	  {
	  Vertex *v0n = new Vertex((tEdge0->verts[0]->x + tEdge0->verts[1]->x)/2 , 
	                           (tEdge0->verts[0]->y + tEdge0->verts[1]->y)/2 ,
	                           (tEdge0->verts[0]->z + tEdge0->verts[1]->z)/2);
	  v0n->ntris = 0;
	  for(int zz=0; zz<vCount ; zz++)
	  {
	  Vertex *currentV = vlist1[zz];
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
	  if(NULL == v0)
	  {
	  vlist1[vCount++] = v0n;
	  v0 = v0n;
	  nverts++;
	  }
	  if(NULL == v3)
	  {
	  vlist1[vCount++] = v3n;
	  v3 = v3n;
	  nverts++;
	  }
	  if(NULL == v4)
	  {
	  vlist1[vCount++] = v4n;
	  v4 = v4n;
	  nverts++;
	  }
	  if(NULL == v5)
	  {
	  vlist1[vCount++] = v5n;
	  v5 = v5n;
	  nverts++;
	  }
	  }
	  
	 
  	   this->tlist[tCount] = new Triangle;
       this->tlist[tCount]->nverts = 3;
       this->tlist[tCount]->verts[0] = v3;
       this->tlist[tCount]->verts[1] = v4;
       this->tlist[tCount]->verts[2] = v0;
	   tCount++;

	   this->tlist[tCount] = new Triangle;
       this->tlist[tCount]->nverts = 3;
       this->tlist[tCount]->verts[0] = v3;
       this->tlist[tCount]->verts[1] = v0;
       this->tlist[tCount]->verts[2] = v5;
	   tCount++;
	  }

	  if(2 == smallerEdgeCount)
	  {
	  // Subdivide the current Triangle into three part
	  Vertex *v3n = new Vertex(v0nnn->x , v0nnn->y , v0nnn->z);
	  Vertex *v4n = new Vertex(v1nnn->x , v1nnn->y , v1nnn->z);
	  Vertex *v5n = new Vertex(v2nnn->x , v2nnn->y , v2nnn->z);
	  // Subdivide the current Triangle into two part
	  if(tEdge0 != NULL && tEdge1 != NULL)
	  {
	  Vertex *v0n = new Vertex((tEdge0->verts[0]->x + tEdge0->verts[1]->x)/2 , 
	  (tEdge0->verts[0]->y + tEdge0->verts[1]->y)/2 ,
	  (tEdge0->verts[0]->z + tEdge0->verts[1]->z)/2);
	  
	  Vertex *v1n = new Vertex((tEdge1->verts[0]->x + tEdge1->verts[1]->x)/2 , 
	  (tEdge1->verts[0]->y + tEdge1->verts[1]->y)/2 ,
	  (tEdge1->verts[0]->z + tEdge1->verts[1]->z)/2);
	  
	  v0n->ntris = 0;
	  v1n->ntris = 0;

	  for(int zz=0; zz<vCount ; zz++)
	  {
	  Vertex *currentV = vlist1[zz];
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
	  if(NULL == v0)
	  {
	  vlist1[vCount++] = v0n;
	  v0 = v0n;
	  nverts++;
	  }
	  if(NULL == v1)
	  {
	  vlist1[vCount++] = v1n;
	  v1 = v1n;
	  nverts++;
	  }
	  if(NULL == v3)
	  {
	  vlist1[vCount++] = v3n;
	  v3 = v3n;
	  nverts++;
	  }
	  if(NULL == v4)
	  {
	  vlist1[vCount++] = v4n;
	  v4 = v4n;
	  nverts++;
	  }
	  if(NULL == v5)
	  {
	  vlist1[vCount++] = v5n;
	  v5 = v5n;
	  nverts++;
	  }
	  }


	  this->tlist[tCount] = new Triangle;
	  this->tlist[tCount]->nverts = 3;
	  this->tlist[tCount]->verts[0] = v0;
	  this->tlist[tCount]->verts[1] = v4;
	  this->tlist[tCount]->verts[2] = v1;
	  tCount++;

	  this->tlist[tCount] = new Triangle;
	  this->tlist[tCount]->nverts = 3;
	  this->tlist[tCount]->verts[0] = v0;
	  this->tlist[tCount]->verts[1] = v1;
	  this->tlist[tCount]->verts[2] = v3;
	  tCount++;

	  this->tlist[tCount] = new Triangle;
	  this->tlist[tCount]->nverts = 3;
	  this->tlist[tCount]->verts[0] = v3;
	  this->tlist[tCount]->verts[1] = v1;
	  this->tlist[tCount]->verts[2] = v5;
	  tCount++;
	  }

	  if(3 == smallerEdgeCount)
	  {
	  // Subdivide the current Triangle into four part
	  Vertex *v3n = new Vertex(v0nnn->x , v0nnn->y , v0nnn->z);
	  Vertex *v4n = new Vertex(v1nnn->x , v1nnn->y , v1nnn->z);
	  Vertex *v5n = new Vertex(v2nnn->x , v2nnn->y , v2nnn->z);
	  // Subdivide the current Triangle into two part
	  if(tEdge0 != NULL && tEdge1 != NULL && tEdge2 != NULL)
	  {
	  Vertex *v0n = new Vertex((tEdge0->verts[0]->x + tEdge0->verts[1]->x)/2 , 
	  (tEdge0->verts[0]->y + tEdge0->verts[1]->y)/2 ,
	  (tEdge0->verts[0]->z + tEdge0->verts[1]->z)/2);
	  
	  Vertex *v1n = new Vertex((tEdge1->verts[0]->x + tEdge1->verts[1]->x)/2 , 
	  (tEdge1->verts[0]->y + tEdge1->verts[1]->y)/2 ,
	  (tEdge1->verts[0]->z + tEdge1->verts[1]->z)/2);
	  
	  Vertex *v2n = new Vertex((tEdge2->verts[0]->x + tEdge2->verts[1]->x)/2 , 
	  (tEdge2->verts[0]->y + tEdge2->verts[1]->y)/2 ,
	  (tEdge2->verts[0]->z + tEdge2->verts[1]->z)/2);
	  
	  v0n->ntris = 0;
	  v1n->ntris = 0;
	  v2n->ntris = 0;

	  for(int zz=0; zz<vCount ; zz++)
	  {
	  Vertex *currentV = vlist1[zz];
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
	  if(NULL == v0)
	  {
	  vlist1[vCount++] = v0n;
	  v0 = v0n;
	  nverts++;
	  }
	  if(NULL == v1)
	  {
	  vlist1[vCount++] = v1n;
	  v1 = v1n;
	  nverts++;
	  }
	  if(NULL == v2)
	  {
	  vlist1[vCount++] = v2n;
	  v2 = v2n;
	  nverts++;
	  }
	  if(NULL == v3)
	  {
	  vlist1[vCount++] = v3n;
	  v3 = v3n;
	  nverts++;
	  }
	  if(NULL == v4)
	  {
	  vlist1[vCount++] = v4n;
	  v4 = v4n;
	  nverts++;
	  }
	  if(NULL == v5)
	  {
	  vlist1[vCount++] = v5n;
	  v5 = v5n;
	  nverts++;
	  }
	  }


	  this->tlist[tCount] = new Triangle;
	  this->tlist[tCount]->nverts = 3;
	  this->tlist[tCount]->verts[0] = v3;
	  this->tlist[tCount]->verts[1] = v0;
	  this->tlist[tCount]->verts[2] = v2;
	  tCount++;

	  this->tlist[tCount] = new Triangle;
	  this->tlist[tCount]->nverts = 3;
	  this->tlist[tCount]->verts[0] = v0;
	  this->tlist[tCount]->verts[1] = v4;
	  this->tlist[tCount]->verts[2] = v1;
	  tCount++;

	  this->tlist[tCount] = new Triangle;
	  this->tlist[tCount]->nverts = 3;
	  this->tlist[tCount]->verts[0] = v0;
	  this->tlist[tCount]->verts[1] = v1;
	  this->tlist[tCount]->verts[2] = v2;
	  tCount++;
	  	  
	  this->tlist[tCount] = new Triangle;
	  this->tlist[tCount]->nverts = 3;
	  this->tlist[tCount]->verts[0] = v1;
	  this->tlist[tCount]->verts[1] = v5;
	  this->tlist[tCount]->verts[2] = v2;
	  tCount++;
	  }
	  
	  if(0 == smallerEdgeCount)
	  {
	  Vertex *v3n = new Vertex(v0nnn->x , v0nnn->y , v0nnn->z);
	      Vertex *v4n = new Vertex(v1nnn->x , v1nnn->y , v1nnn->z);
	      Vertex *v5n = new Vertex(v2nnn->x , v2nnn->y , v2nnn->z);

	  for(int zz=0; zz<vCount ; zz++)
	  {
	  Vertex *currentV = vlist1[zz];
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
	  if(NULL == v3)
	  {
	  vlist1[vCount++] = v3n;
	  v3 = v3n;
	  nverts++;
	  }
	  if(NULL == v4)
	  {
	  vlist1[vCount++] = v4n;
	  v4 = v4n;
	  nverts++;
	  }
	  if(NULL == v5)
	  {
	  vlist1[vCount++] = v5n;
	  v5 = v5n;
	  nverts++;
	  }
	  // Just copy the triangle from the previous model to the current one
	  this->tlist[tCount] = new Triangle;
	  this->tlist[tCount]->nverts = 3;
	  this->tlist[tCount]->verts[0] = v3;//CurrentTriangle->verts[0];
	  this->tlist[tCount]->verts[1] = v4;//CurrentTriangle->verts[1];
	  this->tlist[tCount]->verts[2] = v5;//CurrentTriangle->verts[2];
	  tCount++;
	  //tlist[tCount++] = CurrentTriangle;
	  }
  }

  //Vertex** vlist1 = new Vertex *[4*poly->nverts];
  vlist = new Vertex *[nverts];
  for(i = 0; i<nverts; i++)
  {
	  vlist[i] = vlist1[i];
  }
}

Polyhedron::Polyhedron(FILE *, int Regular)
{
	 int i,j;
  int elem_count;
  char *elem_name;

  // here we use the previous polyhedron (polyR) to create a new polyhedron
  this->tlist = new Triangle *[4*poly->ntris];
  ntris = max_tris = 4*poly->ntris;

  Vertex** vlist1 = new Vertex *[4*poly->nverts];
  
  nverts = 0;
  int vCount = 0;
  for(int jj = 0; jj<poly->ntris; jj++)
  {
	  Triangle *tri = poly->tlist[jj];
      Vertex *v0nn = tri->verts[0];
      Vertex *v1nn = tri->verts[1];
      Vertex *v2nn = tri->verts[2];

	  // instead of doing a new we now use the already available vertices to update to the vlist 
	  Vertex *v3n = new Vertex((v0nn->x + v1nn->x)/2 , (v0nn->y + v1nn->y)/2 ,(v0nn->z + v1nn->z)/2);
	  Vertex *v4n = new Vertex((v1nn->x + v2nn->x)/2 , (v1nn->y + v2nn->y)/2 ,(v1nn->z + v2nn->z)/2);
	  Vertex *v5n = new Vertex((v2nn->x + v0nn->x)/2 , (v2nn->y + v0nn->y)/2 ,(v2nn->z + v0nn->z)/2);
	  v3n->ntris = 0;
	  v4n->ntris = 0;
	  v5n->ntris = 0;
	  
	  Vertex *v0n = new Vertex(v0nn->x , v0nn->y , v0nn->z);
	  Vertex *v1n = new Vertex(v1nn->x , v1nn->y , v1nn->z);
	  Vertex *v2n = new Vertex(v2nn->x , v2nn->y , v2nn->z);

	  v0n->ntris = 0;
	  v1n->ntris = 0;
	  v2n->ntris = 0;
  
	  Vertex *v3 = NULL;
	  Vertex *v4 = NULL;
	  Vertex *v5 = NULL;

	  Vertex *v0 = NULL;
      Vertex *v1 = NULL;
      Vertex *v2 = NULL;

      for(int zz=0; zz<vCount ; zz++)
	  {
	  Vertex *currentV = vlist1[zz];
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

	  if(NULL == v3)
	  {
       vlist1[vCount++] = v3n;	
	   nverts++;	   
	  }
	   
     if(NULL == v4)
	  {
	   vlist1[vCount++] = v4n;
	   nverts++;
	  }
	  
	 if(NULL == v5)
	  {
       vlist1[vCount++] = v5n;
	   nverts++;
	  }
  
	  if(NULL == v0)
	  {
	   vlist1[vCount++] = v0n;
	   nverts++;
	  }
	  
	  if(NULL == v1)
	  {
	   vlist1[vCount++] = v1n;
	   nverts++;
	  }
	  
	  if(NULL == v2)
	  {
	   vlist1[vCount++] = v2n;
	   nverts++;
	  }
  }

  this->vlist = new Vertex *[nverts];
  
  vCount = 0;
  for(int jj = 0; jj<poly->ntris; jj++)
  {
	  Triangle *tri = poly->tlist[jj];
      Vertex *v0nn = tri->verts[0];
      Vertex *v1nn = tri->verts[1];
      Vertex *v2nn = tri->verts[2];

	  /*Vertex *v3 = new Vertex((v0->x + v1->x)/2 , (v0->y + v1->y)/2 ,(v0->z + v1->z)/2);
	  Vertex *v4 = new Vertex((v1->x + v2->x)/2 , (v1->y + v2->y)/2 ,(v1->z + v2->z)/2);
	  Vertex *v5 = new Vertex((v2->x + v0->x)/2 , (v2->y + v0->y)/2 ,(v2->z + v0->z)/2);
      */

	  // instead of doing a new we now use the already available vertices to update to the vlist 
	  Vertex *v3n = new Vertex((v0nn->x + v1nn->x)/2 , (v0nn->y + v1nn->y)/2 ,(v0nn->z + v1nn->z)/2);
	  Vertex *v4n = new Vertex((v1nn->x + v2nn->x)/2 , (v1nn->y + v2nn->y)/2 ,(v1nn->z + v2nn->z)/2);
	  Vertex *v5n = new Vertex((v2nn->x + v0nn->x)/2 , (v2nn->y + v0nn->y)/2 ,(v2nn->z + v0nn->z)/2);
	  v3n->ntris = 0;
	  v4n->ntris = 0;
	  v5n->ntris = 0;
	  
	  //v3n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	  //v4n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	  //v5n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	  //v3n->max_tris = 0;
	  //v4n->max_tris = 0;
	  //v5n->max_tris = 0;
	  
	  Vertex *v0n = new Vertex(v0nn->x , v0nn->y , v0nn->z);
	  Vertex *v1n = new Vertex(v1nn->x , v1nn->y , v1nn->z);
	  Vertex *v2n = new Vertex(v2nn->x , v2nn->y , v2nn->z);

	  v0n->ntris = 0;
	  v1n->ntris = 0;
	  v2n->ntris = 0;

	  //v0n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	  //v1n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	  //v2n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);

	  //v0n->max_tris = 0;
	  //v1n->max_tris = 0;
	  //v2n->max_tris = 0;
	  
	  
	  Vertex *v3 = NULL;
	  Vertex *v4 = NULL;
	  Vertex *v5 = NULL;

	  Vertex *v0 = NULL;
      Vertex *v1 = NULL;
      Vertex *v2 = NULL;

      for(int zz=0; zz<vCount ; zz++)
	  {
	  Vertex *currentV = vlist[zz];
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

	  if(NULL == v3)
	  {
       //v3n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	    vlist[vCount++] = v3n;
	//vlist[vCount++] = new Vertex((v0nn->x + v1nn->x)/2 , (v0nn->y + v1nn->y)/2 ,(v0nn->z + v1nn->z)/2);
	   v3 = v3n;
	   //nverts++;
	   /*v3->x = (v0->x + v1->x)/2;
	   v3->y = (v0->y + v1->y)/2;
	   v3->z = (v0->z + v1->z)/2;*/
	  }
	  else
	  {
	  delete v3n;
	  }
	   
     if(NULL == v4)
	  {
	 // v4n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	   vlist[vCount++] = v4n;
	   //vlist[vCount++] = new Vertex((v1nn->x + v2nn->x)/2 , (v1nn->y + v2nn->y)/2 ,(v1nn->z + v2nn->z)/2);
	   v4 = v4n;
	   //nverts++;
	   /*v4 = vlist[vCount++];
	   v4->x = (v1->x + v2->x)/2;
	   v4->y = (v1->y + v2->y)/2;
	   v4->z = (v1->z + v2->z)/2;*/
	  }	  
	 else
	  {
	  delete v4n;
	  }

	 if(NULL == v5)
	  {
	 // v5n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
       vlist[vCount++] = v5n;
	    //vlist[vCount++] = new Vertex((v2nn->x + v0nn->x)/2 , (v2nn->y + v0nn->y)/2 ,(v2nn->z + v0nn->z)/2);
	   v5 = v5n;
	   //nverts++;
	   /*v5 = vlist[vCount++];
	   v5->x = (v2->x + v0->x)/2;
	   v5->y = (v2->y + v0->y)/2;
	   v5->z = (v2->z + v0->z)/2;*/
	  }
      else
	  {
	  delete v5n;
	  }
	  if(NULL == v0)
	  {
	   //v0n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	   vlist[vCount++] = v0n;
	   //vlist[vCount++] = new Vertex(v0nn->x , v0nn->y , v0nn->z);
	   v0 = v0n;
	   //nverts++;
	  }
	  else
	  {
	  delete v0n;
	  }

	  if(NULL == v1)
	  {
	 //v1n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	   vlist[vCount++] = v1n;
	   //vlist[vCount++] = new Vertex(v1nn->x , v1nn->y , v1nn->z);
	   v1 = v1n;
	   //nverts++;
	  }
	  else
	  {
	  delete v1n;
	  }
	  if(NULL == v2)
	  {
	//v2n->tris = (Triangle **)malloc (sizeof (Triangle *) * 3);
	   vlist[vCount++] = v2n;
	   //vlist[vCount++] = new Vertex(v2nn->x , v2nn->y , v2nn->z);
	   v2 = v2n;
	   //nverts++;
	  }
	  else
	  {
	  delete v2n;
	  }
	  this->tlist[4*jj] = new Triangle;
      this->tlist[4*jj]->nverts = 3;
      this->tlist[4*jj]->verts[0] = v0;
      this->tlist[4*jj]->verts[1] = v3;
      this->tlist[4*jj]->verts[2] = v5;

	  this->tlist[4*jj +1] = new Triangle;
      this->tlist[4*jj+1]->nverts = 3;
      this->tlist[4*jj+1]->verts[0] = v3;
      this->tlist[4*jj+1]->verts[1] = v1;
      this->tlist[4*jj+1]->verts[2] = v4;

	  this->tlist[4*jj+2] = new Triangle;
      this->tlist[4*jj+2]->nverts = 3;
      this->tlist[4*jj+2]->verts[0] = v4;
      this->tlist[4*jj+2]->verts[1] = v2;
      this->tlist[4*jj+2]->verts[2] = v5;

	  this->tlist[4*jj+3] = new Triangle;
      this->tlist[4*jj+3]->nverts = 3;
      this->tlist[4*jj+3]->verts[0] = v3;
      this->tlist[4*jj+3]->verts[1] = v4;
      this->tlist[4*jj+3]->verts[2] = v5;
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
	char mm[1024];
	sprintf(mm, "modified by learnply");
//  append_comment_ply (ply, "modified by simvizply %f");
	  append_comment_ply (ply, mm);
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
    this->vlist[i]->max_tris = 0;

  /* first just count all the face pointers needed for each vertex */

  for (i = 0; i < ntris; i++) {
    f = this->tlist[i];
    for (j = 0; j < f->nverts; j++)
      f->verts[j]->max_tris++;
  }

  /* allocate memory for face pointers of vertices */

  for (i = 0; i < nverts; i++) {
    this->vlist[i]->tris = (Triangle **)
	      malloc (sizeof (Triangle *) * vlist[i]->max_tris);
    this->vlist[i]->ntris = 0;
	 
  }

  /* now actually create the face pointers */

  for (i = 0; i < this->ntris; i++) {
    f = this->tlist[i];
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
//	if (i %1000 == 0)
//	fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
    order_vertex_to_tri_ptrs(vlist[i]);
	
	}
  /* index the edges */

  for (i = 0; i < nedges; i++){
//	if (i %1000 == 0)
//	fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
    elist[i]->index = i;
	}

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
  	 icVector3 norm1(1,0,0);
	 icVector3 norm2(0,1,0);
	 icVector3 norm3(0,0,1);
	 icVector3 bboxpoint1 ;
	 icVector3 bboxpoint2 ;
	 icVector3 bboxpoint3 ;
	 icVector3 bboxpoint4 ;
	 icVector3 bboxpoint5 ;
	 icVector3 bboxpoint6 ;
	 icVector3 bboxpoint7 ;
	 icVector3 bboxpoint8 ;
	 
     bboxpoint1 = center + norm1*fabs(max.entry[0]) + norm2*fabs(max.entry[1]) + norm3*fabs(max.entry[2]);
	 bboxpoint2 = center + norm1*fabs(max.entry[0]) + norm2*fabs(max.entry[1]) - norm3*fabs(min.entry[2]);
	 bboxpoint3 = center + norm1*fabs(max.entry[0]) - norm2*fabs(min.entry[1]) - norm3*fabs(min.entry[2]);
	 bboxpoint4 = center + norm1*fabs(max.entry[0]) - norm2*fabs(min.entry[1]) + norm3*fabs(max.entry[2]);
	 bboxpoint5 = center - norm1*fabs(min.entry[0]) + norm2*fabs(max.entry[1]) + norm3*fabs(max.entry[2]);
	 bboxpoint6 = center - norm1*fabs(min.entry[0]) + norm2*fabs(max.entry[1]) - norm3*fabs(min.entry[2]);
	 bboxpoint7 = center - norm1*fabs(min.entry[0]) - norm2*fabs(min.entry[1]) - norm3*fabs(min.entry[2]);
	 bboxpoint8 = center - norm1*fabs(min.entry[0]) - norm2*fabs(min.entry[1]) + norm3*fabs(max.entry[2]);

	 //icVector3 bboxpoint1 = median + norm1*normval1 + norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint2 = median + norm1*normval1 + norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint3 = median + norm1*normval1 - norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint4 = median + norm1*normval1 - norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint5 = median - norm1*normval1 + norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint6 = median - norm1*normval1 + norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint7 = median - norm1*normval1 - norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint8 = median - norm1*normval1 - norm2*normval2 + norm3*normval3;

	 boundingBox1[0][0] = bboxpoint1.x;
	 boundingBox1[1][0] = bboxpoint1.y;
	 boundingBox1[2][0] = bboxpoint1.z;

	 boundingBox1[0][1] = bboxpoint2.x;
	 boundingBox1[1][1] = bboxpoint2.y;
	 boundingBox1[2][1] = bboxpoint2.z;

	 boundingBox1[0][2] = bboxpoint3.x;
	 boundingBox1[1][2] = bboxpoint3.y;
	 boundingBox1[2][2] = bboxpoint3.z;

	 boundingBox1[0][3] = bboxpoint4.x;
	 boundingBox1[1][3] = bboxpoint4.y;
	 boundingBox1[2][3] = bboxpoint4.z;

	 boundingBox1[0][4] = bboxpoint5.x;
	 boundingBox1[1][4] = bboxpoint5.y;
	 boundingBox1[2][4] = bboxpoint5.z;

	 boundingBox1[0][5] = bboxpoint6.x;
	 boundingBox1[1][5] = bboxpoint6.y;
	 boundingBox1[2][5] = bboxpoint6.z;

	 boundingBox1[0][6] = bboxpoint7.x;
	 boundingBox1[1][6] = bboxpoint7.y;
	 boundingBox1[2][6] = bboxpoint7.z;

	 boundingBox1[0][7] = bboxpoint8.x;
	 boundingBox1[1][7] = bboxpoint8.y;
	 boundingBox1[2][7] = bboxpoint8.z;
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

void Polyhedron::CreateCornerList()
{
	int i=0;
	nCorners = 3*ntris;
  clist = new Corner*[nCorners];
  int cListCount = 0;
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
	clist[cListCount] = new Corner;

	clist[cListCount]->tri = tri;
	clist[cListCount]->vert = v0;
	    corner0 = clist[cListCount];
	    
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
	
	
	 if(0 == pointInEdge0)
	 {
	 sEdge = edge0;
	 }

	 
	 if(0 == pointInEdge1)
	 {
	 sEdge = edge1;
	 }
	 
	 if(0 == pointInEdge2)
	 {
	 sEdge = edge2;
	 }

	clist[cListCount]->edge = sEdge;

	cListCount++;

	// Second corner

	pointInEdge0 = 0;
	pointInEdge1 = 0;
	pointInEdge2 = 0;
	clist[cListCount] = new Corner;

	clist[cListCount]->tri = tri;
	clist[cListCount]->vert = v1;
	    corner1 = clist[cListCount];
	    
	
	//if( (edge0->verts[0]->x != v1->x && edge0->verts[0]->y != v1->y && edge0->verts[0]->z != v1->z) &&
	//	(edge0->verts[1]->x != v1->x && edge0->verts[1]->y != v1->y && edge0->verts[1]->z != v1->z) )
	//{
	//	sEdge = edge0;
	//}

	//if( (edge1->verts[0]->x != v1->x && edge1->verts[0]->y != v1->y && edge1->verts[0]->z != v1->z) &&
	//	(edge1->verts[1]->x != v1->x && edge1->verts[1]->y != v1->y && edge1->verts[1]->z != v1->z) )
	//{
	//	sEdge = edge1;
	//}

	//if( (edge2->verts[0]->x != v1->x && edge2->verts[0]->y != v1->y && edge2->verts[0]->z != v1->z) &&
	//	(edge2->verts[1]->x != v1->x && edge2->verts[1]->y != v1->y && edge2->verts[1]->z != v1->z) )
	//{
	//	sEdge = edge2;
	//}

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

	 if(0 == pointInEdge0)
	 {
	 sEdge = edge0;
	 }

	 
	 if(0 == pointInEdge1)
	 {
	 sEdge = edge1;
	 }
	 
	 if(0 == pointInEdge2)
	 {
	 sEdge = edge2;
	 }

	clist[cListCount]->edge = sEdge;
	clist[cListCount]->previous = corner0;
	cListCount++;

	// Third corner
	pointInEdge0 = 0;
	pointInEdge1 = 0;
	pointInEdge2 = 0;
	clist[cListCount] = new Corner;

	clist[cListCount]->tri = tri;
	clist[cListCount]->vert = v2;
	    corner2 = clist[cListCount];
	    
	
	edge0 = tri->edges[0];
	edge1 = tri->edges[1];
	edge2 = tri->edges[2];
	
	//if( (edge0->verts[0]->x != v2->x && edge0->verts[0]->y != v2->y && edge0->verts[0]->z != v2->z) &&
	//	(edge0->verts[1]->x != v2->x && edge0->verts[1]->y != v2->y && edge0->verts[1]->z != v2->z) )
	//{
	//	sEdge = edge0;
	//}

	//if( (edge1->verts[0]->x != v2->x && edge1->verts[0]->y != v2->y && edge1->verts[0]->z != v2->z) &&
	//	(edge1->verts[1]->x != v2->x && edge1->verts[1]->y != v2->y && edge1->verts[1]->z != v2->z) )
	//{
	//	sEdge = edge1;
	//}

	//if( (edge2->verts[0]->x != v2->x && edge2->verts[0]->y != v2->y && edge2->verts[0]->z != v2->z) &&
	//	(edge2->verts[1]->x != v2->x && edge2->verts[1]->y != v2->y && edge2->verts[1]->z != v2->z) )
	//{
	//	sEdge = edge2;
	//}

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

	     if(0 == pointInEdge0)
	 {
	 sEdge = edge0;
	 }

	 
	 if(0 == pointInEdge1)
	 {
	 sEdge = edge1;
	 }
	 
	 if(0 == pointInEdge2)
	 {
	 sEdge = edge2;
	 }

	clist[cListCount]->edge = sEdge;
	clist[cListCount]->previous = corner1;
	clist[cListCount]->next = corner0;

	corner1->next = corner2;

	corner0->next = corner1;
	corner0->previous = corner2;
	cListCount++;
   
   }

   /// Here we will update the opposite of each of the corners
   for(i = 0; i<nCorners; i++)
   {
	   Corner *CurrentCorner = clist[i];

	   for(int j = 0; j<nCorners; j++)
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
	display();
	}
	}
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
	int kk= -5;
	int yy= 3;
	int zz = kk/yy;
	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHT1);
	float rColor = 0.0f;
	float gColor = 0.0f;
	float bColor = 0.0f;

	int RCOLOR = 0;
	int GCOLOR = 0;
	int BCOLOR = 0;

	Vertex *v1;
	Vertex *v2;
	Vertex *v3;

	int TotoalDeficitCount = 0;
	double TotalAngleDeficit = 0.0;
	// icVector3 norm1(1,0,0);
	//icVector3 norm2(0,1,0);
	//icVector3 norm3(0,0,1);
	// icVector3 bboxpoint1 ;
	// icVector3 bboxpoint2 ;
	// icVector3 bboxpoint3 ;
	// icVector3 bboxpoint4 ;
	// icVector3 bboxpoint5 ;
	// icVector3 bboxpoint6 ;
	// icVector3 bboxpoint7 ;
	// icVector3 bboxpoint8 ;
	int nTriangles = this_poly->ntris;
	float cubeRootTriangle = pow(nTriangles,1.f/3.f);
	//float cIncrement = 1.f/cubeRootTriangle; 
	float cIncrement = 1.f/255;

	//for (i=0; i<this_poly->ntris; i++) {
	//	if (mode == GL_SELECT)
	//	glLoadName(i+1);

	//	Triangle *temp_t=this_poly->tlist[i];

	switch (display_mode) {
	case 0: //Default
	//this_poly->calculate_bounding_box_with_moments();
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

	//case 1: // Based on Unique ID
	//	glDisable(GL_LIGHTING);
	//	glDisable(GL_LIGHT0);
	//	glDisable(GL_LIGHT1);

	//	for (i=0; i<this_poly->ntris; i++) {
	//	if (mode == GL_SELECT)
	//	glLoadName(i+1);

	//	Triangle *temp_t=this_poly->tlist[i];
	//	if (i == this_poly->seed) {
	//	mat_diffuse[0] = 0.0;
	//	mat_diffuse[1] = 0.0;
	//	mat_diffuse[2] = 1.0;
	//	mat_diffuse[3] = 1.0;
	//	} else {
	//	mat_diffuse[0] = 1.0;
	//	mat_diffuse[1] = 1.0;
	//	mat_diffuse[2] = 0.0;
	//	mat_diffuse[3] = 1.0;
	//	}
	//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	//	glBegin(GL_POLYGON);
	//	for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	if (i==this_poly->seed)
	//	glColor3f(0.0, 0.0, 1.0);
	//	else
	//	{
	//	//XD8
	//	//(2.a) Based on ID
	//	//Here we need to provide a unique color to each of the polygon
	//	//Now, number of triangle are:
	//	//The value of i implies the triangle ID
	//	//For that we divide each R G and B values into 255 i.e each new color would add 1/255 to give a unique color	
	//	//   RCOLOR = rColor;
	//	//      GCOLOR = gColor;
	//	//      BCOLOR = bColor;
	//	//if(0==j)
	//	//{
	//	//    glColor3f(rColor,0,0);
	//	//}
	//	//else if(1==j)
	//	//{
	//	//	glColor3f(0,gColor,0);
	//	//}
	//	//else if(2==j)
	//	//{
	//	//	glColor3f(0,0,bColor);
	//	//}
	//	RCOLOR = (i/256)/256;
	//	GCOLOR = (i/256)%256;
	//	BCOLOR = i%256;

	//	rColor = RCOLOR*cIncrement;
	//	gColor = GCOLOR*cIncrement;
	//	bColor = BCOLOR*cIncrement;

	//	glColor3f(rColor,gColor,bColor);
	//	//if(gColor < 1.f)
	//	//{
	//	//	if(bColor < 1.f)
	//	//	{
	//	//	bColor = bColor + cIncrement;
	//	//	}
	//	//	else
	//	//	{
	//	//	bColor = 0.f;
	//	//	gColor = gColor + cIncrement;
	//	//	}
	//	//}
	//	//else
	//	//{
	//	//	gColor = 0.f;
	//	//	rColor = rColor + cIncrement;
	//	//}

	//	// end based on ID
	//	}	
	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//	}
	//	glEnd();
	//	}
	//	break;

	//case 2: // Based on Normal
	//	glDisable(GL_LIGHTING);
	//	glDisable(GL_LIGHT0);
	//	glDisable(GL_LIGHT1);
	//	for (i=0; i<this_poly->ntris; i++) {
	//	if (mode == GL_SELECT)
	//	glLoadName(i+1);

	//	Triangle *temp_t=this_poly->tlist[i];
	//	if (i == this_poly->seed) {
	//	mat_diffuse[0] = 0.0;
	//	mat_diffuse[1] = 0.0;
	//	mat_diffuse[2] = 1.0;
	//	mat_diffuse[3] = 1.0;
	//	} else {
	//	mat_diffuse[0] = 1.0;
	//	mat_diffuse[1] = 1.0;
	//	mat_diffuse[2] = 0.0;
	//	mat_diffuse[3] = 1.0;
	//	}
	//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	//	glBegin(GL_POLYGON);
	//	for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	if (i==this_poly->seed)
	//	glColor3f(0.0, 0.0, 1.0);
	//	else
	//	{
	//	// XD8

	//	// (2.b) Based on Normal	
	//	float rColor = temp_t->normal.x;
	//	float gColor = temp_t->normal.y;
	//	float bColor = temp_t->normal.z;
	//	if(rColor < 0.f)
	//	{
	//	rColor = -rColor;
	//	}
	//	if(gColor < 0.f)
	//	{
	//	gColor = -gColor;
	//	}
	//	if(bColor < 0.f)
	//	{
	//	bColor = -bColor;
	//	}

	//	glColor3f(rColor,gColor,bColor);
	//	//End Based on Normal
	//	}
	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//	}
	//	glEnd();
	//	}
	//	break;

	//case 3: // Based Color to vertex (3D Checkboard)
	//	glDisable(GL_LIGHTING);
	//	glDisable(GL_LIGHT0);
	//	glDisable(GL_LIGHT1);
	//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	//	
	//	for (i=0; i<this_poly->ntris; i++) {
	//	if (mode == GL_SELECT)
	//	glLoadName(i+1);

	//	Triangle *temp_t=this_poly->tlist[i];
	//	glBegin(GL_POLYGON);
	//	for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	if (i==this_poly->seed)
	//	glColor3f(0.0, 0.0, 1.0);
	//	else
	//	{
	//	// (2.c) Color to vertex
	//	float vX = temp_v->x;
	//	float vY = temp_v->y;
	//	float vZ = temp_v->z;

	//	int rInt =floor( vX/L);
	//	int gInt =floor( vY/L);
	//	int bInt = floor(vZ/L);

	//	//int rFloat = vX/L;
	//	//int gFloat = vY/L;
	//	//int bFloat = vZ/L;

	//	//if((rInt-rFloat)>=0.5f)
	//	//{
	//	//	if(rInt>0)
	//	//	rInt = rInt+1;
	//	//	else
	//	//	rInt = rInt -1;
	//	//}

	//	//if((gInt-gFloat)>=0.5f)
	//	//{
	//	//	if(gInt>0)
	//	//	gInt = gInt+1;
	//	//	else
	//	//	gInt = gInt -1;
	//	//}

	//	//if((bInt-bFloat)>=0.5f)
	//	//{
	//	//	if(bInt>0)
	//	//	bInt = bInt+1;
	//	//	else
	//	//	bInt = bInt -1;
	//	//}

	//	if(0 == rInt%2)
	//	{
	//	rColor = 1.f;
	//	}
	//	else
	//	{
	//	rColor = 0.f;
	//	}

	//	if(0 == gInt%2)
	//	{
	//	gColor = 1.f;
	//	}
	//	else
	//	{
	//	gColor = 0.f;
	//	}

	//	if(0 == bInt%2)
	//	{
	//	bColor = 1.f;
	//	}
	//	else
	//	{
	//	bColor = 0.f;
	//	}
	//	glColor3f(rColor, gColor, bColor);
	//	}
	//	// End Color to vertex
	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//	}
	//	glEnd();

	//	// (2.c) Color to vertex
	//	// The color has to be based on a real number as enterd by user
	//	// (2.d)
	//	// The problem is that when the user enters a large number every coordinate/L is 0 (which is even)
	//	// Thus all R , G and B for all vertices on a face are vertices are all 1
	//	// This maps the triangle to a black color
	//	// As we keep on entering a higher value, the entire model keep on looking all black
	//	// Thus we see the entire model as black
	//	// A probable solution can be to make the.... 
	//	// (2.e)
	//	// for icosahedron in the 3d checkboard color mode, the structure does not precisely looks sharp at the corner
	//	// we need to make the corners look sharper to get the sense a real icosahedron 

	//	/*	glPointSize(5);   
	//	glBegin(GL_POINTS); 	

	//	for (j=0; j<3; j++) 
	//	{
	//	Vertex *temp_v = temp_t->verts[j];
	//	float vX = temp_v->x;
	//	float vY = temp_v->y;
	//	float vZ = temp_v->z;

	//	int rInt = vX/L;
	//	int gInt = vY/L;
	//	int bInt = vZ/L;

	//	if(0 == rInt%2)
	//	{
	//	rColor = 1.f;
	//	}
	//	else
	//	{
	//	rColor = 0.f;
	//	}

	//	if(0 == gInt%2)
	//	{
	//	gColor = 1.f;
	//	}
	//	else
	//	{
	//	gColor = 0.f;
	//	}

	//	if(0 == bInt%2)
	//	{
	//	bColor = 1.f;
	//	}
	//	else
	//	{
	//	bColor = 0.f;
	//	}
	//	if(1.f == rColor && 1.f == gColor && 1.f == bColor)
	//	{
	//	rColor = 0.f;
	//	gColor = 0.f;
	//	bColor = 0.f;
	//	}
	//	glColor3f(rColor, gColor, bColor);
	//	glVertex3f(temp_v->x, temp_v->y, temp_v->z);	
	//	}

	//	glEnd();      */   


	//	// (2.e)
	//	// Line drawn with vertices to sense the volume of the model like icosahedron
	//	//v1 = temp_t->verts[0];
	//	//v2 = temp_t->verts[1];
	//	//v3 = temp_t->verts[2];

	//	//glColor3f(0.0, 0.0, 0.0);
	//	//glPointSize(2); 

	//	//glBegin(GL_LINES);	
	//	//glVertex3f(v1->x, v1->y, v1->z);	
	//	//glVertex3f(v2->x, v2->y, v2->z);	
	//	//glEnd();

	//	//glBegin(GL_LINES);	
	//	//glVertex3f(v2->x, v2->y, v2->z);	
	//	//glVertex3f(v3->x, v3->y, v3->z);	
	//	//glEnd();

	//	//glBegin(GL_LINES);	
	//	//glVertex3f(v3->x, v3->y, v3->z);	
	//	//glVertex3f(v1->x, v1->y, v1->z);	
	//	//glEnd();

	//	// End Color to vertex
	//	}
	//	break;

	//case 4: // Axis Aglined Bounding Box

	//	mat_diffuse[0] = 1.0;
	//	mat_diffuse[1] = 1.0;
	//	mat_diffuse[2] = 0.0;
	//	mat_diffuse[3] = 1.0;
	//	//   bboxpoint1 = this_poly->center + norm1*(this_poly->radius/2) + norm2*(this_poly->radius/2) + norm3*(this_poly->radius/2);
	//	//bboxpoint2 = this_poly->center + norm1*(this_poly->radius/2) + norm2*(this_poly->radius/2) - norm3*(this_poly->radius/2);
	//	//bboxpoint3 = this_poly->center + norm1*(this_poly->radius/2) - norm2*(this_poly->radius/2) - norm3*(this_poly->radius/2);
	//	//bboxpoint4 = this_poly->center + norm1*(this_poly->radius/2) - norm2*(this_poly->radius/2) + norm3*(this_poly->radius/2);
	//	//bboxpoint5 = this_poly->center - norm1*(this_poly->radius/2) + norm2*(this_poly->radius/2) + norm3*(this_poly->radius/2);
	//	//bboxpoint6 = this_poly->center - norm1*(this_poly->radius/2) + norm2*(this_poly->radius/2) - norm3*(this_poly->radius/2);
	//	//bboxpoint7 = this_poly->center - norm1*(this_poly->radius/2) - norm2*(this_poly->radius/2) - norm3*(this_poly->radius/2);
	//	//bboxpoint8 = this_poly->center - norm1*(this_poly->radius/2) - norm2*(this_poly->radius/2) + norm3*(this_poly->radius/2);

	//	////icVector3 bboxpoint1 = median + norm1*normval1 + norm2*normval2 + norm3*normval3;
	//	////icVector3 bboxpoint2 = median + norm1*normval1 + norm2*normval2 - norm3*normval3;
	//	////icVector3 bboxpoint3 = median + norm1*normval1 - norm2*normval2 - norm3*normval3;
	//	////icVector3 bboxpoint4 = median + norm1*normval1 - norm2*normval2 + norm3*normval3;
	//	////icVector3 bboxpoint5 = median - norm1*normval1 + norm2*normval2 + norm3*normval3;
	//	////icVector3 bboxpoint6 = median - norm1*normval1 + norm2*normval2 - norm3*normval3;
	//	////icVector3 bboxpoint7 = median - norm1*normval1 - norm2*normval2 - norm3*normval3;
	//	////icVector3 bboxpoint8 = median - norm1*normval1 - norm2*normval2 + norm3*normval3;

	//	//boundingBox1[0][0] = bboxpoint1.x;
	//	//boundingBox1[1][0] = bboxpoint1.y;
	//	//boundingBox1[2][0] = bboxpoint1.z;

	//	//boundingBox1[0][1] = bboxpoint2.x;
	//	//boundingBox1[1][1] = bboxpoint2.y;
	//	//boundingBox1[2][1] = bboxpoint2.z;

	//	//boundingBox1[0][2] = bboxpoint3.x;
	//	//boundingBox1[1][2] = bboxpoint3.y;
	//	//boundingBox1[2][2] = bboxpoint3.z;

	//	//boundingBox1[0][3] = bboxpoint4.x;
	//	//boundingBox1[1][3] = bboxpoint4.y;
	//	//boundingBox1[2][3] = bboxpoint4.z;

	//	//boundingBox1[0][4] = bboxpoint5.x;
	//	//boundingBox1[1][4] = bboxpoint5.y;
	//	//boundingBox1[2][4] = bboxpoint5.z;

	//	//boundingBox1[0][5] = bboxpoint6.x;
	//	//boundingBox1[1][5] = bboxpoint6.y;
	//	//boundingBox1[2][5] = bboxpoint6.z;

	//	//boundingBox1[0][6] = bboxpoint7.x;
	//	//boundingBox1[1][6] = bboxpoint7.y;
	//	//boundingBox1[2][6] = bboxpoint7.z;

	//	//boundingBox1[0][7] = bboxpoint8.x;
	//	//boundingBox1[1][7] = bboxpoint8.y;
	//	//boundingBox1[2][7] = bboxpoint8.z;
	//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	//	glDisable(GL_LIGHTING);
	//	glDisable(GL_LIGHT0);
	//	glDisable(GL_LIGHT1);
	//	glColor3f(0.0, 0.0, 0.0);
	//	glLineWidth(5.0);


	//	glBegin(GL_LINES);	
	//	glVertex3d(boundingBox1[0][0],boundingBox1[1][0],boundingBox1[2][0]);
	//	glVertex3d(boundingBox1[0][1],boundingBox1[1][1],boundingBox1[2][1]);	

	//	glVertex3d(boundingBox1[0][1],boundingBox1[1][1],boundingBox1[2][1]);
	//	glVertex3d(boundingBox1[0][2],boundingBox1[1][2],boundingBox1[2][2]);	

	//	glVertex3d(boundingBox1[0][2],boundingBox1[1][2],boundingBox1[2][2]);
	//	glVertex3d(boundingBox1[0][3],boundingBox1[1][3],boundingBox1[2][3]);	

	//	glVertex3d(boundingBox1[0][3],boundingBox1[1][3],boundingBox1[2][3]);
	//	glVertex3d(boundingBox1[0][0],boundingBox1[1][0],boundingBox1[2][0]);	

	//	glVertex3d(boundingBox1[0][0],boundingBox1[1][0],boundingBox1[2][0]);
	//	glVertex3d(boundingBox1[0][4],boundingBox1[1][4],boundingBox1[2][4]);	

	//	glVertex3d(boundingBox1[0][1],boundingBox1[1][1],boundingBox1[2][1]);
	//	glVertex3d(boundingBox1[0][5],boundingBox1[1][5],boundingBox1[2][5]);	

	//	glVertex3d(boundingBox1[0][2],boundingBox1[1][2],boundingBox1[2][2]);
	//	glVertex3d(boundingBox1[0][6],boundingBox1[1][6],boundingBox1[2][6]);	

	//	glVertex3d(boundingBox1[0][3],boundingBox1[1][3],boundingBox1[2][3]);
	//	glVertex3d(boundingBox1[0][7],boundingBox1[1][7],boundingBox1[2][7]);	

	//	glVertex3d(boundingBox1[0][4],boundingBox1[1][4],boundingBox1[2][4]);
	//	glVertex3d(boundingBox1[0][7],boundingBox1[1][7],boundingBox1[2][7]);	

	//	glVertex3d(boundingBox1[0][7],boundingBox1[1][7],boundingBox1[2][7]);
	//	glVertex3d(boundingBox1[0][6],boundingBox1[1][6],boundingBox1[2][6]);	

	//	glVertex3d(boundingBox1[0][6],boundingBox1[1][6],boundingBox1[2][6]);
	//	glVertex3d(boundingBox1[0][5],boundingBox1[1][5],boundingBox1[2][5]);	

	//	glVertex3d(boundingBox1[0][5],boundingBox1[1][5],boundingBox1[2][5]);
	//	glVertex3d(boundingBox1[0][4],boundingBox1[1][4],boundingBox1[2][4]);

	//	glEnd();

	//	for (i=0; i<this_poly->ntris; i++) {
	//	if (mode == GL_SELECT)
	//	glLoadName(i+1);

	//	Triangle *temp_t=this_poly->tlist[i];

	//	glBegin(GL_POLYGON);
	//	for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	if (i==this_poly->seed)
	//	glColor3f(0.0, 0.0, 1.0);
	//	else
	//	{
	//	// (2.b) Based on Normal	
	//	float rColor = temp_t->normal.x;
	//	float gColor = temp_t->normal.y;
	//	float bColor = temp_t->normal.z;
	//	if(rColor < 0.f)
	//	{
	//	rColor = -rColor;
	//	}
	//	if(gColor < 0.f)
	//	{
	//	gColor = -gColor;
	//	}
	//	if(bColor < 0.f)
	//	{
	//	bColor = -bColor;
	//	}

	//	glColor3f(rColor,gColor,bColor);
	//	//End Based on Normal
	//	//glColor3f(1.0, 1.0, 0.0);
	//	}

	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//	}
	//	glEnd();
	//	}
	//	break;

	//case 5: //Bounding box based on center of gravity
	//	//this_poly->calculate_bounding_box_with_moments();

	//	mat_diffuse[0] = 1.0;
	//	mat_diffuse[1] = 1.0;
	//	mat_diffuse[2] = 0.0;
	//	mat_diffuse[3] = 1.0;

	//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	//	glDisable(GL_LIGHTING);
	//	glDisable(GL_LIGHT0);
	//	glDisable(GL_LIGHT1);
	//	glColor3f(0.0, 0.0, 0.0);
	//	glLineWidth(5.0);


	//	glBegin(GL_LINES);	

	//	// Here we give RED, GREEN and BLUE colors to Major, Middle and Minor eigenvetors
	//	if(1 == RGB[2])
	//	{
	//	glColor3f(1.0, 0.0, 0.0);
	//	}
	//	else if(2 == RGB[2])
	//	{
	//	glColor3f(0.0, 1.0, 0.0);
	//	}
	//	else if(3 == RGB[2])
	//	{
	//	glColor3f(0.0, 0.0, 1.0);
	//	}



	//	//glColor3f(0.0, 1.0, 0.0);
	//	glVertex3d(boundingBox[0][3],boundingBox[1][3],boundingBox[2][3]);
	//	glVertex3d(boundingBox[0][7],boundingBox[1][7],boundingBox[2][7]);	

	//	glVertex3d(boundingBox[0][0],boundingBox[1][0],boundingBox[2][0]);
	//	glVertex3d(boundingBox[0][4],boundingBox[1][4],boundingBox[2][4]);	

	//	glVertex3d(boundingBox[0][1],boundingBox[1][1],boundingBox[2][1]);
	//	glVertex3d(boundingBox[0][5],boundingBox[1][5],boundingBox[2][5]);

	//	glVertex3d(boundingBox[0][2],boundingBox[1][2],boundingBox[2][2]);
	//	glVertex3d(boundingBox[0][6],boundingBox[1][6],boundingBox[2][6]);	


	//	if(1 == RGB[1])
	//	{
	//	glColor3f(1.0, 0.0, 0.0);
	//	}
	//	else if(2 == RGB[1])
	//	{
	//	glColor3f(0.0, 1.0, 0.0);
	//	}
	//	else if(3 == RGB[1])
	//	{
	//	glColor3f(0.0, 0.0, 1.0);
	//	}

	//	//glColor3f(0.0, 0.0, 1.0);
	//	glVertex3d(boundingBox[0][4],boundingBox[1][4],boundingBox[2][4]);
	//	glVertex3d(boundingBox[0][7],boundingBox[1][7],boundingBox[2][7]);

	//	glVertex3d(boundingBox[0][3],boundingBox[1][3],boundingBox[2][3]);
	//	glVertex3d(boundingBox[0][0],boundingBox[1][0],boundingBox[2][0]);

	//	glVertex3d(boundingBox[0][1],boundingBox[1][1],boundingBox[2][1]);
	//	glVertex3d(boundingBox[0][2],boundingBox[1][2],boundingBox[2][2]);	

	//	glVertex3d(boundingBox[0][6],boundingBox[1][6],boundingBox[2][6]);
	//	glVertex3d(boundingBox[0][5],boundingBox[1][5],boundingBox[2][5]);


	//	if(1 == RGB[0])
	//	{
	//	glColor3f(1.0, 0.0, 0.0);
	//	}
	//	else if(2 == RGB[0])
	//	{
	//	glColor3f(0.0, 1.0, 0.0);
	//	}
	//	else if(3 == RGB[0])
	//	{
	//	glColor3f(0.0, 0.0, 1.0);
	//	}

	//	//glColor3f(1.0, 0.0, 0.0);
	//	glVertex3d(boundingBox[0][0],boundingBox[1][0],boundingBox[2][0]);
	//	glVertex3d(boundingBox[0][1],boundingBox[1][1],boundingBox[2][1]);	

	//	glVertex3d(boundingBox[0][2],boundingBox[1][2],boundingBox[2][2]);
	//	glVertex3d(boundingBox[0][3],boundingBox[1][3],boundingBox[2][3]);	

	//	glVertex3d(boundingBox[0][7],boundingBox[1][7],boundingBox[2][7]);
	//	glVertex3d(boundingBox[0][6],boundingBox[1][6],boundingBox[2][6]);	

	//	glVertex3d(boundingBox[0][5],boundingBox[1][5],boundingBox[2][5]);
	//	glVertex3d(boundingBox[0][4],boundingBox[1][4],boundingBox[2][4]);

	//	glEnd();


	//	for (i=0; i<this_poly->ntris; i++) {
	//	if (mode == GL_SELECT)
	//	glLoadName(i+1);

	//	Triangle *temp_t=this_poly->tlist[i];
	//	glBegin(GL_POLYGON);
	//	for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	if (i==this_poly->seed)
	//	glColor3f(0.0, 0.0, 1.0);
	//	else
	//	{
	//	// (2.b) Based on Normal	
	//	float rColor = temp_t->normal.x;
	//	float gColor = temp_t->normal.y;
	//	float bColor = temp_t->normal.z;
	//	if(rColor < 0.f)
	//	{
	//	rColor = -rColor;
	//	}
	//	if(gColor < 0.f)
	//	{
	//	gColor = -gColor;
	//	}
	//	if(bColor < 0.f)
	//	{
	//	bColor = -bColor;
	//	}

	//	glColor3f(rColor,gColor,bColor);
	//	//End Based on Normal
	//	//glColor3f(1.0, 1.0, 0.0);
	//	}

	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//	}
	//	glEnd();
	//	}
	//	break;

	//case 6: // Normal based bounding box

	//	mat_diffuse[0] = 1.0;
	//	mat_diffuse[1] = 1.0;
	//	mat_diffuse[2] = 0.0;
	//	mat_diffuse[3] = 1.0;

	//	//   bboxpoint1 = this_poly->center + norm1*(this_poly->radius/2) + norm2*(this_poly->radius/2) + norm3*(this_poly->radius/2);
	//	//bboxpoint2 = this_poly->center + norm1*(this_poly->radius/2) + norm2*(this_poly->radius/2) - norm3*(this_poly->radius/2);
	//	//bboxpoint3 = this_poly->center + norm1*(this_poly->radius/2) - norm2*(this_poly->radius/2) - norm3*(this_poly->radius/2);
	//	//bboxpoint4 = this_poly->center + norm1*(this_poly->radius/2) - norm2*(this_poly->radius/2) + norm3*(this_poly->radius/2);
	//	//bboxpoint5 = this_poly->center - norm1*(this_poly->radius/2) + norm2*(this_poly->radius/2) + norm3*(this_poly->radius/2);
	//	//bboxpoint6 = this_poly->center - norm1*(this_poly->radius/2) + norm2*(this_poly->radius/2) - norm3*(this_poly->radius/2);
	//	//bboxpoint7 = this_poly->center - norm1*(this_poly->radius/2) - norm2*(this_poly->radius/2) - norm3*(this_poly->radius/2);
	//	//bboxpoint8 = this_poly->center - norm1*(this_poly->radius/2) - norm2*(this_poly->radius/2) + norm3*(this_poly->radius/2);

	//	////icVector3 bboxpoint1 = median + norm1*normval1 + norm2*normval2 + norm3*normval3;
	//	////icVector3 bboxpoint2 = median + norm1*normval1 + norm2*normval2 - norm3*normval3;
	//	////icVector3 bboxpoint3 = median + norm1*normval1 - norm2*normval2 - norm3*normval3;
	//	////icVector3 bboxpoint4 = median + norm1*normval1 - norm2*normval2 + norm3*normval3;
	//	////icVector3 bboxpoint5 = median - norm1*normval1 + norm2*normval2 + norm3*normval3;
	//	////icVector3 bboxpoint6 = median - norm1*normval1 + norm2*normval2 - norm3*normval3;
	//	////icVector3 bboxpoint7 = median - norm1*normval1 - norm2*normval2 - norm3*normval3;
	//	////icVector3 bboxpoint8 = median - norm1*normval1 - norm2*normval2 + norm3*normval3;

	//	//boundingBox1[0][0] = bboxpoint1.x;
	//	//boundingBox1[1][0] = bboxpoint1.y;
	//	//boundingBox1[2][0] = bboxpoint1.z;

	//	//boundingBox1[0][1] = bboxpoint2.x;
	//	//boundingBox1[1][1] = bboxpoint2.y;
	//	//boundingBox1[2][1] = bboxpoint2.z;

	//	//boundingBox1[0][2] = bboxpoint3.x;
	//	//boundingBox1[1][2] = bboxpoint3.y;
	//	//boundingBox1[2][2] = bboxpoint3.z;

	//	//boundingBox1[0][3] = bboxpoint4.x;
	//	//boundingBox1[1][3] = bboxpoint4.y;
	//	//boundingBox1[2][3] = bboxpoint4.z;

	//	//boundingBox1[0][4] = bboxpoint5.x;
	//	//boundingBox1[1][4] = bboxpoint5.y;
	//	//boundingBox1[2][4] = bboxpoint5.z;

	//	//boundingBox1[0][5] = bboxpoint6.x;
	//	//boundingBox1[1][5] = bboxpoint6.y;
	//	//boundingBox1[2][5] = bboxpoint6.z;

	//	//boundingBox1[0][6] = bboxpoint7.x;
	//	//boundingBox1[1][6] = bboxpoint7.y;
	//	//boundingBox1[2][6] = bboxpoint7.z;

	//	//boundingBox1[0][7] = bboxpoint8.x;
	//	//boundingBox1[1][7] = bboxpoint8.y;
	//	//boundingBox1[2][7] = bboxpoint8.z;
	//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	//	glDisable(GL_LIGHTING);
	//	glDisable(GL_LIGHT0);
	//	glDisable(GL_LIGHT1);
	//	glColor3f(0.0, 0.0, 0.0);
	//	glLineWidth(5.0);


	//	glBegin(GL_LINES);	
	//	glVertex3d(boundingBox2[0][0],boundingBox2[1][0],boundingBox2[2][0]);
	//	glVertex3d(boundingBox2[0][1],boundingBox2[1][1],boundingBox2[2][1]);	

	//	glVertex3d(boundingBox2[0][1],boundingBox2[1][1],boundingBox2[2][1]);
	//	glVertex3d(boundingBox2[0][2],boundingBox2[1][2],boundingBox2[2][2]);	

	//	glVertex3d(boundingBox2[0][2],boundingBox2[1][2],boundingBox2[2][2]);
	//	glVertex3d(boundingBox2[0][3],boundingBox2[1][3],boundingBox2[2][3]);	

	//	glVertex3d(boundingBox2[0][3],boundingBox2[1][3],boundingBox2[2][3]);
	//	glVertex3d(boundingBox2[0][0],boundingBox2[1][0],boundingBox2[2][0]);	

	//	glVertex3d(boundingBox2[0][0],boundingBox2[1][0],boundingBox2[2][0]);
	//	glVertex3d(boundingBox2[0][4],boundingBox2[1][4],boundingBox2[2][4]);	

	//	glVertex3d(boundingBox2[0][1],boundingBox2[1][1],boundingBox2[2][1]);
	//	glVertex3d(boundingBox2[0][5],boundingBox2[1][5],boundingBox2[2][5]);	

	//	glVertex3d(boundingBox2[0][2],boundingBox2[1][2],boundingBox2[2][2]);
	//	glVertex3d(boundingBox2[0][6],boundingBox2[1][6],boundingBox2[2][6]);	

	//	glVertex3d(boundingBox2[0][3],boundingBox2[1][3],boundingBox2[2][3]);
	//	glVertex3d(boundingBox2[0][7],boundingBox2[1][7],boundingBox2[2][7]);	

	//	glVertex3d(boundingBox2[0][4],boundingBox2[1][4],boundingBox2[2][4]);
	//	glVertex3d(boundingBox2[0][7],boundingBox2[1][7],boundingBox2[2][7]);	

	//	glVertex3d(boundingBox2[0][7],boundingBox2[1][7],boundingBox2[2][7]);
	//	glVertex3d(boundingBox2[0][6],boundingBox2[1][6],boundingBox2[2][6]);	

	//	glVertex3d(boundingBox2[0][6],boundingBox2[1][6],boundingBox2[2][6]);
	//	glVertex3d(boundingBox2[0][5],boundingBox2[1][5],boundingBox2[2][5]);	

	//	glVertex3d(boundingBox2[0][5],boundingBox2[1][5],boundingBox2[2][5]);
	//	glVertex3d(boundingBox2[0][4],boundingBox2[1][4],boundingBox2[2][4]);

	//	glEnd();


	//	for (i=0; i<this_poly->ntris; i++) {
	//	if (mode == GL_SELECT)
	//	glLoadName(i+1);

	//	Triangle *temp_t=this_poly->tlist[i];
	//	glBegin(GL_POLYGON);
	//	for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	if (i==this_poly->seed)
	//	glColor3f(0.0, 0.0, 1.0);
	//	else
	//	{
	//	// (2.b) Based on Normal	
	//	float rColor = temp_t->normal.x;
	//	float gColor = temp_t->normal.y;
	//	float bColor = temp_t->normal.z;
	//	if(rColor < 0.f)
	//	{
	//	rColor = -rColor;
	//	}
	//	if(gColor < 0.f)
	//	{
	//	gColor = -gColor;
	//	}
	//	if(bColor < 0.f)
	//	{
	//	bColor = -bColor;
	//	}

	//	glColor3f(rColor,gColor,bColor);
	//	//End Based on Normal
	//	//glColor3f(1.0, 1.0, 0.0);
	//	}

	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//	}
	//	glEnd();
	//	}
	//	break;
	case 1: // Show Regular Mesh
	//this_poly->calculate_bounding_box_with_moments();

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
	
	// Show the Triangulated data
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
	// Here we will highlight the point based on deficit 
	// Deficit 1  = Red
	// Deficit 2  = Blue
	// Deficit 3  = Green
	// Deficit -1 = R0G1B1
	// Deficit -2 = R1G0B1
	// Deficit -3 = R1G1B0

	//glPointSize(5);
	//glBegin(GL_POINTS);
	//for(int gg=0; gg<this_poly->nverts;gg++)
	//{
	//	int deficit = 6 - this_poly->vlist[gg]->ntris;
	//	if(deficit != 0)
	//	{
	//	if(1 == deficit)
	//	{
	//	glColor3f(1.0, 0.0, 0.0);
	//	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	//	}
	//	if(2 == deficit)
	//	{
	//	glColor3f(0.0, 1.0, 0.0);
	//	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	//	}
	//	if(3 == deficit)
	//	{
	//	glColor3f(0.0, 0.0, 1.0);
	//	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	//	}
	//	if(-1 == deficit)
	//	{
	//	glColor3f(0.0, 1.0, 1.0);
	//	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	//	}
	//	if(-2 == deficit)
	//	{
	//	glColor3f(1.0, 0.0, 1.0);
	//	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	//	}
	//	if(-3 == deficit)
	//	{
	//	glColor3f(1.0, 1.0, 0.0);
	//	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	//	}
	//	if(-3 > deficit)
	//	{
	//	glColor3f(0.0, 0.0, 0.0);
	//	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	//	}
	//	}
	//}

	//glEnd();
	}
	break;

	case 2: // Show Regular Subdivision of Mesh
	//this_poly->calculate_bounding_box_with_moments();

	mat_diffuse[0] = 1.0;
	mat_diffuse[1] = 1.0;
	mat_diffuse[2] = 0.0;
	mat_diffuse[3] = 1.0;

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	for (i=0; i<polyR->ntris; i++) {
	if (mode == GL_SELECT)
	glLoadName(i+1);

	Triangle *temp_t=polyR->tlist[i];
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

	// Show the Triangulated data
	glLineWidth(1);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glColor3f(0.0, 0.0, 0.0);

	//for(int kk = 0; kk<4 ; kk++)
	//{
	//   Triangle* temp_tt=polyR->tlist[4*i + kk];
	//   glBegin(GL_LINE_LOOP);
	//   for (j=0; j<3; j++) 
	//   {
	//	   Vertex *temp_v = temp_tt->verts[j];	
	//	   glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//   }
	//   glEnd();
	//  }	

	glBegin(GL_LINE_LOOP);
	for (j=0; j<3; j++) {

	Vertex *temp_v = temp_t->verts[j];
	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
	}
	break;

	case 3: // Show irregular Subdivision of Mesh
	//this_poly->calculate_bounding_box_with_moments();

	mat_diffuse[0] = 1.0;
	mat_diffuse[1] = 1.0;
	mat_diffuse[2] = 0.0;
	mat_diffuse[3] = 1.0;

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	for (i=0; i<polyIR->ntris; i++) {
	if (mode == GL_SELECT)
	glLoadName(i+1);

	Triangle *temp_t=polyIR->tlist[i];
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

	// Show the Triangulated data
	glLineWidth(1);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glColor3f(0.0, 0.0, 0.0);
	/*	glBegin(GL_LINE_LOOP);           
	for(int kk = 0; kk<4 ; kk++)
	{
	Triangle* temp_tt=polyIR->tlist[4*i + kk];
	glBegin(GL_LINE_LOOP);
	for (j=0; j<3; j++) 
	{
	Vertex *temp_v = temp_tt->verts[j];	
	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
	}	*/	

	//for(int kk = 0; kk<polyIR->ntris;kk++)
	//{
	//   Triangle* temp_tt=polyIR->tlist[kk];
	//   glBegin(GL_LINE_LOOP);
	//   for (j=0; j<3; j++) 
	//   {
	//	   Vertex *temp_v = temp_tt->verts[j];	
	//	   glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//   }
	//   glEnd();
	//}
	glBegin(GL_LINE_LOOP);
	for (j=0; j<3; j++) {

	Vertex *temp_v = temp_t->verts[j];
	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
	}
	break;
case 4: // Show only 3D ChekerBoard Color for normal mesh
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	if(4 == lastcolorscheme)
	{
	 L = L+0.03;
	}
	else
	{
	lastcolorscheme = 4;
	//L = 0.003;
	}
	
	for (i=0; i<this_poly->ntris; i++) {
	if (mode == GL_SELECT)
	glLoadName(i+1);

	Triangle *temp_t=this_poly->tlist[i];
	glBegin(GL_POLYGON);
	for (j=0; j<3; j++) {

	Vertex *temp_v = temp_t->verts[j];
	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	if (i==this_poly->seed)
	glColor3f(0.0, 0.0, 1.0);
	else
	{
	// (2.c) Color to vertex
	float vX = temp_v->x;
	float vY = temp_v->y;
	float vZ = temp_v->z;

	int rInt =floor( vX/L);
	int gInt =floor( vY/L);
	int bInt = floor(vZ/L);

	if(0 == rInt%2)
	{
	rColor = 1.f;
	}
	else
	{
	rColor = 0.f;
	}

	if(0 == gInt%2)
	{
	gColor = 1.f;
	}
	else
	{
	gColor = 0.f;
	}

	if(0 == bInt%2)
	{
	bColor = 1.f;
	}
	else
	{
	bColor = 0.f;
	}
	glColor3f(rColor, gColor, bColor);
	}
	// End Color to vertex
	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();

	//// Show the Triangulated data
	//glLineWidth(1);
	//glDisable(GL_LIGHTING);
	//glDisable(GL_LIGHT0);
	//glDisable(GL_LIGHT1);
	//glColor3f(0.0, 0.0, 0.0);

	//glBegin(GL_LINE_LOOP);
	//for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//}
	//glEnd();
	}
	break;
 case 5: // Show 3D ChekerBoard Color for regular subdivision
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	if(5 == lastcolorscheme)
	{
	L = L+0.03;
	}
	else
	{
	lastcolorscheme = 5;
	//L = 0.003;
	}
	
	for (i=0; i<polyR->ntris; i++) {
	if (mode == GL_SELECT)
	glLoadName(i+1);

	Triangle *temp_t=polyR->tlist[i];
	glBegin(GL_POLYGON);
	for (j=0; j<3; j++) {

	Vertex *temp_v = temp_t->verts[j];
	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	if (i==this_poly->seed)
	glColor3f(0.0, 0.0, 1.0);
	else
	{
	// (2.c) Color to vertex
	float vX = temp_v->x;
	float vY = temp_v->y;
	float vZ = temp_v->z;

	int rInt =floor( vX/L);
	int gInt =floor( vY/L);
	int bInt = floor(vZ/L);

	if(0 == rInt%2)
	{
	rColor = 1.f;
	}
	else
	{
	rColor = 0.f;
	}

	if(0 == gInt%2)
	{
	gColor = 1.f;
	}
	else
	{
	gColor = 0.f;
	}

	if(0 == bInt%2)
	{
	bColor = 1.f;
	}
	else
	{
	bColor = 0.f;
	}
	glColor3f(rColor, gColor, bColor);
	}
	// End Color to vertex
	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();

	//// Show the Triangulated data
	//glLineWidth(1);
	//glDisable(GL_LIGHTING);
	//glDisable(GL_LIGHT0);
	//glDisable(GL_LIGHT1);
	//glColor3f(0.0, 0.0, 0.0);

	//glBegin(GL_LINE_LOOP);
	//for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//}
	//glEnd();
	}
	break;
  case 6: // Show 3dCheckboard for irregular subdivision
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	if(6 == lastcolorscheme)
	{
	L = L+0.03;
	}
	else
	{
	lastcolorscheme = 6;
	//L = 0.003;
	}
	for (i=0; i<polyIR->ntris; i++) {
	if (mode == GL_SELECT)
	glLoadName(i+1);

	Triangle *temp_t=polyIR->tlist[i];
	glBegin(GL_POLYGON);
	for (j=0; j<3; j++) {

	Vertex *temp_v = temp_t->verts[j];
	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	if (i==this_poly->seed)
	glColor3f(0.0, 0.0, 1.0);
	else
	{
	// (2.c) Color to vertex
	float vX = temp_v->x;
	float vY = temp_v->y;
	float vZ = temp_v->z;

	int rInt =floor( vX/L);
	int gInt =floor( vY/L);
	int bInt = floor(vZ/L);

	if(0 == rInt%2)
	{
	rColor = 1.f;
	}
	else
	{
	rColor = 0.f;
	}

	if(0 == gInt%2)
	{
	gColor = 1.f;
	}
	else
	{
	gColor = 0.f;
	}

	if(0 == bInt%2)
	{
	bColor = 1.f;
	}
	else
	{
	bColor = 0.f;
	}
	glColor3f(rColor, gColor, bColor);
	}
	// End Color to vertex
	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();

	//// Show the Triangulated data
	//glLineWidth(1);
	//glDisable(GL_LIGHTING);
	//glDisable(GL_LIGHT0);
	//glDisable(GL_LIGHT1);
	//glColor3f(0.0, 0.0, 0.0);

	//glBegin(GL_LINE_LOOP);
	//for (j=0; j<3; j++) {

	//	Vertex *temp_v = temp_t->verts[j];
	//	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	//}
	//glEnd();
	}
	break;

	   case 7: // Show normal mesh and 3D ChekerBoard Color
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
	
	// Show the Triangulated data
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
	// Here we will highlight the point based on deficit 
	// Deficit 1  = Red
	// Deficit 2  = Blue
	// Deficit 3  = Green
	// Deficit -1 = R0G1B1
	// Deficit -2 = R1G0B1
	// Deficit -3 = R1G1B0

	
	}
	TotoalDeficitCount = 0;
	glPointSize(8);
	glBegin(GL_POINTS);
	for(int gg=0; gg<this_poly->nverts;gg++)
	{
	int deficit = 6 - this_poly->vlist[gg]->ntris;
	if(deficit != 0)
	{
	TotoalDeficitCount = TotoalDeficitCount +deficit;
	if(1 == deficit)
	{
	glColor3f(1.0, 0.0, 0.0);
	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	}
	if(2 == deficit)
	{
	glColor3f(0.0, 1.0, 0.0);
	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	}
	if(3 == deficit)
	{
	glColor3f(0.0, 0.0, 1.0);
	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	}
	if(-1 == deficit)
	{
	glColor3f(0.0, 1.0, 1.0);
	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	}
	if(-2 == deficit)
	{
	glColor3f(1.0, 0.0, 1.0);
	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	}
	if(-3 == deficit)
	{
	glColor3f(1.0, 1.0, 0.0);
	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	}
	if(-3 > deficit)
	{
	glColor3f(0.0, 0.0, 0.0);
	glVertex3d(this_poly->vlist[gg]->x, this_poly->vlist[gg]->y,this_poly->vlist[gg]->z);
	}
	}
	// Here we calculate the angle deficit
	float totalAngleOfVertex = 0.f;
	for(int mm = 0; mm<this_poly->vlist[gg]->ntris;mm++)	
	{
	Triangle *currentTriangle = this_poly->vlist[gg]->tris[mm];
	Vertex *v0 = this_poly->vlist[gg];
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
	totalAngleOfVertex =totalAngleOfVertex + angle;
	
	}
	}
	double adeficit = 360 - totalAngleOfVertex;
	TotalAngleDeficit =   TotalAngleDeficit + adeficit;
	}

	glEnd();
	cout <<"Total Valence Deficit for Regular Mesh is: "<<TotoalDeficitCount<<endl;
	cout <<"Total Angle Deficit for Regular Mesh is: "<<TotalAngleDeficit<<endl;
	break;
	   
	  
	   case 8: // Show regular subdivision and 3D ChekerBoard Color
	mat_diffuse[0] = 1.0;
	mat_diffuse[1] = 1.0;
	mat_diffuse[2] = 0.0;
	mat_diffuse[3] = 1.0;

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	for (i=0; i<polyR->ntris; i++) {
	if (mode == GL_SELECT)
	glLoadName(i+1);

	Triangle *temp_t=polyR->tlist[i];
	glBegin(GL_POLYGON);
	for (j=0; j<3; j++) {

	Vertex *temp_v = temp_t->verts[j];
	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	glColor3f(1.0, 1.0, 0.0);
	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
	
	// Show the Triangulated data
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
	for(int gg=0; gg<polyR->nverts;gg++)
	{
	int deficit = 6 - polyR->vlist[gg]->ntris;
	if(deficit != 0)
	{
	TotoalDeficitCount = TotoalDeficitCount + deficit;
	if(1 == deficit)
	{
	glColor3f(1.0, 0.0, 0.0);
	glVertex3d(polyR->vlist[gg]->x, polyR->vlist[gg]->y,polyR->vlist[gg]->z);
	}
	if(2 == deficit)
	{
	glColor3f(0.0, 1.0, 0.0);
	glVertex3d(polyR->vlist[gg]->x, polyR->vlist[gg]->y,polyR->vlist[gg]->z);
	}
	if(3 == deficit)
	{
	glColor3f(0.0, 0.0, 1.0);
	glVertex3d(polyR->vlist[gg]->x, polyR->vlist[gg]->y,polyR->vlist[gg]->z);
	}
	if(-1 == deficit)
	{
	glColor3f(0.0, 1.0, 1.0);
	glVertex3d(polyR->vlist[gg]->x, polyR->vlist[gg]->y,polyR->vlist[gg]->z);
	}
	if(-2 == deficit)
	{
	glColor3f(1.0, 0.0, 1.0);
	glVertex3d(polyR->vlist[gg]->x, polyR->vlist[gg]->y,polyR->vlist[gg]->z);
	}
	if(-3 == deficit)
	{
	glColor3f(1.0, 1.0, 0.0);
	glVertex3d(polyR->vlist[gg]->x, polyR->vlist[gg]->y,polyR->vlist[gg]->z);
	}
	if(-3 > deficit)
	{
	glColor3f(0.0, 0.0, 0.0);
	glVertex3d(polyR->vlist[gg]->x, polyR->vlist[gg]->y,polyR->vlist[gg]->z);
	}
	}
	
	// Here we calculate the angle deficit
	float totalAngleOfVertex = 0.f;
	for(int mm = 0; mm<polyR->vlist[gg]->ntris;mm++)	
	{
	Triangle *currentTriangle = polyR->vlist[gg]->tris[mm];
	Vertex *v0 = polyR->vlist[gg];
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
	totalAngleOfVertex =totalAngleOfVertex + angle;
	
	}
	}
	double adeficit = 360 - totalAngleOfVertex;
	TotalAngleDeficit =   TotalAngleDeficit + adeficit;
	}

	glEnd();
	cout <<"Total Valence Deficit for Regular Division Mesh is: "<<TotoalDeficitCount<<endl;
	cout <<"Total Deficit for Regular Division Mesh is: "<<TotalAngleDeficit<<endl;
	break;
	  
	  
	
	   case 9: // Show Irregular subdivision and 3dCheckboard
	mat_diffuse[0] = 1.0;
	mat_diffuse[1] = 1.0;
	mat_diffuse[2] = 0.0;
	mat_diffuse[3] = 1.0;

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	for (i=0; i<polyIR->ntris; i++) {
	if (mode == GL_SELECT)
	glLoadName(i+1);

	Triangle *temp_t=polyIR->tlist[i];
	glBegin(GL_POLYGON);
	for (j=0; j<3; j++) {

	Vertex *temp_v = temp_t->verts[j];
	glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
	glColor3f(1.0, 1.0, 0.0);
	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
	
	// Show the Triangulated data
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
	for(int gg=0; gg<polyIR->nverts;gg++)
	{
	int deficit = 6 - polyIR->vlist[gg]->ntris;
	if(deficit != 0)
	{
	TotoalDeficitCount = TotoalDeficitCount + deficit;
	if(1 == deficit)
	{
	glColor3f(1.0, 0.0, 0.0);
	glVertex3d(polyIR->vlist[gg]->x, polyIR->vlist[gg]->y,polyIR->vlist[gg]->z);
	}
	if(2 == deficit)
	{
	glColor3f(0.0, 1.0, 0.0);
	glVertex3d(polyIR->vlist[gg]->x, polyIR->vlist[gg]->y,polyIR->vlist[gg]->z);
	}
	if(3 == deficit)
	{
	glColor3f(0.0, 0.0, 1.0);
	glVertex3d(polyIR->vlist[gg]->x, polyIR->vlist[gg]->y,polyIR->vlist[gg]->z);
	}
	if(-1 == deficit)
	{
	glColor3f(0.0, 1.0, 1.0);
	glVertex3d(polyIR->vlist[gg]->x, polyIR->vlist[gg]->y,polyIR->vlist[gg]->z);
	}
	if(-2 == deficit)
	{
	glColor3f(1.0, 0.0, 1.0);
	glVertex3d(polyIR->vlist[gg]->x, polyIR->vlist[gg]->y,polyIR->vlist[gg]->z);
	}
	if(-3 == deficit)
	{
	glColor3f(1.0, 1.0, 0.0);
	glVertex3d(polyIR->vlist[gg]->x, polyIR->vlist[gg]->y,polyIR->vlist[gg]->z);
	}
	if(-3 > deficit)
	{
	glColor3f(0.0, 0.0, 0.0);
	glVertex3d(polyIR->vlist[gg]->x, polyIR->vlist[gg]->y,polyIR->vlist[gg]->z);
	}
	}
	
	// Here we calculate the angle deficit
	float totalAngleOfVertex = 0.f;
	for(int mm = 0; mm<polyIR->vlist[gg]->ntris;mm++)	
	{
	Triangle *currentTriangle = polyIR->vlist[gg]->tris[mm];
	Vertex *v0 = polyIR->vlist[gg];
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
	totalAngleOfVertex =totalAngleOfVertex + angle;
	
	}
	}
	double adeficit = 360 - totalAngleOfVertex;
	TotalAngleDeficit =   TotalAngleDeficit + adeficit;
	}

	glEnd();
	cout <<"Total Valence Deficit for Regular Division Mesh is: "<<TotoalDeficitCount<<endl;
	cout <<"Total Deficit for Regular Division Mesh is: "<<TotalAngleDeficit<<endl;
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
	switch(ACSIZE){
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


 void Polyhedron::calculate_bounding_box_with_moments()
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

	 int vCount=0;
	 
	 for(vCount=0; vCount<nverts; vCount++)
	 {
	  xSum += vlist[vCount]->x;
	  ySum += vlist[vCount]->y;
	  zSum += vlist[vCount]->z;
	 }

	 icVector3 center1(xSum/nverts,ySum/nverts,zSum/nverts);
	 for(vCount=0; vCount<nverts; vCount++)
	 {
	   float XX = vlist[vCount]->x - center1.x;
	   float YY = vlist[vCount]->y - center1.y;
	   float ZZ = vlist[vCount]->z - center1.z;

	   AA = AA + XX*XX;
	   BB = BB + XX*YY;
	   CC = CC + XX*ZZ;
	   DD = DD + YY*YY;
	   EE = EE + YY*ZZ;
	   FF = FF + ZZ*ZZ;
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

	 float majoreigen[3];
     float middleeigen[3];
     float minorigen[3];

	 if(normval1>=normval2 && normval1>=normval3)
	 {
	 majoreigen[0] = norm1.x;
	 majoreigen[1] = norm1.y;
	 majoreigen[2] = norm1.z;

	 RGB[0] = 1;
	 if(normval2 >= normval3)
	 {
	   middleeigen[0] = norm2.x;
	   	   middleeigen[1] = norm2.y;
	   middleeigen[2] = norm2.z;	

	   minoreigen[0] = norm3.x;
	   	   minoreigen[1] = norm3.y;
	   minoreigen[2] = norm3.z;	

	   RGB[1] = 2;
	   RGB[2] = 3;
	 }
	 else if(normval3 >= normval2)
	 {
	 middleeigen[0] = norm3.x;
	 middleeigen[1] = norm3.y;
	 middleeigen[2] = norm3.z;	

	 minoreigen[0] = norm2.x;
	 minoreigen[1] = norm2.y;
	 minoreigen[2] = norm2.z;	
	 RGB[2] = 2;
	 RGB[1] = 3;
	 }	 
	 }
	 else if(normval2>=normval1 && normval2>=normval3)
	 {
	 majoreigen[0] = norm2.x;
	 majoreigen[1] = norm2.y;
	 majoreigen[2] = norm2.z;

	 RGB[1] = 1;

	 if(normval1 >= normval3)
	 {
	   middleeigen[0] = norm1.x;
	   	   middleeigen[1] = norm1.y;
	   middleeigen[2] = norm1.z;	

	   minoreigen[0] = norm3.x;
	   	   minoreigen[1] = norm3.y;
	   minoreigen[2] = norm3.z;	
	   RGB[0] = 2;
	   RGB[2] = 3;
	 }
	 else if(normval3 >= normval1)
	 {
	 middleeigen[0] = norm3.x;
	 middleeigen[1] = norm3.y;
	 middleeigen[2] = norm3.z;	

	 minoreigen[0] = norm1.x;
	 minoreigen[1] = norm1.y;
	 minoreigen[2] = norm1.z;	
	 RGB[2] = 2;
	 RGB[0] = 3;
	 }	 
	 }
	 else if(normval3>=normval1 && normval3>=normval2)
	 {
	 majoreigen[0] = norm3.x;
	 majoreigen[1] = norm3.y;
	 majoreigen[2] = norm3.z;

	 RGB[0] = 1;

	 if(normval1 >= normval2)
	 {
	   middleeigen[0] = norm1.x;
	   	   middleeigen[1] = norm1.y;
	   middleeigen[2] = norm1.z;	

	   minoreigen[0] = norm2.x;
	   	   minoreigen[1] = norm2.y;
	   minoreigen[2] = norm2.z;	

	   RGB[1] = 2;
	   RGB[2] = 3;
	 }
	 else if(normval2 >= normval1)
	 {
	 middleeigen[0] = norm2.x;
	 middleeigen[1] = norm2.y;
	 middleeigen[2] = norm2.z;	

	 minoreigen[0] = norm1.x;
	 minoreigen[1] = norm1.y;
	 minoreigen[2] = norm1.z;	
	 
	 RGB[1] = 3;
	     RGB[2] = 2;
	 }	 
	 }


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

	 for(vCount=0; vCount<nverts; vCount++)
	 {
	icVector3 pVector((vlist[vCount]->x -median.x),(vlist[vCount]->y -median.y),(vlist[vCount]->z -median.z));
	double dotProductX = dot(pVector,norm1);
	double dotProductY = dot(pVector,norm2);
	double dotProductZ = dot(pVector,norm3);

	// X Extent
	if(dotProductX>maxDotXProduct)
	{
	maxDotXProduct = dotProductX;
	vertexWithMaxDotX = vCount;
	}
	
	if(dotProductX<minDotXProduct)
	{
	minDotXProduct = dotProductX;
	vertexWithMinDotX = vCount;
	}

	// Y Extent

	if(dotProductY>maxDotYProduct)
	{
	maxDotYProduct = dotProductY;
	vertexWithMaxDotY = vCount;
	}
	
	if(dotProductY<minDotYProduct)
	{
	minDotYProduct = dotProductY;
	vertexWithMinDotY = vCount;
	}

	// Z Extent
	if(dotProductZ>maxDotZProduct)
	{
	maxDotZProduct = dotProductZ;
	vertexWithMaxDotZ = vCount;
	}
	
	if(dotProductZ<minDotZProduct)
	{
	minDotZProduct = dotProductZ;
	vertexWithMinDotZ = vCount;
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

	 //icVector3 bboxpoint1 = median + norm1*normval1 + norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint2 = median + norm1*normval1 + norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint3 = median + norm1*normval1 - norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint4 = median + norm1*normval1 - norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint5 = median - norm1*normval1 + norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint6 = median - norm1*normval1 + norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint7 = median - norm1*normval1 - norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint8 = median - norm1*normval1 - norm2*normval2 + norm3*normval3;

	 boundingBox[0][0] = bboxpoint1.x;
	 boundingBox[1][0] = bboxpoint1.y;
	 boundingBox[2][0] = bboxpoint1.z;

	 boundingBox[0][1] = bboxpoint2.x;
	 boundingBox[1][1] = bboxpoint2.y;
	 boundingBox[2][1] = bboxpoint2.z;

	 boundingBox[0][2] = bboxpoint3.x;
	 boundingBox[1][2] = bboxpoint3.y;
	 boundingBox[2][2] = bboxpoint3.z;

	 boundingBox[0][3] = bboxpoint4.x;
	 boundingBox[1][3] = bboxpoint4.y;
	 boundingBox[2][3] = bboxpoint4.z;

	 boundingBox[0][4] = bboxpoint5.x;
	 boundingBox[1][4] = bboxpoint5.y;
	 boundingBox[2][4] = bboxpoint5.z;

	 boundingBox[0][5] = bboxpoint6.x;
	 boundingBox[1][5] = bboxpoint6.y;
	 boundingBox[2][5] = bboxpoint6.z;

	 boundingBox[0][6] = bboxpoint7.x;
	 boundingBox[1][6] = bboxpoint7.y;
	 boundingBox[2][6] = bboxpoint7.z;

	 boundingBox[0][7] = bboxpoint8.x;
	 boundingBox[1][7] = bboxpoint8.y;
	 boundingBox[2][7] = bboxpoint8.z;
 }


 void Polyhedron::calculate_bounding_box_with_normals()
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

	 int vCount=0;
	 
	 for(vCount=0; vCount<nverts; vCount++)
	 {
	  xSum += vlist[vCount]->x;
	  ySum += vlist[vCount]->y;
	  zSum += vlist[vCount]->z;
	 }

	 icVector3 center1(xSum/nverts,ySum/nverts,zSum/nverts);
	 for(vCount=0; vCount<nverts; vCount++)
	 {
	   /*float XX = vlist[vCount]->x - center1.x;
	   float YY = vlist[vCount]->y - center1.y;
	   float ZZ = vlist[vCount]->z - center1.z;*/

	   // these values are from the average normal to the the vertex now
	   float XX = vlist[vCount]->normal.x;
	   float YY = vlist[vCount]->normal.y;
	   float ZZ = vlist[vCount]->normal.z;

	   AA = AA + XX*XX;
	   BB = BB + XX*YY;
	   CC = CC + XX*ZZ;
	   DD = DD + YY*YY;
	   EE = EE + YY*ZZ;
	   FF = FF + ZZ*ZZ;
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

	 for(vCount=0; vCount<nverts; vCount++)
	 {
	icVector3 pVector((vlist[vCount]->x -median.x),(vlist[vCount]->y -median.y),(vlist[vCount]->z -median.z));
	double dotProductX = dot(pVector,norm1);
	double dotProductY = dot(pVector,norm2);
	double dotProductZ = dot(pVector,norm3);

	// X Extent
	if(dotProductX>maxDotXProduct)
	{
	maxDotXProduct = dotProductX;
	vertexWithMaxDotX = vCount;
	}
	
	if(dotProductX<minDotXProduct)
	{
	minDotXProduct = dotProductX;
	vertexWithMinDotX = vCount;
	}

	// Y Extent

	if(dotProductY>maxDotYProduct)
	{
	maxDotYProduct = dotProductY;
	vertexWithMaxDotY = vCount;
	}
	
	if(dotProductY<minDotYProduct)
	{
	minDotYProduct = dotProductY;
	vertexWithMinDotY = vCount;
	}

	// Z Extent
	if(dotProductZ>maxDotZProduct)
	{
	maxDotZProduct = dotProductZ;
	vertexWithMaxDotZ = vCount;
	}
	
	if(dotProductZ<minDotZProduct)
	{
	minDotZProduct = dotProductZ;
	vertexWithMinDotZ = vCount;
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

	 //icVector3 bboxpoint1 = median + norm1*normval1 + norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint2 = median + norm1*normval1 + norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint3 = median + norm1*normval1 - norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint4 = median + norm1*normval1 - norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint5 = median - norm1*normval1 + norm2*normval2 + norm3*normval3;
	 //icVector3 bboxpoint6 = median - norm1*normval1 + norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint7 = median - norm1*normval1 - norm2*normval2 - norm3*normval3;
	 //icVector3 bboxpoint8 = median - norm1*normval1 - norm2*normval2 + norm3*normval3;

	 boundingBox2[0][0] = bboxpoint1.x;
	 boundingBox2[1][0] = bboxpoint1.y;
	 boundingBox2[2][0] = bboxpoint1.z;

	 boundingBox2[0][1] = bboxpoint2.x;
	 boundingBox2[1][1] = bboxpoint2.y;
	 boundingBox2[2][1] = bboxpoint2.z;

	 boundingBox2[0][2] = bboxpoint3.x;
	 boundingBox2[1][2] = bboxpoint3.y;
	 boundingBox2[2][2] = bboxpoint3.z;

	 boundingBox2[0][3] = bboxpoint4.x;
	 boundingBox2[1][3] = bboxpoint4.y;
	 boundingBox2[2][3] = bboxpoint4.z;

	 boundingBox2[0][4] = bboxpoint5.x;
	 boundingBox2[1][4] = bboxpoint5.y;
	 boundingBox2[2][4] = bboxpoint5.z;

	 boundingBox2[0][5] = bboxpoint6.x;
	 boundingBox2[1][5] = bboxpoint6.y;
	 boundingBox2[2][5] = bboxpoint6.z;

	 boundingBox2[0][6] = bboxpoint7.x;
	 boundingBox2[1][6] = bboxpoint7.y;
	 boundingBox2[2][6] = bboxpoint7.z;

	 boundingBox2[0][7] = bboxpoint8.x;
	 boundingBox2[1][7] = bboxpoint8.y;
	 boundingBox2[2][7] = bboxpoint8.z;
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
	ROTATE(a,j,ip,j,iq)
	}
	for (j=ip+1;j<=iq-1;j++) {
	ROTATE(a,ip,j,j,iq)
	}
	for (j=iq+1;j<n;j++) { 
	ROTATE(a,ip,j,iq,j)
	}
	for (j=0;j<n;j++) {
	ROTATE(v,j,ip,j,iq)
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

void calculate_eigens1(float a[][3], int n, float d[], float v[][3])
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	float aDash[4][4];
	for(int mm =0; mm<3 ; mm++)
	{
	for(int  hh=0; hh<3 ; hh++)
	{
	aDash[mm+1][hh+1] =a[mm][hh]; 
	}
	}
	float dDash[4];
	float vDash[4][4];

	float bDash[4];
	float zDash[4];
	for (ip=1;ip<=n;ip++) { 
	for (iq=1;iq<=n;iq++) vDash[ip][iq]=0.0;
	vDash[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) { 
	bDash[ip]=dDash[ip]=aDash[ip][ip]; 
	zDash[ip]=0.0;
	}

	for (i=1;i<=50;i++) {
	sm=0.0;
	for (ip=1;ip<=n-1;ip++) {
	for (iq=ip+1;iq<=n;iq++)
	sm += fabs(aDash[ip][iq]);
	}
	if (sm == 0.0) {
	for(int dd =0; dd<3 ; dd++)
	{
	d[dd] = dDash[dd+1];
	}

	for(int rr =0; rr<3 ; rr++)
	{
	for(int  pp=0; pp<3 ; pp++)
	{
	v[rr][pp] = vDash[rr+1][pp+1];
	}
	}
	return;
	}
	if (i < 4)
	tresh=0.2*sm/(n*n); 
	else
	tresh=0.0;
	for (ip=1;ip<=n-1;ip++) {
	for (iq=ip+1;iq<=n;iq++) {
	g=100.0*fabs(aDash[ip][iq]);
	if (i > 4 && (float)(fabs(dDash[ip])+g) == (float)fabs(dDash[ip])
	&& (float)(fabs(dDash[iq])+g) == (float)fabs(dDash[iq]))
	aDash[ip][iq]=0.0;
	else if (fabs(aDash[ip][iq]) > tresh) {
	h=dDash[iq]-dDash[ip];
	if ((float)(fabs(h)+g) == (float)fabs(h))
	t=(aDash[ip][iq])/h; 
	else {
	theta=0.5*h/(aDash[ip][iq]);
	t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	if (theta < 0.0) t = -t;
	}
	c=1.0/sqrt(1+t*t);
	s=t*c;
	tau=s/(1.0+c);
	h=t*aDash[ip][iq];
	zDash[ip] -= h;
	zDash[iq] += h;
	dDash[ip] -= h;
	dDash[iq] += h;
	aDash[ip][iq]=0.0;
	for (j=1;j<=ip-1;j++) { 
	ROTATE(aDash,j,ip,j,iq)
	}
	for (j=ip+1;j<=iq-1;j++) {
	ROTATE(aDash,ip,j,j,iq)
	}
	for (j=iq+1;j<=n;j++) { 
	ROTATE(aDash,ip,j,iq,j)
	}
	for (j=1;j<=n;j++) {
	ROTATE(vDash,j,ip,j,iq)
	}

	}
	}
	}
	for (ip=1;ip<=n;ip++) {
	bDash[ip] += zDash[ip];
	dDash[ip]=bDash[ip]; 
	zDash[ip]=0.0; 
	}
	}

	for(int nn =0; nn<3 ; nn++)
	{
	d[nn] = dDash[nn+1];
	}

	for(int kk =0; kk<3 ; kk++)
	{
	for(int zz =0; zz<3 ; zz++)
	{
	v[kk][zz] = vDash[kk+1][zz+1];
	}
	}
}

