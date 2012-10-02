/**
 * @file project.cpp
 * @brief Geometry project
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#include "geometry/project.hpp"
#include "application/opengl.hpp"

/*
   A namespace declaration. All project files use this namespace.
   Add this declaration (and its closing) to all source/headers you create.
   Note that all #includes should be BEFORE the namespace declaration.
 */
namespace _462 {

// definitions of functions for the GeometryProject class

// constructor, invoked when object is allocated
GeometryProject::GeometryProject() { }

// destructor, invoked when object is de-allocated
GeometryProject::~GeometryProject() { }


void GeometryProject::checkMatch(heEdge* edge1, heEdge* edge2)
{
	if(edge1->startVertex->point.operator == (edge2->nextEdge->startVertex->point)){
		if(edge2->startVertex->point.operator == (edge1->nextEdge->startVertex->point)){
		//	printf("%i FOUND MATCH (%f,%f,%f) -> (%f,%f,%f) == (%f,%f,%f) <- (%f,%f,%f)\n",count,edge1->startVertex->point.x,edge1->startVertex->point.y,edge1->startVertex->point.z,edge1->nextEdge->startVertex->point.x,edge1->nextEdge->startVertex->point.y,edge1->nextEdge->startVertex->point.z,edge2->nextEdge->startVertex->point.x,edge2->nextEdge->startVertex->point.y,edge2->nextEdge->startVertex->point.z,edge2->startVertex->point.x,edge2->startVertex->point.y,edge2->startVertex->point.z);
			count++;
			edge1->pairEdge = edge2;
			edge2->pairEdge = edge1;
		}
	}
}

void GeometryProject::findMatches(heEdge* curEdge, heEdge* compEdge)
{
	//Edge 1 vs. 1,2,3
	heEdge* tmpEdge = compEdge;
	checkMatch(curEdge, tmpEdge);
	tmpEdge = tmpEdge->nextEdge;
	checkMatch(curEdge, tmpEdge);
	tmpEdge = tmpEdge->nextEdge;
	checkMatch(curEdge, tmpEdge);
	//Edge 2 vs. 1,2,3
	curEdge = curEdge->nextEdge;
	tmpEdge = compEdge;
	checkMatch(curEdge, tmpEdge);
	tmpEdge = tmpEdge->nextEdge;
	checkMatch(curEdge, tmpEdge);
	tmpEdge = tmpEdge->nextEdge;
	checkMatch(curEdge, tmpEdge);
	//Edge 3 vs. 1,2,3
	curEdge = curEdge->nextEdge;
	tmpEdge = compEdge;
	checkMatch(curEdge, tmpEdge);
	tmpEdge = tmpEdge->nextEdge;
	checkMatch(curEdge, tmpEdge);
	tmpEdge = tmpEdge->nextEdge;
	checkMatch(curEdge, tmpEdge);
}
/**
 * Initialize the project, doing any necessary opengl initialization.
 * @param camera An already-initialized camera.
 * @param mesh The mesh to be rendered and subdivided.
 * @param texture_filename The filename of the texture to use with the mesh.
 *  Is null if there is no texture data with the mesh or no texture filename
 *  was passed in the arguments, in which case textures should not be used.
 * @return true on success, false on error.
 */
bool GeometryProject::initialize( const Camera* camera, const MeshData* mesh, const char* texture_filename )
{

	count = 0;
    this->mesh = *mesh;
	int i, j, indA, indB, indC;
	heEdge* curEdge;
	heEdge* tmpEdge;

	GLfloat light_position[] = { 10.0, 7.0, 15.0, 1.0 };
	GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };



    glClearColor (0.0, 0.0, 0.0, 0.0);
   	glShadeModel (GL_SMOOTH);
   
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
 	glLightfv(GL_LIGHT0, GL_POSITION, light_position);


   	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);

	//
	//load our half edge data structure from mesh
	//

	this->numFaces = this->mesh.num_triangles;
	faceList = (heFace*) calloc(numFaces,sizeof(heFace));

	for(i=0; i<numFaces; i++){

		//malloc vertices
		heVertex* vertA = (heVertex*) malloc(sizeof(heVertex));
		heVertex* vertB = (heVertex*) malloc(sizeof(heVertex));
		heVertex* vertC = (heVertex*) malloc(sizeof(heVertex));

		//find index into vertex array
		indA = this->mesh.triangles[i].vertices[0];
		indB = this->mesh.triangles[i].vertices[1];
		indC = this->mesh.triangles[i].vertices[2];

		//copy position
		vertA->point = this->mesh.vertices[indA].position;
		vertB->point = this->mesh.vertices[indB].position;
		vertC->point = this->mesh.vertices[indC].position;

		//malloc edges
		heEdge* eAB = (heEdge*) malloc(sizeof(heEdge));
		heEdge* eBC = (heEdge*) malloc(sizeof(heEdge));
		heEdge* eCA = (heEdge*) malloc(sizeof(heEdge));

		//set start vertices
		eAB->startVertex = vertA;
		eBC->startVertex = vertB;
		eCA->startVertex = vertC;

		//set edges for vertices
		vertA->edge = eAB;
		vertB->edge = eBC;
		vertC->edge = eCA;

		//set next edges
		eAB->nextEdge = eBC;
		eBC->nextEdge = eCA;
		eCA->nextEdge = eAB;

		//set face edge
		faceList[i].edge = eAB;

		//set edges face
		eAB->leftFace = &faceList[i];
		eBC->leftFace = &faceList[i];
		eCA->leftFace = &faceList[i];
	}

	//find paired edges
	for(i=0; i<numFaces; i++){
		for(j=0; j<numFaces; j++){
			findMatches(faceList[i].edge, faceList[j].edge);
		}
	}




    return true;
}



/**
 * Clean up the project. Free any memory, etc.
 */
void GeometryProject::destroy()
{
  free(faceList);
}

void GeometryProject::myDraw(){

	Vector3 normA, normB, normC, vertA, vertB, vertC;
	int i, indA, indB, indC, numFaces;

	glBegin( GL_LINE_STRIP); 

	for(i=0;i<this->numFaces;i++){
		heEdge* curEdge = faceList[i].edge;
		Vector3 curVertex = curEdge->startVertex->point;
		glVertex3f( curVertex.x, curVertex.y, curVertex.z );
		curEdge = curEdge->nextEdge;
		curVertex = curEdge->startVertex->point;
		glVertex3f( curVertex.x, curVertex.y, curVertex.z );
		curEdge = curEdge->nextEdge;
		curVertex = curEdge->startVertex->point;
		glVertex3f( curVertex.x, curVertex.y, curVertex.z );
	}

	
	glEnd(); 
	glFlush();
}




/**
 * Clear the screen, then render the mesh using the given camera.
 * @param camera The logical camera to use.
 * @see scene/camera.hpp
 */
void GeometryProject::render( const Camera* camera )
{

	Vector3 vPos, vAt, vCen, vUp;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode( GL_PROJECTION ); // set current matrix 
	glLoadIdentity(); // Clears the matrix 

	//set perspective and camera
	gluPerspective(camera->get_fov_degrees(),camera->get_aspect_ratio(),camera->get_near_clip(),camera->get_far_clip());
	vPos = camera->get_position();
	vAt = camera->get_direction();
	vCen = vPos + vAt;
	vUp = camera->get_up();
	gluLookAt(vPos.x, vPos.y, vPos.z, vCen.x, vCen.y, vCen.z, vUp.x, vUp.y,  vUp.z);

    glMatrixMode(GL_MODELVIEW);
   	glLoadIdentity();
	myDraw();


  	
}


/**
 * Subdivide the mesh that we are rendering using Loop subdivison.
 */
void GeometryProject::subdivide()
{
	int i,j;
	heEdge* curEdge;
	Vector3 ab, bc, ca, aPrime, bPrime, cPrime;
	heFace* tmpFaceList = (heFace*) malloc(sizeof(heFace) * 4* numFaces);
	for(i=0; i<numFaces; i++){
		//
		//calculate odd vertices
		//

		curEdge = faceList[i].edge;
		if(curEdge->pairEdge != NULL){
			ab = splitInteriorEdge(curEdge);
			aPrime = computeEvenVertex(curEdge);
		}
		else{
			printf("FOUND A BOUNDARY VERTEX!!\n");
		}
		curEdge = curEdge->nextEdge;
		if(curEdge->pairEdge != NULL){
			bc = splitInteriorEdge(curEdge);
			bPrime = computeEvenVertex(curEdge);
		}
		else{
			printf("FOUND A BOUNDARY VERTEX!!\n");
		}
		curEdge = curEdge->nextEdge;
		if(curEdge->pairEdge != NULL){
			ca = splitInteriorEdge(curEdge);
			cPrime = computeEvenVertex(curEdge);
		}
		else{
			printf("FOUND A BOUNDARY VERTEX!!\n");
		}
		
		updateFaceList(tmpFaceList, i, aPrime, bPrime, cPrime, ab, bc, ca);

	}

	freeFaceList(faceList);
	faceList = tmpFaceList;
	numFaces*=4;

	//find paired edges
	for(i=0; i<numFaces; i++){
		for(j=0; j<numFaces; j++){
			findMatches(faceList[i].edge, faceList[j].edge);
		}
	}
	
}

void GeometryProject::freeFaceList(heFace* faceList){
	
	int i;
	heEdge* bc;
	heEdge* ca;

	for(i=0; i<numFaces; i++){
		free(faceList[i].edge->startVertex);
		bc = faceList[i].edge->nextEdge;
		free(bc->startVertex);
		ca = bc->nextEdge;
		free(ca->startVertex);
		free(bc);
		free(ca);
		free(faceList[i].edge);
	}

	free(faceList);

}

void GeometryProject::createFace(heFace* tmpFaceList, Vector3 a, Vector3 b, Vector3 c){

	//malloc vertices
	heVertex* vertA = (heVertex*) malloc(sizeof(heVertex));
	heVertex* vertB = (heVertex*) malloc(sizeof(heVertex));
	heVertex* vertC = (heVertex*) malloc(sizeof(heVertex));

	//copy position
	vertA->point = a;
	vertB->point = b;
	vertC->point = c;

	//malloc edges
	heEdge* eAB = (heEdge*) malloc(sizeof(heEdge));
	heEdge* eBC = (heEdge*) malloc(sizeof(heEdge));
	heEdge* eCA = (heEdge*) malloc(sizeof(heEdge));

	//set start vertices
	eAB->startVertex = vertA;
	eBC->startVertex = vertB;
	eCA->startVertex = vertC;

	//set edges for vertices
	vertA->edge = eAB;
	vertB->edge = eBC;
	vertC->edge = eCA;

	//set next edges
	eAB->nextEdge = eBC;
	eBC->nextEdge = eCA;
	eCA->nextEdge = eAB;

	//set face edge
	tmpFaceList->edge = eAB;

	//set edges face
	eAB->leftFace = tmpFaceList;
	eBC->leftFace = tmpFaceList;
	eCA->leftFace = tmpFaceList;

}

void GeometryProject::updateFaceList(heFace* tmpFaceList, int i, Vector3 aPrime, Vector3 bPrime, Vector3 cPrime, Vector3 ab, Vector3 bc, Vector3 ca){

	createFace(&tmpFaceList[4*i], aPrime, ab, ca);
	createFace(&tmpFaceList[4*i + 1], ca, bc, cPrime);
	createFace(&tmpFaceList[4*i + 2], ab, bc, ca);
	createFace(&tmpFaceList[4*i + 3], ab, bPrime, bc);

}

Vector3 GeometryProject::computeEvenVertex(heEdge* startEdge){
	int numNeighbors = 0;
	Vector3 neighbors = Vector3::Zero; 
	float math = 0;
	float pi = std::atan(1.0) * 4;
	heEdge* curEdge;
	//
	//calculate even vertices
	//

	//sum neighbors
	curEdge = startEdge->nextEdge;
	neighbors += curEdge->startVertex->point;
	numNeighbors++;
	curEdge = curEdge->nextEdge;
	neighbors += curEdge->startVertex->point;
	numNeighbors++;
	if(curEdge->pairEdge != NULL){
		curEdge = curEdge->pairEdge->nextEdge->nextEdge;
		while(curEdge->nextEdge->startVertex->point.operator ==(startEdge->startVertex->point)){
			if(curEdge->pairEdge != NULL){
				if(curEdge->pairEdge != startEdge){
					neighbors += curEdge->startVertex->point;
					numNeighbors++;
					curEdge = curEdge->pairEdge->nextEdge->nextEdge;
				}
				else{
					break;
				}
			}
			else{
				neighbors += curEdge->startVertex->point;
				numNeighbors++;
				break;
			}
		}
	}
	math = cos(2 * pi / numNeighbors);
	math/=4;
	math+=3.0/8.0;
	math = (5.0/8.0) - pow(math,2);
	math/=numNeighbors;
	neighbors*=math;
	math*=numNeighbors;
	math = 1.0 - math;
	neighbors+=math * startEdge->startVertex->point;
	return neighbors;
}

Vector3 GeometryProject::splitInteriorEdge(heEdge* edge){
	Vector3 a,b,c,d,newVertex;

	a = edge->startVertex->point;
	b = edge->nextEdge->startVertex->point;
	c = edge->nextEdge->nextEdge->startVertex->point;
	d = edge->pairEdge->nextEdge->nextEdge->startVertex->point;

	a = 3 * a / 8;
	b = 3 * b / 8;
	c = c / 8;
	d = d / 8;

	newVertex = a + b + c + d;

	return newVertex;
}



} /* _462 */

