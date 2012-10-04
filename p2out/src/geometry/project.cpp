/**
 * @file project.cpp
 * @brief Geometry project
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#include "geometry/project.hpp"
#include "application/opengl.hpp"
#include "application/imageio.hpp"

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

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_NORMALIZE);

	GLfloat light0_position[] = { 2.0, 2.0, 0.0, 1.0 };
	GLfloat light1_position[] = { -10.0, -15.0, 10.0, 1.0 };
	GLfloat lmodel_ambient[] = { 0.5, 0.5, 0.5, 1.0 };



    glClearColor (0.0, 0.0, 0.0, 0.0);
   	glShadeModel (GL_SMOOTH);
   
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
 	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	 glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

	 hasTexture = false;

	 if(texture_filename != NULL){
		 printf("FOUND FILENAME %s\n",texture_filename);
		 int height;
		 int width;
		 texture = _462::imageio_load_image(texture_filename,&width,&height);
		 if((width != -1) && (height != -1)){
			 printf("Loaded Image!\n");
			 hasTexture = true;
			 glGenTextures(1,&hTex);
			 glBindTexture(GL_TEXTURE_2D, hTex);
			 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
			 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
			 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, 
                   GL_NEAREST);
			 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, 
                   GL_NEAREST);
			 glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, 
                height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 
                texture);
		 }
	 }


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

		//copy normal
		vertA->normal = this->mesh.vertices[indA].normal;
		vertB->normal = this->mesh.vertices[indB].normal;
		vertC->normal = this->mesh.vertices[indC].normal;
		
		if(hasTexture){
		//copy texture
			vertA->texture = this->mesh.vertices[indA].texture_coord;
			vertB->texture = this->mesh.vertices[indB].texture_coord;
			vertC->texture = this->mesh.vertices[indC].texture_coord;
		}

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

		//set edge pairs to null
		eAB->pairEdge = NULL;
		eBC->pairEdge = NULL;
		eCA->pairEdge = NULL;
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

	Vector3 curVertex, curNormal;
	Vector2 curTexture;
	int i, indA, indB, indC, numFaces;

	glEnable(GL_TEXTURE_2D);
	glBegin( GL_TRIANGLES); 

	for(i=0;i<this->numFaces;i++){
		heEdge* curEdge = faceList[i].edge;
		curVertex = curEdge->startVertex->point;
		curNormal = curEdge->startVertex->normal;
		if(hasTexture){
			curTexture = curEdge->startVertex->texture;
			glTexCoord2f(curTexture.x,curTexture.y);
		}
		glNormal3f( curNormal.x, curNormal.y, curNormal.z);
		glVertex3f( curVertex.x, curVertex.y, curVertex.z );
		curEdge = curEdge->nextEdge;
		curVertex = curEdge->startVertex->point;
		curNormal = curEdge->startVertex->normal;
		if(hasTexture){
			curTexture = curEdge->startVertex->texture;
			glTexCoord2f(curTexture.x,curTexture.y);
		}
		glNormal3f( curNormal.x, curNormal.y, curNormal.z);
		glVertex3f( curVertex.x, curVertex.y, curVertex.z );
		curEdge = curEdge->nextEdge;
		curVertex = curEdge->startVertex->point;
		curNormal = curEdge->startVertex->normal;
		if(hasTexture){
			curTexture = curEdge->startVertex->texture;
			glTexCoord2f(curTexture.x,curTexture.y);
		}
		glNormal3f( curNormal.x, curNormal.y, curNormal.z);
		glVertex3f( curVertex.x, curVertex.y, curVertex.z );
	}

	
	glEnd(); 
	glDisable(GL_TEXTURE_2D);
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
	heEdge* tempEdge;
	heEdge* compEdge;
	heFace* tmpFaceList = (heFace*) malloc(sizeof(heFace) * 4* numFaces);
	nvtPair ab, bc, ca, aPrime, bPrime, cPrime;

	for(i=0; i<numFaces; i++){
		//
		//calculate odd vertices
		//

		curEdge = faceList[i].edge;
		if(curEdge->pairEdge != NULL){
			// if we have an interior edge, simply split
			ab = splitInteriorEdge(curEdge);
			tempEdge = curEdge->pairEdge->nextEdge;
			//check if our startVertex is at a boundary
			while(tempEdge->pairEdge != NULL){
				// if we return to our starting place, not a boundary
				if(tempEdge == curEdge){
					aPrime = computeEvenInteriorVertex(curEdge);
					break;
				}
				tempEdge = tempEdge->pairEdge->nextEdge;
			}
			if(tempEdge->pairEdge == NULL){
			aPrime = computeEvenBoundaryVertex(tempEdge);
			}
		}
		else{
			ab = splitBoundaryEdge(curEdge);
			aPrime = computeEvenBoundaryVertex(curEdge);
		}
		curEdge = curEdge->nextEdge;
			if(curEdge->pairEdge != NULL){
			// if we have an interior edge, simply split
			bc = splitInteriorEdge(curEdge);
			tempEdge = curEdge->pairEdge->nextEdge;
			//check if our startVertex is at a boundary
			while(tempEdge->pairEdge != NULL){
				// if we return to our starting place, not a boundary
				if(tempEdge == curEdge){
					bPrime = computeEvenInteriorVertex(curEdge);
					break;
				}
				tempEdge = tempEdge->pairEdge->nextEdge;
			}
			if(tempEdge->pairEdge == NULL){
			bPrime = computeEvenBoundaryVertex(tempEdge);
			}
		}
		else{
			bc = splitBoundaryEdge(curEdge);
			bPrime = computeEvenBoundaryVertex(curEdge);
		}
		curEdge = curEdge->nextEdge;
		if(curEdge->pairEdge != NULL){
			// if we have an interior edge, simply split
			ca = splitInteriorEdge(curEdge);
			tempEdge = curEdge->pairEdge->nextEdge;
			//check if our startVertex is at a boundary
			while(tempEdge->pairEdge != NULL){
				// if we return to our starting place, not a boundary
				if(tempEdge == curEdge){
					cPrime = computeEvenInteriorVertex(curEdge);
					break;
				}
				tempEdge = tempEdge->pairEdge->nextEdge;
			}
			if(tempEdge->pairEdge == NULL){
			cPrime = computeEvenBoundaryVertex(tempEdge);
			}
		}
		else{
			ca = splitBoundaryEdge(curEdge);
			cPrime = computeEvenBoundaryVertex(curEdge);
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

Vector3 GeometryProject::round(Vector3 vector){

	Vector3 rounded;

	rounded.x = ceil(vector.x * 10000.0) / 10000.0;
	rounded.y = ceil(vector.y * 10000.0) / 10000.0;
	rounded.z = ceil(vector.z * 10000.0) / 10000.0;

	return rounded;

}

void GeometryProject::createFace(heFace* tmpFaceList, nvtPair a, nvtPair b, nvtPair c){

	//malloc vertices
	heVertex* vertA = (heVertex*) malloc(sizeof(heVertex));
	heVertex* vertB = (heVertex*) malloc(sizeof(heVertex));
	heVertex* vertC = (heVertex*) malloc(sizeof(heVertex));

	//copy position with slight rounding
	// otherwise paired edges are slightly off
	// and holes open up in the surface
	vertA->point = round(a.vertex);
	vertB->point = round(b.vertex);
	vertC->point = round(c.vertex);

	//copy normal
	vertA->normal = a.normal;
	vertB->normal = b.normal;
	vertC->normal = c.normal;

	//copy textures
	if(hasTexture){
		vertA->texture = a.texture;
		vertB->texture = b.texture;
		vertC->texture = c.texture;
	}

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

	//set edge pairs to null
	eAB->pairEdge = NULL;
	eBC->pairEdge = NULL;
	eCA->pairEdge = NULL;

}

void GeometryProject::updateFaceList(heFace* tmpFaceList, int i, nvtPair aPrime, nvtPair bPrime, nvtPair cPrime, nvtPair ab, nvtPair bc, nvtPair ca){

	/*createFace(&tmpFaceList[4*i], ca, aPrime, ab);
	createFace(&tmpFaceList[4*i + 1], ca, bc, cPrime);
	createFace(&tmpFaceList[4*i + 2], ca, ab, bc);
	createFace(&tmpFaceList[4*i + 3], ab, bPrime, bc);*/

	createFace(&tmpFaceList[4*i], aPrime, ab, ca);
	createFace(&tmpFaceList[4*i + 1], ca, bc, cPrime);
	createFace(&tmpFaceList[4*i + 2], ab, bc, ca);
	createFace(&tmpFaceList[4*i + 3], ab, bPrime, bc);

}

nvtPair GeometryProject::computeEvenBoundaryVertex(heEdge* startEdge){

	Vector3 a,b,newVertex, aNorm, bNorm, newNormVertex;
	Vector2 aTex, bTex, newTex;
	heEdge* temp;
	nvtPair result;

	a = startEdge->nextEdge->startVertex->point;
	aNorm = startEdge->nextEdge->startVertex->normal;
	if(hasTexture){
		aTex = startEdge->nextEdge->startVertex->texture;
	}
	temp = startEdge->nextEdge->nextEdge;
	while(temp->pairEdge != NULL){
		temp = temp->pairEdge->nextEdge->nextEdge;
	}
	b = temp->startVertex->point;
	bNorm = temp->startVertex->normal;
	if(hasTexture){
		bTex = temp->startVertex->texture;
	}

	newVertex = startEdge->startVertex->point;
	newVertex = 3.0 * newVertex / 4.0;
	a = 1.0 * a / 8.0;
	b = 1.0 * b / 8.0;
	newVertex = newVertex + a + b;
	
	newNormVertex = startEdge->startVertex->normal;
	newNormVertex = 3.0 * newNormVertex / 4.0;
	aNorm = 1.0 * aNorm / 8.0;
	bNorm = 1.0 * bNorm / 8.0;
	newNormVertex = newNormVertex + aNorm + bNorm;

	if(hasTexture){
		newTex = startEdge->startVertex->texture;
		newTex = 3.0 * newTex / 4.0;
		aTex = 1.0 * aTex / 8.0;
		bTex = 1.0 * bTex / 8.0;
		newTex = aTex + bTex + newTex;
		result.texture = newTex;
	}

	result.normal = newNormVertex;
	result.vertex = newVertex;
	
	return result;
}

nvtPair GeometryProject::computeEvenInteriorVertex(heEdge* startEdge){
	int numNeighbors = 0;
	Vector3 neighbors = Vector3::Zero; 
	Vector3 normNeighbors = Vector3::Zero;
	Vector2 newTex = Vector2::Zero;
	float math = 0;
	float normMath = 0;
	float texMath = 0;
	float pi = std::atan(1.0) * 4;
	heEdge* curEdge;
	nvtPair result;

	//sum neighbors
	curEdge = startEdge->nextEdge;
	neighbors += curEdge->startVertex->point;
	normNeighbors += curEdge->startVertex->normal;
	if(hasTexture){
		newTex += curEdge->startVertex->texture;
	}
	numNeighbors++;
	curEdge = curEdge->nextEdge;
	neighbors += curEdge->startVertex->point;
	normNeighbors += curEdge->startVertex->normal;
	if(hasTexture){
		newTex += curEdge->startVertex->texture;
	}
	numNeighbors++;
	curEdge = curEdge->pairEdge->nextEdge->nextEdge;

		while(curEdge->nextEdge->startVertex->point.operator ==(startEdge->startVertex->point)){
			if(curEdge->pairEdge != startEdge){
				neighbors += curEdge->startVertex->point;
				normNeighbors += curEdge->startVertex->normal;
				if(hasTexture){
					newTex += curEdge->startVertex->texture;
				}
				numNeighbors++;
				curEdge = curEdge->pairEdge->nextEdge->nextEdge;
			}
			else{
				break;
			}
		}
	
	math = cos(2 * pi / numNeighbors);
	math/=4.0;
	math+=3.0/8.0;
	math = (5.0/8.0) - pow(math,2);
	math/=numNeighbors;
	normMath = math;
	if(hasTexture){
		texMath = math;
	}
	neighbors*=math;
	normNeighbors*=normMath;
	math*=numNeighbors;
	normMath*=numNeighbors;
	math = 1.0 - math;
	normMath = 1.0 - normMath;
	neighbors+=math * startEdge->startVertex->point;
	normNeighbors+=normMath * startEdge->startVertex->normal;

	if(hasTexture){
		newTex*=texMath;
		texMath*=numNeighbors;
		texMath = 1.0 - texMath;
		newTex+=texMath * startEdge->startVertex->texture;
		result.texture = newTex;
	}

	result.normal = normNeighbors;
	result.vertex = neighbors;

	return result;
}

nvtPair GeometryProject::splitBoundaryEdge(heEdge* edge){
	Vector3 a,b,newVertex, aNorm, bNorm, newNormVertex;
	Vector2 newTex, aTex, bTex;
	nvtPair result;

	a = edge->startVertex->point;
	b = edge->nextEdge->startVertex->point;

	aNorm = edge->startVertex->normal;
	bNorm = edge->nextEdge->startVertex->normal;

	a = 1.0 * a / 2.0;
	b = 1.0 * b / 2.0;

	aNorm = 1.0 * aNorm / 2.0;
	bNorm = 1.0 * bNorm / 2.0;

	newVertex = a + b;
	newNormVertex = aNorm + bNorm;

	if(hasTexture){
		aTex = edge->startVertex->texture;
		bTex = edge->startVertex->texture;
		aTex = 1.0 * aTex / 2.0;
		bTex = 1.0 * bTex / 2.0;
		newTex = aTex + bTex;
		result.texture = newTex;
	}

	result.vertex = newVertex;
	result.normal = newNormVertex;

	return result;

}

nvtPair GeometryProject::splitInteriorEdge(heEdge* edge){
	Vector3 a,b,c,d,newVertex, aNorm, bNorm, cNorm, dNorm, newNormVertex;
	Vector2 aTex, bTex, cTex, dTex, newTex;
	nvtPair result;

	a = edge->startVertex->point;
	b = edge->nextEdge->startVertex->point;
	c = edge->nextEdge->nextEdge->startVertex->point;
	d = edge->pairEdge->nextEdge->nextEdge->startVertex->point;

	aNorm = edge->startVertex->normal;
	bNorm = edge->nextEdge->startVertex->normal;
	cNorm = edge->nextEdge->nextEdge->startVertex->normal;
	dNorm = edge->pairEdge->nextEdge->nextEdge->startVertex->normal;

	a = 3.0 * a / 8.0;
	b = 3.0 * b / 8.0;
	c = 1.0 * c / 8.0;
	d = 1.0 * d / 8.0;
	
	aNorm = 3.0 * aNorm / 8.0;
	bNorm = 3.0 * bNorm / 8.0;
	cNorm = 1.0 * cNorm / 8.0;
	dNorm = 1.0 * dNorm / 8.0;

	if(hasTexture){
		aTex = edge->startVertex->texture;
		bTex = edge->nextEdge->startVertex->texture;
		cTex = edge->nextEdge->nextEdge->startVertex->texture;
		dTex = edge->pairEdge->nextEdge->nextEdge->startVertex->texture;

		aTex = 3.0 * aTex / 8.0;
		bTex = 3.0 * bTex / 8.0;
		cTex = 1.0 * cTex / 8.0;
		dTex = 1.0 * dTex / 8.0;

		newTex = aTex + bTex + cTex + dTex;
		result.texture = newTex;
	}

	newVertex = a + b + c + d;
	newNormVertex = aNorm + bNorm + cNorm + dNorm;

	result.normal = newNormVertex;
	result.vertex = newVertex;

	return result;
}



} /* _462 */

