//Library to create and calculate Local Matrix

//Function that returns the node based on some parameters
node selectNode(int i, element e,mesh &m){
	node n;
	switch(i){
		case 1: n = m.getNode(e.getNode1()-1); break;
		case 2: n = m.getNode(e.getNode2()-1); break;
		case 3: n = m.getNode(e.getNode3()-1); break;
        case 4: n = m.getNode(e.getNode4()-1); break;
	}
	return n;
}

//Function that returns the coord based on some parameters
float selectCoord(int c, node n){
	float v;
	switch(c){
		case EQUIS: v = n.getX(); break;
		case YE: v = n.getY(); break;
        case ZETA: v = n.getZ(); break;
	}
	return v;
}

//Function that help us calculate operations like X2-X1 represented as X21
float calcularTenedor(element e, int coord, int i, int j,mesh &m){
	node n1=selectNode(i,e,m),n2=selectNode(j,e,m);

	return selectCoord(coord,n1) - selectCoord(coord,n2);
}

//Function that help us calculate terms like Y3Z4 that translates to Y31Z41 - Y41Z31
float RestaTenedor(float a, float b, float c, float d){
    return (a*b)-(c*d);
}

//Function that returns the magnitude between two coords
float calculateMagnitude(float v1, float v2){
    return sqrt(pow(v1,2)+pow(v2,2));
}

//Functions that returns the Volumme of a local element
float calculateLocalVolume(int i,mesh m){
    float V,A,s,a,b,c;
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    a = calculateMagnitude(n2.getX()-n1.getX(),n2.getY()-n1.getY());
    b = calculateMagnitude(n3.getX()-n2.getX(),n3.getY()-n2.getY());
    c = calculateMagnitude(n3.getX()-n1.getX(),n3.getY()-n1.getY());
    s = (a+b+c)/2;

    A = sqrt(s*(s-a)*(s-b)*(s-c));
    V = (sqrt(2)*A*A*A)/12;
    return V;
}

//CALCULATE ALL LOCALS MATRIX************************************************************************************************

