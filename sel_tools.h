//Library with tools that helps defining each component.

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
float OperarRestaTenedor(element e, int coord1, int coord2,  float value_a, float value_b, mesh &m){
    float a, b, c, d;

    a = calcularTenedor(e, coord1, value_a, 1 ,m);
    b = calcularTenedor(e, coord2, value_b, 1 ,m);
    c = calcularTenedor(e, coord1, value_b, 1 ,m);
    d = calcularTenedor(e, coord2, value_a, 1 ,m);

    return (a*b)-(c*d);
}

//Function Ubicate Sub Matrix on the local Matrix M
void ubicarSubMatriz(Matrix &K,int fi,int ff,int ci,int cf,Matrix M){
    int n = 0, m= 0;
    for(int i=fi;i<=ff;i++){
        for(int j=ci;j<=cf;j++){
            K.at(i).at(j) = M.at(n).at(m);
            m++;
        }
        n++; m = 0;
    }
}