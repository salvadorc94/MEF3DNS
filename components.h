//All the local components needed to create the local Matrix M and local B
//Determinant
float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}


void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

//Jacobiano
float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

//AKA B
void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}

//Special Matrix because of the project
void calculateGamma(int i,Matrix &Gamma,mesh m){
    zeroes(Gamma,12,3);
    element e = m.getElement(i);
    float cord1 = selectCoord(EQUIS,selectNode(1,e,m));
    float cord2 = selectCoord(EQUIS,selectNode(2,e,m));
    float cord3 = selectCoord(EQUIS,selectNode(3,e,m));
    float cord4 = selectCoord(EQUIS,selectNode(4,e,m));

    Gamma.at(0).at(0) = (2*cord1) + cord2 + cord3 + cord4;
    Gamma.at(0).at(1) = 0;
    Gamma.at(0).at(2) = 0;
    Gamma.at(1).at(0) = cord1 + (2*cord2) + cord3 + cord4;
    Gamma.at(1).at(1) = 0;
    Gamma.at(1).at(2) = 0;
    Gamma.at(2).at(0) = cord1 + cord2 + (2*cord3) + cord4;
    Gamma.at(2).at(1) = 0;
    Gamma.at(2).at(2) = 0;
    Gamma.at(3).at(0) = cord1 + cord2 + cord3 + (2*cord4);
    Gamma.at(3).at(1) = 0;
    Gamma.at(3).at(2) = 0;

    Gamma.at(4).at(0) = 0;
    Gamma.at(4).at(1) = (2*cord1) + cord2 + cord3 + cord4;
    Gamma.at(4).at(2) = 0;
    Gamma.at(5).at(0) = 0;
    Gamma.at(5).at(1) = cord1 + (2*cord2) + cord3 + cord4;
    Gamma.at(5).at(2) = 0;
    Gamma.at(6).at(0) = 0;
    Gamma.at(6).at(1) = cord1 + cord2 + (2*cord3) + cord4;
    Gamma.at(6).at(2) = 0;
    Gamma.at(7).at(0) = 0;
    Gamma.at(7).at(1) = cord1 + cord2 + cord3 + (2*cord4);
    Gamma.at(7).at(2) = 0;

    Gamma.at(8).at(0) = 0;
    Gamma.at(8).at(1) = 0;
    Gamma.at(8).at(2) = (2*cord1) + cord2 + cord3 + cord4;
    Gamma.at(9).at(0) = 0;
    Gamma.at(9).at(1) = 0;
    Gamma.at(9).at(2) = cord1 + (2*cord2) + cord3 + cord4;
    Gamma.at(10).at(0) = 0;
    Gamma.at(10).at(1) = 0;
    Gamma.at(10).at(2) = cord1 + cord2 + (2*cord3) + cord4;
    Gamma.at(11).at(0) = 0;
    Gamma.at(11).at(1) = 0;
    Gamma.at(11).at(2) = cord1 + cord2 + cord3 + (2*cord4);
}

void calculateOmmega(int i,Matrix &Omega,mesh m){
    zeroes(Omega,12,3);
    element e = m.getElement(i);
    float cord1 = selectCoord(YE,selectNode(1,e,m));
    float cord2 = selectCoord(YE,selectNode(2,e,m));
    float cord3 = selectCoord(YE,selectNode(3,e,m));
    float cord4 = selectCoord(YE,selectNode(4,e,m));

    Omega.at(0).at(0) = 3*pow(cord1,2) + 2*cord1*(cord2+cord3+cord4) + pow(cord2,2) + cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Omega.at(0).at(1) = 0;
    Omega.at(0).at(2) = 0;
    Omega.at(1).at(0) = pow(cord1,2) + cord1*(2*cord2+cord3+cord4) + 3*pow(cord2,2) + 2*cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Omega.at(1).at(1) = 0;
    Omega.at(1).at(2) = 0;
    Omega.at(2).at(0) = pow(cord1,2) + cord1*(cord2+2*cord3+cord4) + pow(cord2,2) + cord2*(2*cord3+cord4) + 3*pow(cord3,2) + 2*cord3*cord4 + pow(cord4,2);
    Omega.at(2).at(1) = 0;
    Omega.at(2).at(2) = 0;
    Omega.at(3).at(0) = pow(cord1,2) + cord1*(cord2+cord3+2*cord4) + pow(cord2,2) + cord2*(cord3+2*cord4) + pow(cord3,2) + 2*cord3*cord4 + 3*pow(cord4,2);
    Omega.at(3).at(1) = 0;
    Omega.at(3).at(2) = 0;

    Omega.at(4).at(0) = 0;
    Omega.at(4).at(1) = 3*pow(cord1,2) + 2*cord1*(cord2+cord3+cord4) + pow(cord2,2) + cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Omega.at(4).at(2) = 0;
    Omega.at(5).at(0) = 0;
    Omega.at(5).at(1) = pow(cord1,2) + cord1*(2*cord2+cord3+cord4) + 3*pow(cord2,2) + 2*cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Omega.at(5).at(2) = 0;
    Omega.at(6).at(0) = 0;
    Omega.at(6).at(1) = pow(cord1,2) + cord1*(cord2+2*cord3+cord4) + pow(cord2,2) + cord2*(2*cord3+cord4) + 3*pow(cord3,2) + 2*cord3*cord4 + pow(cord4,2);
    Omega.at(6).at(2) = 0;
    Omega.at(7).at(0) = 0;
    Omega.at(7).at(1) = pow(cord1,2) + cord1*(cord2+cord3+2*cord4) + pow(cord2,2) + cord2*(cord3+2*cord4) + pow(cord3,2) + 2*cord3*cord4 + 3*pow(cord4,2);
    Omega.at(7).at(2) = 0;

    Omega.at(8).at(0) = 0;
    Omega.at(8).at(1) = 0;
    Omega.at(8).at(2) = 3*pow(cord1,2) + 2*cord1*(cord2+cord3+cord4) + pow(cord2,2) + cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Omega.at(9).at(0) = 0;
    Omega.at(9).at(1) = 0;
    Omega.at(9).at(2) = pow(cord1,2) + cord1*(2*cord2+cord3+cord4) + 3*pow(cord2,2) + 2*cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Omega.at(10).at(0) = 0;
    Omega.at(10).at(1) = 0;
    Omega.at(10).at(2) = pow(cord1,2) + cord1*(cord2+2*cord3+cord4) + pow(cord2,2) + cord2*(2*cord3+cord4) + 3*pow(cord3,2) + 2*cord3*cord4 + pow(cord4,2);
    Omega.at(11).at(0) = 0;
    Omega.at(11).at(1) = 0;
    Omega.at(11).at(2) = pow(cord1,2) + cord1*(cord2+cord3+2*cord4) + pow(cord2,2) + cord2*(cord3+2*cord4) + pow(cord3,2) + 2*cord3*cord4 + 3*pow(cord4,2);

}

void calculateSigma(int i,Matrix &Sigma,mesh m){
    zeroes(Sigma,3,12);
    element e = m.getElement(i);
    float cord1 = selectCoord(EQUIS,selectNode(1,e,m));
    float cord2 = selectCoord(EQUIS,selectNode(2,e,m));
    float cord3 = selectCoord(EQUIS,selectNode(3,e,m));
    float cord4 = selectCoord(EQUIS,selectNode(4,e,m));

    Sigma.at(0).at(0) = 3*pow(cord1,2) + 2*cord1*(cord2+cord3+cord4) + pow(cord2,2) + cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Sigma.at(0).at(1) = pow(cord1,2) + cord1*(2*cord2+cord3+cord4) + 3*pow(cord2,2) + 2*cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Sigma.at(0).at(2) = pow(cord1,2) + cord1*(cord2+2*cord3+cord4) + pow(cord2,2) + cord2*(2*cord3+cord4) + 3*pow(cord3,2) + 2*cord3*cord4 + pow(cord4,2);
    Sigma.at(0).at(3) = pow(cord1,2) + cord1*(cord2+cord3+2*cord4) + pow(cord2,2) + cord2*(cord3+2*cord4) + pow(cord3,2) + 2*cord3*cord4 + 3*pow(cord4,2);
    Sigma.at(0).at(4) = 0;
    Sigma.at(0).at(5) = 0;
    Sigma.at(0).at(6) = 0;
    Sigma.at(0).at(7) = 0;
    Sigma.at(0).at(8) = 0;
    Sigma.at(0).at(9) = 0;
    Sigma.at(0).at(10) = 0;
    Sigma.at(0).at(11) = 0;

    Sigma.at(1).at(0) = 0;
    Sigma.at(1).at(1) = 0;
    Sigma.at(1).at(2) = 0;
    Sigma.at(1).at(3) = 0;
    Sigma.at(1).at(4) = 3*pow(cord1,2) + 2*cord1*(cord2+cord3+cord4) + pow(cord2,2) + cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Sigma.at(1).at(5) = pow(cord1,2) + cord1*(2*cord2+cord3+cord4) + 3*pow(cord2,2) + 2*cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Sigma.at(1).at(6) = pow(cord1,2) + cord1*(cord2+2*cord3+cord4) + pow(cord2,2) + cord2*(2*cord3+cord4) + 3*pow(cord3,2) + 2*cord3*cord4 + pow(cord4,2);
    Sigma.at(1).at(7) = pow(cord1,2) + cord1*(cord2+cord3+2*cord4) + pow(cord2,2) + cord2*(cord3+2*cord4) + pow(cord3,2) + 2*cord3*cord4 + 3*pow(cord4,2);
    Sigma.at(1).at(8) = 0;
    Sigma.at(1).at(9) = 0;
    Sigma.at(1).at(10) = 0;
    Sigma.at(1).at(11) = 0;

    Sigma.at(2).at(0) = 0;
    Sigma.at(2).at(1) = 0;
    Sigma.at(2).at(2) = 0;
    Sigma.at(2).at(3) = 0;
    Sigma.at(2).at(4) = 0;
    Sigma.at(2).at(5) = 0;
    Sigma.at(2).at(6) = 0;
    Sigma.at(2).at(7) = 0;
    Sigma.at(2).at(8) = 3*pow(cord1,2) + 2*cord1*(cord2+cord3+cord4) + pow(cord2,2) + cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Sigma.at(2).at(9) = pow(cord1,2) + cord1*(2*cord2+cord3+cord4) + 3*pow(cord2,2) + 2*cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Sigma.at(2).at(10) = pow(cord1,2) + cord1*(cord2+2*cord3+cord4) + pow(cord2,2) + cord2*(2*cord3+cord4) + 3*pow(cord3,2) + 2*cord3*cord4 + pow(cord4,2);
    Sigma.at(2).at(11) = pow(cord1,2) + cord1*(cord2+cord3+2*cord4) + pow(cord2,2) + cord2*(cord3+2*cord4) + pow(cord3,2) + 2*cord3*cord4 + 3*pow(cord4,2);

}

void calculateLambda(int i,Matrix &Lambda,mesh m){
    zeroes(Lambda,12,3);
    element e = m.getElement(i);
    float cord1 = selectCoord(ZETA,selectNode(1,e,m));
    float cord2 = selectCoord(ZETA,selectNode(2,e,m));
    float cord3 = selectCoord(ZETA,selectNode(3,e,m));
    float cord4 = selectCoord(ZETA,selectNode(4,e,m));

    Lambda.at(0).at(0) = (2*cord1) + cord2 + cord3 + cord4;
    Lambda.at(0).at(1) = 0;
    Lambda.at(0).at(2) = 0;
    Lambda.at(1).at(0) = cord1 + (2*cord2) + cord3 + cord4;
    Lambda.at(1).at(1) = 0;
    Lambda.at(1).at(2) = 0;
    Lambda.at(2).at(0) = cord1 + cord2 + (2*cord3) + cord4;
    Lambda.at(2).at(1) = 0;
    Lambda.at(2).at(2) = 0;
    Lambda.at(3).at(0) = cord1 + cord2 + cord3 + (2*cord4);
    Lambda.at(3).at(1) = 0;
    Lambda.at(3).at(2) = 0;

    Lambda.at(4).at(0) = 0;
    Lambda.at(4).at(1) = (2*cord1) + cord2 + cord3 + cord4;
    Lambda.at(4).at(2) = 0;
    Lambda.at(5).at(0) = 0;
    Lambda.at(5).at(1) = cord1 + (2*cord2) + cord3 + cord4;
    Lambda.at(5).at(2) = 0;
    Lambda.at(6).at(0) = 0;
    Lambda.at(6).at(1) = cord1 + cord2 + (2*cord3) + cord4;
    Lambda.at(6).at(2) = 0;
    Lambda.at(7).at(0) = 0;
    Lambda.at(7).at(1) = cord1 + cord2 + cord3 + (2*cord4);
    Lambda.at(7).at(2) = 0;

    Lambda.at(8).at(0) = 0;
    Lambda.at(8).at(1) = 0;
    Lambda.at(8).at(2) = (2*cord1) + cord2 + cord3 + cord4;
    Lambda.at(9).at(0) = 0;
    Lambda.at(9).at(1) = 0;
    Lambda.at(9).at(2) = cord1 + (2*cord2) + cord3 + cord4;
    Lambda.at(10).at(0) = 0;
    Lambda.at(10).at(1) = 0;
    Lambda.at(10).at(2) = cord1 + cord2 + (2*cord3) + cord4;
    Lambda.at(11).at(0) = 0;
    Lambda.at(11).at(1) = 0;
    Lambda.at(11).at(2) = cord1 + cord2 + cord3 + (2*cord4);

}

void calculateTau(int i,Vector &Tau,mesh m){
    zeroes(Tau,4);
    element e = m.getElement(i);
    float cord1 = selectCoord(YE,selectNode(1,e,m));
    float cord2 = selectCoord(YE,selectNode(2,e,m));
    float cord3 = selectCoord(YE,selectNode(3,e,m));
    float cord4 = selectCoord(YE,selectNode(4,e,m));

    Tau.at(0) = 3*pow(cord1,2) + 2*cord1*(cord2+cord3+cord4) + pow(cord2,2) + cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Tau.at(1) = pow(cord1,2) + cord1*(2*cord2+cord3+cord4) + 3*pow(cord2,2) + 2*cord2*(cord3+cord4) + pow(cord3,2) + cord3*cord4 + pow(cord4,2);
    Tau.at(2) = pow(cord1,2) + cord1*(cord2+2*cord3+cord4) + pow(cord2,2) + cord2*(2*cord3+cord4) + 3*pow(cord3,2) + 2*cord3*cord4 + pow(cord4,2);
    Tau.at(3) = pow(cord1,2) + cord1*(cord2+cord3+2*cord4) + pow(cord2,2) + cord2*(cord3+2*cord4) + pow(cord3,2) + 2*cord3*cord4 + 3*pow(cord4,2);
}


/*
Matrix createLocalM(int e,mesh &m){
    Matrix matrixA,matrixI,matrixL,matrixG,matrixD;
    float J,Determinant;
    
    // [ A+I  -L+G ]
    // [  D   0 ]
    

    //Matrix A
    Matrix GammaMatrix, Alpha, Beta;

    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(120*Determinant);
    calculateGamma(e,GammaMatrix,m);
    calculateGamma(g_matrix);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    productRealMatrix(real_a, productMatrixMatrix(GammaMatrix,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixA);


    //Matrix K
    Matrix Alpha_t,Beta_t;

    nu = m.getParameter(DYNAMIC_VISCOSITY);
    Ve = calculateLocalVolume(e,m);
    
    float real_k = (float) (nu*Ve)/(Determinant*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    productRealMatrix(real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixK);

    
    //Matrix G
    Matrix Omega;
    
    rho = m.getParameter(DENSITY);
    float real_g = (float) (J/(24*rho*Determinant));

    calculateOmega(Omega);
    productRealMatrix(real_g,productMatrixMatrix(g_matrix,productMatrixMatrix(Alpha,Omega,3,3,4),12,3,4),matrixG);

    //Matrix D
    Matrix g_matrix_t,Omega_t;
    float real_d = (float)(J/(24*Determinant));

    transpose(Omega, Omega_t);
    transpose(g_matrix,g_matrix_t);
    productRealMatrix(real_d,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,g_matrix_t,3,3,12),4,3,12),matrixD);

    //Matrix M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixA,matrixK,12,12));
    ubicarSubMatriz(M,0,11,12,15,matrixG);
    ubicarSubMatriz(M,12,15,0,11,matrixD);

    return M;
}

void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,f;
    Matrix g_matrix;

    calculateF(f, m);

    calculateGamma(g_matrix);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,16);
    productMatrixVector(g_matrix,f,b_aux);
    productRealVector(J/24,b_aux,b);
    
    return b;
}*/