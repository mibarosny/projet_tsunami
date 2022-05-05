# include "tsunami.h"

int nElem, nEdge, nNode, size;
double *nodes;
int *elems;
int *edges;
double *E;
// *U = E + size
// *V = E + 2*size
double *FE;
// *FU = FE + size
// *FV = FE + 2*size

double timestep;
double  xLoc[3],yLoc[3],phi[3],dphidx[3],dphidy[3];
double  xsi,eta,weight,jac;
double  y,e,u,v,x,h;
int     i,j,k,el,elem,mapElem[3],mapEdge[2][2];

void tsunamiComputeInitialize() {
    for (el = 0; el < nElem; el++)
        for (i = 0; i < 3; i++) {
            E[el*3+i] = tsunamiInitialConditionOkada(nodes[3*elems[el*3+i]],
                                                              nodes[3*elems[el*3+i]+1]);
            E[size +   el*3+i] = 0.;
            E[size*2 + el*3+i] = 0.;
        }
}

void femShallowTriangleMap(int index, int map[3]) {
    for (i = 0; i < 3; i++)  map[i] = index*3 + i;
}

void femShallowEdgeMap(int index, int map[2][2]) {
    for (j=0; j < 2; ++j) {
        int node = edges[index*4 + j];
        for (k=0; k < 2; k++) {
            int elem = edges[index*4 + k + 2];
            map[k][j] = (nElem)*3;
            if (elem >= 0) {
                for (i=0; i < 3; i++) {
                    if (elems[elem*3 + i] == node) {
                        map[k][j] = elem*3 + i;  }}}}}
}

void femDiscretePhi(int n, double xsi, double eta, double * phi) {
    if (n == 2) {
        phi[0] = 1 - xsi - eta; phi[1] = xsi; phi[2] = eta;
    } else {
        phi[0] = (1.0 - xsi)/2.0; phi[1] = (1.0 + xsi)/2.0;  
    }
}

double interpolate(double *phi, double *U, int *map, int e, int p, int n) {
    double u = 0.0;
    for (i=0; i < n; i++)
        u += phi[i]*U[e*map[i]+p];
    return u;
}

void femShallowAddIntegralsElements() {
    double *BE = FE;
    double *BU = FE + size;
    double *BV = FE + 2*size;
    double *U = E + size;
    double *V = E + 2*size;    

    for (elem=0; elem < nElem; elem++) {
        femShallowTriangleMap(elem, mapElem);
        int *mapCoord = &(elems[elem*3]);
        for (j=0; j < 3; ++j) {
        	  xLoc[j] = nodes[3*mapCoord[j]];
        	  yLoc[j] = nodes[3*mapCoord[j]+1]; }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        dphidx[0] = (yLoc[1] - yLoc[2])/jac;
        dphidx[1] = (yLoc[2] - yLoc[0])/jac;
        dphidx[2] = (yLoc[0] - yLoc[1])/jac;
        dphidy[0] = (xLoc[2] - xLoc[1])/jac;
        dphidy[1] = (xLoc[0] - xLoc[2])/jac;
        dphidy[2] = (xLoc[1] - xLoc[0])/jac;        
        for (k=0; k < 3; k++) {
            xsi = gaussTriangleXsi[k];
            eta = gaussTriangleEta[k];
            weight = gaussTriangleWeight[k];     
            femDiscretePhi(2,xsi,eta,phi);               
            x = interpolate(phi, nodes, mapCoord, 3, 0, 3); 
            y = interpolate(phi, nodes, mapCoord, 3, 1, 3); 
            h = interpolate(phi, nodes, mapCoord, 3, 2, 3); 
            e = interpolate(phi,E,mapElem,1,0,3);
            u = interpolate(phi,U,mapElem,1,0,3);
            v = interpolate(phi,V,mapElem,1,0,3);
            
            
            double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
  	        double lat = asin(z3d/R)*180/PI;
            double coriolis = 2.*Omega*sin(lat);     
            for (i=0; i < 3; i++) {
                BE[mapElem[i]] += ( (dphidx[i]*h*u + dphidy[i]*h*v)*((4.*R*R+x*x+y*y)/(4.*R*R))
                                + phi[i]*h*((x*u+y*v)/(R*R)) )*jac*weight;
                BU[mapElem[i]] += ( phi[i]*(coriolis*v - Gamma*u) + dphidx[i]*g*e*((4.*R*R+x*x+y*y)/(4.*R*R))
                                + phi[i]*(g*x*e)/(2*R*R))*jac*weight;
                BV[mapElem[i]] += ( phi[i]*(-coriolis*u - Gamma*v) + dphidy[i]*g*e*((4.*R*R+x*x+y*y)/(4.*R*R))
                                + phi[i]*(g*y*e)/(2*R*R))*jac*weight;
            }
        }
    
    }        

}


void femShallowAddIntegralsEdges() {
    double *BE = FE;
    double *BU = FE + size;
    double *BV = FE + 2*size;
    double *U = E + size;
    double *V = E + 2*size;   
    
    double  eL,eR,uL,uR,vL,vR,unL,unR;
    double  qe,qu,qv,x,y,h,factor;

    for (int edge=0; edge < nEdge; edge++) {
        femShallowEdgeMap(edge,mapEdge);
        for (j=0; j < 2; ++j) {
        	  int node = edges[edge*4 + j];
        	  xLoc[j] = nodes[node*3];
        	  yLoc[j] = nodes[node*3 + 1]; }
        int *mapCoord = &(edges[edge*4]);
        int boundary = (mapEdge[1][0] == size-1);
        
        double dxdxsi = (xLoc[1] - xLoc[0]);
        double dydxsi = (yLoc[1] - yLoc[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        jac = norm / 2.0;
        for (k=0; k < 2; k++) {
            xsi = gaussEdgeXsi[k];
            weight = gaussEdgeWeight[k];     
            femDiscretePhi(1,xsi,xsi,phi);           
            
            eL = interpolate(phi,E,mapEdge[0],1,0,2);
            eR = boundary ? eL : interpolate(phi,E,mapEdge[1],1,0,2);
            uL = interpolate(phi,U,mapEdge[0],1,0,2);
            uR = interpolate(phi,U,mapEdge[1],1,0,2);
            vL = interpolate(phi,V,mapEdge[0],1,0,2);
            vR = interpolate(phi,V,mapEdge[1],1,0,2);
            y = interpolate(phi,nodes,mapCoord,3,1,2);
            x = interpolate(phi,nodes,mapCoord,3,0,2);
            h = interpolate(phi,nodes,mapCoord,3,2,2);
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            factor = (4.*R*R+x*x+y*y)/(4.*R*R);

            qe =  0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) )*factor;
            qu =  0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) )*factor;
            qv =  0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) )*factor;        
            for (i=0; i < 2; i++) {
                BE[mapEdge[0][i]] -= qe*phi[i]*jac*weight; 
                BU[mapEdge[0][i]] -= qu*phi[i]*jac*weight; 
                BV[mapEdge[0][i]] -= qv*phi[i]*jac*weight; 
                BE[mapEdge[1][i]] += qe*phi[i]*jac*weight;
                BU[mapEdge[1][i]] += qu*phi[i]*jac*weight;
                BV[mapEdge[1][i]] += qv*phi[i]*jac*weight; }}}

}

void femShallowMultiplyInverseMatrix() {
    double *BE = FE;
    double *BU = FE + size;
    double *BV = FE + 2*size;
         
    double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}};
    double BEloc[3],BUloc[3],BVloc[3];
    
    for (elem=0; elem < nElem; elem++) {
        femShallowTriangleMap(elem, mapElem);
        int *mapCoord = &(elems[elem*3]);
        for (j=0; j < 3; ++j) {
        	  xLoc[j] = nodes[3*mapCoord[j]];
        	  yLoc[j] = nodes[3*mapCoord[j]+1]; }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        for (i=0; i < 3; i++) {
            BEloc[i] = BE[mapElem[i]];
            BUloc[i] = BU[mapElem[i]];
            BVloc[i] = BV[mapElem[i]];
            BE[mapElem[i]] = 0.0; 
            BU[mapElem[i]] = 0.0; 
            BV[mapElem[i]] = 0.0; }
        for (i=0; i < 3; i++) { 
            for (j=0; j < 3; j++) {
                BE[mapElem[i]] += invA[i][j] * BEloc[j] / jac; 
                BU[mapElem[i]] += invA[i][j] * BUloc[j] / jac; 
                BV[mapElem[i]] += invA[i][j] * BVloc[j] / jac; }}}
}


void tsunamiUpdateRungeKutta() {
    int i, j;
    double * E_ = E;
    double * Eold        = malloc(sizeof(double)*size*3);
    double * Epredictor  = malloc(sizeof(double)*size*3);
    double  *FE_ = FE;
    for (i=0; i < size*3; i++) {
        FE_[i] = 0.0;
        Eold[i] = E_[i];  
    }

    const double beta[2]  = {0.0,   1.0  }; 
    const double gamma[2] = {1./2., 1./2.};
    
    E = Epredictor;
    for (j = 0; j < 2; j++) {
        for (i=0; i < size*3; i++) 
            Epredictor[i] = Eold[i] + timestep * beta[j] * FE_[i];

        femShallowAddIntegralsElements();
        femShallowAddIntegralsEdges();
        femShallowMultiplyInverseMatrix();

        for (i=0; i < size*3; i++) 
            E_[i] += timestep * gamma[j] * FE_[i]; 
    }
    E = E_;
    
    free(Eold);
    free(Epredictor);
}

// void euler() {
//     for (i = 0; i < size * 3; i++) {
//         FE[i] = 0.0;
//     }

//     femShallowAddIntegralsElements();
//     femShallowAddIntegralsEdges();
//     femShallowMultiplyInverseMatrix();
//     for (j=0; j < size * 3; j++) {
//         E[j] += timestep * FE[j];
//     }
// }


void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 
    timestep = dt;
   /** lecture du fichier */
    int trash;

    FILE* file = fopen(meshFileName,"r");
    fscanf(file, "Number of nodes %d \n",&nNode);   

    nodes = malloc(sizeof(double)*nNode*3);
    for (i = 0; i < nNode; i++) 
        fscanf(file,"%d : %le %le %le\n",&trash,&nodes[i*3],&nodes[i*3+1],&nodes[i*3+2]); 
        
    fscanf(file, "Number of triangles %d \n",&nElem); 
    elems = malloc(sizeof(int)*3*(nElem));
    size = 3*nElem + 1;
    for (i = 0; i < nElem; i++) 
        fscanf(file,"%d : %d %d %d \n",&trash,&elems[i*3],&elems[i*3+1],&elems[i*3+2]);   

    fscanf(file, "Number of edges %d \n",&nEdge); 
    edges = malloc(sizeof(int)*4*(nEdge));
    for (i = 0; i < nEdge; i++) 
        fscanf(file,"%d : %d %d : %d %d \n",&trash,&edges[i*4],&edges[i*4+1],&edges[i*4+2],&edges[i*4+3]); 

    fclose(file); 

    E = malloc(sizeof(double)*(size)*6);  

    FE = E + 3*size;

    tsunamiComputeInitialize();
    tsunamiWriteFile(baseResultName,0,E+size,E+2*size,E,nElem,3);

    for (int it = 0; it < nmax; it++) {
        
        tsunamiUpdateRungeKutta();
        
        if ((it+1)%sub == 0)
            tsunamiWriteFile(baseResultName,(it+1),E+size,E+2*size,E,nElem,3); 
    }

    free(E);
    free(nodes);
    free(elems);
    free(edges);
}
