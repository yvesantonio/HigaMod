/*----------------------------------------------------------------------------------------------------
		Elliptic Equation with Robin boundary conditions
----------------------------------------------------------------------------------------------------*/

/*	Generating Mesh	*/
 mesh Th=square(20,20);
 
/*	Finite Element Spaces */
 fespace Vh(Th,P1);     // P1 FE space
 
/*	Defining the variational problem	*/
 func f=cos(x)*sin(y);		//  right hand side function 
 func c=min(x,y);             
 func eps =0.1; // Robin boundary condition function

 macro grad(u) [dx(u),dy(u)] //eom <- IMPORTANT : macro ends with a comments
 
varf bilinear(u,v) =  int2d(Th)( grad(u)'*grad(v) + c*u*v) //variational formulation
	+ int1d(Th,1,2,3,4,5)(eps*u*v) ;
varf linear(unused,v) =  int2d(Th)( f*v ) ;

/*	Solving the problem	*/
matrix A = bilinear(Vh,Vh,solver = CG, factorize=0); // assemble the FE matrix
real[int] b = linear(0,Vh); // assemble the right hand side vector
Vh u;
u[] = A^-1*b; //solve the problem


/* Visualize the solution	*/ 
 plot(u,ps="ellipticRobin.eps",value=1, wait=1, fill=1, cmm="Solution u in \Omega");