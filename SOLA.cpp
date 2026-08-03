//SOLA: A SIMPLIFIED PROCEDURE FOR SOLVING THE INCOMPRESSIBLE NAVIER-SOKES EQUATION
//Refer to: "SOLA - A numerical solution algorithm for transient fluid flows" by Hirt et al. 1975, Report No. - LA-5852

#include<iostream>
#include<fstream>
using namespace std;
#define EPS 1e-6	//convergence tolerance limit

const double Re=100.0;	//Reynolds number
const double L=1.0,H=1.0;	//domain size
const int I=80,J=80,COUNT=10000000;	//uniform mesh size and no of iterations (continue up to steady state)
const double dx=L/I,dy=H/J,dt=0.001,om=1.8;	//spatial and temporal resolution, relaxation factor

double P[J+1][I+1],U[J+2][I+1],V[J+1][I+2];	//flow variables at current time step
double P_prev[J+1][I+1],U_prev[J+2][I+1],V_prev[J+1][I+2];	//flow variables at previous time step
double Adv_x[J+1][I],Adv_y[J][I+1],Diff_x[J+1][I],Diff_y[J][I+1];	//advection and diffusion terms of the x and y momentum equations

void U_BC();	//boundary conditions for u
void V_BC();	//boundary conditions for v
void updt_prev();	//update velocity field at the previous time step
void Adv_Diff_x();	//calculate advection and diffusion terms for the x momentum equation
void Adv_Diff_y();	//calculate advection and diffusion terms for the y momentum equation
void mom();	//solve the momentum equations
void press_corr();	//pressure correction procedure
void solve(int file);	//solution algorithm
double continuity();	//check divergence of the flow
double steady(int t);	//check whether steady state is reached or not

void write_mesh();	//write mesh file
void write_UVP(int t);	//write flow variables at specified time intervals

int main()
{
	U_BC(); V_BC();	/impose the boundary conditions
	write_mesh();	//export the mesh
	solve(1000);	//start the solver
	return 0;
}

void U_BC()
{
	for(int i=1;i<I;i++)
	{
		U[0][i]=-U[1][i];	//no slip at bottom boundary
		U[J+1][i]=2.0-U[J][i];	//lid velocity at top boundary
	}
	for(int j=1;j<=J;j++)
	{
		U[j][0]=U[j][I]=0.0;	//no penetration at left and right boundary
	}
}

void V_BC()
{
	for(int j=1;j<J;j++)
	{
		V[j][0]=-V[j][1];	//no slip at left boundary
		V[j][I+1]=-V[j][I];	//no slip at right boundary
	}
	for(int i=1;i<=I;i++)
	{
		V[0][i]=V[J][i]=0.0;	//no penetration at top and bottom boundary
	}
}

void updt_prev()
{
	for(int j=1;j<=J;j++)	//update P
		for(int i=1;i<=I;i++)
			P_prev[j][i]=P[j][i];
	for(int j=1;j<=J;j++)	//update u
		for(int i=1;i<I;i++)
			U_prev[j][i]=U[j][i];
	for(int j=1;j<J;j++)	//update v
		for(int i=1;i<=I;i++)
			V_prev[j][i]=V[j][i];
}

void Adv_Diff_x()
{
	double N,E,W,S;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<I;i++)
		{
			E=((U[j][i]+U[j][i+1])>0.0) ? U[j][i] : U[j][i+1];	//first order upwind scheme
			W=((U[j][i]+U[j][i-1])>0.0) ? U[j][i-1] : U[j][i];
			N=((V[j][i]+V[j][i+1])>0.0) ? U[j][i] : U[j+1][i];
			S=((V[j-1][i]+V[j-1][i+1])>0.0) ? U[j-1][i] : U[j][i];
			Adv_x[j][i]=0.5*(E*(U[j][i]+U[j][i+1])-W*(U[j][i]+U[j][i-1]))/dx	//advection terms
						+0.5*(N*(V[j][i]+V[j][i+1])-S*(V[j-1][i]+V[j-1][i+1]))/dy;
			Diff_x[j][i]=((U[j][i+1]-2.0*U[j][i]+U[j][i-1])/(dx*dx)+(U[j+1][i]-2.0*U[j][i]+U[j-1][i])/(dy*dy))/Re;	//diffusion terms
		}
	}
}

void Adv_Diff_y()
{
	double N,E,W,S;
	for(int j=1;j<J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			E=((U[j][i]+U[j+1][i])>0.0) ? V[j][i] : V[j][i+1];	//first order upwind scheme
			W=((U[j][i-1]+U[j+1][i-1])>0.0) ? V[j][i-1] : V[j][i];
			N=((V[j][i]+V[j+1][i])>0.0) ? V[j][i] : V[j+1][i];
			S=((V[j][i]+V[j-1][i])>0.0) ? V[j-1][i] : V[j][i];
			Adv_y[j][i]=0.5*(E*(U[j][i]+U[j+1][i])-W*(U[j][i-1]+U[j+1][i-1]))/dx	//advection terms
						+0.5*(N*(V[j][i]+V[j+1][i])-S*(V[j][i]+V[j-1][i]))/dy;
			Diff_y[j][i]=((V[j][i+1]-2.0*V[j][i]+V[j][i-1])/(dx*dx)+(V[j+1][i]-2.0*V[j][i]+V[j-1][i])/(dy*dy))/Re;	//diffusion terms
		}
	}
}

void mom()
{
	Adv_Diff_x();	//calculate the advection and diffusion terms
	Adv_Diff_y();
	for(int j=1;j<=J;j++)	//x momentum equation
		for(int i=1;i<I;i++)
			U[j][i]+=dt*(Diff_x[j][i]-Adv_x[j][i]-(P[j][i+1]-P[j][i])/dx);
	for(int j=1;j<J;j++)	//y momentum equation
		for(int i=1;i<=I;i++)
			V[j][i]+=dt*(Diff_y[j][i]-Adv_y[j][i]-(P[j+1][i]-P[j][i])/dy);
}

void press_corr()
{
	double DP,DxDt=1.0/(2.0*dt*(1.0/(dx*dx)+1.0/(dy*dy)));
	int CNT=0;
	do
	{
		for(int j=1;j<=J;j++)
		{
			for(int i=1;i<=I;i++)
			{
				DP=DxDt*om*((U[j][i-1]-U[j][i])/dx+(V[j-1][i]-V[j][i])/dy);	//pressure correction
				if(i!=I) U[j][i]+=dt/dx*DP;	//velocity update
				if((i-1)!=0) U[j][i-1]-=dt/dx*DP;
				if(j!=J) V[j][i]+=dt/dy*DP;
				if((j-1)!=0) V[j-1][i]-=dt/dy*DP;
				P[j][i]+=DP;	//pressure update
			}
		}
		CNT++;
	}
	while((continuity())&&(CNT<=1000));	//check divergence
	//cout<<"CNT = "<<CNT<<endl;
}

void solve(int file)
{
	int CNT=0;
	do
	{
		updt_prev();
		mom();
		press_corr();
		U_BC();	//update the velocity boundary conditions
		V_BC();
		CNT++;
		if((CNT%file)==0)	//file export at regular intervals
			write_UVP(CNT);
	}
	while((steady(CNT))&&(CNT<=COUNT));	//continue until steady state is reached
}

double continuity()
{
	double RES=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			RES+=abs((U[j][i]-U[j][i-1])/dx+(V[j][i]-V[j-1][i])/dy);	//calculate divergence
	//cout<<"RES = "<<abs(RES)<<endl;
	if(abs(RES)<=EPS) return 0;
	else return 1;
}

double steady(int CNT)
{
	double RES_U=0.0,RES_V=0.0,RES_P=0.0;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			RES_U+=abs(U[j][i]-U_prev[j][i]);
			RES_V+=abs(V[j][i]-V_prev[j][i]);
			RES_P+=abs(P[j][i]-P_prev[j][i]);
		}
	}
	RES_U+=RES_V+RES_P;
	if(CNT%100==0) cout<<"RES_steady = "<<abs(RES_U)<<endl;
	if(abs(RES_U)<=10*EPS) return 0;
	else return 1;
}

void write_mesh()
{
	double Xm[I+1],Ym[J+1];
	ofstream p_out("mesh.dat");
	p_out<<"TITLE = \"MESH\""<<endl;
	p_out<<"FILETYPE = GRID"<<endl;
	p_out<<"VARIABLES = \"X\",\"Y\""<<endl;
	p_out<<"ZONE I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK"<<endl;
	Xm[0]=0.0;
	Ym[0]=0.0;
	for(int i=1;i<=I;i++)	//generation of mesh
		Xm[i]=Xm[i-1]+dx;
	for(int j=1;j<=J;j++)
		Ym[j]=Ym[j-1]+dy;
	for(int j=0;j<=J;j++)	//print X co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Xm[i];
		p_out<<endl;
	}
	p_out<<endl<<endl;
	for(int j=0;j<=J;j++)	//print Y co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Ym[j];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MESH WRITE SUCCESSFULL"<<endl;
}

void write_UVP(int t)
{
	string fname="uvp_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"FLOW AND PRESSURE FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"u\",\"v\",\"P\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<0.5*(U[j][i]+U[j][i-1]);	//linear interpolation of the staggered variables
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<0.5*(V[j][i]+V[j-1][i]);	//linear interpolation of the staggered variables
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<P[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
