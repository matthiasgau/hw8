#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

void step(double*, double*, const double, const int);
void Hamilton(double*, double*, double&);

int main(){
  const double pi=3.14159265359;
  const double e=0.6;
  const int dim=2; //dimension of the problem
  double p[dim], q[dim];
  const double Tend=20*pi;
  const double Tstart=0;
  double H;
  
  const double dt=0.0005;
  const int N= int((Tend-Tstart)/dt);
  double t=Tstart;
  
  q[0]=1-e;
  q[1]=0;
  p[0]=0;
  p[1]=sqrt((1+e)/(1-e));
  Hamilton(p,q,H);
  
  ofstream out("data.txt");
  out << t << "\t" << q[0] << "\t" << q[1] << "\t" << H << endl;
  
  for(int i=0; i<N; i++){
    step(p,q,dt,dim);
    Hamilton(p,q,H);
    t+=dt;
    out << t << "\t" << q[0] << "\t" << q[1] << "\t" << H << endl;
  }
  
  
  out.close();
  
  return 0;
}

void step(double* p, double* q, const double dt, const int dim){
  double r=q[0]*q[0]+q[1]*q[1];
  for(int i=0; i<dim;i++){
    p[i]-= dt*q[i]/pow(r,3.0/2.0);
    q[i]+=dt*p[i];
  }
}
void Hamilton(double* p, double* q, double& H){
  H=0.5*(p[0]*p[0]+p[1]*p[1])-1.0/sqrt(q[0]*q[0]+q[1]*q[1]);
}