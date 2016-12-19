#define ARMA_NO_DEBUG
#include <armadillo>
#include <iostream>
#include <vector>
using namespace arma;
using namespace std;

double pdim_arma(const mat & a,const vec & b,const vec & c,const vec & eq)
{
    mat A = a;
    vec B = b;
    vec C = c;

    uword n = a.n_cols;
    vec EQ = eq;

    double max_C = arma::max(arma::max(abs(A)));
    max_C = std::max(max_C,arma::max(abs(C)));
    max_C = std::max(max_C,arma::max(abs(B)));

    for(uword i=0;i<eq.n_elem;i++)
    {
        if(B(i) < 0)
        {
            B(i) *= -1.0;
            A.row(i) *= -1.0;
            EQ(i) *= -1.0;
        }
        if(EQ(i)<0)
        {
            A.insert_cols(A.n_cols,1,true);
            C.insert_rows(C.n_elem,1,true);
            A(i,A.n_cols-1) = 1.0;
        }
        else if(EQ(i)>0)
        {
            A.insert_cols(A.n_cols,1,true);
            C.insert_rows(C.n_elem,1,true);
            A(i,A.n_cols-1) = -1.0;
        }
        else
        {
            A.insert_cols(A.n_cols,1,true);
            C.insert_rows(C.n_elem,1,true);
            A(i,A.n_cols-1) = 1.0;
            C(C.n_elem-1) = 100*max_C;
        }
    }
    
    int N = A.n_cols;
    int M = A.n_rows;

    double sig = 0.99;
    
    vec x = 100.0*max_C * ones<vec>(N);
    vec y = zeros<vec>(M);
    vec s = 100.0*max_C * ones<vec>(N);
    
    vec r_x = A*x -B;
    vec r_y = A.t()*y + s - C;

    double m_xs = mean(x%s);

    vec r_s = x%s - sig*m_xs;

    double b_c = 1+std::max(norm(b),norm(c));

    uword t = 0;
    while(t<100)
    {
	t++;
        r_s = x%s;
        m_xs = mean(x%s);

        double rel_res = sqrt(dot(r_s,r_s)+dot(r_x,r_x)+dot(r_y,r_y))/b_c;
        if(rel_res < 1e-5 && m_xs < 1e-5)
            break;

        r_s = x%s - std::min(0.1,100.*m_xs)*m_xs;
        
        mat d = diagmat(arma::min(x/s,5e15*ones<vec>(N)));
        mat M = A * d * A.t();
        
        vec aux_1 = x % r_y - r_s;
        vec aux_2 = -r_x-A*(aux_1/s);
        vec dy = solve(M,aux_2);
        vec dx = (x%(A.t()*dy)+aux_1)/s;
        vec ds = -(s%dx+r_s)/x;
        
        double tau = std::max(0.9995,1-m_xs);
        double a_p = -1/std::min(arma::min(dx/x),-1.);
        double a_d = -1/std::min(arma::min(ds/s),-1.);
        
        x += tau*a_p*dx;
        y += tau*a_d*dy;
        s += tau*a_d*ds;
        
        r_x = A*x - B;
        r_y = A.t()*y + s - C;
        m_xs = mean(x%s);
        r_s = x%s - sig*m_xs;
    }
    return dot(c,x.rows(0,n-1));
}

double pdim(vector<vector<double> > & a,vector<double> & b,vector<double> & c,vector<double> & eq)
{
	mat A  = zeros<mat>(a.size(),a[0].size());
	vec B  = zeros<vec>(b.size());
	vec C  = zeros<vec>(c.size());
	vec EQ = zeros<vec>(eq.size());
	for(uword i=0;i<A.n_rows;i++)
		for(uword j=0;j<A.n_cols;j++)
			A(i,j) = a[i][j];
	for(uword i=0;i<B.n_elem;i++)
		B(i) = b[i];
	for(uword i=0;i<C.n_elem;i++)
		C(i) = c[i];
	for(uword i=0;i<EQ.n_elem;i++)
		EQ(i) = eq[i];
	return pdim_arma(A,B,C,EQ);
}
