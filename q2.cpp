#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
//#include <iterator>
//#include <iomanip>
//#include <functional>

using std::vector;
using namespace std;

double hermite(double t);

double speed(double t);

double Li(double t, int ii);

double Li_prime(double t, int ii);

double Li_doubleprime(double t, int ii);

void sum_coefficeints();

void product2_coefficients();

void product3_coefficients();

int N(0); // number of grid points
vector<double> x(N+1, 0.0); // the value of x
vector<double> y(N+1, 0.0); // for the value f(x)
vector<double> y_prime(N+1, 0.0); // for the value of f'(x)
vector<double> sum(N+1, 0.0); // sum of coefficients
vector<double> product4(N+1, 0.0); // product of all coefficients
vector<double> product2(N+1, 0.0); // product of pair-wise coefficients
vector<double> product3(N+1, 0.0); // product of triple coefficients
vector<double> denominator(N+1, 0.0); // denominator of the lagrange polynomial

int main()
{
	std::cout<<"Enter the number of grid points"<<std::endl;
	std::cin>>N;
	if(N<=1)
	{
		std::cout<<"Invalid input"<<std::endl;
		return 0;
	}
	
	for(int ii=0;ii<N; ii++)
	{
		std::cout<<"Enter the "<<ii+1<<"th grid point"<<std::endl;
		std::cin>>x[ii];
		std::cout<<"Enter the value of f(x) at "<<ii+1<<"th grid point"<<std::endl;
		std::cin>>y[ii];
		std::cout<<"Enter the value of f'(x) at "<<ii+1<<"th grid point"<<std::endl;
		std::cin>>y_prime[ii];
	}
	sum_coefficeints();
	product2_coefficients();
	product3_coefficients();
	
	
	/* for(int jj=0; jj<N;jj++)
	{
		std::cout<<x[jj]<<" "<<y[jj]<<" "<<y_prime[jj]<<" "<<sum[jj]<<" "<<product2[jj]<<" "<<product3[jj]<<" "<<product4[jj]<<" "<<denominator[jj]<<std::endl;
	} */
	
	double t; // point of evaluation
	std::cout<<"Enter the evaluation point"<<std::endl;
	std::cin>>t;
	//std::cout<<"t ="<<t<<std::endl;
	std::cout<<"The value of H_"<<N-1<<"("<<t<<")="<<hermite(t)<<std::endl;
	double time(0.0); // for recording the time
	double grid(13.0/1000); // it is the delta_x
	double spe(0.0); // records the speed of the car.
	for(int alp=0; alp<1000;alp++)
	{
		if(speed(grid*alp)>=80.667)
		{
			std::cout<<"The car crosses 55mi/h for the first time at t="<<grid*alp<<std::endl;
			time = alp;
			break;
		}
	}
	double temp(0.0);// temperary variable for recording speed in order to save computation time.
	for(int al=int(time);al<N;al++)
	{
		temp = speed(grid*al);
		if(temp>spe)
		{
			spe = temp;
			time = grid*al;
		}
	}
	std::cout<<"The maximum speed of the car is "<<spe<<" at t ="<<time<<std::endl;
	return 0;
}

double speed(double t)
{
	double sp(0.0);
	for(int ii=0;ii<N;ii++)
	{
		sp += y[ii]*(-2.0*Li_prime(x[ii],ii)*Li(t, ii)*Li(t, ii)+2*(1-2*Li_prime(x[ii], ii)*(t-x[ii]))*Li(t, ii)*Li_prime(t, ii))+y_prime[ii]*Li(t, ii)*(Li(t, ii)+2*(t-x[ii])*Li_prime(t, ii));
	}
	return sp;
}

double hermite(double t)
{
	//std::cout<<"Welcome 4"<<std::endl;
	double h_t(0.0);
	for(int ii=0; ii<N; ii++)
	{
		h_t += Li(t, ii)*Li(t, ii)*(y[ii]*(1.0-2.0*Li_prime(x[ii], ii)*(t-x[ii]))+y_prime[ii]*(t-x[ii]));
	}
	//std::cout<<"Thank you 4"<<std::endl;
	return h_t;
}

double Li(double t, int ii)
{
	double l_i(0.0);
	l_i = pow(t,4)-sum[ii]*pow(t,3)+product2[ii]*pow(t,2)-product3[ii]*t+product4[ii];
	l_i = l_i/denominator[ii];
	return l_i;
}

double Li_prime(double t, int ii)
{
	double l_i_prime(0.0);
	l_i_prime = 4*pow(t,3)-3*sum[ii]*pow(t,2)+2*product2[ii]*t-product3[ii];
	l_i_prime= l_i_prime/denominator[ii];
	return l_i_prime;
}

double Li_doubleprime(double t, int ii)
{
	double l_i_doubleprime(0.0);
	l_i_doubleprime = 12*pow(t,2)-6*sum[ii]*t+2*product2[ii];
	l_i_doubleprime = l_i_doubleprime/denominator[ii];
	
	return l_i_doubleprime;
}

void sum_coefficeints()
{
	double sum1(0.0);
	double prod(1.0);
	double deno(1.0);
	//std::cout<<"Welcome 1"<<std::endl;
	for(int ii=0;ii<N;ii++)
	{
		sum1 = 0.0;
		prod = 1.0;
		deno = 1.0;
		for(int jj=0;jj<N;jj++)
		{
			if(jj!=ii)
			{
				sum1 += x[jj];
				prod *= x[jj];
				deno *= ( x[ii]-x[jj]);
			}
		}
		sum[ii]=sum1;
		product4[ii]=prod;
		denominator[ii]=deno;
	}
	//std::cout<<"Thank you 1"<<std::endl;
}

void product2_coefficients()
{
	double sum1(0.0);
	//std::cout<<"Welcome 2"<<std::endl;
	for(int ii=0; ii<N;ii++)
	{
		sum1 = 0.0;
		for(int jj=0;jj<(N-1);jj++)
		{
			if(jj!=ii)
			{
				for(int kk=jj+1;kk<N;kk++)
				{
					if(kk!=ii)
					{
						sum1 += x[kk]*x[jj];
					}
				}
			}
		}
		product2[ii]=sum1;
	}
	//std::cout<<"Thank you 2"<<std::endl;
}

void product3_coefficients()
{
	//std::cout<<"Welcome 3"<<std::endl;
	double sum1(0.0);
	for(int ii=0; ii<N;ii++)
	{
		sum1 = 0.0;
		for(int jj=0;jj<(N-2);jj++)
		{
			if(jj!=ii)
			{
				for(int kk=jj+1;kk<(N-1);kk++)
				{
					if(kk!=ii)
					{
						for(int ll=kk+1;ll<N;ll++)
						{
							if(ll!=ii)
								sum1 += x[jj]*x[kk]*x[ll];
						}
					}
				}
			}
		}
		product3[ii]=sum1;
	}
	//std::cout<<"Thank you 3"<<std::endl;
}
