#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <limits>

struct data{
int row;
int col;
double** matris;
};

data init(std::string in);
double** initialize(int rows, int cols);

data init(std::string in)
{
	data temp;

	std::istringstream iss;
	iss.str(in);
	iss >> temp.row >> temp.col;
	
	temp.matris = initialize(temp.row,temp.col);

	for(int i=0;i<temp.row;++i)
	{
		for(int j=0;j<temp.col;++j)
		iss >> temp.matris[i][j];
	}
	return temp;
}

double** initialize(int rows, int cols)
{
	double** temp;
	temp = (double**)calloc(rows , sizeof(double *));
	for(int i=0 ; i< rows ; ++i)
		temp[i] = (double*)calloc(cols , sizeof(double));
	return temp;
}

int main(int argc, char **argv)
{
	int T;
	int M;
	int N;
	
	int MaxItter = 100;
	int itter = 0;
	double OldLogProb = -std::numeric_limits<double>::max();
	
	// Read the file
	std::vector<std::string> board;
	for (std::string line; std::getline(std::cin, line);)
		board.push_back(line);

	//initialize
	data temp = init(board[0]);
	double** A = temp.matris;
	N = temp.row;
	
	temp = init(board[1]);
	double** B = temp.matris;
	M = temp.col;
	
	temp = init(board[2]);
	double* q = temp.matris[0];
	
	std::vector<int> seq;
	
	int index;
	std::istringstream iss;
	iss.str(board[3]);
	iss >> T;
	for(int i=0;i<T;++i)
	{
		iss >> index;
		seq.push_back(index);
	}
	
	//variabler som går att återanvända
	double* c = (double*) calloc(T,sizeof(double)); //c[time]
	double** alpha = initialize(T, N);  //alpha[time][index]
	double** beta = initialize(T, N);  //beta[time][index]
	double** Gamma = initialize(T,N);
	double*** diGamma = (double***)calloc(T,sizeof(double**)); //diGamma[time][index1][index2]
	for(int t=0;t<T;++t)
		diGamma[t] = initialize(N,N);
	
	/*Felkoll
	std::cout << "A" << std::endl;
	for(int i=0;i<N;++i)
	{
		for(int j=0;j<N;++j)
			std::cout << A[i][j] << " ";
		std::cout<<std::endl;
	}
	
	std::cout << "\nB" << std::endl;
	for(int i=0;i<N;++i)
	{
		for(int j=0;j<M;++j)
			std::cout << B[i][j] << " ";
		std::cout<<std::endl;
	}
	
	std::cout << "\nq" << std::endl;
	for(int i=0;i<N;++i)
		std::cout << q[i] << " ";
	*/
	
	bool GO = true;
	while((itter < MaxItter) && GO)
	{
		
		/**----------Alpha-pass------------------------*/
		
		//alpha-0
		c[0] = 0;
		
		for(int i=0;i<N;++i)
		{
			alpha[0][i] = q[i]*B[i][seq[0]];
			c[0] += alpha[0][i];
		}
		//scale
		c[0] = 1.0/c[0];
		for(int i=0;i<N;++i)
		{
			alpha[0][i] *= c[0];
		}
		
		//alpha-t
		for(int t=1;t<T;++t)
		{
			c[t] = 0;
			for(int i=0;i<N;++i)
			{
				alpha[t][i] = 0;
				for(int j=0;j<N;++j)
				{
					alpha[t][i] += alpha[t-1][j]*A[i][j];
				}
				alpha[t][i] *= B[i][seq[t]];
				c[t] += alpha[t][i];
			}
			
			//scale
			c[t] = 1.0/c[t];
			for(int i=0;i<N;++i)
			{
				alpha[t][i] *= c[t];
			}
		}
		
		/*for(int i=0;i<N;++i)
			std::cout << alpha[T-1][i] << " ";*/
			
		/**----------Beta-pass------------------------*/
		
		//beta T-1
		for(int i=0;i<N;++i)
			beta[T-1][i] = c[T-1];
			
		for(int t=T-2; t>=0;--t)
		{
			for(int i=0;i<N;++i)
			{
				beta[t][i] = 0;
				for(int j=0;j<N;++j)
				{
					beta[t][i] += A[i][j] * B[j][seq[t]] * beta[t+1][j];
				}
				//scale
				beta[t][i] *= c[t];
			}
		}
		
		
		/**----------diGamma & Gamma----------------*/
		
		for(int t=0;t<(T-1);++t)
		{
			double denom = 0;
			for(int i=0;i<N;++i)
				for(int j=0;j<N;++j)
					denom += alpha[t][i] * A[i][j] * B[j][seq[t+1]] * beta[t+1][j];
			
			for(int i=0;i<N;++i)
			{
				Gamma[t][i] = 0;
				for(int j=0;j<N;++j)
				{
					diGamma[t][i][j] = (alpha[t][i] * A[i][j] * B[j][seq[t+1]] * beta[t+1][j]) / denom;
					Gamma[t][i] += diGamma[t][i][j];
				}
			}
	
		}
		
		/**----------Estimate A,b,q------------------*/
		
		/*q*/
		for(int i=0;i<N;++i)
			q[i] = Gamma[0][i];
		
		/*A*/
		for(int i=0;i<N;++i)
		{
			for(int j=0;j<N;++j)
			{
				double numer = 0;
				double denom = 0;
				for(int t=0;t<(T-1);++t)
				{
					numer += diGamma[t][i][j];
					denom += Gamma[t][i];
				}
			A[i][j] = numer/denom;
			}
		}

		/*B*/
		for(int i=0;i<N;++i)
		{
			for(int j=0;j<M;++j)
			{
				double numer = 0;
				double denom = 0;
				for(int t=0;t<(T-1);++t)
				{
					if(seq[t] == j)
						numer += Gamma[t][i];
					denom += Gamma[t][i];
				}
			B[i][j] = numer/denom;
			}
		}
		
		double logProb = 0;
		for(int t=0;t<T;++t)
			logProb += log(c[t]);
		logProb *= -1;
		
		if(logProb < OldLogProb)
		{
			GO = false;
			break;
		}
		
		OldLogProb = logProb;
		itter++;
	}
	
	std::cout << "Iterations " << itter << std::endl;
	
	//Print matrix
	std::cout << N << " " << N;
	for(int i=0;i<N;++i)
		for(int j=0;j<N;++j)
			std::cout << " " << A[i][j];
	
	std::cout << "\n" << N << " " << M;
	for(int i=0;i<N;++i)
		for(int j=0;j<M;++j)
			std::cout << " " << B[i][j];
	std::cout << std::endl;
	
	return 0;
}
