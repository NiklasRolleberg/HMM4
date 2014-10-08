#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>

struct matrix
{
	int row;
	int col;
	double** m;

	//constructor
	matrix(double** mat, int rows, int cols)
	{
		row = rows;
		col = cols;
		m = mat;
	}

	//constructor 2
	matrix(std::string in)
	{
		std::istringstream iss;
		iss.str(in);
		iss >> row >> col;

        m = initialize(row,col);

		for(int i=0;i<row;++i)
		{
			for(int j=0;j<col;++j)
			iss >> m[i][j];
		}
	}

    //creted the 2d matrix
	static double** initialize(int rows, int cols)
	{
	    double** temp;
	    temp = (double**)calloc(rows , sizeof(double *));
		for(int i=0 ; i< rows ; ++i)
			temp[i] = (double*)calloc(cols , sizeof(double));
        return temp;
	}

	double get(int r, int c)
	{
		return m[r][c];
	}

	void print()
	{
		for (int i=0;i<row;++i)
		{
			for(int j=0;j<col;++j)
				std::cout << m[i][j] <<" ";
			std::cout << std::endl;
		}
	}

	void writeMatrix()
	{
		std::cout << row << " " << col;
		for(int i=0;i<row;++i) {
			for(int j=0; j<col;++j)
				std::cout << " " << m[i][j];
		}
	}

	matrix transpose()
	{
		int i,j;
		double** temp = initialize(col,row);
		for(i=0;i<row;++i)
        {
            for(j=0;j<col;++j)
                temp[j][i] = m[i][j];
        }
        return matrix(temp,col,row);
	}

	matrix getCol(int index)
	{
		int newRows = row;
		int newCols = 1;

		double** temp = initialize(newRows,newCols);

		for(int i=0;i<newRows;++i)
			temp[i][0] = m[i][index];
		return matrix(temp,newRows,newCols);
	}

	matrix getRow(int index)
	{
		int newRows = 1;
		int newCols = col;

		double** temp = initialize(newRows,newCols);

		for(int i=0;i<newCols;++i)
			temp[0][i] = m[index][i];
		return matrix(temp,newRows,newCols);
	}

	matrix operator*(matrix H)
	{
		int newRows = row;
		int newCols = H.col;

		if(col != H.row)
			std::cerr << "Fel. matrix 1: " << row << " " << col << ", matrix 2: " << H.row << " " << H.col << std::endl;

		double** temp;
		temp = initialize(newRows,newCols);

		//RÄKNA!
		double t;
		int i,j,k;
		for (i=0;i<newRows;++i)
		{
			for(j=0;j<newCols;++j)
			{
				t = 0;
				//row in first matrix X col in second matrix;
				for(k=0;k<col;++k)
					t+= (m[i][k] * H.get(k,j));
				temp[i][j] = t;
			}
		}
		return matrix(temp,newRows,newCols);
	}

	matrix operator&(matrix H) // element multiplication
	{
		int newRows = row;
		int newCols = col;

		if(row!=H.row || col != H.col)
			std::cout << "FEL " << row << " " << H.row << " , " << col <<" " <<H.col <<std::endl;

		double** temp = initialize(newRows,newCols);

		//RÄKNA!
		int i,j;
		for (i=0;i<newRows;++i)
		{
			for(j=0;j<newCols;++j)
				temp[i][j] = m[i][j]*H.get(i,j);
		}
		return matrix(temp,newRows,newCols);
	}

};

int main(int argc, char **argv)
{
	// Read the file
	std::vector<std::string> board;
	for (std::string line; std::getline(std::cin, line);)
		board.push_back(line);

	matrix A = matrix(board[0]);
	matrix B = matrix(board[1]);
	matrix q = matrix(board[2]);
	std::vector<int> seq;

	//Felkoll

    /*
    std::cout << "\nA" <<std::endl;
    A.print();
    std::cout << "\nB" <<std::endl;
    B.print();
    std::cout << "\nq" <<std::endl;
    q.print();

	//matrix X = A*B*A*B*A*B*A*B*(q*A).transpose(); //ger rätt svar
    //X.print();
    */

	//find out the statesequence
	int Nstates;
	int index;
	std::istringstream iss;
	iss.str(board[3]);
	iss >> Nstates;
	for(int i=0;i<Nstates;++i)
	{
		iss >> index;
		seq.push_back(index);
	}
	
	
	//forward
	
	//matrix with all alpha values
	matrix alpha = matrix(A.initialize(B.row,Nstates),B.row,Nstates);
	
	matrix L = q&B.getCol(seq[0]).transpose();
	
	for(int j=0;j<B.row;++j)
			alpha.m[j][0] = L.get(0,j);
			
	for(int i=1;i<Nstates;++i)
	{
		L = (L*A)&(B.getCol(seq[i]).transpose());
		
		for(int j=0;j<B.row;++j)
			alpha.m[j][i] = L.get(0,j);
	}

	//sum the elements in L
	double sumA = 0;
	for(int i=0;i<L.col;++i)
            sumA += L.get(0,i);
	std::cout << "Forward: " << sumA << std::endl;
	
	//Backward
	
	//matrx with all beta values
	matrix beta = matrix(A.initialize(A.row,Nstates),A.row,Nstates);
	
	double** temp = A.initialize(1,A.col);
	for(int i=0;i<A.col;++i)
		temp[0][i] = 1.0;
		
	matrix b = matrix(temp,1,A.col);
	for(int j=0;j<B.row;++j)
			beta.m[j][Nstates-1] = b.get(0,j);
			
	for(int k=Nstates-2;k>=0;--k)
	{
		//b = (A*B.getCol(seq[k])).transpose()&b; //antagligen fel
		
		for(int i=0;i<A.col;++i)
		{
			double temp = 0;
			for(int j=0;j<A.row;++j)
				temp += A.get(i,j) * B.get(j,seq[k]) * b.get(0,j);
			
			b.m[0][i] = temp;
		}
		
		for(int j=0;j<B.row;++j)
			beta.m[j][k] = b.get(0,j);		
	}
	std::cout << "Backward: ";
	b.print();
	
	std::cout << "\nresten\n" <<std::endl;
	alpha.print();
	std::cout <<"\n"<<std::endl;
	beta.print();
	
	//Di-gamma
	matrix diGamma(A.initialize(A.row,A.col),A.row,A.col);
	
	int t = 0;
	
	for(int i=0;i<A.row;++i)
	{
		for(int j=0;j<A.col;++j)
			diGamma.m[i][j] = (alpha.get(i,t) * A.get(i,j) * B.get(j,t+1))/sumA;
	}
	
	std::cout << "\ndi-gamma\n" << std::endl;
	diGamma.print();
	return 0;
}

