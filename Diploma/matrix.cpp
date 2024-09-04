#include "pch.h"

matrix::matrix(size_t a, size_t b) { // создание матрицы из 0 размера а*b
	matr.assign(a, std::vector<double>(b));
	rows = a;
	cols = b;

	for (size_t i = 0; i < a; i++) {
		for (size_t j = 0; j < b; j++)
		{
			matr[i][j] = 0;
		}
	}
}


matrix::~matrix() {}

std::vector<double> operator*(matrix A, std::vector<double> b) { // оператор умножени¤ дл¤ перемножени¤ матрицы на вектор
	std::vector<double> res(b.size());
	for (int i = 0; i < b.size(); i++) {
		for (int j = 0; j < b.size(); j++) {
			res[i] += A.matr[i][j] * b[j];
		}
	}
	return res;
}

matrix operator+(matrix A, matrix B) // оператор сложени¤
{
	matrix res(A.rows, A.cols);
	for (size_t i = 0; i < A.rows; i++)
	{
		for (size_t j = 0; j < A.cols; j++)
		{
			res.matr[i][j] = A.matr[i][j] + B.matr[i][j];
		}
	}
	return res;
}

matrix operator-(matrix A, matrix B) // оператор вычитани¤
{
	matrix res(A.rows, A.cols);
	for (size_t i = 0; i < A.rows; i++)
	{
		for (size_t j = 0; j < A.cols; j++)
		{
			res.matr[i][j] = A.matr[i][j] - B.matr[i][j];
		}
	}
	return res;
}

matrix operator*(double num, matrix A) // оператор умножени¤ на число
{
	matrix res(A.rows, A.cols);
	for (size_t i = 0; i < A.rows; i++)
	{
		for (size_t j = 0; j < A.cols; j++)
		{
			res.matr[i][j] = num * A.matr[i][j];
		}
	}
	return res;
}

//функци¤ вычеркивани¤ строки и столбца
void matrix::Get_matr(matrix& temp_matr, int indRow, int indCol)
{
	int ki = 0;
	for (int i = 0; i < rows; i++)
	{
		if (i != indRow) {
			for (int j = 0, kj = 0; j < cols; j++) {
				if (j != indCol)
				{
					temp_matr.matr[ki][kj] = matr[i][j];
					kj++;
				}
			}
			ki++;
		}
	}
}

double matrix::Det()
{
	if (rows == cols)
	{
		int n = rows;
		double temp = 0;   //временна¤ переменна¤ дл¤ хранени¤ определител¤
		int k = 1;      //степень
		if (n == 1)
			temp = matr[0][0];
		else if (n == 2)
			temp = matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1];
		else
		{
			for (int i = 0; i < n; i++)
			{
				int m = n - 1;
				matrix temp_matr(m, m);
				Get_matr(temp_matr, 0, i);
				temp = temp + k * matr[0][i] * temp_matr.Det();
				k = -k;
			}
		}
		return temp;
	}
	else
	{
		std::cout << "“.к. матрица не квадратна¤, то определител¤ нет" << '\n';
		return 0;
	}
}

matrix matrix::Transpose()
{
	matrix A_trans(rows, cols);
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = i; j < cols; j++)
		{
			A_trans.matr[i][j] = this->matr[j][i];
			A_trans.matr[j][i] = this->matr[i][j];
		}
	}
	return A_trans;
}

matrix matrix::Obr_matr()
{
	matrix A_obr(rows, cols);
	if (rows == cols)
	{
		int n = rows;
		double det = Det();
		if (det)
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++) {
					int m = n - 1;
					matrix temp_matr(m, m);
					Get_matr(temp_matr, j, i);
					A_obr.matr[i][j] = pow(-1.0, i + j + 2) * temp_matr.Det() / det;
				}
			}
		}
		else
			std::cout << "“.к. определитель матрицы = 0,\nто матрица вырожденна¤ и обратной не имеет!!!" << '\n';
	}
	else
		std::cout << "“.к. матрица не квадратна¤, то обратной нет" << '\n';
	return A_obr;
}
