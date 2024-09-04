#pragma once
class matrix
{
public:
	matrix(size_t a, size_t b); // создание матрицы из 0 размера а*b
	friend std::vector <double> operator*(matrix A, std::vector <double> b); // оператор умножения для перемножения матрицы на вектор
	friend matrix operator+(matrix A, matrix B); // оператор сложения
	friend matrix operator-(matrix A, matrix B); // оператор вычитания
	friend matrix operator*(double num, matrix A); // оператор умножения на число
	~matrix();

	void Get_matr(matrix& temp_matr, int indRow, int indCol); //функция вычеркивания строки и столбца
	double Det(); // Детерминант
	matrix Transpose(); // транспонирование
	matrix Obr_matr(); // обратная матрица

	std::vector<std::vector<double>> matr;
	size_t rows;
	size_t cols;
private:

};