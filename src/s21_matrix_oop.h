#include <stdlib.h>
#include <time.h>

#include <cmath>
#include <cstring>
#include <iostream>

class S21Matrix {
 public:
  // Конструкторы и деструктор
  S21Matrix();
  S21Matrix(unsigned int rows, unsigned int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  // Геттеры и Сеттеры
  unsigned int GetRows_() const;
  unsigned int GetCols_() const;
  void SetRows_(unsigned int rows);
  void SetCols_(unsigned int cols);

  // Функции для задания
  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  // Перегрузка операторов
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(const double num) const;
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(const double num);
  bool operator==(const S21Matrix& other) const;
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  S21Matrix& operator=(const S21Matrix& other);
  double& operator()(int i, int j);

  // Вспомогательные функции
  void FillMatrix();
  void SimpleFill();
  void CopyMatrix(S21Matrix other);
  void FillMyVal(double arr[3][3]);

 private:
  unsigned int rows_;
  unsigned int cols_;
  double** matrix_;

  // Вспомогательные функции
  double RecurseDeterminant(S21Matrix& other);
  S21Matrix MinorMatrix(unsigned int rm_row, unsigned int rm_col,
                        S21Matrix& other);
  void FreeObj();
  void RowsColsSetter(unsigned int rows, unsigned int cols, S21Matrix other);
};

class Matrix_error : public std::exception {
 public:
  Matrix_error(std::string message);
  const char* what() const noexcept override;

 private:
  std::string message;
};