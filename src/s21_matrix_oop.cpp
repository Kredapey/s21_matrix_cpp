#include "s21_matrix_oop.h"

// Конструкторы и деструктор
// ______________________________________________________________________

S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

S21Matrix::S21Matrix(unsigned int rows, unsigned int cols) {
  srand(time(NULL));
  rows_ = rows;
  cols_ = cols;
  matrix_ = new double *[rows_];
  for (unsigned i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : S21Matrix(other.rows_, other.cols_) {
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept
    : rows_(0), cols_(0), matrix_(nullptr) {
  std::swap(rows_, other.rows_);
  std::swap(cols_, other.cols_);
  std::swap(matrix_, other.matrix_);
}

S21Matrix::~S21Matrix() { FreeObj(); }

// ___________________________________________________________________________________________

// Геттеры и Сеттеры
// ___________________________________________________________________________________________

void S21Matrix::SetRows_(unsigned int rows) {
  S21Matrix temp(*this);
  FreeObj();
  rows_ = rows;
  cols_ = temp.cols_;
  matrix_ = new double *[rows_];
  for (unsigned int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
  if (rows_ > temp.rows_) {
    RowsColsSetter(temp.rows_, cols_, temp);
    for (unsigned int i = temp.rows_; i < rows_; i++) {
      for (unsigned int j = 0; j < cols_; j++) {
        matrix_[i][j] = 0;
      }
    }
  } else {
    RowsColsSetter(rows_, cols_, temp);
  }
}

void S21Matrix::SetCols_(unsigned int cols) {
  S21Matrix temp(*this);
  FreeObj();
  cols_ = cols;
  rows_ = temp.rows_;
  matrix_ = new double *[rows_];
  for (unsigned int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
  if (cols_ > temp.cols_) {
    for (unsigned int i = 0; i < rows_; i++) {
      for (unsigned int j = 0; j < temp.cols_; j++) {
        matrix_[i][j] = temp.matrix_[i][j];
      }
      for (unsigned int j = temp.cols_; j < cols_; j++) {
        matrix_[i][j] = 0;
      }
    }
  } else {
    RowsColsSetter(rows_, cols_, temp);
  }
}

unsigned int S21Matrix::GetRows_() const { return rows_; }

unsigned int S21Matrix::GetCols_() const { return cols_; }

// _____________________________________________________________________________________________

// Функции для задания
// _______________________________________________________________________________________________

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  bool result = true;
  constexpr double exp = 1e-07;
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (unsigned int i = 0; i < rows_; i++) {
      for (unsigned int j = 0; j < cols_; j++) {
        if (fabsl(matrix_[i][j] - other.matrix_[i][j]) > exp) {
          result = false;
        }
      }
    }
  } else {
    result = false;
  }
  return result;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (unsigned int i = 0; i < rows_; i++) {
      for (unsigned int j = 0; j < cols_; j++) {
        matrix_[i][j] += other.matrix_[i][j];
      }
    }
  } else {
    throw Matrix_error("The sizes of the matrices do not match");
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (unsigned int i = 0; i < rows_; i++) {
      for (unsigned int j = 0; j < cols_; j++) {
        matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  } else {
    throw Matrix_error("The sizes of the matrices do not match");
  }
}

void S21Matrix::MulNumber(const double num) {
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (rows_ == other.cols_) {
    S21Matrix temp(*this);
    SetCols_(other.cols_);
    for (unsigned int i = 0; i < temp.rows_; i++) {
      for (unsigned int j = 0; j < other.cols_; j++) {
        matrix_[i][j] = 0;
        for (unsigned int k = 0; k < other.rows_; k++) {
          matrix_[i][j] += temp.matrix_[i][k] * other.matrix_[k][j];
        }
      }
    }
  } else {
    throw Matrix_error(
        "The number of columns of the first matrix does not match the number "
        "of rows of the second");
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  double res = 0.0;
  if (rows_ == cols_) {
    res = RecurseDeterminant(*this);
    return res;
  } else {
    throw Matrix_error("The matrix is not square");
  }
}

S21Matrix S21Matrix::CalcComplements() {
  double temp = 0.0;
  S21Matrix result(rows_, cols_);
  if (rows_ == cols_) {
    if (rows_ == 1) {
      result.matrix_[0][0] = Determinant() * pow(-1, 0);
    } else {
      for (unsigned int i = 0; i < rows_; i++) {
        for (unsigned int j = 0; j < cols_; j++) {
          S21Matrix minor = MinorMatrix(i, j, *this);
          temp = minor.Determinant();
          result.matrix_[i][j] = temp * pow(-1, i + j);
        }
      }
    }
    return result;
  } else {
    throw Matrix_error("The matrix is not square");
  }
}

S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix result(1, 1);
  double det = Determinant();
  if (det != 0) {
    if (rows_ == 1) {
      if (matrix_[0][0] != 0) {
        result.matrix_[0][0] = 1.0 / (double)matrix_[0][0];
      } else {
        throw Matrix_error(
            "The only element of the matrix is 0. It cannot be divided by 0");
      }
    } else {
      S21Matrix temp_calc = CalcComplements();
      S21Matrix temp_trans = temp_calc.Transpose();
      temp_trans.MulNumber(1.0 / det);
      result.CopyMatrix(temp_trans);
    }
    return result;
  } else {
    throw Matrix_error("The determinant is 0");
  }
}

// _________________________________________________________________________________

// Вспомогательные функции
// __________________________________________________________________________________

void S21Matrix::FillMatrix() {
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] = ((double)rand() / RAND_MAX) * (200) + (-100);
    }
  }
}

void S21Matrix::SimpleFill() {
  double a = 1;
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] = a++;
    }
  }
}

S21Matrix S21Matrix::MinorMatrix(unsigned int rm_row, unsigned int rm_col,
                                 S21Matrix &other) {
  unsigned int dest_rows = 0, dest_columns = 0;
  S21Matrix res(other.rows_ - 1, other.cols_ - 1);
  for (unsigned int i = 0; i < other.rows_; i++) {
    if (i == rm_row) {
      continue;
    }
    dest_columns = 0;
    for (unsigned int j = 0; j < other.cols_; j++) {
      if (j == rm_col) {
        continue;
      }
      res.matrix_[dest_rows][dest_columns] = other.matrix_[i][j];
      dest_columns++;
    }
    dest_rows++;
  }
  return res;
}

double S21Matrix::RecurseDeterminant(S21Matrix &other) {
  double result = 0.0;
  if (other.rows_ == 1) {
    result = other.matrix_[0][0];
  } else {
    for (unsigned int i = 0; i < other.cols_; i++) {
      S21Matrix temp = MinorMatrix(0, i, other);
      result += (pow(-1, i) * other.matrix_[0][i] * RecurseDeterminant(temp));
    }
  }
  return result;
}

void S21Matrix::CopyMatrix(S21Matrix other) {
  SetRows_(other.rows_);
  SetCols_(other.cols_);
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

void S21Matrix::FillMyVal(double arr[3][3]) {
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] = arr[i][j];
    }
  }
}

void S21Matrix::FreeObj() {
  if (matrix_) {
    for (unsigned int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }
  matrix_ = nullptr;
  cols_ = 0;
  rows_ = 0;
}

void S21Matrix::RowsColsSetter(unsigned int rows, unsigned int cols,
                               S21Matrix other) {
  for (unsigned int i = 0; i < rows; i++) {
    for (unsigned int j = 0; j < cols; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}
// __________________________________________________________________________________

// __________________________________________________________________________________
// Обработка исключений

Matrix_error::Matrix_error(std::string message) : message(message) {}

const char *Matrix_error::what() const noexcept { return message.c_str(); }

// ___________________________________________________________________________________

// ___________________________________________________________________________________
// Перегрузка операторов

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
  S21Matrix res(*this);
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    res.SumMatrix(other);
    return res;
  } else {
    throw Matrix_error("The sizes of the matrices do not match");
  }
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
  S21Matrix res(*this);
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    res.SubMatrix(other);
    return res;
  } else {
    throw Matrix_error("The sizes of the matrices do not match");
  }
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
  S21Matrix res(*this);
  if (rows_ == other.cols_) {
    res.MulMatrix(other);
    return res;
  } else {
    throw Matrix_error(
        "The number of columns of the first matrix does not match the number "
        "of rows of the second");
  }
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix res(*this);
  res.MulNumber(num);
  return res;
}

void S21Matrix::operator+=(const S21Matrix &other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    SumMatrix(other);
  } else {
    throw Matrix_error("The sizes of the matrices do not match");
  }
}

void S21Matrix::operator-=(const S21Matrix &other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    SubMatrix(other);
  } else {
    throw Matrix_error("The sizes of the matrices do not match");
  }
}

void S21Matrix::operator*=(const S21Matrix &other) {
  if (rows_ == other.cols_) {
    MulMatrix(other);
  } else {
    throw Matrix_error(
        "The number of columns of the first matrix does not match the number "
        "of rows of the second");
  }
}

void S21Matrix::operator*=(const double num) { MulNumber(num); }

bool S21Matrix::operator==(const S21Matrix &other) const {
  bool res = EqMatrix(other);
  return res;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  FreeObj();
  std::swap(rows_, other.rows_);
  std::swap(cols_, other.cols_);
  std::swap(matrix_, other.matrix_);
  return *this;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  CopyMatrix(other);
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  if ((i < (int)rows_ && j < (int)cols_) && (i >= 0 && j >= 0)) {
    return matrix_[i][j];
  } else {
    throw Matrix_error("Out of bound exception");
  }
}