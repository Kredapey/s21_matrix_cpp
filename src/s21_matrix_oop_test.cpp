#include "s21_matrix_oop.h"

#include <gtest/gtest.h>

using namespace std;

TEST(construct_tests, default_const) {
  S21Matrix res;
  ASSERT_EQ(res.GetRows_(), 3);
  ASSERT_EQ(res.GetCols_(), 3);
}

TEST(construct_tests, custom_const) {
  S21Matrix res(4, 4);
  ASSERT_EQ(res.GetRows_(), 4);
  ASSERT_EQ(res.GetCols_(), 4);
}

TEST(construct_tests, copy_const) {
  S21Matrix mtrx;
  mtrx.SimpleFill();
  S21Matrix mtrx_copy(mtrx);
  ASSERT_TRUE(mtrx == mtrx_copy);
}

TEST(construct_tests, move_const) {
  S21Matrix mtrx;
  mtrx.SimpleFill();
  S21Matrix mtrx_move(move(mtrx));
  ASSERT_EQ(mtrx.GetRows_(), 0);
  ASSERT_EQ(mtrx.GetCols_(), 0);
  ASSERT_EQ(mtrx_move.GetRows_(), 3);
  ASSERT_EQ(mtrx_move.GetCols_(), 3);
  ASSERT_DOUBLE_EQ(mtrx_move(0, 0), 1);
  ASSERT_DOUBLE_EQ(mtrx_move(2, 2), 9);
}

TEST(set_tests, rows_test) {
  S21Matrix mtrx(1, 1);
  mtrx.SetRows_(2);
  ASSERT_EQ(mtrx.GetRows_(), 2);
}

TEST(set_tests, cols_test) {
  S21Matrix mtrx(1, 1);
  mtrx.SetCols_(2);
  ASSERT_EQ(mtrx.GetCols_(), 2);
}

TEST(eq_test, true_eq) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_TRUE(mtrx1.EqMatrix(mtrx2));
}

TEST(eq_test, false_eq) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.FillMatrix();
  ASSERT_FALSE(mtrx1.EqMatrix(mtrx2));
}

TEST(calculation_tests, sum_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  double arr[3][3] = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
  res.FillMyVal(arr);
  mtrx1.SumMatrix(mtrx2);
  ASSERT_TRUE(mtrx1 == res);
}

TEST(calculation_tests, sum_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1.SumMatrix(mtrx2), Matrix_error);
}

TEST(calculation_tests, sub_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  double arr[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  res.FillMyVal(arr);
  mtrx1.SubMatrix(mtrx2);
  ASSERT_TRUE(mtrx1 == res);
}

TEST(calculation_tests, sub_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1.SubMatrix(mtrx2), Matrix_error);
}

TEST(calculation_tests, mul_num_test) {
  S21Matrix mtrx1;
  mtrx1.SimpleFill();
  S21Matrix res;
  double arr[3][3] = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
  res.FillMyVal(arr);
  mtrx1.MulNumber(2);
  ASSERT_TRUE(mtrx1 == res);
}

TEST(calculation_tests, mul_matrix_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  double arr[3][3] = {{30, 36, 42}, {66, 81, 96}, {102, 126, 150}};
  res.FillMyVal(arr);
  mtrx1.MulMatrix(mtrx2);
  ASSERT_TRUE(mtrx1 == res);
}

TEST(calculation_tests, mul_matrix_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1.MulMatrix(mtrx2), Matrix_error);
}

TEST(matrix_ops, transpose_test) {
  S21Matrix mtrx1;
  S21Matrix myres;
  S21Matrix result;
  mtrx1.SimpleFill();
  myres = mtrx1.Transpose();
  double arr[3][3] = {{1, 4, 7}, {2, 5, 8}, {3, 6, 9}};
  result.FillMyVal(arr);
  ASSERT_TRUE(myres == result);
}

TEST(matrix_ops, determinant_test) {
  S21Matrix mtrx1;
  mtrx1.SimpleFill();
  double res = mtrx1.Determinant();
  ASSERT_DOUBLE_EQ(res, 0);
}

TEST(matrix_ops, determinant_exc_test) {
  S21Matrix mtrx1(3, 4);
  mtrx1.SimpleFill();
  ASSERT_THROW(mtrx1.Determinant();, Matrix_error);
}

TEST(matrix_ops, calc_comp_many_rows_test) {
  S21Matrix mtrx;
  double arr[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};
  mtrx.FillMyVal(arr);
  S21Matrix result;
  double arr_res[3][3] = {{0, 10, -20}, {4, -14, 8}, {-8, -2, 4}};
  result.FillMyVal(arr_res);
  S21Matrix my_res;
  my_res = mtrx.CalcComplements();
  ASSERT_TRUE(my_res == result);
}

TEST(matrix_ops, calc_comp_many_rows_exc_test) {
  S21Matrix mtrx(3, 4);
  mtrx.SimpleFill();
  ASSERT_THROW(mtrx.CalcComplements(), Matrix_error);
}

TEST(matrix_ops, calc_comp_one_row_test) {
  S21Matrix mtrx(1, 1);
  mtrx.SimpleFill();
  S21Matrix my_res;
  my_res = mtrx.CalcComplements();
  ASSERT_TRUE(my_res(0, 0) == 1);
}

TEST(matrix_ops, inv_mtrx_many_rows_exc_test) {
  S21Matrix mtrx;
  mtrx.SimpleFill();
  ASSERT_THROW(mtrx.InverseMatrix(), Matrix_error);
}

TEST(matrix_ops, inv_mtrx_many_rows_test) {
  S21Matrix mtrx(3, 3);
  double arr[3][3] = {{2, 5, 7}, {6, 3, 4}, {5, -2, -3}};
  mtrx.FillMyVal(arr);
  S21Matrix result;
  double arr_res[3][3] = {{1, -1, 1}, {-38, 41, -34}, {27, -29, 24}};
  result.FillMyVal(arr_res);
  S21Matrix my_res;
  my_res = mtrx.InverseMatrix();
  ASSERT_TRUE(my_res == result);
}

TEST(matrix_ops, inv_mtrx_one_row_test) {
  S21Matrix mtrx(1, 1);
  mtrx(0, 0) = 10.0;
  S21Matrix my_res;
  my_res = mtrx.InverseMatrix();
  ASSERT_TRUE(my_res(0, 0) == 0.1);
}

TEST(matrix_ops, inv_mtrx_one_row_exc_test) {
  S21Matrix mtrx(1, 1);
  mtrx(0, 0) = 0.0;
  ASSERT_THROW(mtrx.InverseMatrix(), Matrix_error);
}

TEST(operators_overload, plus_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  S21Matrix my_res;
  double arr[3][3] = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
  res.FillMyVal(arr);
  my_res = mtrx1 + mtrx2;
  ASSERT_TRUE(my_res == res);
}

TEST(operators_overload, plus_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1 + mtrx2, Matrix_error);
}

TEST(operators_overload, minus_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  S21Matrix my_res;
  double arr[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  res.FillMyVal(arr);
  my_res = mtrx1 - mtrx2;
  ASSERT_TRUE(my_res == res);
}

TEST(operators_overload, minus_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1 - mtrx2, Matrix_error);
}

TEST(operators_overload, mul_numb_test) {
  S21Matrix mtrx1;
  mtrx1.SimpleFill();
  S21Matrix res;
  S21Matrix my_res;
  double arr[3][3] = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
  res.FillMyVal(arr);
  my_res = mtrx1 * 2;
  ASSERT_TRUE(my_res == res);
}

TEST(operators_overload, mul_mtrx_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  S21Matrix my_res;
  double arr[3][3] = {{30, 36, 42}, {66, 81, 96}, {102, 126, 150}};
  res.FillMyVal(arr);
  my_res = mtrx1 * mtrx2;
  ASSERT_TRUE(my_res == res);
}

TEST(operators_overload, mul_mtrx_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1 * mtrx2, Matrix_error);
}

TEST(operators_overload, mul_eq_mtrx_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  double arr[3][3] = {{30, 36, 42}, {66, 81, 96}, {102, 126, 150}};
  res.FillMyVal(arr);
  mtrx1 *= mtrx2;
  ASSERT_TRUE(mtrx1 == res);
}

TEST(operators_overload, mul_eq_mtrx_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1 *= mtrx2, Matrix_error);
}

TEST(operators_overload, mul_eq_numb_test) {
  S21Matrix mtrx1;
  mtrx1.SimpleFill();
  S21Matrix res;
  double arr[3][3] = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
  res.FillMyVal(arr);
  mtrx1 *= 2;
  ASSERT_TRUE(mtrx1 == res);
}

TEST(operators_overload, minus_eq_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  double arr[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  res.FillMyVal(arr);
  mtrx1 -= mtrx2;
  ASSERT_TRUE(mtrx1 == res);
}

TEST(operators_overload, minus_eq_mtrx_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1 -= mtrx2, Matrix_error);
}

TEST(operators_overload, plus_eq_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  S21Matrix res;
  double arr[3][3] = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
  res.FillMyVal(arr);
  mtrx1 += mtrx2;
  ASSERT_TRUE(mtrx1 == res);
}

TEST(operators_overload, plus_eq_mtrx_exc_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2(4, 4);
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_THROW(mtrx1 += mtrx2, Matrix_error);
}

TEST(operators_overload, equal_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx1.SimpleFill();
  mtrx2.SimpleFill();
  ASSERT_TRUE(mtrx1 == mtrx2);
}

TEST(operators_overload, copy_eq_test) {
  S21Matrix mtrx1;
  mtrx1.SimpleFill();
  S21Matrix mtrx2;
  mtrx2 = mtrx1;
  ASSERT_TRUE(mtrx1 == mtrx2);
}

TEST(operator_overload, move_eq_test) {
  S21Matrix mtrx1;
  S21Matrix mtrx2;
  mtrx2.SimpleFill();
  S21Matrix res;
  res.SimpleFill();
  mtrx1 = std::move(mtrx2);
  ASSERT_TRUE(mtrx1 == res);
  ASSERT_FALSE(mtrx1 == mtrx2);
}

TEST(operators_overload, index_test) {
  S21Matrix mtrx;
  mtrx.SimpleFill();
  ASSERT_DOUBLE_EQ(mtrx(2, 2), 9);
}

TEST(operators_overload, index_exc_test) {
  S21Matrix mtrx;
  mtrx.SimpleFill();
  ASSERT_THROW(mtrx(3, 3), Matrix_error);
}

TEST(exceptions, what_func_test) {
  try {
    S21Matrix mtrx(3, 3);
    S21Matrix mtrx1(3, 4);
    mtrx.SimpleFill();
    mtrx1.SimpleFill();
    S21Matrix res;
    res = mtrx + mtrx1;
  } catch (Matrix_error& err) {
    ASSERT_STREQ("The sizes of the matrices do not match", err.what());
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}