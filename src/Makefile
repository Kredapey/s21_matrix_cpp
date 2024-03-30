CC = g++
CC_FLAGS = -Wall -Werror -Wextra -std=c++17
OS = $(shell uname -s)
ifeq ($(OS), Linux)
	TEST_LIBS =  -lcheck -lm -lpthread -lrt -lsubunit -lgtest
else
	TEST_LIBS = -lcheck -lm -lpthread -lgtest
endif

.PHONY: all clean test s21_matrix_oop.a gcov_report

all: s21_matrix_oop.a test gcov_report

s21_matrix_oop.o: s21_matrix_oop.cpp
	$(CC) $(CC_FLAGS) -c s21_matrix_oop.cpp -g $(TEST_LIBS)

style: s21_matrix_oop.a
	clang-format -n *.cpp *.h

s21_matrix_oop.a: s21_matrix_oop.o
	$(CC) -c s21_matrix_oop.cpp $(TEST_LIBS)
	ar rcs s21_matrix_oop.a s21_matrix_oop.o

test: s21_matrix_oop.o s21_matrix_oop.a s21_matrix_oop_test.o
	$(CC) s21_matrix_oop.o s21_matrix_oop_test.o -o test s21_matrix_oop.a $(TEST_LIBS)
	./test

s21_matrix_oop_test.o: s21_matrix_oop.cpp
	$(CC) $(CC_FLAGS) -c s21_matrix_oop_test.cpp

gcov_report: s21_matrix_oop.a
	$(CC) $(CC_FLAGS) --coverage s21_matrix_oop_test.cpp s21_matrix_oop.cpp -o gcov_test s21_matrix_oop.a $(TEST_LIBS)
	chmod +x *
	./gcov_test
	lcov -t "gcov_test" -o gcov_test.info --no-external -c -d .
	genhtml -o report/ gcov_test.info
	open ./report/index.html

valgrind: s21_matrix_oop.a test
	CK_FORK=no valgrind --track-origins=yes --leak-check=full --leak-resolution=med --show-leak-kinds=all -s ./test

clean:
	rm -rf *.a
	rm -rf *.o
	rm -rf *.gcno
	rm -rf *.gcda
	rm -rf gcov_test
	rm -rf test
	rm -rf *.info
	rm -rf *.out
	rm -rf *.dSYM
	rm -R report


