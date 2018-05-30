#include "stdafx.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestPatternDecorator {		

	TEST_CLASS(MatrixDecoratorTest) {
	public:
		MatrixDecoratorTest() {
			matrix = new Matrix(size, drawer);
			for (int i = 0; i < size; i++)
				for (int j = 0; j < size; j++)
					matrix->setValue(i, j, matrix_expected[i + j * size]);
		}

		~MatrixDecoratorTest() {
			if(matrix != nullptr)
				delete matrix;
			if (drawer != nullptr)
				delete drawer;
		}

		TEST_METHOD(test_TransposeLowerTriangularMinor) {
			const std::vector<int> minor_idexs = { 0, 2 };

			IMatrix *m = new TransposeMatrix(matrix, drawer);
			m = new LowerTriangularMatrix(m, drawer);
			m = new MinorMatrix(m, drawer, minor_idexs);

			Assert::AreEqual(matrix_expected[0 + 2 * size], m->getValue(1, 0), L"Error! Invalid value");

			delete m;
		}

		TEST_METHOD(test_TransposeMinor) {
			const std::vector<int> minor_idexs = { 1, 2 };

			IMatrix *m = new TransposeMatrix(matrix, drawer);
			m = new MinorMatrix(m, drawer, minor_idexs);

			Assert::AreEqual(matrix_expected[2 + 1 * size], m->getValue(0, 1), L"Error! Invalid value");

			delete m;
		}

		TEST_METHOD(test_LowerTriangularMinor) {
			const std::vector<int> minor_idexs = { 1, 2 };

			IMatrix *m = new LowerTriangularMatrix(matrix, drawer);
			m = new MinorMatrix(m, drawer, minor_idexs);

			double a = m->getValue(1, 1);
			a++;

			Assert::AreEqual(matrix_expected[2 + 2 * size], m->getValue(1, 1), L"Error! Invalid value");

			delete m;
		}

		TEST_METHOD(test_LowerTriangularTranspose) {
			IMatrix *m = new LowerTriangularMatrix(matrix, drawer);
			m = new TransposeMatrix(m, drawer);

			Assert::AreEqual(matrix_expected[2 + 1 * size], m->getValue(1, 2), L"Error! Invalid value");

			delete m;
		}

		TEST_METHOD(test_MinorLowerTriangular) {
			const std::vector<int> minor_idexs = { 0, 2 };

			IMatrix *m = new MinorMatrix(matrix, drawer, minor_idexs);
			m = new LowerTriangularMatrix(m, drawer);

			Assert::AreEqual(matrix_expected[2 + 0 * size], m->getValue(1, 0), L"Error! Invalid value");

			delete m;
		}

	private:
		IMatrix *matrix = nullptr;
		IMatrixDrawer *drawer = nullptr;
		std::vector<double> matrix_expected = {	1, 2, 3, 
												4, 5, 6,
												7, 8, 9 };
		int size = 3;
	};
}