#include "stdafx.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestPatternBridge {
	
	TEST_CLASS(MatrixBridgeTest) {
	public:
	MatrixBridgeTest() {}

	~MatrixBridgeTest() {
		if (matrix != nullptr)
			delete matrix;
	}

	TEST_METHOD(test_DrawMatrixInConsole) {
		IMatrixDrawAPI *drawAPI = new ConsoleMatrixDrawAPI();
		IMatrixDrawer *drawer = new MatrixDrawer(drawAPI);
		matrix = new Matrix(size, drawer);
		str_expected = "Matrix Console API: Draw Outer Boundary;Console API: Draw Inner Boundary;Console API: Draw Elements";

		std::string str_test = matrix->draw();

		Assert::AreEqual(str_test, str_expected);
		Assert::IsNotNull(matrix, L"Error! Pointer on matrix is null value");
	}

	TEST_METHOD(test_DrawMatrixInGraphical) {
		IMatrixDrawAPI *drawAPI = new GraphicalMatrixDrawAPI();
		IMatrixDrawer *drawer = new MatrixDrawer(drawAPI);
		matrix = new Matrix(size, drawer);

		str_expected = "Matrix Graphical API: Draw Outer Boundary;Graphical API: Draw Inner Boundary;Graphical API: Draw Elements";

		std::string str_test = matrix->draw();

		Assert::AreEqual(str_test, str_expected);
		Assert::IsNotNull(matrix, L"Error! Pointer on matrix is null value");
	}

	TEST_METHOD(test_DrawMatrixInHTML) {
		IMatrixDrawAPI *drawAPI = new HTMLMatrixDrawAPI();
		IMatrixDrawer *drawer = new MatrixDrawer(drawAPI);
		matrix = new Matrix(size, drawer);

		str_expected = "Matrix HTML API: Draw Outer Boundary;HTML API: Draw Inner Boundary;HTML API: Draw Elements";

		std::string str_test = matrix->draw();

		Assert::AreEqual(str_test, str_expected);
		Assert::IsNotNull(matrix, L"Error! Pointer on matrix is null value");
	}

	TEST_METHOD(test_DrawSparseMatrixInConsole) {
		IMatrixDrawAPI *drawAPI = new ConsoleMatrixDrawAPI();
		IMatrixDrawer *drawer = new SparseMatrixDrawer(drawAPI);
		matrix = new SparseMatrix(size, drawer);

		str_expected = "SparseMatrix Console API: Draw Outer Boundary;Console API: Draw Non Zero Elements";

		std::string str_test = matrix->draw();

		Assert::AreEqual(str_test, str_expected);
		Assert::IsNotNull(matrix, L"Error! Pointer on matrix is null value");
	}

	TEST_METHOD(test_DrawSparseMatrixInGraphical) {
		IMatrixDrawAPI *drawAPI = new GraphicalMatrixDrawAPI();
		IMatrixDrawer *drawer = new SparseMatrixDrawer(drawAPI);
		matrix = new SparseMatrix(size, drawer);

		str_expected = "SparseMatrix Graphical API: Draw Outer Boundary;Graphical API: Draw Non Zero Elements";

		std::string str_test = matrix->draw();

		Assert::AreEqual(str_test, str_expected);
		Assert::IsNotNull(matrix, L"Error! Pointer on matrix is null value");
	}

	TEST_METHOD(test_DrawSparseMatrixInHTML) {
		IMatrixDrawAPI *drawAPI = new HTMLMatrixDrawAPI();
		IMatrixDrawer *drawer = new SparseMatrixDrawer(drawAPI);
		matrix = new SparseMatrix(size, drawer);

		str_expected = "SparseMatrix HTML API: Draw Outer Boundary;HTML API: Draw Non Zero Elements";

		std::string str_test = matrix->draw();

		Assert::AreEqual(str_test, str_expected);
		Assert::IsNotNull(matrix, L"Error! Pointer on matrix is null value");
	}

	private:
		IMatrix * matrix = nullptr;
		const int size = 10;
		std::string str_expected;
	};
}