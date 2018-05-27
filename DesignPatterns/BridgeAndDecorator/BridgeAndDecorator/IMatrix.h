#ifndef __IMATRIX_H__
#define __IMATRIX_H__

#include "IMatrixDrawer.h"

#include <assert.h>
#include <vector>
#include <algorithm> 

class IMatrix {
public:
	IMatrix(IMatrixDrawer *drawer) { this->drawer = drawer; }
	virtual size_t getSize() { return 0;  }
	virtual double getValue(size_t row, size_t col) { return 0; }
	virtual void setValue(size_t row, size_t col, double value) {}
	virtual IMatrix getComponent() { return *this; }
	std::string draw() { return drawer->drawMatrix(this); }
	virtual ~IMatrix() {}
private:
	IMatrixDrawer *drawer;
};

// ============================= Matrix ============================= //

class Matrix: public IMatrix {
public:
	Matrix(size_t size, IMatrixDrawer *drawer): IMatrix(drawer) {
		assert(size > 0 && "Invalid size of matrix");
		this->values.resize(size * size);
		this->size = size;
	}

	size_t getSize() override { return this->size; }

	double getValue(size_t row, size_t col) override {
		assert((row < size || row >= 0 || col < size || col >= 0) && "Invalid index");

		return this->values[col + row * size];
	} 

	void setValue(size_t row, size_t col, double value) override  {
		assert((row < size || row >= 0 || col < size || col >= 0) && "Invalid index");

		this->values[col + row * size] = value;
	} 
	
	IMatrix getComponent() override {
		return *this;
	}

	~Matrix() {}

private: 
	size_t size;
	std::vector<double> values;
};

// ============================= SparseMatrix ============================= //
typedef struct _SparseMatrix {
	size_t size;					
	std::vector<double> values;
	std::vector<int> columns;
	std::vector<int> rows;
}SMatrix;

class SparseMatrix : public IMatrix {
public:
	SparseMatrix(size_t size, IMatrixDrawer *drawer) : IMatrix(drawer) {
		assert(size > 0 && "Invalid size of matrix");
		this->matrix.size = size;
		this->matrix.values.resize(0);
		this->matrix.rows.resize(0);
		this->matrix.columns.resize(0);
	}

	size_t getSize() override {
		return matrix.size;
	}

	double getValue(size_t row, size_t col) override {
		assert((row < matrix.size || row >= 0 || col < matrix.size || col >= 0) && "Invalid index");

		for (size_t i = 0; i < matrix.rows.size(); ++i)
			if (matrix.rows[i] == row && matrix.columns[i] == col)
				return matrix.values[i];

		return 0;
	}

	void setValue(size_t row, size_t col, double value) override {
		assert((row < matrix.size || row >= 0 || col < matrix.size || col >= 0) && "Invalid index");
		for (size_t i = 0; i < matrix.rows.size(); ++i) {
			if (matrix.rows[i] == row && matrix.columns[i] == col) {
				matrix.values[i] = value;
				return;
			}
		}
		matrix.rows.push_back(static_cast<const int>(row));
		matrix.columns.push_back(static_cast<const int>(col));
		matrix.values.push_back(value);
	}

	IMatrix getComponent() override {
		return *this;
	}

	~SparseMatrix() {}

private:
	SMatrix matrix;
};

// ============================= LowerTriangularMatrix ============================= //

class LowerTriangularMatrix :public IMatrix{
public:
	LowerTriangularMatrix(IMatrix *component, IMatrixDrawer *drawer) : IMatrix(drawer) {
		this->component = component;
	}

	size_t getSize() override {
		return this->component->getSize();
	}

	double getValue(size_t row, size_t col) override {
		if (row < col)
			return 0;
		return component->getValue(row, col);
	}

	void setValue(size_t row, size_t col, double value) override {
		assert((row >= col) && "Set is available only for lower triangular part");
		this->component->setValue(row, col, value);
	}

	IMatrix getComponent() override {
		return *this->component;
	}

	~LowerTriangularMatrix() {}

private:
	IMatrix *component;
};

// ============================= MinorMatrix ============================= //

class MinorMatrix : public IMatrix {
public:
	MinorMatrix(IMatrix *component, IMatrixDrawer *drawer, const std::vector<int> &minor_idexs) : IMatrix(drawer) {
		this->items.resize(minor_idexs.size());
		for (size_t i = 0; i < minor_idexs.size(); ++i)
			items[i] = minor_idexs[i];

		std::sort(items.begin(), items.end(), [](int i, int j) { return (i < j); });

		this->component = component;
	}

	size_t getSize() override {
		return items.size();
	}
	
	double getValue(size_t row, size_t col) override {
		assert((row >= 0 || row < items.size() || col >= 0 || col < items.size()) && "Invalid index");
		return component->getValue(items[row], items[col]);
	}

	void setValue(size_t row, size_t col, double value) override {
		assert((row >= 0 || row < items.size() || col >= 0 || col < items.size()) && "Invalid index");
		component->setValue(items[row], items[col], value);
	}

	IMatrix getComponent() override {
		return *component;
	}

	~MinorMatrix() {}

private: 
	IMatrix *component;
	std::vector<int> items;
};

// ============================= TransposeMatrix ============================= /

class TransposeMatrix : public IMatrix {
public:
	TransposeMatrix(IMatrix *component, IMatrixDrawer *drawer) : IMatrix(drawer) {
		this->component = component;
	}

	size_t getSize() override {
		return component->getSize();
	}

	double getValue(size_t row, size_t col) override {
		return component->getValue(col, row);
	}

	void setValue(size_t row, size_t col, double value) override {
		component->setValue(col, row, value);
	}

	IMatrix getComponent() override {
		return *component;
	}

	~TransposeMatrix() {}
	
private:
	IMatrix *component;
};

#endif // __IMATRIX_H__