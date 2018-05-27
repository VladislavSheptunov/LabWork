#ifndef __IMATRIX_DRAWER_H__
#define __IMATRIX_DRAWER_H__

#include "IMatrixDrawAPI.h"

class IMatrix; // Predestination

class IMatrixDrawer {
public:
	IMatrixDrawer(IMatrixDrawAPI *drawAPI) { this->drawAPI = drawAPI;}
	virtual std::string drawMatrix(IMatrix *matrix) { return ""; }
	virtual ~IMatrixDrawer() {}
protected: 
	IMatrixDrawAPI *drawAPI;
};

// ============================= MatrixDrawer ============================= //

class MatrixDrawer : public  IMatrixDrawer {
public:
	MatrixDrawer(IMatrixDrawAPI *drawAPI) : IMatrixDrawer(drawAPI) {}

	std::string drawMatrix(IMatrix *matrix) override {
		std::string str = "Matrix " + drawAPI->drawOuterBoundary() + ";";
		str += drawAPI->drawInnerBoundary() + ";";
		str += drawAPI->drawElements();
		
		return str;
	}

};

// ============================= SparseMatrixDrawer ============================= //

class SparseMatrixDrawer : public  IMatrixDrawer {
public:
	SparseMatrixDrawer(IMatrixDrawAPI *drawAPI) : IMatrixDrawer(drawAPI) {}

	std::string drawMatrix(IMatrix *matrix) override {
		std::string str = "SparseMatrix " + drawAPI->drawOuterBoundary() + ";";
		str += drawAPI->drawNonZeroElements();

		return str;
	}
};

#endif // __IMATRIX_DRAWER_H__