#ifndef __IMATRIX_DRAW_API_H__
#define __IMATRIX_DRAW_API_H__

#include <iostream>
#include <string>

class IMatrixDrawAPI {
	public:
		virtual std::string drawOuterBoundary() = 0;
		virtual std::string drawInnerBoundary() = 0;
		virtual std::string drawElements() = 0;
		virtual std::string drawNonZeroElements() = 0;
	};

// ============================= ConsoleDrawMatrixAPI ============================= //

class ConsoleMatrixDrawAPI : public IMatrixDrawAPI {
public:
	std::string drawOuterBoundary() override {
		return "Console API: Draw Outer Boundary";
	}

	std::string drawInnerBoundary() override {
		return "Console API: Draw Inner Boundary";
	}

	std::string drawElements() {
		return "Console API: Draw Elements";
	}

	std::string drawNonZeroElements() {
		return "Console API: Draw Non Zero Elements";
	}
};

// ============================= GraphicalMatrixDrawAPI ============================= //

class GraphicalMatrixDrawAPI : public IMatrixDrawAPI {
public:
	std::string drawOuterBoundary() override {
		return "Graphical API: Draw Outer Boundary";
	}

	std::string drawInnerBoundary() override {
		return "Graphical API: Draw Inner Boundary";
	}

	std::string drawElements() {
		return "Graphical API: Draw Elements";
	}

	std::string drawNonZeroElements() {
		return "Graphical API: Draw Non Zero Elements";
	}
};

// ============================= HTMLMatrixDrawAPI ============================= //

class HTMLMatrixDrawAPI : public IMatrixDrawAPI {
public:
	std::string drawOuterBoundary() override {
		return "HTML API: Draw Outer Boundary";
	}

	std::string drawInnerBoundary() override {
		return "HTML API: Draw Inner Boundary";
	}

	std::string drawElements() {
		return "HTML API: Draw Elements";
	}

	std::string drawNonZeroElements() {
		return "HTML API: Draw Non Zero Elements";
	}
};

#endif // __IMATRIX_DRAW_API_HPP__