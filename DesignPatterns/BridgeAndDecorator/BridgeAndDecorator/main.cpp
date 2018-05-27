#include "IMatrix.h"

#include <conio.h>
#include <stdio.h>

#define MSG(msg) printf("%s\n", ##msg)

int main(int argc, char **argv) {
	MSG("For Debug");

	MSG("To end, press any key");
	_getch();

	return EXIT_SUCCESS;
}