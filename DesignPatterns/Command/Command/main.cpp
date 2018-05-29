#include "CommandManager.h"

#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#define MSG(msg) printf("%s\n", ##msg)

int main(int argc, char **argv) {
	MSG("For Debug");

	MSG("To end, press any key");
	_getch();

	return EXIT_SUCCESS;
}