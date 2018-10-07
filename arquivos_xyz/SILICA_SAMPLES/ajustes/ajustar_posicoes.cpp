#include <iostream>
#include <fstream>
#include <limits>

int main(int argc, char *argv[]) {
	if (argc != 2) {
		std::cerr << "<$1> = filename.";
		return -1;
	}

	std::ifstream infile(argv[1], std::ios::binary);
	if (!infile) {
		std::cerr << "Não foi possível abrir o arquivo informado.\n";
		return -1;
	}

	double box = 17.832150;

	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // pula cabeçalho
	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // pula cabeçalho
	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // pula cabeçalho
	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // pula cabeçalho

	double x, y, z;
	while (infile >> x >> y >> z) {
		x = x + box;
		y = y + box;
		z = z + box;

		std::cout << x << "\t" << y << "\t" << z << "\n";
	}

	return 0;
}