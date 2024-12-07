#include <iostream>
#include "sno-de.h"

int main(int argc, char** argv)
{
    if (argc != 1 && argc != 13) {
	cerr << "The number of parameters is incorrect." << endl;
	exit(1);
    }
	
    if (argc == 1) { // default
	sno search;
	search.run();
    }
    else if (argc == 13){ // ./main_sno 1 1000 30 10 10 10 10 10 0.2 0.2 10 1
	sno search(atoi(argv[1]),
		   atoi(argv[2]),
		   atoi(argv[3]),
		   atoi(argv[4]),
		   atoi(argv[5]),
		   atoi(argv[6]),
		   atoi(argv[7]),
		   atoi(argv[8]),
		   atof(argv[9]),
		   atof(argv[10]),
		   atoi(argv[11]),
		   atoi(argv[12])
		   );
	search.run();
    }

    return 0;
}
