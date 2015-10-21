#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cstdlib>
using namespace std;

int atoi_min(char const *value, int min, char const *parameter, char const *name)
{
    std::istringstream iss(value);
    int i;
    char c;
    if (!(iss >> i) || iss.get(c))
    {
        cerr << "error: argument of " << parameter << " must be of type <int>, and greater than or equal to " << min << endl
             << "Check " << name << " --help' for more information." << endl;
        std::exit(1);
    }

    if (i < min)
    {
        cerr << "error: argument of " << parameter << " must be greater than or equal to " << min << endl
             << "Check `" << name << " --help' for more information." << endl;
        std::exit(1);
    }
    return i;
}


int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cerr << "usage: " << argv[0] << " [matrix width] [matrix height]" << endl;
        return 1;
    }

    unsigned seed = atoi_min(argv[1], 0, "seed", argv[0]);
    unsigned width = atoi_min(argv[2], 4, "matrix width", argv[0]);
    unsigned height = atoi_min(argv[3], 4, "matrix height", argv[0]);

    srand(time(NULL)+seed);
    while (height-- > 0)
    {
        for (unsigned i = 0; i < width; ++i)
            cout << (int)rand()%2;
        cout << endl;
    }
    return 0;
}
