

#include <iostream>
    using namespace std;

int main ()
{
    double x[4];
//    int array_length = end(arr) - begin(arr);
int *a;
std::cout << "Length of array = " << (sizeof(a)/sizeof(*a)) << std::endl;
int *p = new int[7];
std::cout << "Length of array = " << (sizeof(p)/sizeof(*p)) << std::endl;

return 0; }

