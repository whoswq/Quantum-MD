#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <windows.h>

//using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

void print_vector_2(double (&vec)[3][3])
{
    cout << "[";
    for (int i = 0; i < 3; i++)
    {
        cout << "[";
        for (int j = 0; j < 3; j++)
        {
            cout << vec[i][j] << ",";
        }
        cout << "], ";
    }
    cout << "]" << endl;
}

double plus(double a, double b)
{
    return a + b;
}

void print(double (&s)[3][3], void (*print_vector_2)(double (&)[3][3]))
{
    print_vector_2(s);
}

int main()
{
    double M_therm[3][3] = {{1, 2, 3},
                            {4, 5, 6},
                            {7, 8, 9}};
    const double x = 1.01;
    cout << plus(x, x) << endl;
    print(M_therm, print_vector_2);

    int i; // 让程序不要运行完就退出
    cin >> i;
    return 0;
}