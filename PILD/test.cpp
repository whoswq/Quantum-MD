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
#include <stdio.h>
#include <time.h>


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


void print(double (&s)[3][3], void (*print_vector_2)(double (&)[3][3]))
{
    print_vector_2(s);
}

int main()
{
    time_t t_start, t_end;
    t_start = time(NULL);
    double M_therm[3][3] = {{1, 2, 3},
                            {4, 5, 6},
                            {7, 8, 9}};
    double vector[3] = {0, 0, 0};
    print_vector_2(M_therm);
    for (int k = 0; k < 3; k++){
        M_therm[2][k] = vector[k];
    }
    print_vector_2(M_therm);
    vector[0] = 100;
    print_vector_2(M_therm);
    t_end = time(NULL);
    printf("time: %.2f s， %d\n", difftime(t_end, t_start), 10);
    int i; // 让程序不要运行完就退出
    cin >> i;
    return 0;
}