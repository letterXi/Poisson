#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main() {
    int N = 20;
    int d = 0;
    int N1, N2;
    int h = 4;
    if (h % 2 == 0) {
        d = h / 2 - 1;
    } else if (h % 2 != 0)
        d = h / 2;
    if (N % 2 != 0) {
        N2 = N - N / 2 + 1 + d;
        N1 = N - N2 + h + 1;
    } else if (N % 2 == 0) {
        N1 = N / 2 + 1 + d;
        N2 = N - N1 + h + 1;
    }
    
    double h_ = 1.0/(N-1);
    double x01 = 0.0;
    double y01 = 0.0;
    double x02 = (N1 - h - 1.0)*h_;
    double y02 = 0.0;

    std::ofstream file1("N1.txt");
    std::ofstream file2("N2.txt");
    file1 << "x,y,z" << std::endl;
    file2 << "x,y,z" << std::endl;

    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N1; i++)
        {
            double x = x01 + h_*i;
            double y = y01 + h_*j;
            file1 << x << ',' << y << ',' << 0 << std::endl;
        }
        for(int i = 0; i < N2; i++)
        {
            double x = x02 + h_*i;
            double y = y02 + h_*j;
            file2 << x << ',' << y << ',' << 0 << std::endl;
        }
    }

    file1.close();
    file2.close();
    return 0;
}
