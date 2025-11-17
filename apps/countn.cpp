#include <string>
#include <iostream>
#include <iomanip>

int main()
{
    int N = 7;
    int d = 0;
    int N1, N2;
    for(int h = 1; h < N-2; h++)
    {
        if (h % 2 == 0)
        {
            d = h / 2 - 1;
        }
        else if (h % 2 != 0)
            d = h / 2;
        if (N % 2 != 0)
        {
            N2 = N - N/2 + 1 + d;
            N1 = N - N2 + h + 1;
        }
        else if (N % 2 == 0)
        {
            N1 = N/2 + 1 + d;
            N2 = N - N1 + h + 1;
        }
        std::cout << std::string(N1, '-') << std::endl;
        std::cout << std::setw(N) << std::setfill(' ') << std::string(N2, '-') << std::endl << std::endl;
        if (N % 2 == 0)
        {
            if (h % 2 != 0)
                std::cout << (N1 == N2) << std::endl;
        } 
        if (N % 2 != 0)
        {
            if (h % 2 == 0)
                std::cout << (N1 == N2) << std::endl;
        }
    }
    return 0;
}
