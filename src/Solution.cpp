#include "../header/Solution.h"

Solution::Solution(const Solution& s) {
    y = new double[n = s.n];
    x = s.x;
    for (int i = 0; i < n; ++i) {
            y[i] = s.y[i];
    }
}
Solution::Solution(int dimension) {
    y = new double[dimension];
    for (int i = 0; i < dimension; ++i) y[i] = 0;
    x = 0;
    n = dimension;
}
Solution& Solution::operator=(const Solution& s) {
    if (this != &s) {
        delete[] y;
        n = s.n;
        x = s.x;
        y = new double[n];
        for (int i = 0; i < n; ++i) y[i] = s.y[i];
    }
    return *this;
}
Solution::~Solution() {
    delete[] y;
}

