#pragma once
/**
 * @brief Simple class.
 * @details A class that contatins solution of ODE at point x and its dimension.
 */
class Solution
{
public:
    Solution(const Solution&);
    Solution(int dimension);
    ~Solution(void);
    Solution& operator=(const Solution& s);
    double x;
    double* y;
    int n;
};