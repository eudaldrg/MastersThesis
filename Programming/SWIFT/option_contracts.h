#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>
#include "my_math.h"

class OptionContract
{
public:

    [[nodiscard]] virtual double GetPayoff(double K, std::size_t m, int k, int k1, int k2, std::size_t iota, bool is_call) const = 0;
};

class EuropeanOptionContract : public OptionContract
{
public:
//    [[nodiscard]] double GetExplicitPayoff
    [[nodiscard]] double GetPayoff(double K, std::size_t m, int k, int k1, int k2, std::size_t iota, bool is_call) const final;

private:
    [[nodiscard]] double GetC_j(std::size_t j, std::size_t J) const;

    [[nodiscard]] double GetI1(double a, double b, std::size_t j, std::size_t J, std::size_t m, int k) const;

    [[nodiscard]] double GetI2(double a, double b, std::size_t j, std::size_t J, std::size_t m, int k) const;
};

class CashOrNothingContract : public OptionContract
{
public:
    [[nodiscard]] double GetPayoff(double K, std::size_t m, int k, int k1, int k2, std::size_t iota, bool is_call) const final;
};
