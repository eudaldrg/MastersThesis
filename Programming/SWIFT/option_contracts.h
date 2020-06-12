#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>
#include "my_math.h"

namespace Swift {
class SwiftParameters;
}

class OptionContract
{
public:

    [[nodiscard]] virtual double GetPayoffNonKComponent(int k, Swift::SwiftParameters const& params, bool is_call) const = 0;
    [[nodiscard]] virtual double GetPayoffKComponent(double K) const = 0;
    [[nodiscard]] double GetPayoff(double K, int k, Swift::SwiftParameters const& params, bool is_call) const
    {
        return GetPayoffKComponent(K) * GetPayoffNonKComponent(k, params, is_call);
    }
};

class EuropeanOptionContract : public OptionContract
{
public:
//    [[nodiscard]] double GetExplicitPayoff
    [[nodiscard]] double GetPayoffNonKComponent(int k, Swift::SwiftParameters const& params, bool is_call) const final;
    [[nodiscard]] double GetPayoffKComponent(double K) const final;

private:
    // TODO: Document why
    [[nodiscard]] double GetBoundedFrom(double from, bool is_call) const
    {
        return is_call ? std::max(from, 0.0) : from;
    }
    [[nodiscard]] double GetBoundedTo(double to, bool is_call) const
    {
        return is_call ? to : std::min(to, 0.0);
    }

    [[nodiscard]] double GetPayoffNonKComponentNewPaper(int k, Swift::SwiftParameters const& params, bool is_call) const;
    [[nodiscard]] double GetPayoffNonKComponentOldPaper(int k, Swift::SwiftParameters const& params, bool is_call) const;

    [[nodiscard]] double V(int k, Swift::SwiftParameters const& params, bool is_call) const;
    [[nodiscard]] double I1(double from, double to, int k, int j, Swift::SwiftParameters const& params) const;
    [[nodiscard]] double I2(double from, double to, int k, int j, Swift::SwiftParameters const& params) const;

};

class CashOrNothingContract : public OptionContract
{
public:
    [[nodiscard]] double GetPayoffNonKComponent(int k, Swift::SwiftParameters const& params, bool is_call) const final;
    [[nodiscard]] double GetPayoffKComponent(double K) const final;
};
