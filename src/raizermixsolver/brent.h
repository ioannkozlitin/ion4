#include <functional>

namespace brent
{

double zero(double a, double b, double t, const std::function<double(double)> &f);

} // end namespace brent
