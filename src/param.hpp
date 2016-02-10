#include <vector>

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &vec) {
  typename std::vector<T>::const_iterator it;
  for (it = vec.begin(); it != vec.end(); ++it) os << *it << " ";
  return os;
}
