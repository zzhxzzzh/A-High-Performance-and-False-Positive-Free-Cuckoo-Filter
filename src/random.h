#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>
#include <stdexcept>

::std::vector<::std::uint32_t> GenerateRandom32(::std::size_t count) {
  if (count > static_cast<::std::size_t>(std::numeric_limits<::std::uint32_t>::max())) {
    throw ::std::length_error("Requested count exceeds the range of 32-bit unsigned integers");
  }

  ::std::vector<::std::uint32_t> numbers;
  numbers.reserve(count);

  ::std::random_device random;
  ::std::mt19937 generator(random());

  for (::std::uint32_t i = 0; i < count; ++i) {
    numbers.push_back(i);
  }

  ::std::shuffle(numbers.begin(), numbers.end(), generator);

  return numbers;
}
