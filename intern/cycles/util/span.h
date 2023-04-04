#ifndef SPAN_H
#define SPAN_H

#include "log.h"

namespace p_std {
namespace span_internal {
// Wrappers for access to container data pointers.
template<typename C> inline constexpr auto get_data_impl(C &c, char) noexcept -> decltype(c.data())
{
  return c.data();
}

template<typename C> inline constexpr auto get_data(C &c) noexcept -> decltype(get_data_impl(c, 0))
{
  return get_data_impl(c, 0);
}

// Detection idioms for size() and data().
template<typename C>
using HasSize = std::is_integral<typename std::decay_t<decltype(std::declval<C &>().size())>>;

// We want to enable conversion from vector<T*> to span<const T* const> but
// disable conversion from vector<Derived> to span<Base>. Here we use
// the fact that U** is convertible to Q* const* if and only if Q is the same
// type or a more cv-qualified version of U.  We also decay the result type of
// data() to avoid problems with classes which have a member function data()
// which returns a reference.
template<typename T, typename C>
using HasData =
    std::is_convertible<typename std::decay_t<decltype(get_data(std::declval<C &>()))> *,
                        T *const *>;

}  // namespace span_internal

inline constexpr std::size_t dynamic_extent = -1;

// span implementation partially based on absl::span from Google's Abseil library.
template<typename T> class span {
 public:
  // Used to determine whether a span can be constructed from a container of
  // type C.
  template<typename C>
  using EnableIfConvertibleFrom = typename std::enable_if_t<span_internal::HasData<T, C>::value &&
                                                            span_internal::HasSize<C>::value>;

  // Used to SFINAE-enable a function when the slice elements are const.
  template<typename U>
  using EnableIfConstView = typename std::enable_if_t<std::is_const<T>::value, U>;

  // Used to SFINAE-enable a function when the slice elements are mutable.
  template<typename U>
  using EnableIfMutableView = typename std::enable_if_t<!std::is_const<T>::value, U>;

  using value_type = typename std::remove_cv_t<T>;
  using iterator = T *;
  using const_iterator = const T *;

  span() : ptr(nullptr), n(0)
  {
  }

  span(T *ptr, size_t n) : ptr(ptr), n(n)
  {
  }

  template<size_t N> span(T (&a)[N]) : span(a, N)
  {
  }

  span(std::initializer_list<value_type> v) : span(v.begin(), v.size())
  {
  }

  // Explicit reference constructor for a mutable `span<T>` type. Can be
  // replaced with make_span() to infer the type parameter.
  template<typename V,
           typename X = EnableIfConvertibleFrom<V>,
           typename Y = EnableIfMutableView<V>>
  explicit span(V &v) noexcept : span(v.data(), v.size())
  {
  }

  // Hack: explicit constructors for std::vector to work around warnings
  // about calling a host function (e.g. vector::size()) form a
  // host+device function (the regular span constructor.)
  template<typename V> span(std::vector<V> &v) noexcept : span(v.data(), v.size())
  {
  }
  template<typename V> span(const std::vector<V> &v) noexcept : span(v.data(), v.size())
  {
  }

  // Implicit reference constructor for a read-only `span<const T>` type
  template<typename V, typename X = EnableIfConvertibleFrom<V>, typename Y = EnableIfConstView<V>>
  constexpr span(const V &v) noexcept : span(v.data(), v.size())
  {
  }

  iterator begin()
  {
    return ptr;
  }

  iterator end()
  {
    return ptr + n;
  }

  const_iterator begin() const
  {
    return ptr;
  }

  const_iterator end() const
  {
    return ptr + n;
  }

  T &operator[](size_t i)
  {
    return ptr[i];
  }

  const T &operator[](size_t i) const
  {
    return ptr[i];
  }

  size_t size() const
  {
    return n;
  };

  bool empty() const
  {
    return size() == 0;
  }

  T *data()
  {
    return ptr;
  }

  const T *data() const
  {
    return ptr;
  }

  T front() const
  {
    return ptr[0];
  }

  T back() const
  {
    return ptr[n - 1];
  }

  void remove_prefix(size_t count)
  {
    // assert(size() >= count);
    ptr += count;
    n -= count;
  }

  void remove_suffix(size_t count)
  {
    // assert(size() > = count);
    n -= count;
  }

  span subspan(size_t pos, size_t count = dynamic_extent)
  {
    size_t np = count < (size() - pos) ? count : (size() - pos);
    return span(ptr + pos, np);
  }

 private:
  T *ptr;
  size_t n;
};

template<int &...ExplicitArgumentBarrier, typename T>
inline constexpr span<T> make_span(T *ptr, size_t size) noexcept
{
  return span<T>(ptr, size);
}

template<int &...ExplicitArgumentBarrier, typename T>
inline span<T> make_span(T *begin, T *end) noexcept
{
  return span<T>(begin, end - begin);
}

template<int &...ExplicitArgumentBarrier, typename T>
inline span<T> make_span(std::vector<T> &v) noexcept
{
  return span<T>(v.data(), v.size());
}

template<int &...ExplicitArgumentBarrier, typename C>
inline constexpr auto make_span(C &c) noexcept
    -> decltype(make_span(span_internal::get_data(c), c.size()))
{
  return make_span(span_internal::get_data(c), c.size());
}

template<int &...ExplicitArgumentBarrier, typename T, size_t N>
inline constexpr span<T> make_span(T (&array)[N]) noexcept
{
  return span<T>(array, N);
}

template<int &...ExplicitArgumentBarrier, typename T>
inline constexpr span<const T> make_const_span(T *ptr, size_t size) noexcept
{
  return span<const T>(ptr, size);
}

template<int &...ExplicitArgumentBarrier, typename T>
inline span<const T> make_const_span(T *begin, T *end) noexcept
{
  return span<const T>(begin, end - begin);
}

template<int &...ExplicitArgumentBarrier, typename T>
inline span<const T> make_const_span(const std::vector<T> &v) noexcept
{
  return span<const T>(v.data(), v.size());
}

template<int &...ExplicitArgumentBarrier, typename C>
inline constexpr auto make_const_span(const C &c) noexcept -> decltype(make_span(c))
{
  return make_span(c);
}

template<int &...ExplicitArgumentBarrier, typename T, size_t N>
inline constexpr span<const T> make_const_span(const T (&array)[N]) noexcept
{
  return span<const T>(array, N);
}

}  // namespace p_std

#endif  // SPAN_H
