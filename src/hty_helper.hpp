/**
 * @file hty_helper.hpp
 * @author Sarutch Supaibulpipat (Pokpong) {ssupaibu@cmkl.ac.th}
 * @brief
 *  Helper functions for HTY file format programs.
 * @version 0.1
 * @date 2024-11-29
 */
#include <algorithm>
#include <list>
#include <optional>
#include <string>
#include <variant>

/**
 * @brief A bidirectional map implementation.
 *
 * @tparam T1 Type of the first key.
 * @tparam T2 Type of the second key.
 */
template <typename T1, typename T2> class BidirectMap {
  private:
    std::map<T1, T2> map1to2;
    std::map<T2, T1> map2to1;

  public:
    BidirectMap(std::vector<std::pair<T1, T2>> pairs)
    {
        for (std::pair<T1, T2> pair : pairs) {
            map1to2.insert({pair.first, pair.second});
            map2to1.insert({pair.second, pair.first});
        }
    }

    /**
     * @brief Retrieves the value mapped to the given key.
     * @param val The key of type `T1`.
     * @return The value of type `T2` corresponding to the key.
     * @throws std::out_of_range if the key is not found.
     */
    T2 get(T1 val) { return map1to2.at(val); }

    /**
     * @brief Retrieves the key mapped to the given value.
     * @param val The value of type `T2`.
     * @return The key of type `T1` corresponding to the value.
     * @throws std::out_of_range if the value is not found.
     */
    T1 get_reverse(T2 val) { return map2to1.at(val); }
};

enum class Num32Type { I, F };

static const std::vector<std::pair<Num32Type, std::string>>
    num32_type_and_str_pairs = {{Num32Type::I, "int"}, {Num32Type::F, "float"}};

BidirectMap num32_type_str_bimap(num32_type_and_str_pairs);

typedef std::variant<std::int32_t, float> Num32;

typedef std::vector<std::string> StrVec_t;
typedef std::vector<std::uint32_t> UInt32Vec_t;

#define get_num32(num, type)                                                   \
    ((type == Num32Type::I) ? std::get<std::int32_t>(num)                      \
                            : std::get<float>(num))
/**
 * @brief Slices a portion of a vector based on the specified range.
 *
 * @tparam T Type of the vector elements.
 * @param v The input vector.
 * @param start The start index (inclusive).
 * @param end The end index (exclusive).
 * @return A new vector containing elements from `start` to `end`.
 */
template <typename T>
std::vector<T> slice_vec(std::vector<T> v, std::size_t start, std::size_t end)
{
    std::vector<T> ref(v.begin() + start, v.begin() + end);
    return ref;
}

/**
 * @brief Retrieves a subset of a vector based on specified indices.
 *
 * @tparam T Type of the vector elements.
 * @tparam I Type of the indices.
 * @param v The input vector.
 * @param idxs A vector of indices specifying the subset.
 * @return A new vector containing the elements at the specified indices.
 */
template <typename T, typename I>
std::vector<T> get_vec_subset(std::vector<T> v, std::vector<I> idxs)
{
    std::vector<T> subset(idxs.size());
    for (I i = 0; i < idxs.size(); i++) {
        subset[i] = v[idxs[i]];
    }
    return subset;
}

/**
 * @brief Trims leading whitespace from a string.
 *
 * @param s The string to be trimmed.
 */
void ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](char ch) {
                return std::isspace(ch) == 0;
            }));
}

/**
 * @brief Trims trailing whitespace from a string.
 *
 * @param s The string to be trimmed.
 */
void rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [](char ch) { return std::isspace(ch) == 0; })
                .base(),
            s.end());
}

/**
 * @brief Trims both leading and trailing whitespace from a string.
 *
 * @param s The string to be trimmed.
 */
void trim(std::string &s)
{
    ltrim(s);
    rtrim(s);
}

/**
 * @brief Trims whitespace from all strings in a vector.
 *
 * @param v The vector of strings to be trimmed.
 */
void trim_all(StrVec_t &v)
{
    std::for_each(v.begin(), v.end(), [](std::string &s) { trim(s); });
}

/**
 * @brief Splits a string into substrings based on the specified delimiter.
 *
 * @param s The string to split.
 * @param delim delimeter to split at.
 * @return A vector of substrings.
 */
StrVec_t split(std::string &s, std::string delim)
{
    std::size_t curr_idx = 0;
    std::list<std::string> str_list;
    while (curr_idx < s.length()) {
        std::size_t next_idx = std::min(s.find(delim, curr_idx), s.length());
        std::string substr = s.substr(curr_idx, next_idx - curr_idx);
        str_list.insert(str_list.end(), substr);
        curr_idx = next_idx + 1;
    }
    StrVec_t substrs(std::begin(str_list), std::end(str_list));
    return substrs;
}

/**
 * @brief Checks if a string contains at least one alphabetic character.
 *
 * @param s The string to check.
 * @return `true` if the string contains an alphabetic character, otherwise
 * `false`.
 */
bool has_alpha(std::string &s)
{
    bool has_alphabet = false;
    for (char c : s) {
        if (std::isalpha(c) != 0) {
            has_alphabet = true;
        }
    }
    return has_alphabet;
}

/**
 * @brief Checks if all elements in a vector are unique.
 * 
 * @tparam T Type of the vector elements.
 * @param v The vector to check.
 * @return `true` if all elements are unique, otherwise `false`.
 */
template <typename T> bool is_unique(std::vector<T> v)
{
    std::sort(v.begin(), v.end());
    return std::adjacent_find(v.begin(), v.end()) == v.end();
}

/**
 * @brief Checks if a string represents a valid integer.
 * 
 * @param s The string to check.
 * @return `true` if the string represents an integer, otherwise `false`.
 */
bool str_is_int(std::string &s)
{
    std::string::const_iterator i = s.begin();
    if (s.empty()) {
        return false;
    } else if (s[0] == '-' || s[0] == '+') {
        i++;
        if (i == s.end()) {
            return false;
        }
    }
    for (; i != s.end(); i++) {
        if (!std::isdigit(*i)) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Checks if a string represents a valid floating-point number.
 * @param s The string to check.
 * @return `true` if the string represents a float, otherwise `false`.
 */
bool str_is_float(std::string &s)
{
    std::string::const_iterator i = s.begin();
    if (s.empty()) {
        return false;
    } else if (s[0] == '-' || s[0] == '+') {
        i++;
        if (i == s.end()) {
            return false;
        }
    }
    std::size_t pt_count = 0;
    for (; i != s.end(); i++) {
        if (*i == '.') {
            pt_count++;
        } else if (!std::isdigit(*i)) {
            return false;
        }
    }
    if (pt_count > 1) {
        return false;
    }
    return true;
}
