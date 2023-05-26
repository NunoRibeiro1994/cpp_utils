#ifndef _STORAGE_NUMERIC_HPP_
#define _STORAGE_NUMERIC_HPP_
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <random>
#include <algorithm>
#include <utility>
#include "storage/matrix.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// templated header-only functions for numeric operations
// provides std::matrix bulk and element-wise arithmetic operator overloads (+-/*) and matrix presets (zero,ones,ident,random,etc);
// provides std::volume bulk and element-wise arithmetic operator overloads (+-/*) and volume presets (zero,ones,ident,random,etc);
// namespace std::numeric includes sorting algorithms, scalar utility functions (norm,deg2rad,rad2deg,floor,clamp,round,etc)

// numerical operators between data types are assumed to be defined (or implicit conversion between types);
// no static assertion is performed, a compilation error will be thrown otherwise anyway;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace std {
namespace storage {
template <typename lCT>
using container_subset = st_subset_base<lCT>;
}  // namespace storage
}  // namespace std
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// subset compound operator overloads
// take rvalue references as arguments in order to 1) be able to modify left operand (vs const&) and 2) allow in-line statements (mat.row()+=5)
template <typename lCT, typename rT>
inline std::storage::st_subset_base<lCT>&& operator*= (std::storage::st_subset_base<lCT>&& _lhs, const rT& _rhs) {
    for (auto& v : _lhs) v*=_rhs;
    return std::move(_lhs);
}
template <typename lCT, typename rT>
inline std::storage::st_subset_base<lCT>&& operator/= (std::storage::st_subset_base<lCT>&& _lhs, const rT& _rhs) {
    for (auto& v : _lhs) v/=_rhs;
    return std::move(_lhs);
}
template <typename lCT, typename rT>
inline std::storage::st_subset_base<lCT>&& operator+= (std::storage::st_subset_base<lCT>&& _lhs, const rT& _rhs) {
    for (auto& v : _lhs) v+=_rhs;
    return std::move(_lhs);
}
template <typename lCT, typename rT>
inline std::storage::st_subset_base<lCT>&& operator-= (std::storage::st_subset_base<lCT>&& _lhs, const rT& _rhs) {
    for (auto& v : _lhs) v-=_rhs;
    return std::move(_lhs);
}
template <typename lCT, typename rT>
inline std::storage::st_subset_base<lCT>&& operator*= (std::storage::st_subset_base<lCT>&& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v*=_rhs[rindex++];
    return std::move(_lhs);
}
template <typename lCT, typename rT>
inline std::storage::st_subset_base<lCT>&& operator/= (std::storage::st_subset_base<lCT>&& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v/=_rhs[rindex++];
    return std::move(_lhs);
}
template <typename lCT, typename rT>
inline std::storage::st_subset_base<lCT>&& operator+= (std::storage::st_subset_base<lCT>&& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v+=_rhs[rindex++];
    return std::move(_lhs);
}
template <typename lCT, typename rT>
inline std::storage::st_subset_base<lCT>&& operator-= (std::storage::st_subset_base<lCT>&& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v-=_rhs[rindex++];
    return std::move(_lhs);
}
template <typename lCT, typename rCT>
inline std::storage::st_subset_base<lCT>&& operator*= (std::storage::st_subset_base<lCT>&& _lhs, std::storage::st_subset_base<rCT>&& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v*=_rhs[rindex++];
    return std::move(_lhs);
}
template <typename lCT, typename rCT>
inline std::storage::st_subset_base<lCT>&& operator/= (std::storage::st_subset_base<lCT>&& _lhs, std::storage::st_subset_base<rCT>&& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v/=_rhs[rindex++];
    return std::move(_lhs);
}
template <typename lCT, typename rCT>
inline std::storage::st_subset_base<lCT>&& operator+= (std::storage::st_subset_base<lCT>&& _lhs, std::storage::st_subset_base<rCT>&& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v+=_rhs[rindex++];
    return std::move(_lhs);
}
template <typename lCT, typename rCT>
inline std::storage::st_subset_base<lCT>&& operator-= (std::storage::st_subset_base<lCT>&& _lhs, std::storage::st_subset_base<rCT>&& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v-=_rhs[rindex++];
    return std::move(_lhs);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// subset non-compound operator overloads
// return std::vector with new values because 1) std_subset_base doesn't hold any data and 2) should not modify input subset container;
// std::vector as return types allows 1) direct subset assignment ("mat.row() = mat.col()+5;"), and 2) matrix/volume initialization with shape ("matrix mat(3,3,mat.all()+7;");
// can't return rvalue reference because vector is instantiated within the function and is goes out of scope
template <typename lCT, typename rT>
inline std::vector<typename lCT::value_type> operator* (const std::storage::st_subset_base<lCT>& _lhs, const rT& _rhs) {
    std::vector<typename lCT::value_type> data(0);
    for (auto& v : _lhs) data.push_back(v*_rhs);
    return data;
}
template <typename lCT, typename rT>
inline std::vector<typename lCT::value_type> operator/ (const std::storage::st_subset_base<lCT>& _lhs, const rT& _rhs) {
    std::vector<typename lCT::value_type> data(0);
    for (auto& v : _lhs) data.push_back(v/_rhs);
    return data;
}
template <typename lCT, typename rT>
inline std::vector<typename lCT::value_type> operator+ (const std::storage::st_subset_base<lCT>& _lhs, const rT& _rhs) {
    std::vector<typename lCT::value_type> data(0);
    for (auto& v : _lhs) data.push_back(v+_rhs);
    return data;
}
template <typename lCT, typename rT>
inline std::vector<typename lCT::value_type> operator- (const std::storage::st_subset_base<lCT>& _lhs, const rT& _rhs) {
    std::vector<typename lCT::value_type> data(0);
    for (auto& v : _lhs) data.push_back(v-_rhs);
    return data;
}
template <typename lCT, typename rT>
inline std::vector<typename lCT::value_type> operator* (const std::storage::st_subset_base<lCT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<typename lCT::value_type> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v*_rhs[rindex++]);
    return data;
}
template <typename lCT, typename rT>
inline std::vector<typename lCT::value_type> operator/ (const std::storage::st_subset_base<lCT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<typename lCT::value_type> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v/_rhs[rindex++]);
    return data;
}
template <typename lCT, typename rT>
inline std::vector<typename lCT::value_type> operator+ (const std::storage::st_subset_base<lCT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<typename lCT::value_type> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v+_rhs[rindex++]);
    return data;
}
template <typename lCT, typename rT>
inline std::vector<typename lCT::value_type> operator- (const std::storage::st_subset_base<lCT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<typename lCT::value_type> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v-_rhs[rindex++]);
    return data;
}
template <typename lCT, typename rCT>
inline std::vector<typename lCT::value_type> operator* (const std::storage::st_subset_base<lCT>& _lhs, const std::storage::st_subset_base<rCT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<typename lCT::value_type> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v*_rhs[rindex++]);
    return data;
}
template <typename lCT, typename rCT>
inline std::vector<typename lCT::value_type> operator/ (const std::storage::st_subset_base<lCT>& _lhs, const std::storage::st_subset_base<rCT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<typename lCT::value_type> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v/_rhs[rindex++]);
    return data;
}
template <typename lCT, typename rCT>
inline std::vector<typename lCT::value_type> operator+ (const std::storage::st_subset_base<lCT>& _lhs, const std::storage::st_subset_base<rCT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<typename lCT::value_type> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v+_rhs[rindex++]);
    return data;
}
template <typename lCT, typename rCT>
inline std::vector<typename lCT::value_type> operator- (const std::storage::st_subset_base<lCT>& _lhs, const std::storage::st_subset_base<rCT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<typename lCT::value_type> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v-_rhs[rindex++]);
    return data;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// std::vector compound operators (to allow chained inline operations)
// std::vector can
template <typename lT, typename rT>
inline std::vector<lT>& operator*= (std::vector<lT>& _lhs, const rT& _rhs) {
    for (auto& v : _lhs) v *= _rhs;
    return _lhs;
}
template <typename lT, typename rT>
inline std::vector<lT> operator/= (std::vector<lT>& _lhs, const rT& _rhs) {
    for (auto& v : _lhs) v /= _rhs;
    return _lhs;
}
template <typename lT, typename rT>
inline std::vector<lT>& operator+= (std::vector<lT>& _lhs, const rT& _rhs) {
    for (auto& v : _lhs) v += _rhs;
    return _lhs;
}
template <typename lT, typename rT>
inline std::vector<lT>& operator-= (std::vector<lT>& _lhs, const rT& _rhs) {
    for (auto& v : _lhs) v -= _rhs;
    return _lhs;
}
template <typename lT, typename rT>
inline std::vector<lT>& operator*= (std::vector<lT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v *= _rhs[rindex++];
    return _lhs;
}
template <typename lT, typename rT>
inline std::vector<lT>& operator/= (std::vector<lT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v /= _rhs[rindex++];
    return _lhs;
}
template <typename lT, typename rT>
inline std::vector<lT>& operator+= (std::vector<lT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v += _rhs[rindex++];
    return _lhs;
}
template <typename lT, typename rT>
inline std::vector<lT>& operator-= (std::vector<lT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v -= _rhs[rindex++];
    return _lhs;
}
template <typename lT, typename rCT>
inline std::vector<lT>& operator*= (std::vector<lT>& _lhs, std::storage::st_subset_base<rCT>&& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v *= _rhs[rindex++];
    return _lhs;
}
template <typename lT, typename rCT>
inline std::vector<lT>& operator/= (std::vector<lT>& _lhs, std::storage::st_subset_base<rCT>&& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v /= _rhs[rindex++];
    return _lhs;
}
template <typename lT, typename rCT>
inline std::vector<lT>& operator+= (std::vector<lT>& _lhs, std::storage::st_subset_base<rCT>&& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v += _rhs[rindex++];
    return _lhs;
}
template <typename lT, typename rCT>
inline std::vector<lT>& operator-= (std::vector<lT>& _lhs, std::storage::st_subset_base<rCT>&& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    int rindex = 0;
    for (auto& v : _lhs) v -= _rhs[rindex++];
    return _lhs;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// std::vector non-compound operators
template <typename lT, typename rT>
inline std::vector<lT> operator* (const std::vector<lT>& _lhs, const rT& _rhs) {
    std::vector<lT> data(0);
    for (auto& v : _lhs) data.push_back(v*_rhs);
    return data;
}
template <typename lT, typename rT>
inline std::vector<lT> operator/ (const std::vector<lT>& _lhs, const rT& _rhs) {
    std::vector<lT> data(0);
    for (auto& v : _lhs) data.push_back(v/_rhs);
    return data;
}
template <typename lT, typename rT>
inline std::vector<lT> operator+ (const std::vector<lT>& _lhs, const rT& _rhs) {
    std::vector<lT> data(0);
    for (auto& v : _lhs) data.push_back(v+_rhs);
    return data;
}
template <typename lT, typename rT>
inline std::vector<lT> operator- (const std::vector<lT>& _lhs, const rT& _rhs) {
    std::vector<lT> data(0);
    for (auto& v : _lhs) data.push_back(v-_rhs);
    return data;
}
template <typename lT, typename rT>
inline std::vector<lT> operator* (const std::vector<lT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<lT> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v*_rhs[rindex++]);
    return data;
}
template <typename lT, typename rT>
inline std::vector<lT> operator/ (const std::vector<lT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<lT> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v/_rhs[rindex++]);
    return data;
}
template <typename lT, typename rT>
inline std::vector<lT> operator+ (const std::vector<lT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<lT> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v+_rhs[rindex++]);
    return data;
}
template <typename lT, typename rT>
inline std::vector<lT> operator- (const std::vector<lT>& _lhs, const std::vector<rT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<lT> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v-_rhs[rindex++]);
    return data;
}
template <typename lT, typename rCT>
inline std::vector<lT> operator* (const std::vector<lT>& _lhs, const std::storage::st_subset_base<rCT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<lT> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v*_rhs[rindex++]);
    return data;
}
template <typename lT, typename rCT>
inline std::vector<lT> operator/ (const std::vector<lT>& _lhs, const std::storage::st_subset_base<rCT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<lT> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v/_rhs[rindex++]);
    return data;
}
template <typename lT, typename rCT>
inline std::vector<lT> operator+ (const std::vector<lT>& _lhs, const std::storage::st_subset_base<rCT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<lT> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v+_rhs[rindex++]);
    return data;
}
template <typename lT, typename rCT>
inline std::vector<lT> operator- (const std::vector<lT>& _lhs, const std::storage::st_subset_base<rCT>& _rhs) {
    assert(_lhs.size() <= _rhs.size());
    std::vector<lT> data(0);
    int rindex = 0;
    for (auto& v : _lhs) data.push_back(v-_rhs[rindex++]);
    return data;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VECTOR NUMERIC PRESETS
namespace std {
namespace vector_presets {
template <typename _type>
vector<_type> linspace(int _s, _type _low, _type _high) {
    static_assert(std::is_arithmetic<_type>::value, "PRESET ONLY AVAILABLE TO NUMERIC TYPES");
    vector<_type> lin(_s);
    _type step = ((_high-_low)/(_s-1));
    for (int i=0; i<_s; i++) lin[i] = static_cast<_type>(step*i);
    return lin;
}
}  // namespace vector_presets
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// NUMERIC TYPEDEFS
typedef matrix<double> matrixDouble;
typedef matrix<int>    matrixInt;
typedef matrix<size_t> matrixUInt;
typedef matrix<float>  matrixFloat;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MATRIX NUMERIC PRESETS
namespace matrix_presets {
template <typename _type>
std::matrix<_type> zeros(int _r, int _c) {
    static_assert(std::is_arithmetic<_type>::value, "PRESET ONLY AVAILABLE TO NUMERIC TYPES");
    return std::matrix<_type>(_r, _c, 0);
}
template <typename _type>
std::matrix<_type> ones(int _r, int _c) {
    static_assert(std::is_arithmetic<_type>::value, "PRESET ONLY AVAILABLE TO NUMERIC TYPES");
    return std::matrix<_type>(_r, _c, 1);
}
template <typename _type>
std::matrix<_type> square(int _s) {
    static_assert(std::is_arithmetic<_type>::value, "PRESET ONLY AVAILABLE TO NUMERIC TYPES");
    return std::matrix<_type> (_s, _s, 0.0);
}
template <typename _type>
std::matrix<_type> list(int _s) {
    static_assert(std::is_arithmetic<_type>::value, "PRESET ONLY AVAILABLE TO NUMERIC TYPES");
    return std::matrix<_type> (_s, 1, static_cast<_type>(0));
}
template <typename _type>
std::matrix<_type> identity(int _s) {
    static_assert(std::is_arithmetic<_type>::value, "PRESET ONLY AVAILABLE TO NUMERIC TYPES");
    std::matrix<_type> ident(_s, _s, 0.0);
    ident[ident.diag()] = std::vector<_type> (_s, static_cast<_type>(1));
    return ident;
}
template <typename _type>
std::matrix<_type> random(int _r, int _c) {
    static_assert(std::is_arithmetic<_type>::value, "PRESET ONLY AVAILABLE TO NUMERIC TYPES");
    std::vector<_type> vec(_r*_c);
    std::mt19937 rng(std::time(nullptr));
    std::normal_distribution<float> normal(0, 50);
    for (auto& val : vec) val = normal(rng);
    return std::matrix<_type> (vec, _r, _c);
}
template <typename _type>
std::matrix<_type> linspace(int _s, _type _low, _type _high) {
    static_assert(std::is_arithmetic<_type>::value, "PRESET ONLY AVAILABLE TO NUMERIC TYPES");
    std::matrix<_type> lin(1, _s);
    _type step = ((_high-_low)/(_s-1));
    for (int i=0; i < _s; i++) lin(0, i) = step*i;
    return lin;
}
}  // namespace matrix_presets
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SORTING METHODS (for std::vector, std::matrix and std::volume)
template <typename _type>
std::vector<size_t> sort_indexes(const std::vector<_type> &vec) {
    // initialize original index locations
    std::vector<size_t> idx(vec.size());
    std::iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in vec
    std::sort(idx.begin(), idx.end(), [&vec](size_t i1, size_t i2) { return vec[i1] < vec[i2]; });
    return idx;
}
// sorts within the whole matrix & volume, returns matrix with order
// for segment sorting, call each segment separately e.g. sort_indexes(mat[mat.row(47)]);
template <typename _type>
std::matrix<size_t> sort_indexes(const std::matrix<_type> &mat) {
    return std::matrix<_type>(sort_indexes(mat.vector()), mat.rows(), mat.cols());
}
// template <typename _type>
// std::matrix<size_t> sort_indexes(const std::volume<_type> &vol) {
//  return std::volume<_type>(sort_indexes(vol.vector()),vol.levels(),vol.rows(),vol.cols());
// }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SCALAR NON-MATRIX NUMERIC STUFF
inline float norm(float x, float y) {
    return sqrt(pow(x, 2) + pow(y, 2));
}
inline float norm(float x, float y, float z) {
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}
inline float deg2rad(float deg) {
    return deg * M_PI / 180.0;
}
inline float rad2deg(float rad) {
    return rad * 180.0 / M_PI;
}
inline float sign(float value) {
    return value < 0 ? -1.0 : 1.0;
}
inline float clamp(float value, float min, float max) {
  if (min > max) {
    assert(0);
    return value;
  }
  return value < min ? min : value > max ? max : value;
}
inline float round(float val, int dig) {
    float vv, mult, rnd;
    mult = 1.0;
    for (int i = 0; i < dig; ++i) mult=mult*10;
    vv = mult*val;
    //
    if (vv < 0) {
        rnd = ceil(vv - 0.5);
    } else {
        rnd = floor(vv + 0.5);
    }
    //
    return (rnd/mult);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace std
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif  // _STORAGE_NUMERIC_HPP_
