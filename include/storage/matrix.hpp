#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <utility>
#include <ctime>
#include <cassert>
#include <iostream>
#include <sstream>  // for std::matrix<T>::load()
#include <fstream>  // for std::matrix<T>::load()
#include <algorithm>
#include "storage/subset.hpp"
#include "storage/type_check.hpp"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace std {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO(joao): replace assert with "if () throw"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper classes
template < typename T > class matrix;  // forward declaration required for 'using'
template < typename dT >
using matrix_subset = storage::st_subset_base< matrix< dT > >;
template < typename dT >
using matrix_subset_const = storage::st_subset_base< const matrix< dT > >;

//------------------------------------------------------------------------------
/// @brief      Class implementaing a STL-like 2D container
///
/// @tparam     T     Underlying _data type.
///
template < typename T >
class matrix {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Regarding std::matrix< T > member structure:
    // A) Constructors and Assignment Operators:
    //   C++ "rule of three": classes should define 1) destructor; 2) copy constructor; and 3) copy assigment operator "=" overload;
    //   C++ "rule of five": classes should define 1) destructor; 2) copy constructor; 3) move constructor; 4) copy assigment operator "=" overload; 5) move assignment operator "=" overload;
    //   In this particular class, since no dynamic memory is being directly managed, there is no need for any of these five methods
    // B) Const methods:
    //   For const correctness, methods that do not modify the contents are "const" keyword, making them available when passing by const reference;
    // C) Move semantics and lvalue/rvalue arguments
    //   It used to be generally recommended best practice to use pass by const ref for all types, except for
    //   1) builtin types (char, int, double, etc.), 2) iterators and 3) function objects (lambdas, classes deriving from std::*_function); [cf. https://stackoverflow.com/a/270435/2340002]
    // D) Subcript indices return type
    //   A specific storage::subset class was defined that holds an index vector on std::matrix< T > _data vector
    //   This allows both read and write acess to multiple elements on the matrix _data (e.g. row(), column(), etc)
    //   Access (i.e. std::vector row = mat[mat.row(47)]) returns a copy of the original _data, because vectors always own their content (can't have multiple vectors mapping the same memory)
    //   [cf https://stackoverflow.com/a/52939707/2340002]
    // E) Error handling:
    //   assert() is used to ensure code behaves as expected, in alternative to exception handling, which should be used for run-time or i/o error
    //   [cf https://stackoverflow.com/a/12887257/2340002]
    // F) Type conversion
    //   Templated members allow different argument and return types when interacting with std::vector and specialized containers, assuming types are casteable;
    //   A cast operator is also provided to allow casting between matrix types;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 protected:
    //--------------------------------------------------------------------------
    /// @brief      Number of matrix rows.
    ///
    size_t _rows;

    //--------------------------------------------------------------------------
    /// @brief      Number of matrix columns.
    ///
    size_t _cols;

    //--------------------------------------------------------------------------
    /// @brief      Underlying _data container.
    ///
    vector< T > _data;

 public:
    //--------------------------------------------------------------------------
    /// @brief      Element type.
    ///
    typedef T value_type;

    //--------------------------------------------------------------------------
    /// @brief      Static 'factory' named constructor that initializes a matrix from a text file
    ///
    /// @param[in]  path  Path to source file.
    /// @param[in]  sep   Separator character to use. Defaults to ',' (.csv format default).
    /// @param[in]  skip  Number of lines to skip @ beginning for file. Useful if file format includes multi-line header information. Defaults to 1 (assumes first line w/ column labels)
    ///
    /// @return     Newly instantiated std::matrix< T > object wherein each row is populated with contents of each line in given file @ *path*
    ///
    /// @note       Temporary/develo0pment implementation. will be replaced by methods within std::io namespace (cf io.hpp)
    ///
    static matrix< T > load(const string& path, char sep = ',', size_t skip = 1);

    // //--------------------------------------------------------------------------
    // /// @brief      Constructs a new instance.
    // ///
    // matrix();

    // //--------------------------------------------------------------------------
    // /// @brief      Constructs a new instance.
    // ///
    // /// @param[in]  rows  Number of rows
    // /// @param[in]  cols  Number of cols
    // ///
    // matrix(size_t rows, size_t cols);

    //--------------------------------------------------------------------------
    /// @brief      Constructs a new instance.
    ///
    /// @param[in]  rows  Number of rows
    /// @param[in]  cols  Number of cols
    /// @param[in]  args  Arguments (parameter pack) to be forward to T constructor.
    ///
    template < typename... Args, typename = typename enable_if< is_constructible< T, Args... >::value >::type >
    matrix(size_t rows = 0, size_t cols = 0, Args&&... args);
    // matrix(size_t _rows, size_t _cols, const T& _val);


    //explicit matrix(const std::vector< T >& in):

    //--------------------------------------------------------------------------
    /// @brief      Constructs a new instance, from a nested vector instance.
    ///
    /// @param[in]  in   Input vector of vectors
    ///
    /// @tparam     iT   Input _data type, which must be convertible to T.
    ///
    // template < typename iT >
    // explicit matrix(const vector< vector< iT > >& in);






    //--------------------------------------------------------------------------
    /// @brief      Constructs a new instance, from a vector of generic container type (e.g. nested vectors).
    ///
    template < typename iT = vector< T >, typename = typename enable_if< is_generic_container< iT >() >::type >
    explicit matrix(const vector< iT >& in);

    //--------------------------------------------------------------------------
    /// @brief      Constructs a new instance, from a vector of generic container type.
    ///
    /// @todo       consider marking as implict for simpler syntax (return expresions and assignment operations)
    /// 
    template < typename iT, typename = typename enable_if< is_nd_container< iT >()>::type >
    explicit matrix(const storage::st_subset_base< iT >& in);    // from multi-dimensional container (e.g. volume, volume subset, vector, etc) with shape info (.shape() and .position())

    template < typename iT, typename = typename enable_if< !is_nd_container< iT >()>::type, typename = void >
    explicit matrix(const storage::st_subset_base< iT >& in);    // from non-nd containers, constructs unidimensional matrix;

    //--------------------------------------------------------------------------
    /// @brief      Constructs a new [*rows* x *cols*] matrix, sourcing elements from a generic container.
    ///
    /// @param[in]  rows       Number of rows
    /// @param[in]  cols       Number of cols
    /// @param[in]  in         Container/data source
    ///
    /// @tparam     iT         Type of input container
    ///
    template < typename iT, typename = typename enable_if< is_generic_container< iT >()>::type >
    matrix(size_t rows, size_t cols, const iT& in);          // from generic container (e.g. volume, volume subset, vector, etc), reshapes!

    //--------------------------------------------------------------------------
    /// @brief      Destroys the object.
    ///
    ~matrix() = default;

    //--------------------------------------------------------------------------
    /// @brief      Dereferencing/indirection operator.
    ///
    /// @return     Returns first element in matrix i.e. @ (0,0)
    ///
    T& operator*();

    //--------------------------------------------------------------------------
    /// @brief      Assignment operator.
    ///
    /// @param[in]  in    Input value to assign to all elements in matrix.
    ///
    /// @tparam     iT    Type of the input argument.
    ///
    template < typename iT >
    matrix< T >& operator=(const iT& in);

    //--------------------------------------------------------------------------
    /// @brief      Random access operator, for element @ given (*row*, *col*) position.
    ///
    /// @param[in]  row   Row index.
    /// @param[in]  col   Column index.
    ///
    /// @return     Reference to element @ (*row*, *col*)
    ///
    T& operator()(size_t row, size_t col);

    //--------------------------------------------------------------------------
    /// @brief      Random access operator, for element @ given (*row*, *col*) position (*const overload*).
    ///
    /// @param[in]  row   Row index.
    /// @param[in]  col   Column index.
    ///
    /// @return     Const reference to element @ (*row*, *col*)
    ///
    const T& operator()(size_t row, size_t col) const;

    //--------------------------------------------------------------------------
    /// @brief      Random access operator, for element @ given *idx* position.
    ///
    /// @param[in]  idx   Element index.
    ///
    /// @return     Reference to element @ *idx* position.
    ///
    /// @todo       Return line subset instead of random accessed element in order to allow [row][col] access.
    ///             Add at(size_t) member for random access.
    ///
    T& operator[](size_t idx);

    //--------------------------------------------------------------------------
    /// @brief      Random access operator, for element @ given (*row*, *col*) position (*const overload*).
    ///
    /// @param[in]  idx   Element index.
    ///
    /// @return     Const reference to element @ *idx* position.
    ///
    /// @todo       Return line subset instead of random accessed element in order to allow [row][col] access.
    ///             Add at(size_t) member for random access.
    ///
    const T& operator[](size_t idx) const;

    //--------------------------------------------------------------------------
    /// @brief      Random access operator, for multiple elements @ given *idx*.
    ///
    /// @param[in]  idx  Vector of element indexes of output matrix.
    ///
    /// @return     Matrix subset instance (as matrix_subset< T >) referencing all elements @ *idx* valid values.
    ///
    matrix_subset< T > operator[](const vector< size_t >& idx);

    //--------------------------------------------------------------------------
    /// @brief      Random access operator, for multiple elements @ given *idx* (*const overload*).
    ///
    /// @param[in]  idx  Vector of element indexes of output matrix.
    ///
    /// @return     Const matrix subset instance (as matrix_subset< T >) referencing all elements @ *idx* valid values.
    ///
    matrix_subset_const< T > operator[](const vector< size_t >& idx) const;

    //--------------------------------------------------------------------------
    /// @brief      Random access operator, for multiple elements described by given matrix subset.
    ///
    /// @param[in]  idx  Vector of element indexes of output matrix.
    ///
    /// @return     Matrix subset instance (as matrix_subset< T >) referencing all elements @ *idx* valid values.
    ///
    matrix_subset< T > operator[](const matrix_subset< T >& sbst);

    //--------------------------------------------------------------------------
    /// @brief      Random access operator, for multiple elements described by given matrix subset (*const overload*).
    ///
    /// @param[in]  idx  Vector of element indexes of output matrix.
    ///
    /// @return     Const matrix subset instance (as matrix_subset< T >) referencing all elements @ *idx* valid values.
    ///
    matrix_subset_const< T > operator[](const matrix_subset< T >& sbst) const;

    //--------------------------------------------------------------------------
    /// @brief      Access underlying STL container (read-only).
    ///
    /// @return     Const reference to underlying data container.
    ///
    const vector< T >& elements() const;

    //--------------------------------------------------------------------------
    /// @brief      Access raw data.
    ///
    /// @return     Const pointer to first element.
    ///
    const T* data() const;

    //--------------------------------------------------------------------------
    /// @brief      Get number of matrix rows.
    ///
    size_t rows() const;

    //--------------------------------------------------------------------------
    /// @brief      Get number of matrix columns
    ///
    size_t cols() const;

    //--------------------------------------------------------------------------
    /// @brief      Get total number of elements in the underlying container.
    ///
    size_t size() const;

    //--------------------------------------------------------------------------
    /// @brief      Gets matrix shape.
    ///
    /// @return     Vector of size 2 matching { rows, cols }.
    ///
    vector< size_t > shape() const;

    //--------------------------------------------------------------------------
    /// @brief      Gets matrix position i.e. (*row*, *column*) of element @ position *idx*.
    ///
    /// @param[in]  idx    Element index.
    ///
    /// @return     Vector of size 2 matching { row, col }.
    ///
    vector< size_t > position(size_t idx) const;

    //--------------------------------------------------------------------------
    /// @brief      Determines if square matrix.
    ///
    /// @return     True if number of rows == number of columns, False otherwise.
    ///
    bool isSquare() const;

    //--------------------------------------------------------------------------
    /// @brief      Determines if empty matrix.
    ///
    /// @return     True if both number of rows and number of columns > 0, False otherwise.
    ///
    bool isEmpty()  const;

    //--------------------------------------------------------------------------
    /// @brief      Reserves capacity for matrix expansion.
    ///
    /// @param[in]  rows  Maximum number of rows.
    /// @param[in]  cols  Maximum number of columns.
    ///
    void reserve(size_t rows, size_t cols);

    //--------------------------------------------------------------------------
    /// @brief      Get indexes of all matrix elements.
    ///
    /// @note       Required in order to construct matrix_subset object referencing all matrix elements.
    ///             matrix_subset<> is not constructible from matrix<>, and can't declare cast operator to vector< size_t >
    ///             as it causes ambiguity when T = size_t.
    ///
    /// @note       Review & debug how allID() is implemented.
    ///
    /// @todo       Declare protected/private - no need to be part of the public interface.
    ///
    vector< size_t > allID() const;

    //--------------------------------------------------------------------------
    /// @brief      Get indexes of all matrix elements in given *row*.
    ///
    /// @param[in]  row    Row index.
    ///
    /// @todo       Declare protected/private - no need to be part of the public interface.
    ///
    vector< size_t > rowID(size_t row) const;

    //--------------------------------------------------------------------------
    /// @brief      Get indexes of all matrix elements in given *col*.
    ///
    /// @param[in]  col    Col index.
    ///
    /// @todo       Declare protected/private - no need to be part of the public interface.
    ///
    vector< size_t > colID(size_t col) const;

    //--------------------------------------------------------------------------
    /// @brief      Get indexes of all matrix elements in matrix diagonal (upper diagonal if matrix is not square).
    ///
    /// @todo       Declare protected/private - no need to be part of the public interface.
    ///
    vector< size_t > diagID() const;

    //--------------------------------------------------------------------------
    /// @brief      Get indexes of all matrix elements in the block of rows from *first_row* to *last_row*
    ///             and columns from *first_col* to *last_col*.
    ///
    /// @todo       Declare protected/private - no need to be part of the public interface.
    ///
    vector< size_t > blockID(size_t first_row, size_t last_row, size_t first_col = 0, size_t last_col = 0) const;

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all matrix elements.
    ///
    matrix_subset< T > all();

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all matrix elements (*const overload*).
    ///
    matrix_subset_const< T > all() const;

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all elements in given *row*.
    ///
    matrix_subset< T > row(size_t row);

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all elements in given *row* (*const overload*).
    ///
    matrix_subset_const< T > row(size_t row) const;

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all elements in given *col*.
    ///
    matrix_subset< T > col(size_t col);

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all elements in given *col* (*const overload*).
    ///
    matrix_subset_const< T > col(size_t col) const;

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all elements in matrix diagonal (upper diagonal if matrix is not square).
    ///
    matrix_subset< T > diag();

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all elements in matrix diagonal (upper diagonal if matrix is not square) (*const overload*).
    ///
    matrix_subset_const< T > diag() const;

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all elements in the block of rows from *first_row* to *last_row*
    ///             and columns from *first_col* to *last_col*.
    ///
    matrix_subset< T > block(size_t first_row, size_t last_row, size_t first_col = 0, size_t last_col = 0);

    //--------------------------------------------------------------------------
    /// @brief      Get subset including all elements in the block of rows from *first_row* to *last_row*
    ///             and columns from *first_col* to *last_col* (*const overload*).
    ///
    matrix_subset_const< T > block(size_t first_row, size_t last_row, size_t first_col = 0, size_t last_col = 0) const;

    //--------------------------------------------------------------------------
    /// @brief      Create new matrix including all elements in the block of rows from *first_row* to *last_row*
    ///             and columns from *first_col* to *last_col* (*const overload*).
    ///
    /// @note       Wraps around block() and constructor matrix(const matrix_subset< >&) for a simpler syntax i.e
    ///             return this->submatrix<>(1, 5, 7, 10) instead of 'std::matrix<>(this->block(1, 5, 7, 10))'       
    ///
    template < typename oT = value_type >
    matrix< oT > submat(size_t first_row, size_t last_row, size_t first_col = 0, size_t last_col = 0) const;

    //--------------------------------------------------------------------------
    /// @brief      Assigns new data to matrix content.
    ///
    /// @param[in]  in    Input object with values to assign,
    ///
    /// @tparam     iT    Type of input container.
    ///
    /// @note       Matrix is not reshapen, input container must have *at least* size() values
    ///
    template < typename iT = vector< T >, typename = typename enable_if< is_generic_container< iT >() >::type >
    void set(const iT& in);

    //--------------------------------------------------------------------------
    /// @brief      Assigns new data to matrix content.
    ///
    /// @param[in]  in    Input object with values to assign.
    ///
    /// @tparam     iT    Type of input container.
    ///
    /// @note       Matrix is not reshapen, input container must have *at least* size() values
    ///
    /// @note       Fallback overload in case input type does not behave like a container type but is implicit convertible to a vector< value_type >.
    ///
    template < typename iT = vector< T >, typename = typename enable_if< !is_generic_container< iT >() >::type,
                                          typename = typename enable_if< is_convertible< iT, vector< value_type > >::value >::type >
    void set(const iT& in);

    //--------------------------------------------------------------------------
    /// @brief      Assigns new data to matrix content.
    ///
    /// @param[in]  in    Input object with values to assign.
    ///
    /// @tparam     iT    Type of input container.
    ///
    /// @note       Matrix is not reshapen, input container must have *at least* size() values
    ///
    /// @note       Fallback overload in case input type does not behave like a container type but is implicit convertible to a vector< value_type >.
    ///
    template < typename iT = vector< T >, typename = typename enable_if< is_convertible< iT, value_type >::value >::type >
    void fill(const iT& in);

    //--------------------------------------------------------------------------
    /// @brief      Adds a new row to the matrix with values from given *in*.
    ///
    /// @param[in]  in    Input object with values to assign to new row.
    ///
    /// @tparam     iT    Type of input container.
    ///
    template < typename iT = vector< T >, typename = typename enable_if< is_generic_container< iT >() >::type >
    void pushRow(const iT& in);

    //--------------------------------------------------------------------------
    /// @brief      Adds a new row to the matrix with values from given *in*.
    ///
    /// @param[in]  in    Input object with values to assign to new row.
    ///
    /// @note       Fallback overload in case input type does not behave like a container type but is implicit convertible to a vector< value_type >.
    ///
    template < typename iT = vector< T >, typename = typename enable_if< !is_generic_container< iT >() >::type,
                                          typename = typename enable_if< is_convertible< iT, vector< value_type > >::value >::type >
    void pushRow(const iT& in);

    //--------------------------------------------------------------------------
    /// @brief      Adds a new column to the matrix with values from given *in*.
    ///
    /// @param[in]  in    Input object with values to assign to new column.
    ///
    /// @tparam     iT    Type of input container.
    ///
    template < typename iT = vector< T >, typename = typename std::enable_if< std::is_generic_container< iT >() >::type >
    void pushCol(const iT& in);

    //--------------------------------------------------------------------------
    /// @brief      Adds a new column to the matrix with values from given *in*.
    ///
    /// @param[in]  in    Input object with values to assign to new column.
    ///
    /// @note       Fallback overload in case input type does not behave like a container type but is implicit convertible to a vector< value_type >.
    ///
    template < typename iT = vector< T >, typename = typename enable_if< !is_generic_container< iT >() >::type,
                                          typename = typename enable_if< is_convertible< iT, vector< value_type > >::value >::type >
    void pushCol(const iT& in);

    //--------------------------------------------------------------------------
    /// @brief      Adds a new empty row to the matrix.
    ///
    /// @todo       Remove in favour of default empty container for templated overload.
    ///
    void pushRow();

    //--------------------------------------------------------------------------
    /// @brief      Adds a new empty column to the matrix.
    ///
    /// @todo       Remove in favour of default empty container for templated overload.
    ///
    void pushCol();

    //--------------------------------------------------------------------------
    /// @brief      Deletes last row of matrix.
    ///
    void popRow();

    //--------------------------------------------------------------------------
    /// @brief      Deletes last column of matrix.
    ///
    void popCol();

    //--------------------------------------------------------------------------
    /// @brief      Deletes given *row*.
    ///
    void deleteRow(size_t row);

    //--------------------------------------------------------------------------
    /// @brief      Deletes given *col*.
    ///
    void deleteCol(size_t col);

    //--------------------------------------------------------------------------
    /// @brief      Reshapes the matrix to given shape *new_rows* x *new_cols* .
    ///
    void reshape(size_t new_rows, size_t new_cols);

    //--------------------------------------------------------------------------
    /// @brief      Reshapes the matrix to given shape *new_rows* x *new_cols* .
    ///
    /// @note       Alternative to reshape(), less efficient, provided for consistency.
    ///             reshape() should be used instead.
    ///
    void resize(size_t new_rows, size_t new_cols);

    //--------------------------------------------------------------------------
    /// @brief      Flips matrix (== transpose())
    ///
    /// @todo       Rename to transpose().
    ///
    void flip();

    //--------------------------------------------------------------------------
    /// @brief      Clears the matrix (elements are destroyed)
    ///
    void clear();

    //--------------------------------------------------------------------------
    /// @brief      Get iterator to the beginning of the matrix.
    ///
    typename vector< T >::iterator begin();

    //--------------------------------------------------------------------------
    /// @brief      Get iterator to the end of the matrix.
    ///
    typename vector< T >::iterator end();

    //--------------------------------------------------------------------------
    /// @brief      Get iterator to the beginning of the matrix (*const overload*).
    ///
    typename vector< T >::const_iterator begin() const;

    //--------------------------------------------------------------------------
    /// @brief      Get iterator to the end of the matrix (*const overload*).
    ///
    typename vector< T >::const_iterator end()   const;

    //--------------------------------------------------------------------------
    /// @brief      Conversion operator to vector object (by reference)
    ///
    template < typename oT >
    operator vector< oT >&();

    //--------------------------------------------------------------------------
    /// @brief      Conversion operator to vector object (by const reference)
    ///
    template < typename oT >
    operator const vector< oT >&() const;

    //--------------------------------------------------------------------------
    /// @brief      Conversion operator to matrix object of different type (by reference)
    ///
    template < typename oT >
    operator matrix< oT >&();

    //--------------------------------------------------------------------------
    /// @brief      Conversion operator to matrix object of different type (by const reference)
    ///
    template < typename oT >
    operator const matrix< oT >&() const;
};

/// @cond

template < typename T >
matrix< T > matrix< T >::load(const string& path, char sep, size_t skip) {
    string line;
    string word;
    stringstream ss;
    // stringstream val;
    ifstream file(path.data());
    if (!file.good()) {
        throw runtime_error(string(__func__) + "(): Unable to load file!");
    }
    file.precision(10);

    /// @todo add size arguments to pre-allocated matrix container (should be faster,  less resizes/reshapes during load)

    matrix< T > out;
    out.reserve(1000, 1000);
    vector< T > line_vals;
    int l = 0;
    // int v = 0;
    double value = 0.0;

    if (file.is_open()) {
        while (getline(file, line)) {
            l++;
            if (l < skip + 1) continue;

            // // std::vector< T > line_vals;
            // size_t start;
            // size_t end = 0;
            // while ((start = line.find_first_not_of(sep, end)) != std::string::npos) {
            //     end = line.find(sep, start);
            //     auto word = line.substr(start, end - start);
            //     if (word.empty()) printf("EMPTY WORD!!!\n");
            //     T val;
            //     stringstream(word) >> val;
            //     line_vals.push_back(val);
            //     // if (single_split && end != std::string::npos) {
            //     //     words.push_back(data.substr(end + 1));
            //     //     break;
            //     // }
            // }
            // std::cout << "start: " << (start == std::string::npos) << std::endl;

            ss.str(line);
            while (getline(ss, word, sep)) {
                // v++;
                if (word.empty()) {
                    value = 0.0;
                } else {
                    // val.str(word);
                    istringstream(word) >> value;
                }
                line_vals.push_back(value);
                // if (v < out.cols()) {
                //     out(out.rows() - 1, v++) = value;
                // } else {
                //     out.reshape(out.rows(), out.cols() + 1);
                // }
                // val.clear();
                // std::cout << v << std::endl;
            }

            // fill untilk expected size
            // required in cases where there is a missing value @ end of line
            while (line_vals.size() < out.cols()) {
                line_vals.push_back(0.0);
            }
            // resize matrix when required!
            // if at least a row has been 
            // if (out.rows() && line_vals.size() > out.cols()) {
            //     out.reshape(out.rows(), line_vals.size());
            // }
            // 
            // std::cout << l << ": " << line_vals.size() << std::endl;
            out.pushRow(line_vals);
            line_vals.resize(0);
            line_vals.reserve(out.cols());  // may be redundant!
            ss.clear();
            // v = 0;
        }
        file.close();
    }

    return out;
}




// template < typename T >
// matrix< T >::matrix() {
//     _data  = vector< T > (0);
//     _rows = 0;
//     _cols = 0;
//     _data.reserve(10*10);
// }


// template < typename T >
// matrix< T >::matrix(size_t _rows, size_t _cols) {
//     assert(_rows > 0 && _cols > 0);
//     _data  = vector< T > ( _rows * _cols, T());
//     _rows = _rows;
//     _cols = _cols;
//     _data.reserve(2* _rows * _cols);
// }


template < typename T >
template < typename... Args, typename  >
matrix< T >::matrix(size_t rows, size_t cols, Args&&... args) {
    // assert(rows > 0 && cols > 0);
    // _data  = vector< T > (rows * cols, _val);
    _rows = rows;
    _cols = cols;
    _data.reserve(_rows * _cols);
    for (int i = 0; i < (rows * cols); i++) {
        _data.emplace_back(std::forward< Args >(args)...);
    }
}


// template < typename T >
// template < typename iT >
// matrix< T >::matrix(const vector< vector< iT > >& in) {
//     _rows = in.size();
//     _cols = 0;
//     _data.reserve(_rows * _cols);
//     for (size_t r = 0; r < _rows; r++) {
//         if (in[r].size() > _cols) {
//             _cols = in[r].size();
//         }
//     }
//     _data = vector< T >(_rows * _cols);
//     for (size_t r = 0; r < _rows; r++) {
//         for (size_t c = 0; c < _cols; ++c) {
//             if (c < in[r].size()) {
//                 (*this)(r, c) = static_cast< T >(in[r][c]);
//             }
//         }
//     }
// }


template < typename T >
template < typename iT, typename >
matrix< T >::matrix(const vector< iT >& in) {
    _rows = in.size();
    _cols = 0;
    _data.reserve(_rows * (_cols + 1));
    for (size_t r = 0; r < _rows; r++) {
        if (in[r].size() > _cols) {
            _cols = in[r].size();
        }
    }
    _data = vector< T >(_rows * _cols);
    for (size_t r = 0; r < _rows; r++) {
        for (size_t c = 0; c < _cols; ++c) {
            if (c < in[r].size()) {
                (*this)(r, c) = static_cast< T >(in[r][c]);
            }
        }
    }
}


template < typename T >
template < typename iT, typename >
matrix< T >::matrix(const storage::st_subset_base< iT >& in) {
    static_assert(is_nd_container< iT >(), "INVALID INPUT CONTAINER!");
    vector< size_t > inDims = in.source()->shape();
    assert(inDims.size() > 0);
    size_t iR = 0;
    size_t iC = 0;
    vector< size_t > rs(in.size(), 0);
    vector< size_t > cs(in.size(), 0);
    if (inDims.size() > 0) {
        iC = inDims.size()-1;
    }
    if (inDims.size() > 1) {
        iR = inDims.size()-2;
    } else {
        iR = 1;
    }
    size_t minR = in.source_size();
    size_t minC = in.source_size();
    size_t maxR = 0;
    size_t maxC = 0;
    // assuming dims and pos return values major->minor; only the last 2 values are used
    for (size_t i = 0; i < in.size(); i++) {
        vector< size_t > pos = in.source()->position(in.index()[i]);
        if (inDims.size() > 1)  {
            rs[i] = pos[iR];
        }
        if (inDims.size() > 0) {
            cs[i] = pos[iC];
        }
        minR = min(minR, rs[i]);
        minC = min(minC, cs[i]);
        maxR = max(maxR, rs[i]);
        maxC = max(maxC, cs[i]);
    }
    _rows = (maxR - minR + 1);
    _cols = (maxC - minC + 1);
    _data = vector< T >(_rows * _cols, T());
    for (size_t i = 0; i < in.size(); i++) {
        _data[(rs[i] - minR) * _cols + (cs[i] - minC)] = static_cast< T >(in[i]);
    }
}


template < typename T >
template < typename iT, typename, typename >  // input type == generic container, requires size() and operator[](size_t)
matrix< T >::matrix(const storage::st_subset_base< iT >& in) {
    static_assert(is_generic_container< iT >(),  "INVALID INPUT CONTAINER!");  // should be generic container
    _rows = 1;
    _cols = in.size();
    _data = vector< T >(0);
    _data.reserve(in.size());
    for (size_t i = 0; i < (_rows * _cols); i++) {
        _data.push_back(in[i]);
    }
}


template < typename T >
template < typename iT, typename >  // input type == generic container, requires size() and operator[](size_t)
matrix< T >::matrix(size_t rows, size_t cols, const iT& in) {
    static_assert(is_generic_container< iT >(),  "INVALID INPUT CONTAINER!");
    assert(rows > 0 && cols > 0);
    assert(in.size() >=  rows * cols);
    _data.reserve(rows * cols);
    for (size_t i = 0; i < rows * cols; ++i) {
        _data.push_back(static_cast< T >(in[i]));
    }
    _rows = rows;
    _cols = cols;
}



template < typename T >
T& matrix< T >::operator*() {
    assert(_data.size() > 0);
    return _data[0];
}


template < typename T >
template < typename iT >
matrix< T >& matrix< T >::operator=(const iT& in) {  // input type == generic container, requires size() and operator[](int)
    static_assert(is_generic_container< iT >(), "");
    assert(in.size() >= _data.size());
    for (size_t i = 0; i < _data.size(); i++) {
        _data[i] = static_cast< T >(in[i]);
    }
    return *this;
}


template < typename T >
T& matrix< T >::operator()(size_t row, size_t col) {
    assert(row < _rows && col < _cols);
    return _data[row * _cols + col];
}


template < typename T >
const T& matrix< T >::operator()(size_t row, size_t col) const {
    assert(row < _rows && col < _cols);
    return _data[row * _cols + col];
}


template < typename T >
T& matrix< T >::operator[](size_t idx) {
    assert(idx < _data.size());
    return _data[idx];
}


template < typename T >
const T& matrix< T >::operator[](size_t idx) const {
    assert(idx < _data.size());
    return _data[idx];
}


template < typename T >
matrix_subset< T > matrix< T >::operator[](const vector< size_t >& idx) {
    return matrix_subset< T >(this, idx);  // storage::st_subset_base constructor asserts index validity
}


template < typename T >
matrix_subset_const< T > matrix< T >::operator[](const vector< size_t >& idx) const {
    return matrix_subset_const< T >(this, idx);  // storage::st_subset_base constructor asserts index validity
}


template < typename T >
matrix_subset< T > matrix< T >::operator[](const matrix_subset< T >& sbst) {
    if (sbst.source() != this && sbst.source_size() > size()) assert(sbst.index().back() < size());
    return matrix_subset< T >(sbst);  // copy constructor of subset
}


template < typename T >
matrix_subset_const< T > matrix< T >::operator[](const matrix_subset< T >& sbst) const {
    if (sbst.source() != this && sbst.source_size() > size()) assert(sbst.index().back() < size());
    return matrix_subset_const< T >(sbst);  // copy constructor of const_subset
}


template < typename T >
const vector< T >& matrix< T >::elements() const {
    return _data;
}


template < typename T >
const T* matrix< T >::data() const {
    return _data.data();
}


template < typename T >
size_t matrix< T >::rows() const {
    return _rows;
}


template < typename T >
size_t matrix< T >::cols() const {
    return _cols;
}


template < typename T >
size_t matrix< T >::size() const {
    return _data.size();
}


template < typename T >
vector< size_t > matrix< T >::shape() const {
    return vector< size_t >({_rows, _cols});
}


template < typename T >
vector< size_t > matrix< T >::position(size_t idx) const {
    // @todo -> review!
    // vector< size_t > dims = shape();
    // vector< size_t > ind(2);
    // size_t comp = idx;
    // for (size_t i = 0; i <= 2; i++) {
    //     size_t dim_sz = 1;
    //     for (size_t j = i+1; j < 2; j++) {
    //         dim_sz *= dims[j];
    //     }
    //     while (comp-dim_sz >= 0) {
    //         ind[i]++;
    //         comp -= dim_sz;
    //     }
    // }
    // return ind;

    assert(idx <= size());

    size_t ii = idx;
    size_t ri = 0;
    while (ii >= cols()) {
        ii -= cols();
        ri++;
    }
    size_t ci = idx - ri * cols();

    return { ri, ci };
}


template < typename T >
bool matrix< T >::isSquare() const {
    return (_rows == _cols);
}


template < typename T >
bool matrix< T >::isEmpty() const {
    return (!_rows && !_cols);
}


template < typename T >
void matrix< T >::reserve(size_t rows, size_t cols) {
    _data.reserve(rows * cols);
}


template < typename T >
vector< size_t > matrix< T >::allID() const {
    size_t i = 0;
    vector< size_t > ids(0);
    ids.reserve(size());
    for (auto& val : _data) {
        ids.push_back(i++);
    }
    return ids;
}


template < typename T >
vector< size_t > matrix< T >::rowID(size_t row_) const {
    assert(row_ < _rows);
    vector< size_t > ids(0);
    ids.reserve(_cols);
    for (size_t c = 0; c < _cols; c++) {
        ids.push_back(row_ * _cols + c);
    }
    return ids;
}


template < typename T >
vector< size_t > matrix< T >::colID(size_t col_) const {
    assert(col_ < _cols);
    vector< size_t > ids(0);
    ids.reserve(_rows);
    for (size_t r = 0; r < _rows; r++) {
        ids.push_back(r * _cols + col_);
    }
    return ids;
}


template < typename T >
vector< size_t > matrix< T >::diagID() const {
    vector< size_t > ids(0);
    ids.reserve(_rows > _cols ? _rows : _cols);
    for (size_t i = 0; i < min(_rows, _cols); i++) {
        ids.push_back(i * (_cols + 1));
    }
    return ids;
}


template < typename T >
vector< size_t > matrix< T >::blockID(size_t first_row, size_t last_row, size_t first_col, size_t last_col) const {
    // NOTE: does not work with 1 col/1 line matrix
    assert(first_row < last_row && last_row <= rows());
    assert(first_col <= last_col && last_col <= cols());
    if (!last_col) {
        last_col = cols();
    }
    vector< size_t > ids(0);
    ids.reserve(size());
    for (size_t r = first_row; r < last_row; ++r) {
        for (size_t c = first_col; c < last_col; ++c) {
            ids.push_back(r * _cols + c);
        }
    }
    return ids;
}


template < typename T >
matrix_subset< T > matrix< T >::all() {
    return matrix_subset< T >(this, allID());
}


template < typename T >
matrix_subset_const< T > matrix< T >::all() const {
    return matrix_subset_const< T >(this, allID());
}


template < typename T >
matrix_subset< T > matrix< T >::row(size_t row) {
    return matrix_subset< T >(this, rowID(row));
}


template < typename T >
matrix_subset_const< T > matrix< T >::row(size_t row) const {
    return matrix_subset_const< T >(this, rowID(row));
}


template < typename T >
matrix_subset< T > matrix< T >::col(size_t col) {
    return matrix_subset< T >(this, colID(col));
}


template < typename T >
matrix_subset_const< T > matrix< T >::col(size_t col) const {
    return matrix_subset_const< T >(this, colID(col));
}


template < typename T >
matrix_subset< T > matrix< T >::diag() {
    return matrix_subset< T >(this, diagID());
}


template < typename T >
matrix_subset_const< T > matrix< T >::diag() const {
    return matrix_subset_const< T >(this, diagID());
}


template < typename T >
matrix_subset< T > matrix< T >::block(size_t first_row, size_t last_row, size_t first_col, size_t last_col) {
    return matrix_subset< T >(this, blockID(first_row, last_row, first_col, last_col));
}


template < typename T >
matrix_subset_const< T > matrix< T >::block(size_t first_row, size_t last_row, size_t first_col, size_t last_col) const {
    return matrix_subset_const< T >(this, blockID(first_row, last_row, first_col, last_col));
}


template < typename T >
template < typename oT >
matrix< oT > matrix< T >::submat(size_t first_row, size_t last_row, size_t first_col, size_t last_col) const {
    return matrix< oT >(block(first_row, last_row, first_col, last_col));
}


template < typename T >
template < typename iT, typename >
void matrix< T >::set(const iT& in) {
    static_assert(is_generic_container< iT >(), "INVALID INPUT CONTAINER!");
    assert(in.size() >= _data.size());
    for (size_t i = 0; i < _data.size(); i++) {
        _data[i] = in[i];
    }
}


template < typename T >
template < typename iT, typename, typename >
void matrix< T >::set(const iT& in) {
    set(vector< value_type >(in));
}


template < typename T >
template < typename iT, typename >
void matrix< T >::fill(const iT& in) {
    for (size_t i = 0; i < _data.size(); i++) {
        _data[i] = in;
    }
}


template < typename T >
template < typename iT, typename >
void matrix< T >::pushRow(const iT& in) {
    static_assert(is_generic_container< iT >(), "INVALID INPUT CONTAINER!");
    assert(in.size() > 0);
    if (_cols) {
        assert(in.size() == _cols);
    }
    for (auto& v : in) {
        _data.push_back(v);
    }
    if (!_cols) {
        _cols = in.size();
    }
    _rows++;
}


template < typename T >
template < typename iT, typename, typename >
void matrix< T >::pushRow(const iT& in) {
    pushRow(vector< value_type >(in));
}


template < typename T >
template < typename iT, typename >
void matrix< T >::pushCol(const iT& in) {
    static_assert(is_generic_container< iT >(), "INVALID INPUT CONTAINER!");
    assert(in.size() > 0);
    if (_rows) {
        assert(in.size() == _rows);
    }
    for (int r = in.size() - 1; r >= 0; r--) {
        _data.insert(_data.begin() + (r * _cols + _cols), in[r]);
    }
    if (!_rows) {
        _rows = in.size();
    }
    _cols++;
}


template < typename T >
template < typename iT, typename, typename >
void matrix< T >::pushCol(const iT& in) {
    pushCol(vector< value_type >(in));
}


template < typename T >
void matrix< T >::pushRow() {
    assert(_rows > 0 && _cols > 0);              // can only add empty row in non-empty matrices (last value is used as placeholder)
    if constexpr (is_convertible< float, T >()) {
        pushRow(vector< T > (_cols, static_cast< T >(0.0)));
        // // alternative #1 (more efficient, no temporary vector)  [needs testing!]
        // for (int c = 0; c < _cols; c++) { _data.push_back(static_cast< T >(0.0)); }
        // // alternative #2 call insert  [needs testing!]
        // _data.insert(_data.end(), _cols, static_cast< T >(0.0));
    } else {
        pushRow(vector< T > (_cols /* elements are default-constructed */));
    }
}


template < typename T >
void matrix< T >::pushCol() {
    assert(_rows > 0 && _cols > 0);             // can only add empty col in non-empty matrices (last value is used as placeholder)
    // pushCol(vector< T > (_rows, _data.back()));
    if constexpr (is_convertible< float, T >()) {
        pushCol(vector< T > (_rows, static_cast< T >(0.0)));
    } else {
        pushCol(vector< T > (_rows /* elements are default-constructed */));
    }
}


template < typename T >
void matrix< T >::popRow() {
    assert(_rows > 0);
    if (_rows == 1) {
        clear();
    } else {
        _data.resize((_rows-1)*_cols);
        _rows--;
    }
}


template < typename T >
void matrix< T >::popCol() {
    assert(_cols > 0);
    if (_cols == 1) {
        clear();
    } else {
        for (int r = _rows-1; r >= 0; r--) {
            _data.erase(_data.begin() + (r * _cols + (_cols - 1)));
        }
        _cols--;
    }
}


template < typename T >
void matrix< T >::deleteRow(size_t _r) {
    assert(_r < _rows);
    if (_rows == 1) {
        clear();
    } else {
        _data.erase(_data.begin() + (_r * _cols), _data.begin() + (_r + 1 * _cols));  // range erase only works for contiguous _data (rows)
        _rows--;
    }
}


template < typename T >
void matrix< T >::deleteCol(size_t _c) {
    assert(_c < _cols);
    if (_cols == 1) {
        clear();
    } else {
        auto old_cols = _cols;
        for (int r = _rows-1; r >= 0; r--) {
            // throw std::runtime_error("std::matrix<T>::" + std::string(__func__) + " is not working at the moment!")
            _data.erase(_data.begin() + (r * old_cols + _c));
        }
        _cols--;
    }
}


template < typename T >
void matrix< T >::reshape(size_t new_rows, size_t new_cols) {
    if (!new_rows || !new_cols) {
        clear();
    }
    if (!_rows) {
        _rows = new_rows;
    }
    // match cols
    while (_cols > new_cols) {
        popCol();
    }
    while (_cols < new_cols) {
        pushCol(vector< T > (_rows, T()));
    }
    // match rows
    while (_rows > new_rows) {
        popRow();
    }
    while (_rows < new_rows) {
        pushRow(vector< T > (_cols, T()));
    }
}


template < typename T >
void matrix< T >::resize(size_t new_rows, size_t new_cols) {  // resize() is much less efficient than reshape()
    if (!new_rows || !new_cols) {
        clear();
    }
    if (new_cols == _cols) {
        // faster && straightforward due to row-major storage
        _data.resize(new_rows * _cols);
    } else {
        // _data needs to be reshaped for the second row onwards
        vector< T > tmp_data = _data;
        _data.resize(new_rows * new_cols);
        for (size_t r = 1; r < new_rows; r++) {
            for (size_t c = 0; c < new_cols; c++) {
                _data[r * new_cols + c] = tmp_data[r * _cols + c];
            }
        }
        tmp_data.clear();
    }
    _rows = new_rows;
    _cols = new_cols;
}


template < typename T >
void matrix< T >::flip() {
    matrix< T > tmp = *this;
    _rows = tmp.cols();
    _cols = tmp.rows();
    for (size_t r = 0; r < _rows; ++r) {
        row(r) = tmp.col(r);
    }
}


template < typename T >
void matrix< T >::clear() {
    _data.resize(0);
    _rows = 0;
    _cols = 0;
}


template < typename T >
typename vector< T >::iterator matrix< T >::begin() {
    return _data.begin();
}


template < typename T >
typename vector< T >::iterator matrix< T >::end() {
    return _data.end();
}


template < typename T >
typename vector< T >::const_iterator matrix< T >::begin() const {
    return _data.cbegin();
}


template < typename T >
typename vector< T >::const_iterator matrix< T >::end() const {
    return _data.cend();
}


template < typename T >
template < typename oT >
matrix< T >::operator vector< oT >&() {
    // return vector< oT >(_data.begin(), _data.end());
    return _data;
}


template < typename T >
template < typename oT >
matrix< T >::operator const vector< oT >&() const {
    // return vector< oT > (_data.begin(), _data.end());
    return _data;
}


template < typename T >
template < typename oT >
matrix< T >::operator matrix< oT >&() {
    return matrix< oT > (_rows, _cols, vector< oT >(_data.begin(), _data.end()));
}


template < typename T >
template < typename oT >
matrix< T >::operator const matrix< oT >&() const {
    return matrix< oT > (_rows, _cols, vector< oT >(_data.begin(), _data.end()));
}

/// @endcond

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// COMPARISON OPERATORS (NON-MEMBERS)
// negation(!) operator is defined by default
// with the exception of the == operator, only matrix size is compared
// use subset boolean operators (e.g. mat.all()==vol.level(int) or mat.col(int)==mat.row(int)) for element-wise comparison
template <typename _ltype, typename _rtype>
inline bool operator==(const std::matrix<_ltype>& _lmat, const std::matrix<_rtype>& _rmat) {
    if (_lmat.size() != _rmat.size()) {
        return false;
    } else {
        for (size_t i=0; i < _lmat.size(); i++) {
            if (_lmat[i] != _rmat[i]) return false;
        }
    }
    return true;
}
template <typename _ltype, typename _rtype>
inline bool operator!=(const std::matrix<_ltype>& _lmat, const std::matrix<_rtype>& _rmat) {
    return !operator==(_lmat, _rmat);
}
template <typename _ltype, typename _rtype>
inline bool operator< (const std::matrix<_ltype>& _lmat, const std::matrix<_rtype>& _rmat) {
    return (_lmat.size() < _rmat.size());
}
template <typename _ltype, typename _rtype>
inline bool operator> (const std::matrix<_ltype>& _lmat, const std::matrix<_rtype>& _rmat) {
    return  operator< (_rmat, _lmat);
}
template <typename _ltype, typename _rtype>
inline bool operator<=(const std::matrix<_ltype>& _lmat, const std::matrix<_rtype>& _rmat) {
    return !operator> (_lmat, _rmat);
}
template <typename _ltype, typename _rtype>
inline bool operator>=(const std::matrix<_ltype>& _lmat, const std::matrix<_rtype>& _rmat) {
    return !operator< (_lmat, _rmat);
}
/////////////////////////////////////////////////////////////
// template <typename _ltype, typename _rtype>
// inline bool operator==(const std::matrix<_ltype>& _lmat, const std::storage::st_subset_base<_rtype>& _rcont){
//  if (_lmat.size() != _rcont.size()) return false;
//  else for (int i=0; i<_lmat.size(); i++) if (_lmat[i]!=_rcont[i]) return false;
//  return true;
// };
// template <typename _ltype, typename _rtype>
// inline bool operator!=(const std::matrix<_ltype>& _lmat, const std::storage::st_subset_base<_rtype>& _rcont){
//  return !operator==(_lmat,_rcont);
// };
// template <typename _ltype, typename _rtype>
// inline bool operator< (const std::matrix<_ltype>& _lmat, const std::storage::st_subset_base<_rtype>& _rcont){
//  return (_lmat.size() < _rcont.size());
// };
// template <typename _ltype, typename _rtype>
// inline bool operator> (const std::matrix<_ltype>& _lmat, const std::storage::st_subset_base<_rtype>& _rcont){
//  return  operator< (_rcont,_lmat);
// };
// template <typename _ltype, typename _rtype>
// inline bool operator<=(const std::matrix<_ltype>& _lmat, const std::storage::st_subset_base<_rtype>& _rcont){
//  return !operator> (_lmat,_rcont);
// };
// template <typename _ltype, typename _rtype>
// inline bool operator>=(const std::matrix<_ltype>& _lmat, const std::storage::st_subset_base<_rtype>& _rcont){
//  return !operator< (_lmat,_rcont);
// };
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace std
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
