#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "InitializationException.h"
#include "OperatorsException.h"
#include "Complex.h"
#include <vector>
#include <iostream>
#include <thread>
#include <mutex>

/** The number of rows in the default matrix */
const int MINIMAL_NUM_ROWS = 1;

/** The number of columns in the default matrix */
const int MINIMAL_NUM_COLS = 1;

/** enum representing the type of exception should be thrown. */
enum EXCEPTION_TYPE
{
    INITIALIZE,
    ADDITION,
    SUBTRACTION,
    MULTIPLICATION,
    BOUNDS,
    SQUARE
};

bool PARALLEL_FLAG = false;

/**
 * Template class representing the matrix object.
 */
template <class T>
class Matrix
{
private:

    /** The number of rows and columns in the matrix. */
    unsigned int _rows, _cols;

    /** A vector holding all the variables of the matrix.*/
    std::vector<T> _matrix;

    /** A void function checking if the given arguments fit the request, for example checking
     * if the matrix is square before calculating it's transpose.
     * @param row - The number of rows of the other matrix (if the operator acts on two matrices).
     * @param col - The number of columns of the other matrix (if the operator acts on two
     *              matrices).
     * @param e - The type of exception that should be thrown if the matrix does not fit.
     */
    void _isRequestLegal(const unsigned int row, const unsigned int col,
                         const EXCEPTION_TYPE& e) const;

    /** The mutex object, used to make the threads update the result matrix in order. */
    std::mutex mtx;

public:

    /**
     * A static method, used to update the static member PARALLEL_FLAG.
     * @param isParallel - updates the flag to this, if the flag has been changed an informative
     * message will be printed.
     */
    static void setParallel(bool isParallel)
    {
        if(isParallel != PARALLEL_FLAG)
        {
            if(isParallel)
            {
                std::cout << "Generic Matrix mode changed to Parallel mode" << std::endl;
            }
            else
            {
                std::cout << "Generic Matrix mode changed to non-Parallel mode" << std::endl;
            }
            PARALLEL_FLAG = isParallel;
        }
    }

    /**
     * Parallel multiplication.
     * @param other - The given matrix.
     * @return - A new matrix.
     */
    Matrix<T> parallelMult(const Matrix<T> &other) const;

    /**
     * Parallel addition.
     * @param other - The given matrix.
     * @return - A new matrix.
     */
    Matrix<T> parallelAdd(const Matrix<T> &other) const;

    /**
     * Calculates a single row of the result matrix in the addition of two matrices.
     * @param row - The number of the row in the result matrix.
     * @param other - the RHS matrix.
     * @param result - the result matrix.
     */
    void _addRow(unsigned int row, const Matrix<T> &other, Matrix<T> &result) const;

    /**
     * Calculates a single row of the result matrix in the multiplication of two matrices.
     * @param row - The number of the row in the result matrix.
     * @param other - the RHS matrix.
     * @param result - the result matrix.
     */
    void _multRow(unsigned int row, const Matrix<T> &other, Matrix<T> &result) const;

    /**
     * Constructor, creating a matrix with default values.
     * @return - Matrix object.
     */
    Matrix();

    /**
     * Constructor, creating a matrix of the given size.
     * @param rows - The number of rows in the matrix.
     * @param cols - The number of colums in the matrix.
     * @return - Matrix object.
     */
    Matrix(unsigned int rows, unsigned int cols);

    /**
     * Constructor, creates a matrix of given size with the given cells in it.
     * @param rows - The number of rows in the matrix.
     * @param cols - The number of colums in the matrix.
     * @param cells - The given cells to enter to the matrix.
     * @return - Matrix object.
     */
    Matrix(unsigned int rows, unsigned int cols, const std::vector<T>& cells);

    /**
     * Copy constructor.
     */
    Matrix(const Matrix<T> &other);

    /**
     * Destructor.
     */
    ~Matrix();

    /**
     * Move constructor
     */
    Matrix(Matrix<T> && other);

    /**
     * Assignment operator "=", changes 'this' matrix to be a deep copy of the given matrix.
     * @param other - The given matrix.
     */
    Matrix<T> &operator=(const Matrix<T> &other);


    /**
     * Addition operator "+", creates a new matrix by the standart addition  operation between
     * matrices.
     * @param other - The given matrix.
     * @return - A new matrix.
     */
    Matrix<T> operator+(const Matrix<T> &other) const;

    /**
     * Subtraction operator "-", creates a new matrix by the standart subtraction operation between
     * matrices.
     * @param other - The given matrix.
     * @return - A new matrix.
     */
    Matrix<T> operator-(const Matrix<T> &other) const;

    /**
     * Multiplication operator "*", creates a new matrix by the standart multiplication operation
     * between matrices.
     * @param other - The given matrix.
     * @return - A new matrix.
     */
    Matrix<T> operator*(const Matrix<T> &other) const;

    /**
     * @return - The number of columns in the matrix.
     */
    unsigned int cols() const;

    /**
     * @return - The number of rows in the matrix.
     */
    unsigned int rows() const;

    /**
     * @return - True iff the matrix has same number of rows and columns.
     */
    bool isSquareMatrix() const;

    /**
     * Comparison operator "==".
     * @param other - The matrix to compare to.
     * @return - True iff they are of the same size and hold the same objects.
     */
    bool operator==(const Matrix<T> &other) const;

    /**
     * Comparison operator "!=".
     * @param other - The matrix to compare to.
     * @return - True iff they are not of the same size or hold different objects.
     */
    bool operator!=(const Matrix<T> &other) const;

  /**
    * Stream print operator <<
    * @return Prints the given matrix into the given stream, in the following format:
    *         Each row is printed in a different row, and after each value a tab is printed.
    */
    template <typename U>
    friend std::ostream& operator<<(std::ostream &os, const Matrix<U> &matrix);

    /**
     * Location operator.
     * @return The variable in the given row and column in the matrix (given as lvalue thus can
     *         be modified).
     */
    T &operator()(unsigned int row, unsigned int col);

    /**
     * Location operator.
     * @return The variable in the given row and column in the matrix.(given as rvalue thus can
     *         not be modified).
     */
    const T &operator()(unsigned int row, unsigned int col) const;

    /**
     * @return The transpose of this matrix
     */
    Matrix<T> trans() const;

    /** Redefine name for convenience purposes.*/
    typedef typename std::vector<T>::const_iterator const_iterator;
    /**
     * @return An iterator that starts at the end of the matrix
     */
    const_iterator end();

    /**
     * @return An iterator that starts at the beginning of the matrix
     */
    const_iterator begin();

};

///------------------------------------ Constructors ----------------------------------------------
template <class T>
Matrix<T>::Matrix()
{
    Matrix(MINIMAL_NUM_ROWS, MINIMAL_NUM_COLS);
}

template <class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols)
try : _rows(rows),
      _cols(cols),
      _matrix(rows * cols, T(0))
{
    _isRequestLegal(rows, cols, INITIALIZE);
}
catch (std::exception& e)
{
    std::cout << e.what() << std::endl;
    throw e;
}


template <class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const std::vector<T> &cells)
try : _rows(rows),
      _cols(cols),
      _matrix(cells)
{
    _isRequestLegal(rows, cols, INITIALIZE);
}
catch (std::exception& e)
{
    std::cout << e.what() << std::endl;
    throw e;
}


template <class T>
Matrix<T>::Matrix(Matrix<T> && other) : _rows(std::move(other.rows())),
                                       _cols(std::move(other.cols())),
                                       _matrix(std::move(other._matrix))
{
}

template <class T>
Matrix<T>::Matrix(const Matrix<T> &other) : _rows(other.rows()),
                                            _cols(other.cols()),
                                            _matrix(other._matrix)
{
}

template <class T>
Matrix<T>::~Matrix()
{
    _matrix.clear();
}

///------------------------------------- Operators ------------------------------------------------
template <class T>
bool Matrix<T>::operator==(const Matrix<T> &other) const
{
    return _matrix == other._matrix;
}

template <class T>
bool Matrix<T>::operator!=(const Matrix<T> &other) const
{
    return !(_matrix == other._matrix);
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const
{
    try
    {
        _isRequestLegal(other.rows(), other.cols(), ADDITION);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        throw e;
    }
    if(PARALLEL_FLAG)
    {
        return parallelAdd(other);
    }
    Matrix<T> res = Matrix(_rows, _cols);
    for(unsigned int row = 0; row < _rows; row++)
    {
        for(unsigned int col = 0; col < _cols; col++)
        {
            res(row, col) = (*this)(row, col) + other(row, col);
        }
    }
    return res;
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const
{
    try
    {
        _isRequestLegal(other.rows(), other.cols(), SUBTRACTION);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        throw e;
    }
    Matrix<T> res = Matrix(_rows, _cols);
    for(unsigned int row = 0; row < _rows; row++)
    {
        for(unsigned int col = 0; col < _cols; col++)
        {
            res(row, col) = (*this)(row, col) - other(row, col);
        }
    }
    return res;
}


template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const
{
    try
    {
        _isRequestLegal(other.rows(), other.cols(), MULTIPLICATION);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        throw e;
    }
    if(PARALLEL_FLAG)
    {
        return parallelMult(other);
    }
    Matrix<T> res = Matrix(_rows, other.cols());
    for(unsigned int thisRow = 0; thisRow < _rows; thisRow++)
    {
        for(unsigned int otherCol = 0; otherCol < other.cols(); otherCol++)
        {
            T sum = T(0);
            for(unsigned int thisCol = 0; thisCol < _cols; thisCol++)
            {
                sum += (*this)(thisRow, thisCol) * other(thisCol, otherCol);
            }
            res(thisRow, otherCol) = sum;
        }
    }
    return res;
}

template <class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &other)
{
    Matrix<T> tmp = Matrix(other);
    _rows = tmp.rows();
    _cols = tmp.cols();
    Matrix<T>& tmpRef = tmp;
    (*this)._matrix.swap(tmpRef._matrix);
    return *this;
}

template <class T>
const T& Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
    try
    {
        _isRequestLegal(row, col, BOUNDS);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        throw e;
    }
    return _matrix.at(row * _cols + col);
}

template <class T>
T& Matrix<T>::operator()(unsigned int row, unsigned int col)
{
    try
    {
        _isRequestLegal(row, col, BOUNDS);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        throw e;
    }
    return _matrix.at(row * _cols + col);
}

template <class T>
Matrix<T> Matrix<T>::trans() const
{
    try
    {
        _isRequestLegal(_rows, _cols, SQUARE);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        throw e;
    }
    Matrix<T> res = Matrix(_cols, _rows);
    for(unsigned int curRow = 0; curRow < _rows; curRow++)
    {
        for(unsigned int curCol = 0; curCol < _cols; curCol++)
        {
            T v = (*this)(curRow, curCol);
            res(curCol, curRow) = v;
        }
    }
    return res;
}

/**
 * A specific implementation of transpose for the Complex type
 * @return The transpose of the matrix.
 */
template <>
Matrix<Complex> Matrix<Complex>::trans() const
{
    try
    {
        _isRequestLegal(_rows, _cols, SQUARE);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        throw e;
    }
    Matrix<Complex> res = Matrix<Complex>(_cols, _rows);
    for(unsigned int i = 0; i < res._matrix.size(); i++)
    {
        int div = i / _rows;
        int mod = i % _rows;
        res._matrix[i] = _matrix[_cols * mod + div].conj();
    }
    return res;
}

template <class U>
std::ostream& operator<<(std::ostream &os, const Matrix<U> &matrix)
{
    unsigned int i, j;
    for(i = 0; i < matrix._rows; i++)
    {
        for(j = 0; j < matrix._cols; j++)
        {
            os << (matrix(i, j)) << "\t";
        }
        os << "\n";
    }
    return os;
}

///-------------------------------------- Iterators  ----------------------------------------------
template <class T>
typename Matrix<T>::const_iterator Matrix<T>::begin()
{
    return _matrix.begin();
}

template <class T>
typename Matrix<T>::const_iterator Matrix<T>::end()
{
    return _matrix.end();
}

///-------------------------------------- Getters -------------------------------------------------
template <class T>
unsigned int Matrix<T>::cols() const
{
    return _cols;
}

template <class T>
unsigned int Matrix<T>::rows() const
{
    return _rows;
}

///------------------------------------ Parallel Methods ------------------------------------------

template <class T>
Matrix<T> Matrix<T>::parallelAdd(const Matrix<T> &other) const
{
    Matrix<T> res = Matrix(_rows, _cols);
    std::vector<std::thread> threads(_rows);
    for(unsigned int i = 0; i < _rows; i++)
    {
        threads[i] = std::thread(&Matrix<T>::_addRow, this, i, std::ref(other), std::ref(res));
    }

    for(unsigned int i = 0; i < _rows; i++)
    {
        threads[i].join();
    }
    return res;
}

template <class T>
Matrix<T> Matrix<T>::parallelMult(const Matrix<T> &other) const
{
    Matrix<T> res = Matrix(_rows, other._cols);
    std::vector<std::thread> threads(_rows);
    for(unsigned int i = 0; i < _rows; i++)
    {
        threads[i] = std::thread(&Matrix<T>::_multRow, this, i, std::ref(other), std::ref(res));
    }
    for(unsigned int i = 0; i < _rows; i++)
    {
        threads[i].join();
    }
    return res;
}

template <class T>
void Matrix<T>::_addRow(unsigned int row, const Matrix<T> &other, Matrix<T> &result) const
{
    //Lock the result matrix letting one thread modification at a time.
    std::lock_guard<std::mutex>(result.mtx);
    for (unsigned int col = 0; col < _cols; col++)
    {
        result(row, col) = (*this)(row, col) + other(row, col);
    }
}

template <class T>
void Matrix<T>::_multRow(unsigned int row, const Matrix &other, Matrix &result) const
{
    //Lock the result matrix letting one thread modification at a time.
    std::lock_guard<std::mutex>(result.mtx);
    for(unsigned int otherCol = 0; otherCol < other.cols(); otherCol++)
    {
        T sum = T(0);
        for(unsigned int thisCol = 0; thisCol < _cols; thisCol++)
        {
            sum += (*this)(row, thisCol) * other(thisCol, otherCol);
        }
        result(row, otherCol) = sum;
    }
}

///--------------------------------------- Others -------------------------------------------------
template <class T>
bool Matrix<T>::isSquareMatrix() const
{
    return _rows == _cols;
}

template <class T>
void Matrix<T>::_isRequestLegal(const unsigned int row, const unsigned int col,
                                const EXCEPTION_TYPE &e) const
{
    if(e == INITIALIZE)
    {
        if ((row == 0 && col > 0) || (row > 0 && col == 0))
        {
            throw IllegalDimensions();
        }
    }
    else if(e == ADDITION)
    {
        if(_rows != row || _cols != col)
        {
            throw IllegalAddition();
        }
    }
    else if(e == SUBTRACTION)
    {
        if(_rows != row || _cols != col)
        {
            throw IllegalSubtraction();
        }
    }
    else if(e == MULTIPLICATION)
    {
        if(_cols != row)
        {
            throw IllegalMultiplication();
        }
    }
    else if(e == BOUNDS)
    {
        if(row < 0 || col < 0)
        {
            throw CoordinatesOutOfBounds();
        }
    }
    else if(e == SQUARE)
    {
        if(!isSquareMatrix())
        {
            throw CantTransposeNonSquare();
        }
    }
}
#endif