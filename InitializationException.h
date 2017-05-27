#ifndef EX3_INITIALIZATIONEXCEPTION_H
#define EX3_INITIALIZATIONEXCEPTION_H

#include <exception>
#define ILLEGAL_DIMENSIONS "The given dimensions are illegal, make sure you didn't try to "\
                           "initialize a matrix with non positive numbers."

#define OUT_OF_BOUNDS "The given row or column number is illegal, make sure both are not "\
                      "negative numbers."

#define NON_SQUARE_MATRIX "Cant transpose a non-square matrix."

/**
 * class representing exception in matrix size.
 */
class IllegalDimensions : public std::exception
{
    virtual const char* what() const throw()
    {
        return ILLEGAL_DIMENSIONS;
    }
};

/**
 * class representing exception in matrix size.
 */
class CoordinatesOutOfBounds : public std::exception
{
    virtual const char* what() const throw()
    {
        return OUT_OF_BOUNDS;
    }
};

/**
 * class representing exception in matrix size.
 */
class CantTransposeNonSquare : public std::exception
{
    virtual const char* what() const throw()
    {
        return NON_SQUARE_MATRIX;
    }
};

#endif
