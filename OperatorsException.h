#ifndef EX3_OPERATORSEXCEPTIONS_H
#define EX3_OPERATORSEXCEPTIONS_H

#include <exception>
#define ILLEGAL_ADDITION "The given addition is impossible, make sure the matrices are of the " \
                         "same dimensions."

#define ILLEGAL_SUBTRACTION "The given subtraction is impossible, make sure the matrices are of " \
                            "the same dimensions."

#define ILLEGAL_MULTIPLICATION "The given multiplication is impossible, make sure that the left" \
                               "matrix has same number of rows as the number of columns in the " \
                               "right one."

/**
 * class representing exception in matrix operations.
 */
class IllegalAddition : public std::exception
{
    virtual const char* what() const throw()
    {
        return ILLEGAL_ADDITION;
    }
};

/**
 * class representing exception in matrix operations.
 */
class IllegalSubtraction : public std::exception
{
    virtual const char* what() const throw()
    {
        return ILLEGAL_SUBTRACTION;
    }
};

/**
 * class representing exception in matrix operations.
 */
class IllegalMultiplication : public std::exception
{
    virtual const char* what() const throw()
    {
        return ILLEGAL_MULTIPLICATION;
    }
};

#endif