module Utils

export finIndxMax, swap!, multiply!, getIndexFortranReading

function finIndxMax(X, iStart, iEnd)
    # Julia, as Fortran, read matrices row by row for each column
    max = abs(X[iStart])
    imax = iStart
    for i in (iStart + 1):iEnd
        xi = abs(X[i])
        if xi > max
            imax = i
            max  = xi
        end
    end
    return imax
end

# mutates X
function swap!(X, k0, k)
    if k != k0
        temp  = X[k0]
        X[k0] = X[k]
        X[k]  = temp
    end
end

# mutates X
function unrolledMultiply!(X, alpha, iStart, iEnd)
    n = iEnd - iStart + 1
    m = mod(n ,5)
    if m != 0
        for i in iStart:(iStart + m - 1)
            X[i] *= alpha
        end    
    end
    # 
    if n >= 5
        for i in (iStart + m):5:iEnd
            X[i]     *= alpha
            X[i + 1] *= alpha
            X[i + 2] *= alpha
            X[i + 3] *= alpha
            X[i + 4] *= alpha
        end
    end
end

function unrolledDivide!(X, alpha, iStart, iEnd)
    n = iEnd - iStart + 1
    m = mod(n ,5)
    if m != 0
        for i in iStart:(iStart + m - 1)
            X[i] /= alpha
        end    
    end
    # 
    if n >= 5
        for i in (iStart + m):5:iEnd
            X[i]     /= alpha
            X[i + 1] /= alpha
            X[i + 2] /= alpha
            X[i + 3] /= alpha
            X[i + 4] /= alpha
        end
    end
end


function divide!(X, divisor, iStart, iEnd, sfmin)
    if abs(divisor) >= sfmin
        mutltiplicador = 1.0 / divisor
        unrolledMultiply!(X, mutltiplicador, iStart, iEnd)
    else
        unrolledDivide!(X, divisor, iStart, iEnd)
    end
end

function getIndexFortranReading(row, col, nRange)
    # Similar as Fortran reading
    k = ((col - 1) * nRange) + row
    return k
end


end
## to understand how A(j,j) is passed by reference in Fortran:

# Program Hello
# implicit none
# Print *, "Hello World"
# integer :: X(3,3)
# X(1,1) = 1
# X(2,1) = 2
# X(3,1) = 3
# X(1,2) = 4
# X(2,2) = 5
# X(3,2) = 6
# X(1,3) = 7
# X(2,3) = 8
# X(3,3) = 9
# call print_a( X(3,1) )


# contains
# subroutine print_a( A )
#   integer, intent(in) :: A(3,*)
#   integer :: i,j
#   do i = 1, 3
#         write(*,*) A(i,1), A(i,2), A(i,3)
#   end do
# end subroutine


# End Program Hello

# output
# 3           6           9
# 4           7           0 <<<< you should handle those 2 numbers so that they are not used by the rest of the program
# 5           8 -1496688640 <<<< I think in Fortran they were somehow avoided or converted to zero.