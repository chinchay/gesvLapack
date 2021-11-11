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
        swap0!(X, k0, k)
    end
end

function swap0!(X, k0, k)
        temp  = X[k0]
        X[k0] = X[k]
        X[k]  = temp
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


function swapRowARowB!(M, rowA, rowB, colStart, colEnd)
    if rowA != rowB
        for col in colStart:colEnd
            temp = M[rowA, col]
            M[rowA, col] = M[rowB, col]
            M[rowB, col] = temp
        end
    end
end

function swapManyRows!(M, listOfPairRows, k1, k2, colStart, colEnd)
    for rowA in k1:k2
        rowB = listOfPairRows[rowA]
        swapRowARowB!(M, rowA, rowB, colStart, colEnd)
    end

end

# modified from http://www.netlib.org/lapack/explore-html/d7/d6b/dlaswp_8f_source.html#l00114
function dlaswp!(A, colMin, colMax, ipiv, k1, k2)
    n = colMax - colMin + 1
    n32 = div(n, 32) * 32
    #
    if n32 != 0
        for i in 1:32:n32
            col_o = colMin + i - 1
            col_f = col_o + 31
            swapManyRows!(A, ipiv, k1, k2, col_o, col_f)
    end
    #
    if n32 != n
        n32 += 1
        col_o = colMin + n32 - 1
        col_f = n
        swapManyRows!(A, ipiv, k1, k2, col_o, col_f)
    end
end

# function initialize(rangeA)
#     global T = zeros(Float32, rangeA, rangeA)
# end

# adapted from:
# http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_ga6a0a7704f4a747562c1bd9487e89795c.html#ga6a0a7704f4a747562c1bd9487e89795c
function trsm_LLNU!(A, ro, co, T, lo, m, n)
    # dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) 
    # side  = 'L'
    # ul    = 'L'
    # tA    = 'N'
    # dA    = 'U'   
    # T was already defined in initialize()
    # global T .= 0.0
    # B .= 0.0 # no need to do that

    # copy A[l2,l2+1,l2+2,...] into B
    cont = 1
    mn   = m * n
    for l in lo:mn
        T[cont] = A[l]
        cont += 1
    end
    

    # ... things
    # lside = true # if lsame(side,'L')
    # nrowa = m # if lsame(side,'L')
    # upper = false  # lsame(uplo,'U')
    # 
    # nounit = false # lsame(diag,'N')    
    
    info = 0
    #
    # *           Form  B := inv( A ) * B.
    for j in 1:n
        for k in 1:m
            if B[k,j] != 0.0
                K = co + k - 1
                for i in (k + 1):m
                    I = ro + i - 1
                    T[i,j] -= ( T[k,j] * A[I,K] )
                end
            end
        end
    end
    # 

    # mutate A:
    cont = 1
    for l in lo:mn
        A[l] = T[cont]
        cont += 1
    end
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