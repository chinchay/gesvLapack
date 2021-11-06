module Utils


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