module GESVfunctions

const dlmach_val = dlmach()

# source (LAPACK is a software package provided by Univ. of Tennessee):
# this is a modification for Julia language, for square matrices!
# http://www.netlib.org/lapack/explore-html/d8/d72/dgesv_8f_source.html
function dgesv!(n, A, B, ipiv)
    # n : range of the square matrix
    # A: dimension matrix (n,n), lds
    # IPIV is INTEGER array, dimension (n)
    # B: dimension matrix (n,n)
    # INFO is INTEGER see http://www.netlib.org/lapack/explore-html/d8/d72/dgesv_8f_source.html
    
    #* Compute the LU factorization of A.
    # CALL dgetrf( n, n, a, lda, ipiv, info )
    dgetrf!(A, n, ipiv, info)
    # A is mutated according to:
    # http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html#ga0019443faea08275ca60a734d0593e60
    # "A is DOUBLE PRECISION array, dimension (LDA,N)
    # On entry, the M-by-N matrix to be factored.
    # On exit, the factors L and U from the factorization
    # A = P*L*U; the unit diagonal elements of L are not stored.    "

    if info == 0
        #* Solve the system A*X = B, overwriting B with X.
        # * A and B are mutated!
        # CALL dgetrs( 'No transpose', n, nrhs, a, lda, ipiv, b, ldb, info )
        dgetrs!("No transpose", n, nrhs, A, lda, ipiv, B, ldb, info)
    end
end
# needs to define dgetrf, dgetrs for dgesv

# needs to define dgemm, dgetrf2, dlaswp, dtrsm, xerbla, ilaenv for dgetrf
# A is mutated
function dgetrf!(A, n, ipiv, info)
    one = 1.0

    # Determine the block size for this environment.
    # nb = ilaenv(1, "DGETRF", " ", n, n, -1, -1)
    # for the above case, it resulted nb=64:
    nb = 64 #
    if (nb <= 1) || (nb >= n) 
        # Use unblocked code.
        dgetrf2(n, n, A, n, ipiv, info)
    else
        #  Use blocked code.
        for j in 1:nb:n # 1, 1+nb, 1+(2*nb), ..., n
            jb = min(n - j + 1, nb)
            
            # Factor diagonal and subdiagonal blocks and test for exact singularity.
            dgetrf2( n - j + 1, jb, A[j,j], n, ipiv[j], iinfo )
            # TODO: it returns iinfo, rearrange it

            #  Adjust INFO and the pivot indices.
            if (info == 0) && (iinfo > 0)
                info = iinfo + j - 1
            end

            for i in j: min(n, j + jb - 1)
                ipiv[i] += j - 1
            end
            
            # Apply interchanges to columns 1:J-1.
            dlaswp( j - 1, A, n, j, j + jb - 1, ipiv, 1 )

            if j + jb <= n
                # Apply interchanges to columns J+JB:N
                dlaswp(n - j - jb + 1, A[1, j + jb], n, j, j + jb - 1, ipiv, 1)

                # Compute block row of U.
                dtrsm("Left", "Lower", "No transpose", "Unit", jb, n - j - jb + 1, one, A[j,j], n, A[j, j + jb], n )

                if j + jb <= m
                    # Update trailing submatrix.
                    dgemm("No transpose", "No transpose", n - j - jb + 1, n - j - jb + 1, jb, -one, A[j + jb, j], lda, A[j, j + jb], n, one, A[j + jb, j + jb], n )
                end
            end
end


#TODO: define ilaenv function
#* I found that ilaenv just returns nb=64 for `ilaenv(1, "DGETRF", " ", n, n, -1, -1)` in dgetrf!(...)
# nb =   ilaenv(1,     "DGETRF", " ", n, n, -1, -1)
# function ilaenv( ISPEC, NAME, OPTS,  N1, N2, N3, N4 )
# end

# dgetrf2(n, m, X, p, ipiv, info)
function dgetrf2( m, n, A, lda, IPIV, info )
    one  = 1.0
    zero = 0.0
    #
    if m == 1
        #* Use unblocked code for one row case
        #* Just need to handle IPIV and INFO
        IPIV[1] = 1
        if A[1, 1] == zero
            info = 1
        end
    elseif n == 1
        #* Use unblocked code for one column case

        #* Compute machine safe minimum
        # sfmin = dlamch('S')
        sfmin = dlmach_val # defined as global variable in the current module

        #* Find pivot and test for singularity
        i = idamax( m, A[1, 1], 1 )
        IPIV[1] = i
        
        if  A[i, 1] != zero
            #* Apply the interchange
            if i != 1
                temp    = A[1, 1]
                A[1, 1] = A[i, 1]
                A[i, 1] = temp
            end

            #* Compute elements 2:M of the column                
            if abs( A[1, 1] ) >= sfmin
                dscal( m - 1, one / A[1, 1], A[2, 1], 1 )
            else
                for i in 1:(m - 1)
                    A[1 + i, 1] = A[1 + i, 1] / A[1, 1]
                end
            end
        else
            info = 1
        end
    else
        min_mn = min(m, n)
        #* Use recursive code
        n1 = min_mn / 2
        n2 = n - n1

        # *               [ A11 ]
        # *        Factor [ --- ]
        # *               [ A21 ]   
        dgetrf2( m, n1, A, lda, IPIV, iinfo )  
        
        if (info == 0) && (iinfo > 0)
            info = iinfo
        end

        # *                              [ A12 ]
        # *        Apply interchanges to [ --- ]
        # *                              [ A22 ]
        n1plus1  = n1 + 1
        dlaswp( n2, A[1, n1plus1], lda, 1, n1, IPIV, 1 )

        # *        Solve A12            
        dtrsm( 'L', 'L', 'N', 'U', n1, n2, one, A, lda, A[1, n1plus1], lda )        

        # *        Update A22
        dgetrf2( m - n1, n2, A[n1plus1, n1plus1], lda, IPIV[n1plus1], iinfo )   
        
        # *        Factor A22
        dgetrf2( m - n1, n2, A[n1plus1, n1plus1], lda, IPIV[n1plus1], iinfo )  
        
        #* Adjust INFO and the pivot indices
        if (info == 0) && (iinfo > 0)
            info = iinfo + n1
        end

        for i in n1plus1:min_mn
            IPIV[i] += n1
        end

        #* Apply interchanges to A21
        dlaswp( n1, A[1, 1], lda, n1plus1, min_mn, IPIV, 1 )
    end
end


# http://www.netlib.org/lapack/explore-html/d5/dd4/dlamch_8f_source.html#l00068
function dlmach() # for dlamch('S')
    # https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/
    # 
    # The distance between 1.0 and the next larger representable floating-point value of Float64.
    # eps() : same as eps(Float64) : 2.220446049250313e-16    
    # eps(0.0) : 5.0e-324 
    # prevfloat(Inf)  : 1.7976931348623157e308
    # nextfloat(-Inf) : -1.7976931348623157e308 

    zero = 0.0
    one  = 1.0
    #
    # eps_ = eps(zero) * 0.5 # it will return 0.0 in REPL julia
    eps_ = eps(zero) # it returns 5.0e-324
    #
    small = one / prevfloat(Inf)
    sfmin = eps(zero) # tiny(zero)
    if small >= sfmin
        #* Use SMALL plus a bit, to avoid the possibility of rounding
        #* causing overflow when computing  1/sfmin.
        sfmin = small * (one + eps_)
    end
    return sfmin
end


end # module