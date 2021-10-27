module GESVfunctions

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
    dgetrf!(n, n, A, n, ipiv, info)
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
function dgetrf!( M, N, A, LDA, IPIV, INFO )
    one = 1.0

    # Determine the block size for this environment.
    nb = ilaenv(1, "DGETRF", " ", n, n, -1, -1)
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






end # module