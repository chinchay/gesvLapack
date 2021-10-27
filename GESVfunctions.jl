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
    dgetrf(n, n, A, n, ipiv, info)

    if info == 0
        #* Solve the system A*X = B, overwriting B with X.
        # * A and B are mutated!
        # CALL dgetrs( 'No transpose', n, nrhs, a, lda, ipiv, b, ldb, info )
        dgetrs!("No transpose", n, nrhs, A, lda, ipiv, B, ldb, info)
    end
end
