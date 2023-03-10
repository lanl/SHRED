module lapack
    implicit none
    
    ! This is the precision that LAPACK "d" routines were compiled with (typically
    ! double precision, unless a special compiler option was used while compiling
    ! LAPACK). This "dp" is only used in lapack.f90
    ! The "d" routines data type is defined as "double precision", so
    ! we make "dp" the same kind as 0.d0 ("double precision"), so
    ! as long as LAPACK and this file were compiled with the same compiler options,
    ! it will be consistent. (If for example all double precision is promoted to
    ! quadruple precision, it will be promoted both in LAPACK and here.)
    integer, parameter:: dp=kind(0.d0)
    
    interface

        SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
        import :: dp
        CHARACTER          TRANS
        INTEGER            LDA, M, N, NRHS
        INTEGER            INCX, INCY
        REAL(dp)           ALPHA, BETA
        REAL(dp)           A( LDA, * ), X( * ), Y( * )
        END SUBROUTINE
    
        SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        import :: dp
        INTEGER            INFO, LDA, LDB, N, NRHS
        INTEGER            IPIV( * )
        REAL(dp)           A( LDA, * ), B( LDB, * )
        END SUBROUTINE
    
        SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                           EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, &
                           IWORK, INFO )
        import :: dp
        CHARACTER          EQUED, FACT, TRANS
        INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
        REAL(dp)           RCOND
        INTEGER            IPIV( * ), IWORK( * )
        REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), &
                           C( * ), FERR( * ), R( * ), WORK( * ), X( LDX, * )
        END SUBROUTINE
    
        SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        import :: dp
        INTEGER            INFO, LDA, LDB, N, NRHS
        INTEGER            IPIV( * )
        COMPLEX(dp)        A( LDA, * ), B( LDB, * )
        END SUBROUTINE
    
        SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                           EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, &
                           WORK, RWORK, INFO )
        import :: dp
        CHARACTER          EQUED, FACT, TRANS
        INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
        REAL(dp)           RCOND
        INTEGER            IPIV( * )
        REAL(dp)           BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * )
        COMPLEX(dp)        A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), &
                           X( LDX, * )
        END SUBROUTINE
    
        SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
        import :: dp
        INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
        INTEGER            IPIV( * )
        REAL(dp)           AB( LDAB, * ), B( LDB, * )
        END SUBROUTINE
    
        SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
        import :: dp
        CHARACTER          UPLO
        INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
        INTEGER            IPIV( * )
        REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, &
                           LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, &
                           IWORK, INFO )
        import :: dp
        CHARACTER          FACT, UPLO
        INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS
        REAL(dp)           RCOND
        INTEGER            IPIV( * ), IWORK( * )
        REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
                           BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
        END SUBROUTINE
    
        SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
                           LIWORK, INFO )
        import :: dp
        CHARACTER          JOBZ, UPLO
        INTEGER            INFO, LDA, LIWORK, LWORK, N
        INTEGER            IWORK( * )
        REAL(dp)           A( LDA, * ), W( * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
                           VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                           LWORK, IWORK, IFAIL, INFO )
        import :: dp
        CHARACTER          JOBZ, RANGE, UPLO
        INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
        REAL(dp)           ABSTOL, VL, VU
        INTEGER            IFAIL( * ), IWORK( * )
        REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), &
                           Z( LDZ, * )
        END SUBROUTINE
    
        SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
                          BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
        import :: dp
        CHARACTER          JOBVL, JOBVR
        INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
        REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
                           B( LDB, * ), BETA( * ), VL( LDVL, * ), &
                           VR( LDVR, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, &
                           ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, IHI, &
                           LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, &
                           LWORK, IWORK, BWORK, INFO )
        import :: dp
        CHARACTER          BALANC, JOBVL, JOBVR, SENSE
        INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
        REAL(dp)           ABNRM, BBNRM
        LOGICAL            BWORK( * )
        INTEGER            IWORK( * )
        REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), &
                           BETA( * ), LSCALE( * ), RCONDE( * ), RCONDV( * ), &
                           RSCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                          LDVR, WORK, LWORK, INFO )
        import :: dp
        CHARACTER          JOBVL, JOBVR
        INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
        REAL(dp)           A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), &
                           WORK( * ), WR( * )
        END SUBROUTINE
    
        SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, &
                           VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, &
                           RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
        import :: dp
        CHARACTER          BALANC, JOBVL, JOBVR, SENSE
        INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
        REAL(dp)           ABNRM
        INTEGER            IWORK( * )
        REAL(dp)           A( LDA, * ), RCONDE( * ), RCONDV( * ), &
                           SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), &
                           WI( * ), WORK( * ), WR( * )
        END SUBROUTINE
    
        SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
                          WORK, LWORK, RWORK, INFO )
        import :: dp
        CHARACTER          JOBVL, JOBVR
        INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
        REAL(dp)           RWORK( * )
        COMPLEX(dp)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
                           WORK( * )
        END SUBROUTINE
    
        SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, &
                           LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, &
                           RCONDV, WORK, LWORK, RWORK, INFO )
        import :: dp
        CHARACTER          BALANC, JOBVL, JOBVR, SENSE
        INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
        REAL(dp)           ABNRM
        REAL(dp)           RCONDE( * ), RCONDV( * ), RWORK( * ), SCALE( * )
        COMPLEX(dp)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
                           WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
                           LWORK, IWORK, LIWORK, INFO )
        import :: dp
        CHARACTER          JOBZ, UPLO
        INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
        INTEGER            IWORK( * )
        REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
        END SUBROUTINE
    
        REAL(dp) FUNCTION DLAMCH( CMACH )
        import :: dp
        CHARACTER          CMACH
        END FUNCTION
    
        INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
        CHARACTER*( * )    NAME, OPTS
        INTEGER            ISPEC, N1, N2, N3, N4
        END FUNCTION
    
        SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
        import :: dp
        INTEGER            INFO, LDA, M, N
        INTEGER            IPIV( * )
        COMPLEX(dp)        A( LDA, * )
        END SUBROUTINE
    
        SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        import :: dp
        CHARACTER          TRANS
        INTEGER            INFO, LDA, LDB, N, NRHS
        INTEGER            IPIV( * )
        COMPLEX(dp)         A( LDA, * ), B( LDB, * )
        END SUBROUTINE
    
        SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        import :: dp
        INTEGER            INFO, LDA, LWORK, N
        INTEGER            IPIV( * )
        COMPLEX(dp)        A( LDA, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
        import :: dp
        INTEGER            INFO, LDA, M, N
        INTEGER            IPIV( * )
        REAL(dp)           A( LDA, * )
        END SUBROUTINE
    
        SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        import :: dp
        INTEGER            INFO, LDA, LWORK, N
        INTEGER            IPIV( * )
        REAL(dp)           A( LDA, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
        import :: dp
        CHARACTER          JOBZ, UPLO
        INTEGER            INFO, LDA, LWORK, N
        REAL(dp)           RWORK( * ), W( * )
        COMPLEX(dp)        A( LDA, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, &
                           LRWORK, IWORK, LIWORK, INFO )
        import :: dp
        CHARACTER          JOBZ, UPLO
        INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
        INTEGER            IWORK( * )
        REAL(dp)           RWORK( * ), W( * )
        COMPLEX(dp)        A( LDA, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE ZHEGVD( ITYPE,  JOBZ,  UPLO,  N,  A,  LDA,  B, LDB, W, &
                           WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &
                           INFO )
        import :: dp
        CHARACTER          JOBZ, UPLO
        INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N
        INTEGER            IWORK( * )
        REAL(dp)           RWORK( * ), W( * )
        COMPLEX(dp)        A( LDA, * ), B( LDB, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
                           WORK, LWORK, INFO )
        import :: dp
        INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
        REAL(dp)           RCOND
        INTEGER            JPVT( * )
        REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
                           WORK, LWORK, RWORK, INFO )
        import :: dp
        INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
        REAL(dp)           RCOND
        INTEGER            JPVT( * )
        REAL(dp)           RWORK( * )
        COMPLEX(dp)        A( LDA, * ), B( LDB, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, &
                           LDVT, WORK, LWORK, INFO )
        import :: dp
        CHARACTER          JOBU, JOBVT
        INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
        REAL(dp)           A( LDA, * ), S( * ),  U( LDU,  * ), VT( LDVT, * ), &
                           WORK( * )
        END SUBROUTINE
    
        SUBROUTINE ZGESVD( JOBU, JOBVT,  M,  N,  A,  LDA, S, U, LDU, VT, LDVT, &
                           WORK, LWORK, RWORK, INFO )
        import :: dp
        CHARACTER          JOBU, JOBVT
        INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
        REAL(dp)           RWORK( * ), S( * )
        COMPLEX(dp)        A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
        END SUBROUTINE
    
        SUBROUTINE DSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, &
                           LIWORK, INFO )
        import :: dp
        CHARACTER          JOBZ
        INTEGER            INFO, LDZ, LIWORK, LWORK, N
        INTEGER            IWORK( * )
        REAL(dp)           D( * ), E( * ), WORK( * ), Z( LDZ, * )
        END SUBROUTINE
    
        SUBROUTINE XERBLA( SRNAME, INFO )
        CHARACTER*(*)      SRNAME
        INTEGER            INFO
        END SUBROUTINE

        SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
          import :: dp
                 REAL(dp)   ALPHA,BETA
                 INTEGER K,LDA,LDC,N
                 CHARACTER TRANS,UPLO
                 REAL(dp)   A(LDA,*),C(LDC,*)
        END SUBROUTINE

    
    ! BLAS
    
        SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
        import :: dp
        INTEGER INCX,INCY,N
        COMPLEX(dp) ZX(*),ZY(*)
        END SUBROUTINE
    
        SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
        import :: dp
        integer :: INCX, INCY, N
        real(dp) :: DA, DX(*), DY(*)
        END SUBROUTINE
    
        SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        import :: dp
        DOUBLE PRECISION ALPHA,BETA
        INTEGER K,LDA,LDB,LDC,M,N
        CHARACTER TRANSA,TRANSB
        REAL(dp) A(LDA,*),B(LDB,*),C(LDC,*)
        END SUBROUTINE
    
        SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        import :: dp
        complex(dp) :: ALPHA,BETA
        INTEGER K,LDA,LDB,LDC,M,N
        CHARACTER TRANSA,TRANSB
        complex(dp) :: A(LDA,*),B(LDB,*),C(LDC,*)
        END SUBROUTINE
    
        real(dp) FUNCTION DNRM2(N,X,INCX)
        import :: dp
        integer :: INCX, N
        real(dp) :: X(*)
        END FUNCTION
    
        SUBROUTINE DSCAL(N,DA,DX,INCX)
        import :: dp
        real(dp) :: DA, DX(*)
        integer :: INCX, N
        END SUBROUTINE
    
        SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        import :: dp
        REAL(dp) ALPHA,BETA
        INTEGER LDA,LDB,LDC,M,N
        CHARACTER SIDE,UPLO
        REAL(dp) A(LDA,*),B(LDB,*),C(LDC,*)
        END SUBROUTINE
    
        SUBROUTINE DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
        import :: dp
        INTEGER  INFO, LDA, LWORK, M, N
        REAL(dp) A(LDA, *), TAU(*), WORK(*)
        END SUBROUTINE
    
        SUBROUTINE DORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
        import :: dp
        INTEGER  INFO, K, LDA, LWORK, M, N
        REAL(dp) A(LDA,*), TAU(*), WORK(*)
        END SUBROUTINE
    
        SUBROUTINE DPOTRF(UPLO, N, A, LDA, INFO)
        import :: dp
        CHARACTER UPLO
        INTEGER  INFO, LDA, N
        REAL(dp) A(LDA,*)
        END SUBROUTINE
    
        SUBROUTINE ZPOTRF(UPLO, N, A, LDA, INFO)
        import :: dp
        CHARACTER UPLO
        INTEGER  INFO, LDA, N
        COMPLEX(dp) A(LDA,*)
        END SUBROUTINE
    
        SUBROUTINE ZTRTRI(UPLO, DIAG, N, A, LDA, INFO)
        import :: dp
        CHARACTER UPLO, DIAG
            INTEGER  INFO, LDA, N
    COMPLEX(dp) A(LDA,*)
            END SUBROUTINE
    
    SUBROUTINE ZGELQF(M, N, A, LDA, TAU, WORK, LWORK,  INFO)
            import :: dp
            INTEGER  INFO, LDA, N, M, LWORK
    COMPLEX(dp) A(LDA,*), TAU(*), WORK(*)
            END SUBROUTINE
    
    SUBROUTINE ZUNGQL(M, N, K, A, LDA, TAU, WORK, LWORK,  INFO)
            import :: dp
            INTEGER  INFO, LDA, K, N, M, LWORK
    COMPLEX(dp) A(LDA,*), TAU(*), WORK(*)
            END SUBROUTINE
    
    SUBROUTINE ZGELQ2(M, N, A, LDA, TAU, WORK, INFO)
            import :: dp
            INTEGER  INFO, LDA, N, M 
        COMPLEX(dp) A(LDA,*), TAU(*), WORK(*)
        END SUBROUTINE
    
        SUBROUTINE ZUNML2(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO)
        import :: dp
        CHARACTER SIDE, TRANS
        INTEGER  INFO, LDA, N, M, K, LDC
        COMPLEX(dp) A(LDA,*), TAU(*), WORK(*), C(LDC,*)
        END SUBROUTINE
    
        SUBROUTINE ZTRCON(NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, RWORK, INFO)
        import :: dp
        CHARACTER NORM, UPLO, DIAG
        INTEGER  INFO, LDA, N
        REAL(dp) RWORK(N), RCOND
        COMPLEX(dp) A(LDA,*), WORK(2*N)
        END SUBROUTINE
    
          SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        import :: dp
          COMPLEX(dp) ALPHA
          INTEGER LDA,LDB,M,N
          CHARACTER DIAG,SIDE,TRANSA,UPLO
          COMPLEX(dp) A(LDA,*),B(LDB,*)
        END SUBROUTINE
          SUBROUTINE ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        import :: dp
          complex(dp) ALPHA,BETA
          INTEGER K,LDA,LDC,N
          CHARACTER TRANS,UPLO
          complex(dp) A(LDA,*),C(LDC,*)
        END SUBROUTINE
          SUBROUTINE ZTRTTP( UPLO, N, A, LDA, AP, INFO )
        import :: dp
             CHARACTER          UPLO
             INTEGER            INFO, N, LDA
            COMPLEX(dp)         A( LDA, * ), AP( * )
        END SUBROUTINE
          SUBROUTINE ZTPTTR( UPLO, N, AP, A, LDA, INFO )
        import :: dp
             CHARACTER          UPLO
             INTEGER            INFO, N, LDA
            COMPLEX(dp)         A( LDA, * ), AP( * )
        END SUBROUTINE
          SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        import :: dp
          COMPLEX(dp) ALPHA,BETA
          INTEGER INCX,INCY,LDA,M,N
          CHARACTER TRANS
          COMPLEX(dp) A(LDA,*),X(*),Y(*)
        END SUBROUTINE
          SUBROUTINE ZHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, &
                             LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, &
                             IFAILR, INFO )
          import :: dp
          CHARACTER          EIGSRC, INITV, SIDE
          INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N
          LOGICAL            SELECT( * )
          INTEGER            IFAILL( * ), IFAILR( * )
          DOUBLE PRECISION   RWORK( * )
          COMPLEX(dp)         H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), &
                              W( * ), WORK( * )
        END SUBROUTINE

        SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ, &
                             WORK, LWORK, INFO )
           import :: dp
           INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
           CHARACTER          COMPZ, JOB
           COMPLEX(dp)         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
          END SUBROUTINE

    end interface
    
    contains
    
end module
    