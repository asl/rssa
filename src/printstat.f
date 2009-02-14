c
c     (C) Rasmus Munk Larsen, Stanford University, 2000
c
      subroutine clearstat
      implicit none
      include 'stat.h'
      nopx = 0
      nreorth = 0
      ndot = 0
      nitref = 0
      nbsvd = 0
      nrestart = 0
      tmvopx = 0
      tgetu0 = 0
      tupdmu = 0
      tupdnu = 0
      tintv = 0
      tlanbpro = 0
      treorth = 0
      treorthu = 0
      treorthv = 0
      telru = 0
      telrv = 0
      tbsvd = 0
      tnorm2 = 0
      tdot = 0
      tlansvd = 0
      nlandim = 0
      nsing = 0
      tritzvec = 0
      trestart = 0
      end

      subroutine printchar(label)
      implicit none
      character*(*) label
      integer nc
      external printchar0
      nc = len(label)
      call printchar0(label, nc)
      end

      subroutine printint(label, data)
      implicit none
      character*(*) label
      integer nc, data
      external printint0
      nc = len(label)
      call printint0(label, nc, data)
      end

      subroutine printdbl(label, data)
      implicit none
      character*(*) label
      integer nc
      real data
      double precision outd
      external printdbl0
      nc = len(label)
      outd = data
      call printdbl0(label, nc, outd)
      end

      subroutine printstat
      implicit none

      include 'stat.h'

      call printchar(' +-------------------------------------------------
     &-----------+')
      call printint(' Dimension of Lanczos basis                  = ',
     &     nlandim)
      call printint(' Number of singular values requested         = ',
     &     nsing)
      call printint(' Number of restarts                          = ',
     &     nrestart)
      call printint(' Number of matrix-vector multiplications     = ',
     &     nopx)
      call printint(' Number of reorthogonalizations              = ',
     &     nreorth)
      call printint(' Number of inner products in reorth.         = ',
     &     ndot)
c      print *,'Number of iterative refinement steps        = ',nitref
      call printint(' Number of bidiagonal SVDs calculated        = ',
     &     nbsvd)
      call printchar('')
      call printchar('')

      call printdbl('  Time spent doing matrix-vector multiply    = ',
     &     tmvopx)
      call printdbl('  Time spent generating starting vectors     = ',
     &     tgetu0)
      call printdbl('    Time spent reorthogonalizing U_{j+1}     = ',
     &     treorthu)
      call printdbl('    Time spent reorthogonalizing V_{j}       = ',
     &     treorthv)
      call printdbl('  Time spent reorthogonalizing               = ',
     &     treorth)
      call printdbl(' Total Time spent in LANBPRO                 = ',
     &     tlanbpro)

c      print *
c      print *,'Time spent updating mu-recurrence           = ',tupdmu
c      print *,'Time spent updating nu-recurrence           = ',tupdnu
c      print *,'Time spent on local reorth. on U_{j+1}      = ',telru
c      print *,'Time spent on local reorth. on V_{j+1}      = ',telrv
c      print *,'Time spent in PDNORM2                       = ',tnorm2
c      print *,'Time spent in PDDOT                         = ',tdot
      call printchar('')
      call printdbl('  Time spent in LANBPRO                      = ',
     &     tlanbpro)
      call printdbl('  Time spent computing bidiagonal SVDs       = ',
     &     tbsvd)
      call printdbl('  Time spent doing implicit restarts         = ',
     &     trestart)
      call printdbl('  Time spent computing Ritz vectors          = ',
     &     tritzvec)
      call printchar('')
      call printdbl(' Total Time spent in LANSVD                  = ',
     &     tlansvd)
      call printchar(' +-------------------------------------------------
     &----------+')
      end
