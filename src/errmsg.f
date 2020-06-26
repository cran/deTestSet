c ===================================================================================
c print R-messages
c ===================================================================================

      subroutine intpr_k(label, nchar, data, ndata)
      integer nchar, ndata
      character*(*) label
      integer data(ndata)

       if (ndata .eq. 1) then
         call rprintfi1(label//char(0), data(1))
       else if (ndata .eq. 2) then
         call rprintfi2(label//char(0), data(1), data(2))
       else if (ndata .eq. 3) then
         call rprintfi3(label//char(0), data(1), data(2), data(3))
       else if (ndata .ge. 4) then
         call rprintfi4(label//char(0),data(1),data(2),data(3),data(4))
       endif 
      end subroutine


      subroutine dblepr_k(label, nchar, data, ndata)
      integer nchar, ndata
      character*(*) label
      double precision data(ndata)

       if (ndata .eq. 1) then
         call rprintfd1(label//char(0), data(1))
       else if (ndata .eq. 2) then
         call rprintfd2(label//char(0), data(1), data(2))
       else if (ndata .eq. 3) then
         call rprintfd3(label//char(0), data(1), data(2), data(3))
       else if (ndata .ge. 4) then
         call rprintfd4(label//char(0),data(1),data(2),data(3),data(4))
       endif 
      end subroutine

C just a string
      subroutine rprint(msg)
      character (len=*) msg

      call rprintf(msg//char(0))

      end subroutine 


C printing with one integer and a double
      subroutine rprintid(msg, i1, d1)
      character (len=*) msg
      double precision d1
      integer i1
         call rprintfd1(msg//char(0), d1)
         call rprintfi1(' '//char(0), i1)
      end subroutine 


C printing with doubles
C
      subroutine rprintd1(msg, d1)
      character (len=*) msg
      double precision d1
        call rprintfd1(msg//char(0), d1)
      end subroutine 
C
      subroutine rprintd2(msg, d1, d2)
      character (len=*) msg
      double precision d1, d2
        call rprintfd2(msg//char(0), d1, d2)
      end subroutine 
C
      subroutine rprintd3(msg, d1, d2, d3)
      character (len=*) msg
      double precision d1, d2, d3
        call rprintfd3(msg//char(0), d1, d2, d3)
      end subroutine 
C
      subroutine rprintd4(msg, d1, d2, d3, d4)
      character (len=*) msg
      double precision d1, d2, d3, d4
        call rprintfd4(msg//char(0), d1, d2, d3, d4)
      end subroutine 


C printing with integers
      subroutine rprinti1(msg, i1)
      character (len=*) msg
      integer i1
        call rprintfi1(msg//char(0), i1)
      end subroutine 
C
      subroutine rprinti2(msg, i1, i2)
      character (len=*) msg
      INTEGER i1, i2
        call rprintfi2(msg//char(0), i1, i2)
      end subroutine 
C
      subroutine rprinti3(msg, i1, i2, i3)
      character (len=*) msg
      INTEGER i1, i2, i3
        call rprintfi3(msg//char(0), i1, i2, i3)
      end subroutine
