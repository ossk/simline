      subroutine printerror(status)

C     Print out the FITSIO error messages to the user

      integer status
      character errtext*30,errmessage*80

C     check if status is OK (no error); if so, simply return
      if (status .le. 0)return

C     get the text string which describes the error
 1    call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

C     read and print out all the error messages on the FITSIO stack
 2    call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end

      subroutine deletefile(filename,status)

C     A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

C     simply return if status is greater than zero
      if (status .gt. 0)return

C     Get an unused Logical Unit Number to use to open the FITS file
 1    call ftgiou(unit,status)

C     try to open the file, to see if it exists
 2    call ftopen(unit,filename,1,blocksize,status)
      if (status .eq. 0)then
C         file was opened;  so now delete it 
 3        call ftdelt(unit,status)
      else if (status .eq. 103)then
C         file doesn't exist, so just reset status to zero and clear errors
          status=0
 4        call ftcmsg
      else
C         there was some other error opening the file; delete the file anyway
          status=0
 5        call ftcmsg
          call ftdelt(unit,status)
      end if

C     free the unit number for later reuse
 6    call ftfiou(unit, status)
      end

