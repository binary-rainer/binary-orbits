;
; lmatrix.pro
; Created:     Wed Jun 20 16:25:32 2007 by Rkoehler@lx40
; Last change: Wed Jul  2 15:27:51 2014
;
; PURPOSE:
;	run lmorbitw on many orbits in matrix from svd_hr_T0
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
PRO lmatrix, inmat, matrix, materr, ObsDat=obsdat, Fitel=fitel, FILE=file, STEP=step, DCHIMAX=dchimax

  if not keyword_set(fitel) then fitel= replicate(1,7)	;; fit all 7 elements
  if not keyword_set(file)  then file="lmatrix.sav"
  if not keyword_set(step)  then step=10
  if not keyword_set(dchimax) then dchimax=9.

  sz = long( size(inmat,/dim) / step)

  ;; new lmorbitw returns a bigger struct (incl. No. of iterations)
  ;;matrix = replicate(inmat[0], sz[0], sz[1])
  ;;materr = matrix

  t0 = systime(/seconds)
  allpnts= sz[0]*sz[1]

  chimin = min(inmat.error)
  alloced= 0
  LMcount= 0

  print,"Min error",chimin
  print,"=> max error to start lmfit:",chimin+dchimax
  print,"Fit elements:",fitel
  junk='' & read,prompt="Hit return",junk

  for i1=0,sz[0]-1 do begin
      for i2=0,sz[1]-1 do begin
          print,"=====",i1,i2," =================================="
          el = inmat[ i1*step, i2*step ]

          if el.error le chimin+dchimax then begin
             lmorbitw, el, del, ObsDat=obsdat, Fitel=fitel
             LMcount++

             if alloced eq 0 then begin
                matrix = replicate( el, sz[0], sz[1])
                materr = replicate(del, sz[0], sz[1])
                alloced= 1
             endif
             matrix[i1,i2]= el
             materr[i1,i2]= del
          endif
          donepts= i1*sz[1]+i2+1
          print,"Finished ",donepts," of ",allpnts," points"
          t= (systime(/seconds) - t0) * (allpnts-donepts)/donepts/60
          if t lt 100 then print,"Done in",t," minutes" $
            else           print,"Done in",t/60," hours"
      endfor
      print,"Saving results in ",file
      save,matrix,materr,obsdat,file=file
  endfor
  print,"Ran LM-fit",LMcount," times"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
