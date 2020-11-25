;
; lmfinderr.pro
; Created:     Tue Jun 19 13:57:30 2007 by Rkoehler@lx40
; Last change: Fri Jun 12 16:52:40 2015
;
; PURPOSE:
;	find errors of orbital elements
;
; 27-May-2015: added option scaleChiSq to lmfinderr (adapted from TWA5/lmorbit.pro)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Function lmfinderr1, prefix, indvar, depvar, deperr, chisq0, inelem, dinpar, ifix,$
                     MASS=mass, PARARR=sigarr, CHIARR=chiarr, CNVARR=cnvarr

  ;print,"Chisq^2_min =",chisq0
  sigarr= dblarr(100)
  chiarr= dblarr(100)
  cnvarr= intarr(100)	;; convergenve flag

  chiarr[0]= chisq0
  cnvarr[0]= 1		;; assume the main fit converged

  fita= replicate(1,7)
  fita[ifix]= 0		;; don't fit param ifix

  ;if abs(dinpar) gt 5 then sigstp= round(dinpar/10.) else sigstp= dinpar
  sigstp= dinpar

  ;print,"sigma step: ",sigstp

  for try=0,9 do begin
     param= double([ inelem.periJD, $
                     inelem.period, $
                     inelem.axis,   $
                     inelem.excent, $
                     inelem.periast,$
                     inelem.ascnode,$
                     inelem.inclin  ])

     sigarr[0]= 0.	;;param[ifix]

     ;print,'LMfinderr1:',param

     isig= 0
     ;;print,"Try",try,", isig",isig
     repeat begin
        mass = (param[2]*param[2]*param[2]) / (param[1]*param[1])
        ;;print,"Iteration",isig,", mass",mass
        ;; the returned mass should be the last mass *before* chisq goes to heaven
        isig++
        ;;param= double([ inelem.periJD, $
        ;;                inelem.period, $
        ;;                inelem.axis,   $
        ;;                inelem.excent, $
        ;;                inelem.periast,$
        ;;                inelem.ascnode,$
        ;;                inelem.inclin  ])
        ;;param[ifix] += double(isig*sigstp)  ;;nsig*dinpar
        param[ifix] += double(sigstp)
        ;;sigarr[isig]= param[ifix]
        sigarr[isig]= sigarr[isig-1]+double(sigstp)

        res = LMFit( indvar, depvar, param, MEASURE_Errors=deperr, FITA=fita,$
                     ITMAX=444, Function_name='LMorbwfunc', /double,$
                     CHISQ=chisq, CONVERGENCE=convflg, ITER=iter, SIGMA=sigma, COVAR=covar )

        print,Format='(A,": Sigma #",I3,", chi^2 =",F7.2,", Convergence:",I2," (",I4," iterations)",A)',$
              prefix,isig,chisq,convflg,iter,string(13B)
        chiarr[isig]= chisq
        cnvarr[isig]= convflg
        ;;print,sigarr[0:isig<4]
        plot, sigarr[0:isig], chiarr[0:isig], yrange=[chisq0,chisq0+2]
        noconv= where( cnvarr[0:isig] ne 1, cnt)
        if cnt gt 0 then oplot,sigarr[noconv],chiarr[noconv],psym=4
        wait,1./50.

        ;;;isig++
        ;;;param[ifix] += double(sigstp)
        ;;;sigarr[isig]= sigarr[isig-1]+double(sigstp)

     endrep until (chisq gt chisq0+1) or isig ge 99
     if isig lt 3 then begin
        sigstp /= 10.
        ;;print,''
        print,"Reducing step to",sigstp
     endif else break
  endfor

  sigarr= sigarr[0:isig]
  chiarr= chiarr[0:isig]
  cnvarr= cnvarr[0:isig]

  print,Format='(A,": Sigma #",I3,", chi^2=",F7.2,", Convergence:",I2," (",I3," iter.), confidence int.=",G10.4,$)',$
        prefix,isig,chisq,convflg,iter,abs(sigstp*isig)

  ;junk='' & read,prompt=' ',junk
  print,'' & wait,1.

  return,{ periJD:  param[0],$
           period:  param[1],$
           axis:    param[2],$
           excent:  param[3],$
           periast: param[4],$
           ascnode: param[5],$
           inclin:  param[6],$
           Mass:    param[2]*param[2]*param[2] / param[1]/param[1],$
           error:   chisq    }
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Function lmfindMasserr1, prefix, indvar, depvar, deperr, chisq0, inelem, dinpar, ifix,$
                     MASS=mass, PARARR=sigarr, CHIARR=chiarr

  if ifix ne 2 then print,"Mass should be fixed if you call lmfindMasserr1!"

  ;print,"Chisq^2_min =",chisq0
  sigarr= dblarr(100)
  chiarr= dblarr(100)
  cnvarr= intarr(100)	;; convergenve flag

  chiarr[0]= chisq0
  cnvarr[0]= 1		;; assume the main fit converged

  fita= replicate(1,7)
  fita[ifix]= 0		;; don't fit param ifix

  ;if abs(dinpar) gt 5 then sigstp= round(dinpar/10.) else sigstp= dinpar
  sigstp= dinpar

  print,"LMasserr: P=",inelem.period,", a=",inelem.axis,", M=",inelem.mass,", step=",sigstp
  ;print,"sigma step: ",sigstp

  param= double([ inelem.periJD, $
                  inelem.period, $
                  inelem.Mass,   $	;; hope there is one ;-)
                  inelem.excent, $
                  inelem.periast,$
                  inelem.ascnode,$
                  inelem.inclin  ])

  sigarr[0]= param[ifix]

  for try=0,9 do begin
     isig= 0
     ;;print,"Try",try,", isig",isig
     repeat begin
        mass = (param[2]*param[2]*param[2]) / (param[1]*param[1])
        ;;print,"Iteration",isig,", mass",mass
        ;; the returned mass should be the last mass *before* chisq goes to heaven
        isig++
        param= double([ inelem.periJD, $
                        inelem.period, $
                        inelem.Mass,   $
                        inelem.excent, $
                        inelem.periast,$
                        inelem.ascnode,$
                        inelem.inclin  ])
        param[ifix] += double(isig*sigstp)  ;;nsig*dinpar
        ;;param[ifix] += double(sigstp)
        sigarr[isig]= param[ifix]

        res = LMFit( indvar, depvar, param, MEASURE_Errors=deperr, FITA=fita,$
                     ITMAX=444, Function_name='LMorbMassfunc', /double,$
                     CHISQ=chisq, CONVERGENCE=convflg, ITER=iter, SIGMA=sigma, COVAR=covar )

        ;;a = param[2]*0.1476
        ;;print,"axis:",param[2],", period:",param[1],", mass:",a*a*a/(param[1]*param[1])
        print,Format='(A,": Sigma #",I3,", chi^2 =",F7.2,", Convergence:",I2," (",I4," iterations)",A,$)',$
              prefix,isig,chisq,convflg,iter,string(13B)
        chiarr[isig]= chisq
        cnvarr[isig]= convflg
        ;;print,sigarr[0:isig<4]
        plot, sigarr[0:isig], chiarr[0:isig], yrange=[chisq0,chisq0+2]
        noconv= where( cnvarr[0:isig] ne 1, cnt)
        if cnt gt 0 then oplot,sigarr[noconv],chiarr[noconv],psym=4
        wait,1./50.
     endrep until (chisq gt chisq0+1) or isig ge 99
     if isig lt 3 then begin
        sigstp /= 10.
        print,''
        print,"Reducing step to",sigstp
     endif else break
  endfor

  sigarr= sigarr[0:isig]
  chiarr= chiarr[0:isig]

  print,Format='(A,": Sigma #",I3,", chi^2=",F7.2,", Convergence:",I2," (",I3," iter.), confidence int.=",G10.4,$)',$
        prefix,isig,chisq,convflg,iter,abs(sigstp*isig)

  ;junk='' & read,prompt=' ',junk
  print,'' & wait,1.

  return,{ periJD:  param[0],$
           period:  param[1],$
           axis:   (param[2]*param[1]*param[1])^(1./3.),$
           excent:  param[3],$
           periast: param[4],$
           ascnode: param[5],$
           inclin:  param[6],$
           Mass:    param[2],$
           error:   chisq    }
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_chiscan, parname, xn, chin, xp, chip, cnvn, cnvp

  set_plot,"PS"
  x  = [reverse(xn),xp]
  chi= [reverse(chin),chip]
  plot, x, chi, xtitle=parname, ytitle='chi!u2!n', yst=16

  if N_elements(cnvn) gt 0 and N_elements(cnvp) gt 0 then begin
     cnv= [reverse(cnvn),cnvp]
     noconv= where( cnv ne 1, cnt)
     if cnt gt 0 then oplot, x[noconv], chi[noconv],psym=4
  endif
  set_plot,"X"
end
;;
;; find number of significant digits after decimal point (at most 9)
;;
Function sigdig, err_in

  ndig= 0
  err = abs(err_in)	;; don't mess with caller's variable

  while(err lt 1.1 and ndig lt 9) do begin
     ndig++
     err *= 10.
     ;;print,"sigdig:",ndig,err
  endwhile
  return,ndig
end

PRO print_valplusminus, OUT, txt, val, pos, neg

  ndig= sigdig(abs(pos)<abs(neg))	;; use smaller of the two errors
  if ndig gt 0 then begin
     vfm = string(Format='("F7.",I1)',ndig)
     pfm = string(Format='("F",I1,".",I1)',(ndig+2)<9,ndig)
     nfm = string(Format='("F",I1,".",I1)',(ndig+3)<9,ndig)
  endif else begin
     vfm = "I7"		;; print integers
     pfm = "I4"
     nfm = "I4"
     val = round(val)	;; round correctly (print with I-code does not round!)
     pos = round(pos)
     neg = round(neg)
  endelse
  printf,OUT, txt,val, pos, neg, Format='(A-48,"& $",'+vfm+',"$ & $\,^{+",'+pfm+',"}_{",'+nfm+',"}$\\")'
  ;printf,OUT, "%% ",txt,val, pos, neg, ndig	;; debug smart subroutines...
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Distance and Disterr in kpc
;;
PRO lmfinderr, obs, el, del, Distance, Disterr, OUTFILE=outfile, Comment=comment, SCALECHISQ=scaleChiSq

  if not keyword_set(outfile) then outfile='orbit'

  pnames= ["PeriJD ","Period ","Axis   ","Excent ","Periast","Ascnode","Inclin "]

  if el.ascnode ge 180. then begin
     el.ascnode -= 180.
     el.periast -= 180.
  endif

  print,"Input Elements:" ;;& help, el, /str
  tagnms= tag_names(el)
  for i=0,N_tags(del)-1 do print,format='(A7,":",F10.3," +/-",F7.3)',tagnms[i],el.(i),del.(i)

  openw,OUT,outfile+'.tex',/get_lun
  printf,OUT,'%'
  printf,OUT,'% ',outfile,'.tex'
  printf,OUT,'% created by lmfinderr on ',systime()
  printf,OUT,'%'
  if keyword_set(comment) then printf,OUT,'% '+comment
  printf,OUT,'%'
  printf,OUT,'\begin{table}'
  printf,OUT,'\caption{Parameters of the best orbital solution.}'
  printf,OUT,'\label{',outfile,'Tab}'
  printf,OUT,"%% stretch a bit so that sub- and superscripts don't overlap"
  printf,OUT,'\renewcommand{\arraystretch}{1.3}'
  printf,OUT,'\begin{center}'
  printf,OUT,'\begin{tabular}{lr@{}l}'
  printf,OUT,'\noalign{\vskip1pt\hrule\vskip1pt}'
  printf,OUT,'Orbital Element				& \multicolumn{2}{c}{Value} \\'
  printf,OUT,'\noalign{\vskip1pt\hrule\vskip1pt}'

  N_obs= N_elements(obs)

  print,N_obs," observations => dof =",(2*N_obs-7)
  chisq0= el.error
  chisq_unscaled = chisq0	;; keep this for printed table
  print,"Chisq^2_min =",chisq0

  errscale = 1.		;; do not scale chi^2
  if keyword_set(scaleChiSq) then begin
     errscale = sqrt(chisq0/(2*N_obs-7))
     chisq0 /= errscale*errscale	;; this should be 2*N_obs-7
  endif
  print,"Scaling errors by",errscale

  eSep_scaled = obs.eSep * errscale
  ePA_scaled  = obs.ePA  * errscale

  indvar = [ obs.MJD, -obs.MJD ]
  depvar = [ obs.Sep,  obs.PA  ]
  deperr = [ eSep_scaled, ePA_scaled ]
  ;;
  ;; check chi^2
  ;;
  EAnom = compute_EAnom(obs.MJD, el)
  Elements_to_ThieleInnes, el, X,Y
  mk_orbit, el, X,Y, EAnom, xorb,yorb, xnode,ynode, xexpect,yexpect

  sepexp = sqrt(xexpect*xexpect + yexpect*yexpect)
  PA_exp = atan(xexpect,yexpect) * 180./!pi
  dr = abs(obs.Sep- sepexp) / eSep_scaled
  dp = abs(obs.PA - PA_exp) / ePA_scaled

  mychi= total( dr*dr + dp*dp)
  print,"My Chi^2:",mychi," (with scaled errors - should be dof)"
  print,"My red.Chi^2:",mychi/(2*N_obs-7)," (with scaled errors - should be 1.0)"

  param= double([ el.periJD,$
                  el.period,$
                  el.axis,  $
                  el.excent,$
                  el.periast,$
                  el.ascnode,$
                  el.inclin  ])


  print,"Parameters:",param
  res = LMFit( indvar, depvar, param, MEASURE_Errors=deperr, $
               ITMAX=444, Function_name='LMorbwfunc', /double,$
               CHISQ=chisq, CONVERGENCE=convflg, ITER=iter, SIGMA=sigma, COVAR=covar )

  print,"New Chisq^2 =",chisq
  print,"Parameters:",param

  delarr= replicate({ periJD:  0.d,$
                      period:  0.d,$
                      axis:    0.d,$
                      excent:  0.d,$
                      periast: 0.d,$
                      ascnode: 0.d,$
                      inclin:  0.d,$
                      Mass:    0.d,$
                      error:   0.d  }, 16)

  set_plot,"PS"
  device,file=outfile+'_chiscan.ps'
  set_plot,"X"

  ;;
  ;; PERIASTRON
  ;;
  if el.excent gt 0.001 then begin
      step = ceil(del.periJD/10.)
      print,"Periastron: (LMsigma = ",del.periJD,", step = ",step,")"
      delarr[0]= lmfinderr1("earlier",indvar,depvar,deperr,chisq0, el,-step, 0, Pararr=JDn, Chiarr=Chin, CNVARR=Cnvn)
      delarr[1]= lmfinderr1("  later",indvar,depvar,deperr,chisq0, el, step, 0, Pararr=JDp, Chiarr=Chip, CNVARR=Cnvp)
      plot_chiscan, "Time of Periastron [days]", JDn, Chin, JDp, Chip, Cnvn,Cnvp

      neg= delarr[0].periJD-el.periJD
      pos= delarr[1].periJD-el.periJD

      JD = el.periJD+2400000.5d	;; orbit-fit uses MJD
      printf,OUT, round(JD), pos, neg,$
        FORMAT='("Date of periastron $T_0$			& $",I7,"$ & $\,^{+",I4,"}_{",I5,"}$\\")'
      printf,OUT,JD,Format="('						& (',C(CYI,' ',CMoA,' ',CDI),')\span\\')"
  endif
  ;;
  ;; PERIOD
  ;;
  step = del.period
  if step eq 0. then step= el.period/1000.
  print,''
  print,"Period: (LMsigma = ",del.period,", step =",step,")"
  delarr[2]= lmfinderr1(" faster",indvar,depvar,deperr,chisq0, el,-step, 1,MASS=nmass, Pararr=Pn, Chiarr=Chin, CNVARR=Cnvn)
  delarr[3]= lmfinderr1(" slower",indvar,depvar,deperr,chisq0, el, step, 1,MASS=pmass, Pararr=Pp, Chiarr=Chip, CNVARR=Cnvp)
  plot_chiscan, "Period [years]", Pn, Chin, Pp, Chip, Cnvn,Cnvp

  neg= delarr[2].period-el.period
  pos= delarr[3].period-el.period

  ;ndig = sigdig(pos<neg)	;; use smaller of the two errors
  ;vfmt = string(Format='("F7.",I1)',ndig)
  ;pfmt = string(Format='("F",I1,".",I1)',(ndig+2)<9,ndig)
  ;nfmt = string(Format='("F",I1,".",I1)',(ndig+3)<9,ndig)
  ;printf,OUT, el.period, pos, neg,$
  ;  Format='("Period $P$ (years)				&  & $",'+vfmt+',"$ & $\,^{",'+pfmt+',"}_{",'+nfmt+',"}$\\")'

  print_valplusminus, OUT, "Period $P$ (years)", el.period, pos, neg
  if el.period lt 3 then $
    printf,OUT, el.period*365.25, pos*365.25, neg*365.25,$
    Format='("Period $P$ (days)				& $",F7.2,"$ & $\,^{+",F5.2,"}_{",F6.2,"}$\\")'

  print,nmass,pmass
  ndmass1= (nmass<pmass)-el.mass & print,"neg. Delta Mass =", ndmass1*Distance^3
  pdmass1= (nmass>pmass)-el.mass & print,"pos. Delta Mass =", pdmass1*Distance^3
  ;;
  ;; SEMI-MAJOR AXIS
  ;;
  step = fix(del.axis)
  if step eq 0. then step= del.axis
  print,''
  print,"Axis: (LMsigma = ",del.axis,", step = ",step,")"
  delarr[4]= lmfinderr1("smaller",indvar,depvar,deperr,chisq0, el,-step,2,MASS=nmass, Pararr=an, Chiarr=Chin, CNVARR=Cnvn)
  delarr[5]= lmfinderr1(" larger",indvar,depvar,deperr,chisq0, el, step,2,MASS=pmass, Pararr=ap, Chiarr=Chip, CNVARR=Cnvp)
  plot_chiscan, "Semi-major axis [mas]", an, Chin, ap, Chip, Cnvn,Cnvp

  neg= delarr[4].axis-el.axis
  pos= delarr[5].axis-el.axis
  print,"Semi-major axis:",el.axis," +",pos,neg

  ndig= sigdig(pos<neg)		;; use smaller of the two errors
  vfm = string(Format='("F7.",I1)',ndig)
  pfm = string(Format='("F",I1,".",I1)',(ndig+2)<9,ndig)
  nfm = string(Format='("F",I1,".",I1)',(ndig+3)<9,ndig)

  print_valplusminus, OUT, "Semi-major axis $a$ (mas)", el.axis, pos, neg
  ;;Format='("Semi-major axis $a$ (mas)			& $",'+vfm+',"$ & $\,^{+",'+pfm+',"}_{",'+nfm+',"}$\\")'

  pos=  sqrt( (el.axis*Disterr*el.axis*Disterr) + (pos*Distance*pos*Distance))
  neg= -sqrt( (el.axis*Disterr*el.axis*Disterr) + (neg*Distance*neg*Distance))
  print_valplusminus, OUT, "Semi-major axis $a$ (AU)", el.axis*Distance, pos,neg
  ;;Format='("Semi-major axis $a$ (AU)			& $",F7.2,"$ & $\,^{+",F5.2,"}_{",F6.2,"}$\\")'

  ndmass2= (nmass<pmass)-el.mass & print,"neg. Delta Mass =", ndmass2*Distance^3
  pdmass2= (nmass>pmass)-el.mass & print,"pos. Delta Mass =", pdmass2*Distance^3
  ;;
  ;; ECCENTRICITY
  ;;
  print,''
  print,"Eccentricity: (LMsigma = ",del.excent,", step = 0.01)"
  delarr[6]= lmfinderr1("less",indvar,depvar,deperr,chisq0, el,-0.01, 3, Pararr=en, Chiarr=Chin, CNVARR=Cnvn)
  delarr[7]= lmfinderr1("more",indvar,depvar,deperr,chisq0, el, 0.01, 3, Pararr=ep, Chiarr=Chip, CNVARR=Cnvp)
  plot_chiscan, "Eccentricity", en, Chin, ep, Chip, Cnvn,Cnvp

  neg= delarr[6].excent-el.excent
  pos= delarr[7].excent-el.excent
  print,"Eccentricity:",el.excent," +",pos,neg
  print_valplusminus, OUT, "Eccentricity $e$", el.excent, pos, neg
  ;Format='("Eccentricity $e$				& $",F7.3,"$ & $\,^{+",F5.3,"}_{",F6.3,"}$\\")'
  ;;
  ;; PERIASTRON
  ;;
  if el.excent gt 0.001 then begin
      print,''
      print,"Periastron: (LMsigma = step = ",del.periast,")"
      delarr[8]= lmfinderr1("smaller",indvar,depvar,deperr,chisq0, el,-del.periast, 4, Pararr=on, Chiarr=Chin, CNVARR=Cnvn)
      delarr[9]= lmfinderr1(" larger",indvar,depvar,deperr,chisq0, el, del.periast, 4, Pararr=op, Chiarr=Chip, CNVARR=Cnvp)
      plot_chiscan, "Argument of Periastron [deg]", on, Chin, op, Chip, Cnvn,Cnvp

      neg= delarr[8].periast-el.periast
      pos= delarr[9].periast-el.periast
      print,"Argument of Periastron:",el.periast," +",pos,neg
      print_valplusminus, OUT, "Argument of periastron $\omega$ ($^\circ$)", el.periast, pos,neg
      ;Format='("Argument of periastron $\omega$ ($^\circ$)	& $",F7.2,"$ & $\,^{+",F5.2,"}_{",F6.2,"}$\\")'
  endif
  ;;
  ;; NODE
  ;;
  print,''
  print,"Node: (LMsigma = ",del.ascnode,", step = ",del.ascnode/10.,")"
  delarr[10]= lmfinderr1("smaller",indvar,depvar,deperr,chisq0, el,-del.ascnode/10., 5, Pararr=On, Chiarr=Chin, CNVARR=Cnvn)
  delarr[11]= lmfinderr1(" larger",indvar,depvar,deperr,chisq0, el, del.ascnode/10., 5, Pararr=Op, Chiarr=Chip, CNVARR=Cnvp)
  plot_chiscan, "Position angle of line of nodes [deg]", On, Chin, Op, Chip, Cnvn,Cnvp

  neg= delarr[10].ascnode-el.ascnode
  pos= delarr[11].ascnode-el.ascnode
  print,"Ascending Node:",el.ascnode," +",pos,neg
  print_valplusminus, OUT, "P.A. of ascending node $\Omega$ ($^\circ$)", el.ascnode, pos,neg
  ;Format='("P.A. of ascending node $\Omega$ ($^\circ$)	& $",F7.2,"$ & $\,^{+",F6.3,"}_{",F7.3,"}$\strut\\")'
  ;;
  ;; INCLINATION
  ;;
  print,''
  print,"Inclination: (LMsigma = ",del.inclin,", step = ",del.inclin,")"
  delarr[12]= lmfinderr1("smaller",indvar,depvar,deperr,chisq0, el,-del.inclin, 6, Pararr=in, Chiarr=Chin, CNVARR=Cnvn)
  delarr[13]= lmfinderr1(" larger",indvar,depvar,deperr,chisq0, el, del.inclin, 6, Pararr=ip, Chiarr=Chip, CNVARR=Cnvp)
  plot_chiscan, "Inclination [deg]", in, Chin, ip, Chip, Cnvn,Cnvp

  neg= delarr[12].inclin-el.inclin
  pos= delarr[13].inclin-el.inclin
  print,"Inclination:",el.inclin," +",pos,neg
  print_valplusminus, OUT, "Inclination $i$ ($^\circ$)", el.inclin, pos,neg
  ;Format='("Inclination $i$ ($^\circ$)			& $",F7.2,"$ & $\,^{+",F5.2,"}_{",F6.2,"}$\\")'
  ;;
  ;; MASS
  ;;
  printf,OUT, el.mass,$
    Format='("% System mass $M_S$ ($\rm mas^3/year^2$)		& $",G10.4,"$ & \\")'

  printf,OUT, pdmass1, ndmass1,$
    Format='("% Mass error from scanning $P$			& & $\,^{+",G8.3,"}_{",G9.3,"}$\\")'
  printf,OUT, pdmass2, ndmass2,$
    Format='("% Mass error from scanning $a$			& & $\,^{+",G8.3,"}_{",G9.3,"}$\\")'

  ndmass = (ndmass1 < ndmass2)	;; the minimum is the larger absolute error
  pdmass = (pdmass1 > pdmass2)

  print,''
  print,"Mass: guess from a&P: +",pdmass," -",ndmass
  if pdmass le 0. then begin
     pdmass= el.mass/100.
     print,"      (going to use",pdmass," for positive step)"
  endif
  delarr[14]= lmfindMasserr1("lighter",indvar,depvar,deperr,chisq0, el, ndmass/10., 2, Pararr=mn, Chiarr=Chin)
  delarr[15]= lmfindMasserr1("heavier",indvar,depvar,deperr,chisq0, el, pdmass/10., 2, Pararr=mp, Chiarr=Chip)
  plot_chiscan, "Mass [separation^3/years^2]", mn, Chin, mp, Chip, Cnvn,Cnvp

  ;printf,OUT, el.mass, delarr[15].Mass-el.mass, delarr[14].Mass-el.mass,$
  ;  Format='("$M_S$ as independent par.\ ($\rm mas^3/year^2$) & $",G10.4,"$ & $\,^{+",G8.3,"}_{",G9.3,"}$\\")'

  ;; ignore the errors from P and a, this is the error I believe now
  pos = delarr[15].Mass-el.mass
  neg = delarr[14].Mass-el.mass
  print_valplusminus,OUT,"System mass $M_S$ ($\rm mas^3/year^2$)", el.mass, pos, neg

  MDerr = 3.*el.mass * Distance*Distance * Disterr	;; mass error caused by D-error
  print,"Mass error because of Distance error:",MDerr
  printf,OUT, MDerr,$
    Format='("Mass error from distance error ($M_\odot$)	& $	$ & $\pm ",F7.3,"$\\")'

  neg *= Distance^3
  pos *= Distance^3

  pdmass=  sqrt(pos*pos + MDerr*MDerr) & print,"pos mass err",pdmass
  ndmass= -sqrt(neg*neg + MDerr*MDerr) & print,"neg mass err",ndmass

  Mass = el.mass * Distance^3	;; in Msun
  print_valplusminus, OUT, "System mass $M_S$ ($M_\odot$)", Mass, pdmass, ndmass
  ;Format='("System mass $M_S$ ($M_\odot$)			& $",F7.4,"$ & $\,^{+",F6.4,"}_{",F7.4,"}$\\")'

  printf,OUT,Format='("%%$\chi^2$	       				& $",F7.3,"$	\\")',chisq_unscaled
  printf,OUT,Format='("reduced $\chi^2$				& $",F7.1,"$\\")',chisq_unscaled/(2*N_obs-7)
  ;;3.9$ &     		&  & $    1.4$ &	\\

  printf,OUT,'\noalign{\vskip1pt\hrule\vskip1pt}'
  printf,OUT,'\end{tabular}'
  printf,OUT,'\end{center}'
  printf,OUT,'\end{table}'

  EAnom = compute_EAnom(obs.MJD, el)
  Elements_to_ThieleInnes, el, X,Y
  mk_orbit, el, X,Y, EAnom, xorb,yorb, xnode,ynode, xexpect,yexpect

  sepexp = sqrt(xexpect*xexpect + yexpect*yexpect)
  PA_exp = atan(xexpect,yexpect) * 180./!pi
  dr = abs(obs.Sep- sepexp) / obs.eSep
  dp = abs(obs.PA - PA_exp) / obs.ePA
  ;; we want to print the "real" deviations, not those with scaled errors!

  Epoch = MJD2Epoch(obs.MJD)

  printf,OUT,"% Julian date  Observed______________  Expected______________  Sigma_____ Ref"
  for i=0,N_obs-1 do $
    printf,OUT, FORMAT='(%"%% %10.3f:  %6.2f mas %7.2f deg, %6.2f mas %7.2f deg, %4.2f %5.2f %s")',$
           Epoch[i], obs[i].Sep, obs[i].PA, sepexp[i], PA_exp[i], dr[i], dp[i], obs[i].Ref

  close,OUT & free_lun,OUT

  set_plot,'PS' & device,/close & set_plot,'X'

  save,el,delarr,Mass,pdmass,ndmass,chisq0,obs,file=outfile+'_err.sav'
end
