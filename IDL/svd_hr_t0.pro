;
; svd_HR_T0.pro
;
; Created:     Tue Oct  7 14:46:22 2003 by Koehler@sun47
; Last change: Tue Jan 12 22:35:53 2016
;
; PURPOSE:
;	make everyone happy (with exceptions)
;
; CALLING SEQUENCE:
;	see below
;
;	periJD_start= central MJD0, will scan this +-Period/2.
;	periJD_npnt = number of points in MJD0 per grid step
;	periJD_mstp = minimum step size in MJD0.
;		 The grid search will be repeated with finer grids
;		 until the grid step is smaller than this
;
; INPUT:
;	range of parameters to scan
;
; OUTPUT:
;	orbital parameters that fit almost, but not quite,
;	the observations
;
; COMMENT:
;	This program sucks, you should write a better one
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;
PRO svd_HR_T0, periJD_start, periJD_npnt,periJD_mstp,$
               period_start, period_end, period_step,$
               excent_start, excent_end, excent_step,$
               PLOT=plot, LOG_PERIOD=log_period,$
               OBSDAT=obsdat, PSfile=psfile, RESULT=bestel, MATRIX=matrix

  if N_params() lt 9 then begin
      print,"USAGE: svd_p_and_e, periJD_start, periJD_npnt,periJD_mstp,$"
      print,"                    period_start, period_end, period_step,$"
      print,"                    excent_start, excent_end, excent_step $"
      print,"                 [, PSfile=psfile, /PLOT, Result=result]"
      print,"(periJD = Julian Day - 2400000)"
      return
  endif

  ;; make sure we float
  periJD_start= float(periJD_start)
  periJD_mstp = float(periJD_mstp)
  period_start= float(period_start)
  period_end  = float(period_end)
  period_step = float(period_step)
  excent_start= float(excent_start)
  excent_end  = float(excent_end)
  excent_step = float(excent_step)

  if not keyword_set(obsdat) then measurements,obsdat
  N_obs = N_elements(obsdat)

  Months = [ "Ups", "Jan", "Feb", "Mar", "Apr", "May", "Jun",$
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" ]

  el = { orbelements, $
         periJD : periJD_start, $
	 period : period_start, $
	 axis   : 0., $
         excent : excent_start, $
         periast: 0., $
         ascnode: 0., $
         inclin : 0., $
         mass   : 0., $
         error  : 1e38 }
  ;; make sure max error in plot-block at end of period-loop matches!

  bestel = el

  nJD0 = periJD_npnt	;;was: fix((periJD_end - periJD_start)/periJD_step)+1

  if keyword_set(log_period) then $
    nPer = fix(alog(period_end/period_start)/alog(period_step))+1 $
  else $
    nPer = fix((period_end - period_start)/period_step)+1

  nExc = fix((excent_end - excent_start) / excent_step)+1

  print,nPer," periods,",nJD0," t0s,",nExc," excentricities"

  if nPer * nExc < 1e6 then $
      matrix = replicate(el,nPer,nExc)

  if keyword_set(PLOT) and !D.name eq "X" then begin
      window,0,title="Orbit"
      window,1,title="Matrix",xsize=nPer,ysize=nExc
      wset,0
  endif

  A = fltarr(2,N_obs)
  Ae= fltarr(2,N_obs)

  time0 = systime(/sec)
  ;window,2	;; for chi^2(T0)
  ;wset,0

  for iPer= 0, nPer-1 do begin
      print, "Period",iPer," of",nPer
      if iPer gt 0 then begin
          t= (systime(/seconds) - time0) * (nPer-iPer)/iPer/60
          if t lt 100 then print,"Done in",t," minutes" $
            else           print,"Done in",t/60," hours"
      endif

      if keyword_set(log_period) then $
        period= period_start * period_step ^ iPer $
      else $
        period= period_start + period_step * iPer

      el.period= period

      for iExc= 0, nExc-1 do begin
          excent= excent_start + excent_step * iExc
          el.excent= excent

          periJD_cent= periJD_start		;; center of grid
          ;;periJD_step= period*365.25/periJD_npnt
          periJD_step= period*365.25/10.		;; pretend that 10 grid steps are one period

          MJDs = fltarr(periJD_npnt)
          MJDe = fltarr(periJD_npnt)
          periJD_rcnt=0

          repeat begin
              periJD_step= periJD_step*10./periJD_npnt  ;; scan 10 grid pnts around best

              ;;print,"MJD0: center",periJD_cent,", step =",periJD_step," days
              for iJD0= 0, periJD_npnt-1 do begin
                  el.periJD= periJD_cent + periJD_step * (iJD0-periJD_npnt/2)

                  if excent lt 1. then begin
                      EAnom = elliptic_EAnom(obsdat.MJD, el)

                      A[0,*] = cos(EAnom) - el.excent
                      A[1,*] = sqrt(1. - el.excent*el.excent) * sin(EAnom)
                      ;; Montenbruck p.58, Formeln fuer r*cos(v) und r*sin(v)
                  endif else if excent eq 1.0 then begin
                      ;; EAnom is really true anom, simplify mk_orbit-call later
                      EAnom = parabolic_Anom(obsdat.MJD, el, radius)

                      A[0,*] = radius * cos(EAnom)
                      A[1,*] = radius * sin(EAnom)
                  endif else begin
                      EAnom = hyper_EAnom(obsdat.MJD, el)

                      A[0,*] = el.excent - cosh(EAnom)
                      A[1,*] = sqrt(el.excent*el.excent - 1.) * sinh(EAnom)
                      ;; Montenbruck p.66, Formeln fuer r*cos(v) und r*sin(v)
                  endelse

                  Ae[0,*] = A[0,*] / obsdat.ex
                  Ae[1,*] = A[1,*] / obsdat.ex
                  svdc, Ae, W,U,V ;;& print,"Wx:",W & print,"Ux:",U & print,"Vx:",V
                  X = svsol(U,W,V, obsdat.x/obsdat.ex)

                  Ae[0,*] = A[0,*] / obsdat.ey
                  Ae[1,*] = A[1,*] / obsdat.ey
                  svdc, Ae, W,U,V
                  Y = svsol(U,W,V, obsdat.y/obsdat.ey)

                  mk_orbit, el, X,Y, EAnom, xorb,yorb, xnode,ynode, xexpect,yexpect

                  sepexp = sqrt(xexpect*xexpect + yexpect*yexpect)
                  PA_exp = atan(xexpect,yexpect) * 180./!pi
                  ;;order is correct, we swapped x and y compared to usual conventions

                  dr = (obsdat.Sep- sepexp) / obsdat.eSep
                  dp = (obsdat.PA - PA_exp)
                  ;; 12jan2016: handle both cases where
                  ;; one angle has been shifted by 360 deg
                  idx= where(dp gt  180.,cnt)  &  if cnt gt 0 then dp[idx] -= 360.
                  idx= where(dp lt -180.,cnt)  &  if cnt gt 0 then dp[idx] += 360.
                  dp /= obsdat.ePA
                  el.error = total(dr*dr + dp*dp)

                  MJDs[iJD0] = el.periJD
                  MJDe[iJD0] = el.error

                  ;;dx = (xexpect - obsdat.x) / obsdat.ex
                  ;;dy = (yexpect - obsdat.y) / obsdat.ey
                  ;;el.error = total(dx*dx + dy*dy)

                  ;;print,"JD0=",el.periJD,", P=",el.period,", e=",el.excent,", Err^2: ", el.error
                  ;;print,"dr: ",dr  &  print,"dp: ",dp

                  if el.error lt matrix[iPer,iExc].error then begin
                      ThieleInnes_to_Elements, X,Y, el
                      matrix[iPer,iExc]= el

                      if el.error lt bestel.error then begin
                          ;;print,"X:",X
                          ;;print,"Y:",Y
                          bestel= el
                          bestX = X
                          bestY = Y
                          if keyword_set(plot) then $
                            quick_plot_orbit, obsdat.MJD, obsdat.x,obsdat.y, bestel
                      endif
                  endif
              endfor ;; periJD

              el = matrix[iPer,iExc]	;; best elements at this (p,e)
              periJD_cent= el.periJD

              ;wset,2
              ;if periJD_rcnt eq 0 then begin
              ;    plot,MJDs,MJDe,/ylog,xrange=[50000,60000],yrange=[100,700],xtitle="periJD",ytitle="chi^2"
              ;endif else oplot,MJDs,MJDe,linest=periJD_rcnt
              ;mi = min(matrix[iPer,0:iExc].error, midx)
              ;oplot, [matrix[iPer,midx].periJD], [mi], psym=1
              ;wset,0

              if (iExc mod 10) eq 0 then begin
                  print,Format='("JD0=",F8.1,"+-",F7.1,", P=",G9.4,", e=",F7.4,", Err^2:",F9.4)',$
                    el.periJD, periJD_step, el.period, el.excent, el.error
              endif
              periJD_rcnt++

          endrep until (periJD_step le periJD_mstp)

      endfor ;; excent

      caldat, bestel.periJD+2400000.5, Mon, Day, Year, hour, min, sec
      Mon = Months[Mon]
      print,"BEST:"
      print,Format='(%"	periMJD = %9.1f = %2d-%s-%4d, %2d:%02d:%02d")',$
        bestel.periJD, Day, Mon, Year, hour, min, sec
      print,Format='("	period =",F9.3," years")',bestel.period
      print,Format='("	excent =",F7.4)', bestel.excent
      print,Format='("	mass   =",G14 )', bestel.mass
      print,Format='("	Err^2  =",F14.9)',bestel.error

      if keyword_set(PLOT) then begin
          m = matrix[*,*].error
          mi= min(m[where(m lt 1e38)],max=ma)
          idx= where(m gt ma)
          if idx[0] gt -1 then m[idx] = (mi+ma)/2.
          wset,1
          ;;tvscl,m < ((9*mi+ma)/10.)
          tvscl,m < (mi+2.)
          wset,0
      endif
  endfor ;; period - end of the big loop

  caldat, bestel.periJD+2400000.5, Mon, Day, Year, hour, min, sec
  Mon = Months[Mon]

  print,"Best fit:"
  print,Format='(%"periMJD= %9.1f = %2d-%s-%4d, %2d:%02d:%02d")',$
    bestel.periJD, Day, Mon, Year, hour, min, sec
  print,"period =",bestel.period," yr"
  print,"axis   =",bestel.axis
  print,"excent =",bestel.excent
  print,"periast=",bestel.periast
  print,"ascnode=",bestel.ascnode
  print,"inclin =",bestel.inclin
  print,"Mass   =",bestel.mass," [mas^3/yr^2]"
  print,"error  =",bestel.error

  EAnom = compute_EAnom(obsdat.MJD, bestel)
  Elements_to_ThieleInnes, bestel, X,Y

  mk_orbit, bestel, X,Y, EAnom, xorb,yorb, xnode,ynode, xexpect,yexpect

  if keyword_set(psfile) then begin
      set_plot,"PS"
      device,/landscape,file=psfile
  endif

  quick_plot_orbit, obsdat.MJD, obsdat.x, obsdat.y, bestel

  ;;plot,obsdat.x,obsdat.y,psym=3,xrange=[20,-20],yrange=[-15, 15],ystyle=1
  ;; note reverse x-axis
  oplot, [0],[0], psym=7, symsize=2
  ;
  ;  Plot the error bars (IDL is doof)
  FOR I = 0,N_obs-1 DO begin
      r  = obsdat[i].Sep
      dr = obsdat[i].eSep
      pa = obsdat[i].PA  * !pi/180.
      dpa= obsdat[i].ePA * !pi/180.
      plots, [ (r+dr)*sin(pa), (r-dr)*sin(pa) ], [ (r+dr)*cos(pa), (r-dr)*cos(pa) ]
      plots, [ r *sin(pa+dpa), r *sin(pa-dpa) ], [ r *cos(pa+dpa), r *cos(pa-dpa) ]
  endfor
  ;;oplot,xorb, yorb, linestyle=1
  ;;oplot,xnode,ynode,linestyle=3
  ;;oplot,[0.,xorb[0]], [0.,yorb[0]],linestyle=2	;xorb[0] -> nu=0 -> periastron
  ;;for i=0,N_obs-1 do $
  ;;  oplot, [obsdat[i].x, xexpect[i]], [obsdat[i].y,yexpect[i]]

  sepexp = sqrt(xexpect*xexpect + yexpect*yexpect)
  PA_exp = atan(xexpect,yexpect) * 180./!pi
  dr = abs(obsdat.Sep- sepexp) / obsdat.eSep
  dp = abs(obsdat.PA - PA_exp)
  idx = where(dp ge 360.,cnt)
  if cnt gt 0 then dp[idx] -= 360.
  dp /= obsdat.ePA

  print,"Julian date  Observed____________  Expected____________  Sigma__"
  for i=0,N_obs-1 do $
    print, FORMAT='(%"%10.3f:  %5.1f mas %6.1f deg, %5.1f mas %6.1f deg, %3.1f %4.1f")',$
    obsdat[i].MJD, obsdat[i].Sep, obsdat[i].PA, sepexp[i], PA_exp[i], dr[i], dp[i]

  if keyword_set(psfile) then begin
      device,/close
      set_plot,"X"
  endif
end
