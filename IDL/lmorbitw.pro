;
; lmorbitw.pro
; Created:     Tue Jun 19 13:57:30 2007 by Rkoehler@lx40
; Last change: Wed Jul  2 13:53:03 2014
;
; PURPOSE:
;	save the world (using Levenberg-Marquart)
;
; This version knows derivatives for hyperbolic orbits!
;
;
Function LMorbpoint, MJD, el,X,Y, DERIV=deriv,$
                     DT0=dT0, DPERIOD=DPer, DECC=decc, DPERIAST=dpast, DNODE=dnode, DINC=dinc

  el.excent = abs(el.excent)
  if el.excent eq 1.0 then el.excent= 0.999	;; yes, we cheat

  if el.excent lt 1. then begin
      EAnom= elliptic_EAnom( abs(MJD), el)

      sinE = sin(EAnom)
      cosE = cos(EAnom)
      cEfac= 1 - el.excent*cosE
      sqre2= sqrt(1.-el.excent*el.excent)
      rcnu = cosE - el.excent
      rsnu = sqre2 * sinE
  endif else if el.excent eq 1.0 then begin
      ;; EAnom is really true anom, simplify mk_orbit-call later
      EAnom= parabolic_Anom( abs(MJD), el, radius)
      rcnu = radius * cos(EAnom)
      rsnu = radius * sin(EAnom)
  endif else begin
      EAnom= hyper_EAnom( abs(MJD), el)

      cosE = cosh(EAnom)	;; equivalent to elliptic orbit => same var names
      sinE = sinh(EAnom)
      cEfac= el.excent*cosE - 1.
      sqre2= sqrt(el.excent*el.excent - 1.)
      rcnu = el.excent - cosE
      rsnu = sqre2 * sinE
      ;; Montenbruck p.66, Formeln fuer r*cos(v) und r*sin(v)
  endelse

  xexpect= X[0] * rcnu + X[1] * rsnu
  yexpect= Y[0] * rcnu + Y[1] * rsnu

  if keyword_set(deriv) then begin
     if el.excent eq 1. then $
        print,"Derivatives only supported for elliptic and hyperbolic orbits! I got",el.excent

     comf = 2.*!dpi / cEfac / el.period / 365.25	;; common factor
     drc_dT =  sinE * comf
     drs_dT = -sqre2*cosE * comf
     dx_dT0 = drc_dT * X[0] + drs_dT * X[1]
     dy_dT0 = drc_dT * Y[0] + drs_dT * Y[1]

     M = (abs(MJD) - el.periJD) / (el.period*365.25) ;mean anomaly in revolutions
     M = (M-floor(M))* 2.*!pi                        ;in radian in [0,2pi]
     comf = M / cEfac / el.period
     drc_dP=  sinE * comf
     drs_dP= -sqre2*cosE * comf
     dx_dP = drc_dP * X[0] + drs_dP * X[1]
     dy_dP = drc_dP * Y[0] + drs_dP * Y[1]

     if el.excent lt 1. then begin
        drc_de = -sinE*sinE / cEfac - 1.
        drs_de = sinE * rcnu / sqre2 / cEfac
     endif else begin
        drc_de = 1. + sinE*sinE / cEfac
        drs_de = -sinE * rcnu / sqre2 / cEfac
     endelse
     dx_de = drc_de * X[0] + drs_de * X[1]
     dy_de = drc_de * Y[0] + drs_de * Y[1]

     peri = el.periast * !dtor	;; dtor = !pi/180.
     node = el.ascnode * !dtor
     inc  = el.inclin  * !dtor

     sosO = sin(peri)*sin(node)  &  socO = sin(peri)*cos(node)
     cosO = cos(peri)*sin(node)  &  cocO = cos(peri)*cos(node)
     cosi = cos(inc)
     sini = sin(inc)

     dx_dperi= -el.axis * (rcnu*(sosO-cocO*cosi) + rsnu*(cosO+socO*cosi))
     dy_dperi= -el.axis * (rcnu*(socO+cosO*cosi) + rsnu*(cocO-sosO*cosi))

     dx_dnode=  el.axis * (rcnu*(cocO-sosO*cosi) - rsnu*(socO+cosO*cosi))
     dy_dnode= -el.axis * (rcnu*(cosO+socO*cosi) - rsnu*(sosO-cocO*cosi))

     dx_di= -el.axis * sini * (rcnu*socO + rsnu*cocO)
     dy_di=  el.axis * sini * (rcnu*sosO + rsnu*cosO)

  endif

  res = dblarr(8)
  radius = sqrt(xexpect*xexpect + yexpect*yexpect)

  if MJD gt 0. then begin
      res[0]= radius
      if keyword_set(deriv) then begin
          ;; !radeg = 180./!pi => convert from rad to deg
          res[1]= (xexpect*dx_dT0   + yexpect*dy_dT0  ) / radius
          res[2]= (xexpect*dx_dP    + yexpect*dy_dP   ) / radius
          res[3]= radius / el.axis
          res[4]= (xexpect*dx_de    + yexpect*dy_de   ) / radius
          res[5]= (xexpect*dx_dperi + yexpect*dy_dperi) / radius / !radeg
          res[6]= (xexpect*dx_dnode + yexpect*dy_dnode) / radius / !radeg
          res[7]= (xexpect*dx_di    + yexpect*dy_di   ) / radius / !radeg
      endif
      return, res
  endif else begin
      res[0]= atan(xexpect,yexpect) * 180./!pi
      ;;order is correct, we swapped x and y compared to usual conventions
      if keyword_set(deriv) then begin
         r2= radius*radius
         res[1]= (yexpect*dx_dT0   - xexpect*dy_dT0  ) / r2 * !radeg
         res[2]= (yexpect*dx_dP    - xexpect*dy_dP   ) / r2 * !radeg
         res[3]= 0.
         res[4]= (yexpect*dx_de    - xexpect*dy_de   ) / r2 * !radeg
         res[5]= (yexpect*dx_dperi - xexpect*dy_dperi) / r2
         res[6]= (yexpect*dx_dnode - xexpect*dy_dnode) / r2
         res[7]= (yexpect*dx_di    - xexpect*dy_di   ) / r2
         ;; dPA/di is dimensionless, NOT in radian!
      endif
      return, res
   endelse
end


Function LMorbwfunc, MJD, param
  common MyLM,count

  if N_elements(param) ne 7 then begin
      print,"LMorbwfunc: Unexpected number of elements!!"
      return,0
  endif

  ;print,"Count:",count
  ;print,"param:",param
  ;count++

  el = { periJD : param[0],$
	 period : param[1],$
	 axis   : param[2],$
         excent : param[3],$
         periast: param[4],$
         ascnode: param[5],$
         inclin : param[6],$
         mass   : 0., $
         error  : 1e9 }

  Elements_to_ThieleInnes, el,X,Y

  res = LMorbpoint(MJD, el,X,Y, /Deriv)

  ;;if MJD gt 0. then print,"LMfunc:",MJD,res[5],numperi,res[6],numnode;;res[7],numinc
  return, res
end

;;
;; like above, but use Mass instead of axis as free parameter
;;
Function LMorbMassfunc, MJD, param
  common MyLM,count

  if N_elements(param) ne 7 then begin
      print,"LMorbMassfunc: Unexpected number of elements!!"
      return,0
  endif

  ;print,"Count:",count
  ;print,"param:",param
  ;count++

  el = { periJD : param[0],$
	 period : param[1],$
	 axis   : (param[2]*param[1]*param[1])^(1./3.),$	;; a^3 = M P^2
         excent : param[3],$
         periast: param[4],$
         ascnode: param[5],$
         inclin : param[6],$
         mass   : 0., $
         error  : 1e9 }

  ;;print,"LMass: period ",param[1],", mass ",param[2],", axis =",el.axis

  Elements_to_ThieleInnes, el,X,Y

  res = LMorbpoint(MJD, el,X,Y, /Deriv)


  if MJD gt 0. then begin
     res[3] = res[0] / 3./param[2]	;; r = r0 * a = r0 * (MP^2)^(1/3)
					;; dr/dM = r0 * a/(3M) = r/(3M)
  endif else begin
     res[3] = 0.		;; P.A. does not depend on axis
  endelse


  ;;if MJD gt 0. then print,"LMfunc:",MJD,res[5],numperi,res[6],numnode;;res[7],numinc
  return, res
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Fitel is an array with zero for fixed elem,
;;	and not zero for an element to be fit (passed as FITA to lmfit)

PRO lmorbitw, el, del, Fitel=Fitel, ObsDat=obs

  print,"Input Elements:"
  help, el, /str

  ;if not keyword_set(obs) then obs= observations()	;; new style measurements fct.
  N_obs= N_elements(obs)

  indvar = [ obs.MJD, -obs.MJD ]	;; untested without further parameters :-)
  depvar = [ obs.Sep,  obs.PA  ]
  deperr = [ obs.eSep, obs.ePA ]


  Elements_to_ThieleInnes, el, X, Y

  param= double([ el.periJD,$
                  el.period,$
                  el.axis,  $
                  el.excent,$
                  el.periast,$
                  el.ascnode,$
                  el.inclin  ])

  pnames= ["PeriJD ","Period ","Axis   ","Excent ","Periast","Ascnode","Inclin "]


  bestchisq = 1d99
  fiddle_with_sigmas= 0	;; don't waste time

  ;;for outer=1,2 do begin
  ;;   print
  ;;   print,"===== Outer iteration",outer," ====="
     res = LMFit( indvar, depvar, param, MEASURE_Errors=deperr, FITA=Fitel,$
               ITMAX=111, Function_name='LMorbwfunc', /double,$
               CHISQ=chisq, CONVERGENCE=convflg, ITER=iter, SIGMA=sigma, COVAR=covar )

     for i=0,6 do print,pnames[i],":",param[i]," +/-",sigma[i]

     print,"Chi^2:",chisq
     print,"Convergence:",convflg," (",iter," iterations)"

     ;;if fiddle_with_sigmas then begin
     ;;   bestparam = param
     ;;   bestsigma = sigma
     ;;   bestchisq = chisq
     ;;   bestcovar = bestcovar
     ;;
     ;;   step = sigma / 20.
     ;;   step[0] /= 5.	;; step in periJD
     ;;   ;;step[1] *= 10.	;; step in period
     ;;
     ;;   for pm=0,13 do begin
     ;;      p = pm mod 7
     ;;      newparam= param
     ;;      if pm le 6 then $
     ;;         newparam[p] -= step[p] $
     ;;      else newparam[p] += step[p]
     ;;
     ;;      ;;print,Format='(I2," ",A7,":",F14.5)', pm,pnames[p], newparam[p]
     ;;
     ;;      res = LMFit( indvar, depvar, newparam, MEASURE_Errors=deperr,$
     ;;                   ITMAX=42, Function_name='LMorbwfunc', /double,$
     ;;                   CHISQ=newchisq, CONVERGENCE=newconvflg, ITER=iter, SIGMA=newsigma, COVAR=newcovar )
     ;;
     ;;      print,Format='(I2," ",A7,":",F14.5," +/- ",F10.5,", Chi^2:",F11.5,", Convergence:",I2,",",I4," iteratons")',$
     ;;            pm,pnames[p], newparam[p], newsigma[p], newchisq, newconvflg, iter
     ;;
     ;;      if newchisq lt bestchisq then begin
     ;;         bestparam = newparam
     ;;         bestsigma = newsigma
     ;;         bestchisq = newchisq
     ;;         bestcovar = newcovar
     ;;      endif
     ;;   endfor
     ;;
     ;;   print,"best chisq:",bestchisq
     ;;
     ;;   if bestchisq gt chisq-0.01 then break
     ;;
     ;;   param = bestparam
     ;;   sigma = bestsigma
     ;;endif
  ;;endfor


  ;;;print,"Covar P-a:",covar[1,2]
  ;;
  ;; find the mass error
  ;;
  ;; M = a^3/P^2
  ;;mass = param[2]*param[2]*param[2] / (param[1]*param[1])
  P = param[1]
  a = param[2]
  mass = a*a*a / (P*P)

  print,"Nick's formula:"
  dma = 3.*a*a   * sigma[2] / (p*p)
  dmp = 2.*a*a*a * sigma[1] / (p*p*p)
  adp = a/p
  dmap= 12.* adp*adp*adp*adp*adp * covar[1,2]
  dmass= sqrt( dma*dma + dmp*dmp - dmap )
  ;;print,"dma =",dma,", dmp =",dmp,", dmap =",sqrt(dmap)
  print,"Mass =",mass," +/-",dmass," (mas^3/yr^2)"


  el = { LMorbelements, $
         periJD : param[0],$
	 period : param[1],$
	 axis   : param[2],$
         excent : param[3],$
         periast: param[4],$
         ascnode: param[5],$
         inclin : param[6],$
         mass   : mass, $
         error  : chisq,$
         iter   : iter $
       }

  del= { orbelements, $
         periJD : sigma[0],$
	 period : sigma[1],$
	 axis   : sigma[2],$
         excent : sigma[3],$
         periast: sigma[4],$
         ascnode: sigma[5],$
         inclin : sigma[6],$
         mass   : dmass, $
         error  : 0. }

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
