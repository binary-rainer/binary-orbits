;
; print_residuals.pro
; Created:     Wed Jul  2 15:52:12 2014 by Koehler@Quorra
; Last change: Thu Jul  3 10:37:42 2014
;
; PURPOSE:
;	print residuals of orbit fit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro print_residuals, obs, el

  EAnom = compute_EAnom(obs.MJD, el)
  Elements_to_ThieleInnes, el, Xth,Yth

  mk_orbit, el, Xth,Yth, EAnom, xorb,yorb, xnode,ynode, xexpect,yexpect

  sepexp = sqrt(xexpect*xexpect + yexpect*yexpect)
  PA_exp = atan(xexpect,yexpect) * 180./!pi
  dr = abs(obs.Sep- sepexp) / obs.eSep
  dp = abs(obs.PA - PA_exp) / obs.ePA

  Epoch = MJD2Epoch(obs.MJD)

  print,"Julian date  Observed____________  Expected____________  Sigma____ Ref"
  for i=0,N_elements(obs)-1 do $
    print, FORMAT='(%"%10.3f:  %5.1f mas %6.1f deg, %5.1f mas %6.1f deg, %4.1f %4.1f %s")',$
	Epoch[i], obs[i].Sep, obs[i].PA, sepexp[i], PA_exp[i], dr[i], dp[i], obs[i].Ref

  chisq = total( dr*dr + dp*dp)
  rchisq= chisq/(2*N_elements(obs)-7)
  print,"chi^2:",chisq
  print,"red.chi^2:",rchisq
end

pro print_residuals_xy, obs, el

  EAnom = compute_EAnom(obs.MJD, el)
  Elements_to_ThieleInnes, el, Xth,Yth

  mk_orbit, el, Xth,Yth, EAnom, xorb,yorb, xnode,ynode, xexpect,yexpect

  dx = abs(obs.x - xexpect)
  dy = abs(obs.y - yexpect)

  etng = obs.ePA * !pi/180. * obs.sep
  etot = sqrt(obs.esep*obs.esep + etng*etng)
  dtot = sqrt(dx*dx + dy*dy) / etot

  evec = sqrt(obs.ex*obs.ex + obs.ey*obs.ey)	;; length of error vector
  dvec = sqrt(dx*dx + dy*dy) / evec	;; lenght of (o-c)-vector / length of error vector

  dx /= obs.ex
  dy /= obs.ey

  dqua= sqrt(dx*dx + dy*dy)

  print,"Julian date  Observed___________  Expected___________  Sigma____ Ref"
  for i=0,N_elements(obs)-1 do $
    print, FORMAT='(%"%10.3f:  %5.1f mas %5.1f mas, %5.1f mas %5.1f mas, %4.1f %4.1f %s %4.1f %4.1f")',$
	obs[i].MJD, obs[i].x, obs[i].y, xexpect[i], yexpect[i], dx[i], dy[i], obs[i].Ref, dtot[i], dqua[i]

  Ndof = (2*N_elements(obs)-7)

  chisq = total( dx*dx + dy*dy)
  rchisq= chisq/Ndof
  print,"chi^2:",chisq
  print,"red.chi^2:",rchisq

  print,'Sum( sqrt((dx/ex)^2 + (dy/ey)^2)): ',total(dqua),', reduced:',total(dqua)/Ndof

  print,"Chi^2 from 1D-error:",total(dtot),", reduced:",total(dtot)/Ndof

  print,"Chi^2 from vector:  ",total(dvec),", reduced:",total(dvec)/Ndof
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
