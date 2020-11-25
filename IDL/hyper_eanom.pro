;
; hyper_eanom.pro
; Created:     Wed Jul  2 11:40:15 2014 by Koehler@Quorra
; Last change: Wed Jul  2 11:40:43 2014
;
; PURPOSE:
;	compute eccentric anomaly for hyperbolic orbit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hyper_EAnom, JDobs, el

  ;; definition for the period of a non-periodic orbit:
  ;; M = 2pi/period * (t-t0) = sqrt(GM/a^3) * (t-t0)
  ;; => period = 2pi * sqrt(a^3/G/M)

  ;; in solar units [Msun,yr,AU]:  G = 4pi^2
  ;; => period = sqrt(a^3/M)
  ;; => M = a^3/period^2 (what a surprise :-)

  M = (JDobs - el.periJD) * 2.*!pi/ (el.period*365.25)  ;;pseudo-mean anomaly in radian
  ;;print,"M=",M

  H1= (6.*abs(M))^(1./3.)
  i = where(M lt 0.)
  if i[0] ne -1 then H1[i] = -H1[i]

  for i=0,100 do begin
      H0= H1
      H1= H0 - (el.excent * sinh(H0)-H0-M) / (el.excent * cosh(H0) - 1.)
      if max(abs(H0-H1)) lt 1e-7 then break
  endfor
  return,H1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
