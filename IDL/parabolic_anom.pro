;
; parabolic_anom.pro
; Created:     Wed Jul  2 11:39:16 2014 by Koehler@Quorra
; Last change: Wed Jul  2 11:40:07 2014
;
; PURPOSE:
;	compute anomaly for parabolic orbit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function parabolic_Anom, JDobs, el, r

  ;; Define a period for a parabolic orbit:
  ;; Mean Anom = 2pi/period * (t-t0) = sqrt(GM/2/q^3) (t-t0)
  ;; => period = 2pi * sqrt(2q^3/G/M)

  ;; in solar units [Msun,yr,AU]:  G = 4pi^2
  ;; => M = 2q^3 / period^2

  ;; note that the orbit is not periodic, A goes from -inf to +inf

  A = 1.5 * (JDobs - el.periJD) * 2.*!pi / (el.period*365.25)
  hlp= sqrt(A*A+1.)
  tanVh = (hlp+A)^(1./3.) - (hlp-A)^(1./3.)	;; = tan(nu/2)
  V = 2. * atan(tanVh)
  r = 1. + tanVh*tanVh
  return,V
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
