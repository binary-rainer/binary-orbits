;
; elliptic_eanom.pro
; Created:     Wed Jul  2 11:38:09 2014 by Koehler@Quorra
; Last change: Wed Jul  2 11:38:54 2014
;
; PURPOSE:
;	compute eccentric anomaly for elliptic orbit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function elliptic_EAnom, JDobs, el

  M = (JDobs - el.periJD) / (el.period*365.25) ;mean anomaly in revolutions
  M = (M-floor(M))* 2.*!pi  	;in radian in [0,2pi]
  ;;print,"M=",M

  if el.excent eq 0. then return,M

  E1= M
  for i=0,100 do begin
      E0= E1
      E1= M + el.excent * sin(E0)
      if max(abs(E0-E1)) lt 1e-7 then break
  endfor
  ;;print,"E=",E1
  return,E1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
