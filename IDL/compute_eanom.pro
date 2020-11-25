;
; compute_eanom.pro
; Created:     Thu Feb 27 11:20:47 2014 by Koehler@Quorra
; Last change: Wed Jul  2 11:40:46 2014
;
; PURPOSE:
;	compute eccentric anomaly
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function compute_EAnom, JDobs, el, radius
  ;; radius is used only with parabolic orbs

  if el.excent lt 1. then begin
      return, elliptic_EAnom(JDobs,el)
  endif else if el.excent eq 1. then begin
      return, parabolic_Anom(JDobs,el,radius)
  endif else begin
      return, hyper_EAnom(JDobs,el)
  endelse
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
