;
; thieleinnes_to_elements.pro
; Created:     Thu Feb 27 11:38:40 2014 by Koehler@Quorra
; Last change: Wed Aug 10 13:19:06 2016
;
; PURPOSE:
;	compute elements struct from Thiele-Innes-constants
;
; 10-aug-2016: added math to make sure 0 < el.ascnode < 180deg
;
pro ThieleInnes_to_Elements, X,Y, el

  opO = atan( X[0]-Y[1], X[1]+Y[0] )	 	; omega + Omega
  omO = atan(-X[0]-Y[1],-X[1]+Y[0] )		; omega - Omega

  el.periast = (opO + omO) * 180./2./!pi ; omega = argument of periastron
  el.ascnode = (opO - omO) * 180./2./!pi ; Omega = P.A. of asc. node

  if  el.ascnode lt 0.  then begin
     el.ascnode = el.ascnode + 180.
     el.periast = el.periast - 180.
     ;; swap both by 180 deg, to get 0 < Omega < 180.
     ;; (Hilditch p.52: "convention demands...")
  endif
  if  el.periast lt 0.  then el.periast = el.periast + 360.

  ;;
  ;; avoid singularity for sin(opO) = 0
  ;;
  if abs(sin(opO)*sin(omO)) gt abs(cos(opO)*cos(omO)) then begin

     help= -((X[0]+Y[1])*sin(opO)) / ((X[0]-Y[1]) * sin(omO))
     cosinc= (1.-help) / (1.+help)
     el.axis  = (X[0]-Y[1]) / sin(opO) / (1.+cosinc)

  endif else begin

     help= ((X[1]-Y[0])*cos(opO)) / ((X[1]+Y[0]) * cos(omO))
     cosinc= (1.+help) / (1.-help)
     el.axis  = (X[1]+Y[0]) / cos(opO) / (1.+cosinc)
  endelse

  el.inclin= acos(cosinc) * 180./!pi

  el.mass= el.axis * el.axis * el.axis / (el.period * el.period)
  if el.excent eq 1. then el.mass= el.mass*2.

  ;;print,"axis =",el.axis
  ;;print,"arg of peri =",el.periast
  ;;print,"ascending node =",el.ascnode
  ;;print,"inclination = +-",el.inclin
end
