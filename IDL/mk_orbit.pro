;
; mk_orbit.pro
; Created:     Thu Feb 27 11:40:54 2014 by Koehler@Quorra
; Last change: Thu Mar 28 11:02:16 2019
;
; PURPOSE:
;	compute coordinates for plotting an orbit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function rad_from_cnu, excent, cnu
;; cnu = cos(true Anomaly)

  if excent ne 1. then begin
      return, abs(1. - excent * excent) / (1. + excent * cnu)
      ;; Montenbruck p. 54/63 and 124, a = 1 with Thiele-Innes
  endif else if excent eq 1. then begin
      return, 2. / (1. + cnu)
  endif
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mk_orbit, el, X,Y, E_obs, xorb,yorb, xnode,ynode, xexpect,yexpect

  ;;print,"mkorbit:",X,Y
  ;
  ; we want to plot the shape of the orbit only (we do not care about the times),
  ; so we loop over the true anomaly
  ;
  nu = findgen(361) * !pi/180.

  cnu= cos(nu)
  snu= sin(nu)
  r  = rad_from_cnu(el.excent, cnu)

  good= where(r gt 0.)
  r   = r[good]
  cnu = cnu[good]
  snu = snu[good]

  xorb = r * (X[0] * cnu + X[1] * snu) ;# in direction PA=90
  yorb = r * (Y[0] * cnu + Y[1] * snu) ;# in direction PA=0

  xnode = fltarr(2) & ynode = fltarr(2)

  nu = -el.periast * !pi/180.	; periast = angle (node -> periastron)
  cnu= cos(nu)
  snu= sin(nu)
  r  = abs(rad_from_cnu(el.excent, cnu))
  xnode[0] = r * (X[0] * cnu + X[1] * snu)
  ynode[0] = r * (Y[0] * cnu + Y[1] * snu)

  ; nu + 180deg => cos -> -cos, sin -> -sin
  r  = abs(rad_from_cnu(el.excent,-cnu))
  xnode[1] = -r * (X[0] * cnu + X[1] * snu)
  ynode[1] = -r * (Y[0] * cnu + Y[1] * snu)

  ;cnu= cos(nu_obs)
  ;snu= sin(nu_obs)
  ;r  = (1. - excent * excent) / (1. + excent * cnu)
  ;xexpect= r * (X[0] * cnu + X[1] * snu)
  ;yexpect= r * (Y[0] * cnu + Y[1] * snu)

  ; NOTE: thiele-innes contain axis already
  if el.excent lt 1. then begin
      rcnu = cos(E_obs) - el.excent
      rsnu = sqrt(1.-el.excent*el.excent) * sin(E_obs)
  endif else if el.excent eq 1. then begin
      ;; E_obs is in fact true anomaly here
      r    = 2./(1.+cos(E_obs))
      rcnu = r * cos(E_obs)
      rsnu = r * sin(E_obs)
  endif else begin
      rcnu = el.excent - cosh(E_obs)
      rsnu = sqrt(el.excent*el.excent - 1.) * sinh(E_obs)
  endelse
  xexpect= X[0] * rcnu + X[1] * rsnu
  yexpect= Y[0] * rcnu + Y[1] * rsnu
end
