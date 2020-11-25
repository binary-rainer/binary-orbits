;
; elements_to_thieleinnes.pro
; Created:     Thu Feb 27 11:37:41 2014 by Koehler@Quorra
; Last change: Thu Feb 27 11:38:20 2014
;
; PURPOSE:
;	compute Thiele-Innes-constants from elements struct
;

pro Elements_to_ThieleInnes, el, X,Y

  peri = el.periast * !pi/180.
  node = el.ascnode * !pi/180.
  inc  = el.inclin * !pi/180.

  cp = cos(peri) & sp = sin(peri)
  cn = cos(node) & sn = sin(node)
  ci = cos(inc)
  ;puts "P: $co, $so, K: $cO, $sO, I: $ci"

  X = fltarr(2)
  Y = fltarr(2)
  X[0] = el.axis * ( cp * sn + sp * cn * ci)
  Y[0] = el.axis * ( cp * cn - sp * sn * ci)
  X[1] = el.axis * (-sp * sn + cp * cn * ci)
  Y[1] = el.axis * (-sp * cn - cp * sn * ci)
  ;;print, "X:", X
  ;;print, "Y:", Y
end
