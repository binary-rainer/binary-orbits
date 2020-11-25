;
; plot_orbpnt.pro
; Created:     Tue Mar 18 14:59:28 2014 by Koehler@Quorra
; Last change: Tue Mar 18 15:00:35 2014
;
; PURPOSE:
;	plot one point in an orbit with error ellipse
;	(data in the usual struct with Sep, eSep, PA, ePA)
;
PRO plot_orbpnt, obs, Color=c

  PA = obs.PA * !pi/180.
  daz= obs.Sep * obs.ePA * !pi/180.	;; delta azimuth

  w = findgen(360) * !pi/180.
  xx= (obs.sep + obs.eSep*cos(w)) * sin(PA) - daz*sin(w) * cos(PA)
  yy= (obs.sep + obs.eSep*cos(w)) * cos(PA) + daz*sin(w) * sin(PA)

  if keyword_set(c) then $
    oplot, xx, yy, color=c $
    else oplot, xx, yy
end
