;
; quick_plot_orbit.pro
; Created:     Wed Jul  2 11:24:57 2014 by Koehler@Quorra
; Last change: Wed Jul  2 11:56:58 2014
;
; PURPOSE:
;	plot an orbit real quick
;	used by svd_HR_T0
;
PRO quick_plot_orbit, JDobs, x_obs, y_obs, el

  xmax= max(abs(x_obs))*1.1
  ymax= max(abs(y_obs))*1.1*4./3.
  xmax= xmax>ymax
  ymax= xmax*3./4.

  plot,x_obs,y_obs,psym=1,xrange=[xmax,-xmax],yrange=[-ymax,ymax],/xstyle,/ystyle,/isotrop
  ;; note reverse x-axis

  EAnom = compute_EAnom(JDobs,el)
  Elements_to_ThieleInnes, el, X,Y
  mk_orbit, el,X,Y, EAnom, xorb,yorb, xnode,ynode, xexpect,yexpect
  oplot,xorb, yorb, linestyle=4
  oplot,xnode,ynode,linestyle=3
  oplot,[0.,xorb[0]], [0.,yorb[0]],linestyle=2

  for i=0,N_elements(x_obs)-1 do $
    oplot, [x_obs[i],xexpect[i]], [y_obs[i],yexpect[i]]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
