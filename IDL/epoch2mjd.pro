;
; epoch2mjd.pro
; Created:     Thu Jul  3 10:36:35 2014 by Koehler@Quorra
; Last change: Thu Jul  3 10:37:07 2014
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; function to convert Epoch to MJD,
;; because I can never remember where to find this
;; (e.g. http://scienceworld.wolfram.com/astronomy/JulianEpoch.html,
;;	 http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?epoch)
;;
Function Epoch2MJD, Epoch
  return, (Epoch-2000.0)*365.25 + 51544.5
end
