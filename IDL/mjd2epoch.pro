;
; mjd2epoch.pro
; Created:     Thu Jul  3 10:37:15 2014 by Koehler@Quorra
; Last change: Thu Jul  3 10:37:38 2014
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; function to convert MJD to Epoch,
;; because I can never remember where to find this
;; (e.g. http://scienceworld.wolfram.com/astronomy/JulianEpoch.html,
;;	 http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?epoch)
;;
Function MJD2Epoch, MJD
  return, 2000.0 + (MJD-51544.5) / 365.25
end
