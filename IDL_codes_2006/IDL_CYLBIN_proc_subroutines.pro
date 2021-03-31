;; This function returns a 3-element array of
;;; coefficients of a mean + 1 harmonic fit.
;;; Coefficient 0 is the mean (divergence), 1 is sin(x) or u
;;; component, 2 is cosine or v component.
;;; coeff_error_std is a 3-vector of error estimates in these.

function VAD_1harm, az_in_radians, Vr, input_error_std, $
                    u_guess, v_guess, fitted_output, coeff_output_stderr, $
                    chisq_goodness_of_fit

 coeffs = [0,u_guess,v_guess] ;; initial guess                                                                                                   
 result = svdfit(az_in_radians, Vr, A=coeffs, $
                 FUNCTION_NAME='mean_and_1harm', $
                 MEASURE_ERRORS=input_error_std, SIGMA=coeff_output_stderr, $
                 YFIT = fitted_output, CHISQ=chisq_goodness_of_fit)
 on_error,2
 return, result
end

;;; This is needed by VAD_1harm, it is optimized
function mean_and_1harm, X, M
  return, [ 1.0, sin(x), cos(x) ]
end

pro opt_grid,outputfile,nv,nlon,nlat,nt,v,lon,lat,t,d
 openw,1,outputfile
 printf,1,nv,nlon,nlat,nt
 printf,1,format='(4a20)',v
 printf,1,lon
 printf,1,lat
 printf,1,format='(4(e13.5,1x))',t
 printf,1,format='(4(e13.5,1x))',d
 close,1
 end


pro ipt_grid,inputfile,nv,nlon,nlat,nt,v,lon,lat,t,d

 nt=1.0
 openr,1,inputfile
 readf,1,nv,nlon,nlat,nt
 print,nv,nlon,nlat,nt
 v=strarr(nv)
 lon=fltarr(nlon)
 lat=fltarr(nlat)
 t=fltarr(nt)
 d=fltarr(nv,nlon,nlat,nt)
 readf,1,format='(4a20)',v
 print,v
 readf,1,lon
 print,lon
 readf,1,lat
 print,lat
 readf,1,t
 readf,1,d
 close,1
 
end

pro calc_var,x,f,y_mean,y_var
 n=total(f)
 y_mean=total(f*x)/n
 y_var=(total(f*(x-y_mean)^2.0)/(n-1.0))^0.5
end

pro radar_Z_R,a,b,dbz,Ze,rrate
 
;         assume Ze=a*R^b, Ze in mm^6/m^3, R in mm/hr
 
 Ze=10.0^(dbz/10.0)
 rrate=(Ze/a)^(1.0/b)/3600.0    ; kg/m^2/s
 
end


pro dbz_Vt,dbz,iwc,Vt
;         output:
;           iwc    Kg/m^3
;           Vt     m/s
 
 Ze=10.0^(dbz/10.0)
 Zei=4.68*Ze    ;;;  Churchill and Houze 1984
                         ;;;  (i.e. dbz_i=dbz+6.7) Braun and Houze 1995
 iwc=8.0E-6*Zei^0.61     ;;;  kg/m^3, Herzegh and Hobbs 1980
 Vt=3.29*iwc^0.16        ;;;  m/s Heymsfield and Donner 1990, JAS     
 
end


pro dbz_Vt2,dbz,H,Vt
;         input
;           H  height in km
;         output:
;           Vt     m/s
;         ref: 
;	    Lee, W., Jou, B. J.-D., Chang, P.-L., and Marks, F. D., 2000: 
;            Tropical Cyclone Kinematic Structure Retrieved from Single-Doppler 
;            Radar Observations. Part III: Evolution and Structures of Typhoon 
;            Alex (1987). Mon. Wea. Rev., 128, 39824001.
;         note:
;           In EPIC and TEPPS, bright band locates at 575-625 mb (+3C-0C, or
;           4.25-4.75km)
 
 
 Ze=10.0^(dbz/10.0)
 
;boundary of bright band, H1<H2
 H1=5.0                         ;km
 H2=5.0                         ;km
 
 const=(exp(H/9.58))^0.4
 if(H le H1) then begin
    Vt=2.6*Ze^0.107*const  
 endif 
 if(H ge H2) then begin
    Vt=0.817*Ze^0.063*const
 endif
 if(H gt H1 and H lt H2) then begin
;	    const=(exp(H1/9.58))^0.4
;            Vt1=2.6*Ze^0.107*const  
    
;	    const=(exp(H2/9.58))^0.4
;	    Vt2=0.817*Ze^0.063*const
    
;	    Vt=(Vt1*(H2-H)+Vt2*(H-H1))/(H2-H1)
 endif
 
end


pro radar_microphysics,dbz,Ze,Zei,rwc,iwc,rrate,rratei,Vti
 
 Ze=10.0^(dbz/10.0)
 Zei=4.68*Ze    ;;;  Churchill and Houze 1984
                         ;;;  (i.e. dbz_i=dbz+6.7) Braun and Houze 1995
 rwc   = 5.5E-7*Ze^0.80      ;;;  kg/m^3, Leary and Houze 1979
 iwc   = 8.0E-6*Zei^0.61     ;;;  kg/m^3, Herzegh and Hobbs 1980
 rrate = 3.611E-6*Ze^0.8     ;;;  kg/m^2/s, Hudlow 1979 GATE
 
 calc_Vt,iwc,Vti             ;;;  m/s
 rratei= iwc*Vti             ;;;  kg/m^2/s
 
end


pro calc_Vt,iwc,Vt
;         Heymsfield and Donner 1990, JAS
;         input:  
;           iwc    Kg/m^3
;         output:
;           Vt     m/s
 
 Vt=3.29*iwc^0.16             
 
end

pro box,x,y,w,h,color

x1=x-0.5*w
x2=x+0.5*w
y1=y-0.5*h
y2=y+0.5*h
x3=[x1,x2,x2,x1,x1]
y3=[y1,y1,y2,y2,y1]
polyfill,[x1,x2,x2,x1],[y1,y1,y2,y2],col=color
oplot,x3,y3
end

;+
; NAME:
;   ROUND_ANY
;
;
; PURPOSE:
;
;   Round off the number (or array) N to the nearest Dth decimal and
;     output as a string array.
;
;
; CALLING SEQUENCE:
;
;   nstring = ROUND_ANY(N,[D=d,/ADD_SIGN,SIG=sig])
;
; INPUTS:
;
;   N: The number (or array of numbers) to be rounded.
;
;
; OPTIONAL KEYWORD PARAMETERS:
;
;   D: The number of decimal places. Either a scale or a vector.
;      If D is negative, pads rounds to an integer and pads with zeros.
;      If missing (and no SIG), then just chop off the zeroes after the decimal.
;
;   ADD_SIGN: Add a plus sign (+) to the begining of the string for
;      positive values.
;
;   SIG: The number of significant digits to round off to, rather than
;        decimal places.
;
;
; OUTPUTS:
;
;   nstring: The STRING (or STRARR) of the rounded number(s).
;
;
; MODIFICATION HISTORY:
;	Written by: Chris Torrence.
;	Modified by: Eric Leuliette 9/30/94
;   Modified by: CT July 1996, added SIG keyword.
;-
FUNCTION ROUND_ANY,N,dummy1,dummy2,D=D,ADD_SIGN=add_sign,SIG=sig
	ON_ERROR,2
	IF (N_PARAMS() EQ 3) THEN N=dummy2
	s = SIZE(N)
	type = s(s(0)+1)
	IF (type EQ 0) THEN RETURN,"Undefined"
	IF (type GT 6) THEN RETURN,N  ; structure or string type

	s = N_ELEMENTS(N)
	
	IF (N_ELEMENTS(D) LE 0) THEN D = -999  ; default is zero digits
	IF (N_ELEMENTS(D) EQ 1) THEN D = REPLICATE(([D])(0),s)
	IF (N_ELEMENTS(D) LT s) THEN MESSAGE,"Not enough D values"
	IF (NOT KEYWORD_SET(sig)) THEN sig = 0

	stri=STRARR(s)

	FOR i = 0,(s - 1) DO BEGIN

		d1=D(i)   ; # of places after the decimal
		n1=N(i)   ; copy of the actual number
		IF (n1 EQ FIX(n1)) THEN n1 = FIX(n1)
		sgn = SIGN(n1)

		IF (sig GT 0) THEN BEGIN
			IF (n1 NE 0.) THEN BEGIN  ; sigdigs: round & find d1
				n1 = ABS(n1)   ; make sure to take ABS
				expnt = ALOG10(n1)   ; find Exponent
				expnt = FIX(expnt) - $
					(expnt LT 0.0)*(expnt NE FIX(expnt))   ; correct for n1 < 1
				d1 = (sig - expnt - 1)   ; find conversion exponent
				n1 = LONG(n1*10.^d1 + 0.499999) ; convert n1 & truncate to sig digs
				n1 = sgn*DOUBLE(n1)*10.^(-d1) ; convert back to true# & restore sign
			ENDIF ELSE BEGIN
				d1 = sig
				n1 = 1E-37
				sig = 0
			ENDELSE
		ENDIF


		IF ((sig EQ 0) AND (d1 EQ -999)) THEN BEGIN
			str1 = STRCOMPRESS(n1,/REM)
			decimal = STRPOS(str1,'.')
			IF (decimal GT -1) THEN BEGIN
				WHILE (STRMID(str1,STRLEN(str1)-1,1) EQ '0') DO BEGIN 
					str1 = STRMID(str1,0,STRLEN(str1)-1)
				ENDWHILE
			ENDIF
		ENDIF ELSE BEGIN
			CASE SIGN(d1) OF    ; 0=0.0, -1=negative, +1=positive
			 0: str1=STRCOMPRESS(ROUND(n1),/REM)
			 1: str1=STRCOMPRESS(STRING(n1, $
				FORMAT='(F25.'+STRCOMPRESS(FIX(d1),/REM)+')'),/REM)
			-1: BEGIN
				str1 = STRCOMPRESS(ROUND(n1),/REM)
				sl = STRLEN(str1)+d1
				IF (sl LT 0) THEN str1 = STRING(REPLICATE(48B,-sl))+str1
		    	END
			ENDCASE
		ENDELSE
		
		stri(i) = str1
	ENDFOR

	IF (KEYWORD_SET(add_sign)) THEN BEGIN
		pos = [WHERE(N GT 0)]
		IF pos(0) NE -1 THEN stri(pos) = '!9+!X' + stri(pos)
		pos = [WHERE(N EQ 0)]
		IF pos(0) NE -1 THEN stri(pos) = ' ' + stri(pos)
		neg = [WHERE(N LT 0)]
		IF neg(0) NE -1 THEN stri(neg) = STRING(177b) + STRMID(stri(neg),1,255)
	ENDIF

	IF (s EQ 1) THEN strin=stri(0)
	RETURN,strin
END

