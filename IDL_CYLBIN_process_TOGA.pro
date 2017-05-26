;;; This is the IDL code to process and plot hourly VAD histogram
;;; files from CYLBIN. A code was written by Jialin Lin in extensive
;;; back and forth conceptual design with Brian Mapes over 2002-5, and
;;; extensively edited for clarity by Brian Mapes in August 2006. This
;;; edit is so complex that it is almost new. It started w/basic needs like 
;;; commenting blocks of code, but turned into more. Rightness is evident in the
;;; plotted outputs. Istead of  the orginal structure - acomplex
;;; loop over all experiments - this code just works from an input file
;;; List, which contains a list of the histogram filenames hour by
;;; hour to be processed. Also there are some needed parameters. All
;;; these go in a file called IDL_CYLBIN_params.idl which will be
;;; included via the @IDL_CYLBIN_params.idl line below. 

;;; For output, an NHOURS x NP array is filled with data and sent to ouputs
;;; at the end. In other words, this program just processes the hours available.
;;; It leaves the job of filling a seamless regular grid of hours with
;;; data and MISSING flags as a separate processing step. 

;;; Each hour's results are plotted in a postscript file with 2 pages per hour -
;;; the first page each hour is a bunch of postage-stamp VAD plots (every other
;;; altitude level and 3 different radii), the second page is the
;;; hourly summary of echo (planview and CFAD) plus profiles of wind
;;; and divergence.  Finally, a whole-set mean diagram is made. 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SUBROUTINES TUCKED IN HERE

@IDL_CYLBIN_proc_subroutines.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Main Program

pro IDL_CYLBIN_process_TOGA
  
  close,/all
; cp IDL_CYLBIN_params_EPICstratus.idl IDL_CYLBIN_params.idl 
; DO THIS IN A SHELL SCRIPT THAT CALLS THE PROCESSING 
  
@IDL_CYLBIN_params.pro
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Grab and count files (hours)
;;; INPUT LIST: an ASCII list of lines each with the format
;;; YYYYMMHHDD.header
  
  hourly_filename_list='' ;; initial (to be chopped off at end)
  headerfilename=''
;;; Grab all in list
  openr,2,FileNameList
  while not eof(2) do begin
     readf, 2, headerfilename
     hourly_filename_list = [hourly_filename_list, headerfilename]
  endwhile
  close,2
  
  hourly_filename_list = hourly_filename_list[1:*]
  nt = n_elements(hourly_filename_list)
  print, 'NUMBER OF HOURS TO PROCESS: ', nt
  
 
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; CONTROL FLAGS
;;; Apparently B&W plotting code was prepared for publication
;;; purposes. I cut out a lot so it won't work any more! 
  BW_plots_flag=0 ;; keep it 0 unless you make it work
  BW_plots_string=['','bw']
  
;;; The Doppler unfolding requires a first guess wind, which could
;;; come from various sources. For purity, it should come from the VAD
;;; itself, although with sparse data this can be unstable. This is an
;;; topic for agorithmic experimentation still. For now, set iguess=0
;;; and try to do pure VAD with no recourse to other datasets.  
  
  cguess=['','extguess']
  wind_firstguess=0
;print,'first guess wind? 0-radar only, 1-use sonde or ncep when echo coverage is small'
  
;;; Postscript output is very useful, but takes a while to draw
;;; (device=2). device=1 (screen) is only useful for a single hour, in
;;; a debugging mode. device=0 (no plots) speeds up the work. 
  
  IP_FIRST = 0 ;; normal operation: start from the lowest level
;plot_output_device=2 & debug_unfolding=0  ;;; No debugging, output to postscript. Normal. 
  plot_output_device=2 & debug_unfolding=2  ;;; Long debugging, output to postscript. 
  
;;; DEBUGGING FLAGS 1 = screen device
;plot_output_device=1 & debug_unfolding=1  ;; debug_unfolding = 1 for step 1 (relative folding),
;plot_output_device=1 & debug_unfolding=2  ;; debug_unfolding = 2 for step 2 (absolute folding)
;IP_FIRST = 8;; DEBUGGING: Start from a problem level 3=825, 18 = 75mb in 50mb binning
;;; DEBUGGING FLAGS
  
  if ( debug_unfolding eq 2 ) then begin
;!p.multi=[0,8,8]
     !p.multi=[0,6,5]
  endif
  
;print,'0-no plot, 1-window, 2-ps'
  
;;; Vertical coordinate: p or z? Data are in z, but rebinned to p as
;;; needed for analysis. This is useful since upper levels are often
;;; sparsely sampled and p layers encompass several z layers up there
;;; so data are pooled. Also for divergence mass coordinates are
;;; convenient. 
;;; Without loss of generality, let "p" represent
;;; "pooled" and we always want p coordinates. If z is wanted,
;;; just make p the same as z. 
;vert_coord_string=['pres','hgt']
  p_or_z_coord=0 
;print,'vertical coordinate? 0-pressure, 1-height'
  
;;; This was called "quality control" but it's a series of
;;; experimental processing steps that are best left out. I may excise
;;; it from the code later, for now let's just label is more clearly
  
  qc_string=['','expqc']
  experimental_qc_flag=0
;print,'quality control? 0-n, 1-y'
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; define arrays containing raw data for one column (one hour of data) in 
;;; Z coordinate ;24*12*36 grids. See IDL_CYLBIN_params.idl for these
  
  dbzhist = lonarr(ndbz,naz,nr,nz)
  vrhist = lonarr(nvr,naz,nr,nz)
  widhist = lonarr(nwid,naz,nr,nz)
  sdzhist = lonarr(nsdz,ndbz,naz,nr)
  
  if (hists_are_shortint eq 1) then begin
     dbzhist = fix(dbzhist)
     vrhist = fix(vrhist)
     widhist = fix(widhist)
     sdzhist = fix(sdzhist)
  endif
  vrsum = fltarr(naz,nr,nz)
  vr_guess = fltarr(naz,nr,nz)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Por pooling to p vertical coordinate (presure, or pooled) 
;;; must fill ip2iz and iz2ip conversion arrays (read their names like
;;; functions, with iz and ip as inputs and outputs). 
  
  iz2ip=fltarr(nz)              ; pressure level index array for height levels
  ip2iz=fltarr(np)              ; height level index array for pressure levels
  
  for iz=0,nz-1 do begin
     ip=fix((psond(iz)-p(0))/dp+0.5) ;; discretized by dp
     iz2ip(iz)=ip
  endfor
  
  for ip=0,np-1 do begin
     temp=abs(psond-p(ip))
     ip2iz(ip)=where(temp eq min(temp))
  endfor
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Pre-calculate the sine and cosine of the different angles (zenith
;;; = zen, azimuth = theta)
  
  coszen=fltarr(nr,nz)
  sinzen=fltarr(nr,nz)
  for ir=0,nr-1 do begin
     for k=0,nz-1 do begin
        temp=atan(z(k),r(ir))
        coszen(ir,k)=cos(temp)
        sinzen(ir,k)=sin(temp)
     endfor
  endfor
  
  theta = !pi/2 - (0.5+findgen(naz))*!pi*2.0/naz ;; azimuth angle of bin centers
  sintheta = sin( theta ) #(1+fltarr(nz)) 
  costheta = cos( theta ) #(1+fltarr(nz))
  range_km = (0.5+findgen(nr))*dr*0.001 ;; range array (horizontal range)
  PI_180 = !pi/180.
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; For polar plots, neex x, y, coordinates of box centers.
  theta2 = [theta,theta(0)]
  tt = theta2 # replicate(1,n_elements(range_km)) ;; 2D (theta,r) arrays
  rr = replicate(1,n_elements(theta2)) # range_km
  xx = rr * cos(tt)
  yy = rr * sin(tt)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Define some useful things for rainrate calculations, Z-R
;;;; relations and such
  
  dbz_values = (findgen(ndbz,naz,nr) mod NDBZ)
  radar_microphysics,dbz_values,Ze0,Zei0,rwc0,iwc0,rrate_GATEZR,rratei,Vti
  
;;; Use a,b coefficients for "local" Z-R relationship 
  radar_Z_R,a1,b1,dbz_values,Ze0,rrate_localZR ; rrate_local in kg/m^2/s
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; calculate the area of each spatial cell, as a function of range
  
  range = findgen(nr+1)*dr ;;; meters
  aa = !pi * range^2       ;;; sq meters
  aa = aa - shift(aa,1)
  aa = aa(1:*)/(naz+0.0) 
  
;;; fill 2D area array (square meters) - function of range only
  
  area2D = fltarr(naz,nr)
  for j = 0,nr-1 do $
     area2D(*,j) = aa(j)
  
;;; fill 3D area array (square meters) - function of range only
  
  area = fltarr(ndbz,naz,nr)
  for j = 0,nr-1 do $
     area(*,*,j) = aa(j)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Define and initialize the VAD data ouputs in a rigid output
;;; 'structure': an array called d01 with variable names called v0
;;; Yes, horrible programming style, Jialin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nv0=30
  v01=['cover total','dbz (db)','Ze (mm^6/m^3)', $
       'rwc (kg/m^3)','iwc (kg/m^3)', $
       'rain_gate (mm/h)','rain_local (mm/h)', $
       'spectral width (m/s)','u_7r (m/s)','v_7r (m/s)','wspd_7r (m/s)', $
       'div_1r (1/s)','div_3r (1/s)','div_5r (1/s)', $
       'std_div_1r (1/s)','std_div_3r (1/s)','std_div_5r (1/s)', $
       'chi_fit_1r','chi_fit_3r','chi_fit_5r', $
       'az_obs_1r (deg)','az_obs_3r (deg)','az_obs_5r (deg)', $
       'gapmax_1r (deg)','gapmax_3r (deg)','gapmax_5r (deg)', $
       'std_u_7r (m/s)','std_v_7r (m/s)','chi_fit_7r','stratiform fraction']
  d01=fltarr(nv0,nr,np,nt)-9999.0
  
;;; Define the much smaller output dataset 3 (just these 2 variables)
  nv03=2
  v03=['cover_48km','cover_88km']
  nr03=ndbz-1
  r03=dbz(1:*)
  d03=fltarr(nv03,nr03,np,nt)-9999.0
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; define arrays for storing data for savefile of outputs
  timestrings=strarr(nt)
  cover_0db = fltarr(nt,nr,np)
  cover_15db = fltarr(nt,nr,np)
  cover_30db = fltarr(nt,nr,np)
  rain_plotalt_GATEZR = fltarr(nt,nr)
  rain_plotalt_localZR = fltarr(nt,nr)
  stratiform_fraction = fltarr(nt,nr)
  cfad0 = fltarr(ndbz-1,np,nt)
  uwind = fltarr(nt,np) ;;; VAD wind
  vwind = fltarr(nt,np)
  u_std = fltarr(nt,np) ;;; Wind std error from sampling variability
  v_std = fltarr(nt,np) 
  divergence = fltarr(nt,nr,np) ;; VAD divergence
  div_stdev = fltarr(nt,nr,np)  ;;; std error of div from sampling variability
  chi_uv = fltarr(nt,np)        ;;; goodness of fit for wind
  chi_div = fltarr(nt,nr,np)    ;;; godoness of fit parameter for div
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; open postscript file for output plots   
  
  file_plot=PSFILE_STUB+'.ps'
;  '_'+qc_string(experimental_qc_flag) $
;  +cguess(wind_firstguess)+BW_plots_string(BW_plots_flag) +'.ps'
  
  if(plot_output_device eq 2) then begin
     tall
     device, file=file_plot
  endif
  
  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; At this stage, it is time to read in the list of filenames for 
;;; all the hourly data to be processed: filenames all end in .header,
;;; but the key input data are histograms
;;; of dBZ (called Zhist), SDZ, a boxcar stdev of dbz along the ray 
;;; that feels for sharp gradients (convective vs. stratiform measure),
;;; Vr (Vhist), and Doppler spectral width
;;; (called Whist), along with a Vsum file containing the floating
;;; point sum of all Doppler velocities in each spatial grid cell. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;                                                                               
;;; Here is the main time loop over hours, opening datafiles and
;;; grabbing their contents.
  
  for ihr=0,nt-1 do begin
     
     headerfilename=hourly_filename_list(ihr)
     len =strlen(headerfilename)-7
     filebase =strmid(headerfilename,0,len)
     hour =strmid(filebase,9,/REVERSE_OFFSET)
     timestrings(ihr)=hour
     
     print,ihr, headerfilename
    
     
;;; Unzip all files for the hour. Then, grab the 3 hists and Vrsum. 
;;; TIME POOLING? COULD JUST WRAP THIS IN A LOOP AND SUM THE ARRAYS
;;; NAW, 1 HOUR POOLS ARE LONG ENOUGH 
     
;   spawn,'gunzip '+filebase+'*gz'
     openr,1,filebase+'.Zhist.gz', SWAP_ENDIAN=NEED_TO_SWAP_ENDIAN_FLAG, /compress
     readu,1,dbzhist
     close,1
     openr,1,filebase+'.Vhist.gz', SWAP_ENDIAN=NEED_TO_SWAP_ENDIAN_FLAG, /compress
     readu,1,vrhist
     close,1
     openr,1,filebase+'.Whist.gz', SWAP_ENDIAN=NEED_TO_SWAP_ENDIAN_FLAG, /compress
     readu,1,widhist
     close,1
     openr,1,filebase+'.SDZhist.gz', SWAP_ENDIAN=NEED_TO_SWAP_ENDIAN_FLAG, /compress
     readu,1,sdzhist
     close,1
     openr,1,filebase+'.Vsum.gz', SWAP_ENDIAN=NEED_TO_SWAP_ENDIAN_FLAG, /compress
     readu,1,vrsum
     close,1
     
     
;;; If we are processing batch (to ps), go ahead and rezip them 
;   spawn,'gzip '+filebase+'*hist* ' +filebase+'*sum*'
     
;;; Make floating point arrays. If negative values (which can happen
;;; if short integers are incremented past 32767), then some
;;; cleverness will be needed to decide how many observations there
;;; really were -- shoulda used long integers!!
     
     fdbzhist=float(dbzhist)
     fvrhist=float(vrhist)
     fwidhist=float(widhist)
     fsdzhist=float(sdzhist)
     
     if(min([dbzhist,vrhist,widhist]) lt 0 and HISTS_ARE_SHORTINT) then begin
        print,'Negative histogram values - careful (short int overtopped)'
        fdbzhist = fdbzhist + 32767.*2*(dbzhist lt 0) 
        fvrhist = fvrhist + 32768.*2*(vrhist lt 0) 
        fwidhist = fwidhist + 32768.*2*(widhist lt 0) 
     endif
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; For pressure coordinates, we will pool the histograms and sum the
;;; sum arrays within sets of Z levels (up to 5 0.5km Z levels in a
;;; 50mb pressure level at high altitudes)
     
;;;; These are p-coordinate versions of the histogram and Vr sum
;;;; arrays, declared here for size and to make them 0 each hour
     
     fdbzhist_p=fltarr(ndbz,naz,nr,np)
     fvrhist_p =fltarr(nvr,naz,nr,np)
     widhist_p =fltarr(nwid,naz,nr,np)
     sdzhist_p =fltarr(nwid,ndbz,naz,nr) ;; not a function of height 
     
     Vr_p=fltarr(naz,nr,np)     ;; starts as Vrsum_p, then is divided by nobs, unfolded etc. 
     Vr_std=fltarr(naz,nr,np)+99. ;; a standard deviation estimate for above
     ;; (sampling standard deviation based on Vr histogram)
     Vr_fit = fltarr(naz,nr,np) ;; a harmonic fit based radial velocity
     div_fit = fltarr(nr,np)    ;; ditto for divergence
     div_std = fltarr(nr,np)    ;; std of above
     u_p=fltarr(nr,np)          ;; zonal wind profile, for different radii
     v_p=fltarr(nr,np)          ;; v wind
     nobsVr_p=fltarr(naz,nr,np) ;; number of obs in each p bin
     cfad_p=fltarr(ndbz-1,np)   ;; a contoured frequency diagram in pressure
     nobswid_p=fltarr(naz,nr,np) ;; nobs for spectral width
     wid_p=fltarr(naz,nr,np)     ;; mean spectral width
     nobsSDZ=fltarr(naz,nr)      ;; nobs for SDZ
     SDZ=fltarr(naz,nr)          ;; mean SDZ
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
;;; Pool dbZhist only, and compute echo statistics. 
;;; Vr and width are pooled in Doppler VAD processing loop. 
     
     for iz=0,nz-1 do $
        fdbzhist_p[*,*,*, iz2ip(iz) ] = fdbzhist_p[*,*,*, iz2ip(iz) ] + $
        fdbzhist[*,*,*, iz ]
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      
;;; divide reflectivity historgam by total number of observations to 
;;; get fractional coverage
     
     nobs_p_dbz = total(fdbzhist_p,1)
     for idbz=0,ndbz-1 do $ 
        fdbzhist_p(idbz,*,*,*)=fdbzhist_p(idbz,*,*,*)/(nobs_p_dbz >1)
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate statistics of reflectivity, and
; store in d01, d03. These are coverage, rainrate, etc. 
     
     for k=0,np-1 do begin ;;; loop over levels
        
        cover=fdbzhist_p(*,*,*,k)*area        ;m^2 of detectable echo
        area_j=total(area(0,*,CFAD_RANGE_BINS)) ;m^2 total area
        cfad_p(*,k)=total(total(cover(1:*,*,CFAD_RANGE_BINS),3),2)/area_j
        d03(0,*,k,ihr)=cfad_p(*,k)
        area_j=total(area(0,*,CFAD_RANGE_BINS2))
        cfad_p(*,k)=total(total(cover(1:*,*,CFAD_RANGE_BINS2),3),2)/area_j
        d03(1,*,k,ihr)=cfad_p(*,k)
        
        for j=0,nr-1 do begin        ;;; loop over radii for curcular stats
           area_j=total(area(0,*,0:j)) ;m^2
           cover_j=cover(1:*,*,0:j)    ;m^2
           
           d01(0,j,k,ihr)=total(cover_j)/area_j
           cover_0db(ihr,j,k)=total(cover(1:*,*,0:j))/area_j
           cover_15db(ihr,j,k)=total(cover(16:*,*,0:j))/area_j
           cover_30db(ihr,j,k)=total(cover(31:*,*,0:j))/area_j
           
;;;; For a given Z-R, Z-IWC, etc, the histograms of dbZ translate into 
;;;; area averaged rainrate R, IWC, etc.  
           d01(1,j,k,ihr)=total(dbz_values(1:*,*,0:j)*cover_j)/area_j
           d01(2,j,k,ihr)=total(Ze0(1:*,*,0:j)*cover_j)/area_j
           d01(3,j,k,ihr)=total(rwc0(1:*,*,0:j)*cover_j)/area_j
           d01(4,j,k,ihr)=total(iwc0(1:*,*,0:j)*cover_j)/area_j
;         d01(5,j,k,ihr)=total(rrate_GATEZR (1:*,*,0:j)*cover_j)/area_j  ;; rrate is kg m-2 s-2
;         d01(6,j,k,ihr)=total(rrate_localZR(1:*,*,0:j)*cover_j)/area_j;; so time integral is mm
;;; I doubt the correctness of the preceding2 lines d01(5,6). Reset later. 
           
        endfor
     endfor
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; DOPPLER: wind and divergence profiles
;;;;; 1. unfold Vr values within each spatial bin, to make all the values
;;;;;   in Vrhist be in the same nyquits interval (right or wrong) 
;;;;; 2. unfold the resulting spatial bin mean Vr to physically correctness.
;;;;;    This requires a first guess wind profile u(p), v(p). For
;;;;;    purity, do a VAD processing of the radar data themselves  to
;;;;;    get this first guess, using data from a very broad range pool
;;;;;    of rpool = 3, i.e. 7 range bins centered around bin 7. For
;;;;;    this initial wind profile, try first guess winds of various
;;;;;    speeds and directions, maximizing goodness-of-fit for the VAD
;;;;;    harmonic fit. For vertical continuity, we start 
;;;;;    at the surface and work up, disallowing vertical wind shear > MAX_DU
;;;;;    Temporal continuity could be brought in from previous
;;;;;    hours, etc etc. There is lots to try in difficult cases. In simple
;;;;;    cases (great coverage, weak winds) any approach is OK. 
;;;;; 3. After getting the u(p) and v(p) profiles and unfolding
;;;;;    all Vr values to be consistent with that, do VAD analyses
;;;;;    again at all radii and all range pools (rpool = 2,1,0) and
;;;;;    this time keep the divergence profiles as a key output, but
;;;;;    ignore the wind (since radius-dependent VAD wind doesn't make
;;;;;    much sense scientifically).  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
;;;; 1. UNFOLDING OR DEALIASING
     
;;; Figure out the bins of the Vrhist array
     nVrbins=fix((vnyq*2.0/dvr)+0.999)   
     if(nVrbins lt 2) then begin
        print, 'Too few Vr bins - something weird?'
       ; stop
     endif
;   central_bin=fix(nVrbins/2)
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Entering a big loop here - over range pools. note m index set
;;;; here, used later in d01 stuff I dare not touch.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     nrpool_halfs=[0,1,2,3]
     for irpool=3,0,-1 do begin
        m=11+irpool ;;; CODED value for where in d01 array the results go
        nrpool_half=nrpool_halfs(irpool)
        nrpool=nrpool_half*2+1           ;;; # of range bins being pooled
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Construct range and height (pressure) pools, by looping over iz
;;; and ira. First, reset initial values to zero. 
        
        nobsVr_p(*,*,*)=0.0
        fvrhist_p(*,*,*,*)=0.0
        Vr_p(*,*,*)=0.0
        nobswid_p(*,*,*)=0.0
        wid_p(*,*,*)=0.0
        nobsSDZ(*,*)=0.0
        SDZ(*,*)=0.0
        
        for ira=nrpool_half,nr-1-nrpool_half do begin
           for ira1=ira-nrpool_half,ira+nrpool_half do begin
              for iz=0,nz-1 do begin
                 ip=iz2ip(iz)
;;; Histogram pooling 
                 fvrhist_p(*,*,ira,ip)=fvrhist_p(*,*,ira,ip) $
                                       +fvrhist(*,*,ira1,iz)
                 nobswid_p(*,ira,ip)=nobswid_p(*,ira,ip) $
                                     +total(fwidhist(0:nwid-2,*,ira,iz))
;;; Sum up the Vr sum
                 Vr_p(*,ira,ip)=Vr_p(*,ira,ip)+vrsum(*,ira1,iz)
;;; Sum up spectral widths to build a mean
                 wid_p(*,ira,ip)=wid_p(*,ira,ip) $
                                 +total(wid(0:nwid-2)*fwidhist(0:nwid-2,*,ira,iz))
              endfor            ; iz loop
              
;;; SDZ not a function of height - 0 bin is "missing" data so go 1:*
              nobsSDZ(*,ira)=nobsSDZ(*,ira) $
                             +total(total( fSDZhist(1:*,2:*,*,ira) ,2),1) ;;; loop over iaz weighted-summing the histogram of non-missing values
              for iaz = 0,naz-1 do $
                 SDZ(iaz,ira)=SDZ(iaz,ira) + $
                 total(SDZ_values* total(fSDZhist(1:*,2:*,iaz,ira),2) ,1)
              
           endfor ;;ira1
        endfor    ;;ira
        
;;; Compute nobsVr_p
        nobsVr_p=total(fvrhist_p(1:nVrbins,*,*,*) ,1)
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate mean Vr from Vrsum/(# of observation), likewise for wid
        
        Vr_p = Vr_p/(nobsVr_p > 1) ;;; the >1 trick is to avoid dividing 0 by 0
        wid_p = wid_p/(nobswid_p > 1)
        SDZ = SDZ/(nobsSDZ > 1)
        
;;; Declare MISSING (-9999) if there are < min_obs (just 1). 
        
        temp=where(nobsVr_p lt nobs_min)
        if(temp(0) ge 0) then Vr_p(temp)=-9999.0
        temp=where(nobswid_p lt nobs_min)      
        if(temp(0) ge 0) then wid_p(temp)=-9999.0
        
;;; Declare MISSING wind guesses initially      
        u_p(*,*)=-9999.0
        v_p(*,*)=-9999.0
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Unfolding step 1: unfold values in each spatial bin to a common
; nyquist interval. Try all periodic shifts of the Vr histogram and
; choose the one that minimizes the variance implied by that
; histogram. That scenario implies that a certain number of the values
; are folded with respect to the others, so adjust Vr_p in the way
; implied by the number of folded values. 
;;;;; This is done spatial bin by spatial bin, so loop over p, r, az
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; p loop
;      for ip=0,np-1 do begin
        for ip=IP_FIRST,np-1 do begin ;;; For debugging, use IP_FIRST
           iz=ip2iz(ip)
           
;;; Vterm - terminal Fallspeed at this altitude
           dbz1=findgen(100) ;; a dummy array of dbz values
           dbz_Vt2,dbz1,z(iz),Vti
           
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; r loop - what range of r's?
           ira1=0
           ira2=nr-1
           if(irpool eq 3) then begin ;;; this is the big fat mean-wind VAD
              ira1=7
              ira2=7
           endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; r loop
           for ira=ira1,ira2 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; azimuth loop
              for iaz=0,naz-1 do begin
                 
;print, 'ip, ira, iaz ',ip,ira,iaz
                 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; UNFOLD VR WITHIN THIS SPATIAL BIN
;;; OK, grab the Vr histogram (floating point - fractional values)               
                 
                 Vr_p_temp=Vr_p(iaz,ira,ip)
                 fvrhist_temp =reform(fvrhist_p(1:nVrbins,iaz,ira,ip))
                 nobs_temp=total(fvrhist_temp)
                 
;;; used for unfolding-debugging plots only
                 fvrhist_tempHR00 =rebin(fvrhist_temp,10*nVrbins,/sample)
                 
;;; Correct fvrhist_temp (the local copy of the histogram) for the partial
;;; filling of the last bin, if 2* the nyquits interval does not divide
;;; perfectly into the Vr histogram bin width (partial bin filling)
;;; This is ONLY used for the histogram -> variance estimate step.
;;; One must be careful not to unfold observations that didn't occur.
                 
                 lastbin_underfill_inv = dvr/((2*Vnyq) mod fix(dvr))
                 fvrhist_temp_cor=fvrhist_temp
                 fvrhist_temp_cor(nVrbins-1)=fvrhist_temp(nVrbins-1)*lastbin_underfill_inv
                 
;;; skip if empty
                 if(Vr_p(iaz,ira,ip) lt -9990.0 or nobs_temp lt nobs_min) then begin
                    Vr_p(iaz,ira,ip) = -9999.0
                    Vr_std(iaz,ira,ip) = 9999.0
                    print, 'No good data - ', Vr_p_temp, ' nobs ', nobs_temp
                    goto,jump22
                 endif
                 
;;; Interpolate histogram to high res (HR) to get a better variance estimate
                 fvrhist_tempHR1 =rebin(fvrhist_temp_cor,10*nVrbins,/sample)
                 vrbins=findgen(10*nVrbins)*dvr/10.0 ;; vr bin values
                 
;;; Compute mean and variance (really, stdev), from histogram. Mean
;;; is not used since we have a more accurate mean via Vrsum.
                 calc_var,vrbins,fvrhist_tempHR1,vr_mean,vr_var 
                 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; IF DEBUGGING LEVEL 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; UNFOLDONG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
                 if ( debug_unfolding eq 1 ) then begin
                    print, 'P, IRA, IRAPOOL, IAZ ', p(ip), ira, irpool, iaz
;print, '      Histogram '
;print, fvrhist_temp
;print, ' corrected ...'
;print, fvrhist_temp_cor
;print, ' '
                    !p.multi=0
                    
                   ; stop
;;;;;;;;;;;;;;; As an aside (inspired by AMMA weird noise characteristics), 
;;; Plot the grand pooled stratospheric noise Vr histogram
;;; And subtitle with its variance, correcting for last-bin underfill
                    
                    strat_hist = total(total( vrhist(1:nVrbins,*,*,35) ,2),2)
                    strat_hist(nVrbins-1)=strat_hist(nVrbins-1)*lastbin_underfill_inv
                    strat_histHR =rebin(strat_hist,10*nVrbins,/sample)
                    calc_var,vrbins,strat_histHR,vr_meanstrat,vr_varstrat
                    
                    vrbinvalues = -vnyq + 4*findgen(nVrbins)
                    plot, 2+vrbinvalues, total(total( vrhist(1:nVrbins,*,*,35) ,2),2), psym=10, ysty=0, $
                          tit='Histogram of raw Vr values in 17.5-18 km layer', xtit='m/s', ytit='number',$
                          xra=[-vnyq, -vnyq+4*nVrbins], subtit='mean '+str(vr_meanstrat)+', std '+str(vr_varstrat)
                    ver, vnyq, lines=1
                  ;  stop
                    
;;; OK, back to main debugging task for this level-range-az bin
                    !p.multi=[0,3,4]
                    nhr = n_elements(fvrhist_tempHR1)
                    vr_range = findgen(nhr)/(nhr-1) * 2*vnyq -vnyq ;; just for debug plots
                    plot, vr_range, fvrhist_tempHR1>0.1, $
                          tit=str(iaz)+' ORIGINAL var ' +str(vr_var), /ylog, $
                          yrange=[1,max(fvrhist_tempHR1)*2], subtitle='red dits: uncorrected'
                    oplot, vr_range, fvrhist_tempHR00>0.1, lines=1, color=250
                    ver, Vr_p(iaz,ira,ip), color=250
                 endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; ENDIF DEBUGGING LEVEL 1
                 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Loop over periodic shifts of the histogram, seeking minimum
;;; variance: when all the Vrs are most internally consistent.
                 
                 for central_bin=1,nVrbins-1 do begin
                    fvrhist2=[fvrhist_temp_cor(central_bin:*), $
                              fvrhist_temp_cor(0:central_bin-1)]
                    fvrhist_tempHR=rebin(fvrhist2,10*nVrbins,/sample)
                    calc_var,vrbins,fvrhist_tempHR,vr_mean,vr_var2
                    
;; Set large variance when too little data are present
;; nobs_min_var=nVrbins*2 ;; based loosely on histogram theory?
                    if(nobs_temp lt nVrbins*2) then vr_var2=stdev(vrbins)
                    
                    if ( debug_unfolding eq 1 ) then $
                       plot, fvrhist_tempHR  >0.1, tit=str(central_bin)+' var '+str(vr_var2),/ylog
                    
;;; test for lower variance. 
;;; If lowest yet, GRAB IT and adjust Vr_p_temp accordingly for
;;; partial folding implied by this histogram shift. NOTE only 
;;; adjust the actual obs (fvrhist_temp) not the corrected
;;; pseudo-histogram fvrhist_temp_cor
                    
                    if(vr_var2 lt vr_var) then begin
                       vr_var=vr_var2
                       Vr_p_temp=Vr_p(iaz,ira,ip) $
                                 -2*vnyq*total(fvrhist_temp(central_bin:*))/nobs_temp
                       
                       if ( debug_unfolding eq 1 ) then begin
                          plot, vr_range, fvrhist_tempHR1>0.1, $
                                tit=str(central_bin)+'IS BETTER: var '+str(vr_var), $
                                xtit='Vr (m/s)', ytit='# obs', yra=[0,1.1*max(fvrhist_tempHR1>0.1)]
;                          subtit='vert. line shows consequence'
                          ver, Vr_p(iaz,ira,ip), color=250
                          ver, vr_p_temp, color=150
;;; Vertical lines show the new mean and +/- stdev of the wind
                          ver, vr_p_temp+2*VNYQ, color=200
                          ver, vr_p_temp+2*VNYQ+[-1,1]*vr_var, lines=1, color=200
                          ver, vr_p_temp-2*VNYQ, color=100
                          ver, vr_p_temp-2*VNYQ+[-1,1]*vr_var, lines=1, color=100
                       endif
                       
                    endif
                 endfor         ; loop over central_bins
                 
;;;;; ADJUST ERROR BARS UPWARD FOR NONTRUSTED DATA
; ADJ 1. Set aminimum # of obs allowed for trusting the variance , or else 
; assign a large value (stdev(vrbins)). This step fixes a former 
; problem: regions with few obs got a small sample variance and hence
; were perversely viewed as highly certain in the VAD harmonic fit. 
; USED ABOVE      if(nobs_temp lt nobs_min_var) then vr_var=stdev(vrbins)
;;; Could be inflated further if problem recurs
                 
; ADJ 2. AMMA data has stratospheric (no meas) Vr "speckle" that has a
; nonuniform (U-shaped) Vr histogram, systematically. This makes the
; Vr values from stage-1 unfolding be near +/- Vnyq, rather than near
; zero or random. So I want to inflate the error bars on Vr estimates
; that come from histograms as broad as this speckle. I don't
; understand the origins of this - deep in radar electronics - so it
; will be deployment specific (hopefully constant thru AMMA? will
; check). Use a parameter in the .params.idl file to define this. Some
; previous parameter files may not have this, so check if it exists
; before using it so this code can be used with or without this
; kluge. 

;;;;;;; Commented start: Nov 09, 2015 ;;;;;;;;;;;;;;;;;
;
;                 if ( n_elements(vr_speckle_stdev_floor) GT 0) then $
;                    if( (vr_var gt vr_speckle_stdev_floor) and $
;                        (abs(vr_p_temp) gt LIKELY_VR_FLOOR_SPECKLE) ) then $
;                           vr_var = inflated_speckle_stdev

;;;;;;; Commented end: Nov 09, 2015 ;;;;;;;;;;;;;;;;;

;;;;;;; Replaced the above speckle correction with the new one below: Nov 09, 2015 ;;;;;;;;;;;;;;;;;
;;;;; In some cases, variance estimate is too small, because of random speckle treated as real
;;;; data- to correct for that, we are assigning large error bars, which means these values
;;;;; will not be used for the fitting

mean_speckle_hist=mean(fvrhist_temp(0:nVrbins-2));
std_speckle_hist=stdev(fvrhist_temp(0:nVrbins-2));

if ((std_speckle_hist/mean_speckle_hist) lt flathist_thresh) then $
  vr_var=inflated_speckle_stdev;
  

                 
;;;;;;;;;;;;; DONE WITH STAGE 1 UNFOLDING AND SUBSEQUENT STDEV ADJUSTMENTS
;;; Having found minimum variance, update the main Vr_p mean Vr and
;;; its sampling standard deviation Vr_std ("variance" means std not
;;; std^2 here)
                 Vr_p(iaz,ira,ip)=Vr_p_temp
                 Vr_std(iaz,ira,ip)=vr_var
;print, 'Final stdev: ', vr_var, Vr_std(iaz,ira,ip)
;print, 'FINAL MEAN WIND: ', vr_p_temp, '     at ',iaz, ira, ip
                 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; VTERM correction
;correction for hydrometer falling velocity
                 
                 Vt1=Vti(1+dbz_offset1:ndbz-1+dbz_offset1)
                 fdbzhist1=reform(fdbzhist_p(1:ndbz-1,iaz,ira,ip))
                 Vt1a=0.0
                 fdbzhist1a = total(fdbzhist1) >1e-10
                 if(fdbzhist1a gt 0.0) then Vt1a=total(Vt1*fdbzhist1)/fdbzhist1a
                 Vr_p(iaz,ira,ip) = Vr_p(iaz,ira,ip)+Vt1a*sinzen(ira,iz)
                 
;shift to within +-vnyq of zero. Why? Does it matter?? Try commenting
;it out
                 temp=(0.0-Vr_p(iaz,ira,ip))/(2.0*vnyq)
                 if(abs(temp) ge 0.5) then $
                    Vr_p(iaz,ira,ip) = Vr_p(iaz,ira,ip)+2*vnyq*fix(temp+0.5*sign(temp))
                 
                 jump22: ;;; skip steps where data are missing
              endfor     ;; iaz loop closes. Still in range loop, p loop, and rpool loop. 
              
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;VAD7 range pool 7 VAD derives mean wind profile. Until jump26
              if(irpool eq 3) then begin ;derive first guess wind by VAD
                 
                 Vr1=Vr_p(*,ira,ip) ;; Make a trial copy of Vr at this pool & pressure level
                 
;;; Count only the good ones
                 where_Vr_is_good=where(Vr1 gt -9990.0)
;;; Surface clutter is characterized by too narrow a histogram of vr,
;;; try killing it this way
;if (ip lt 2) then where_Vr_is_good=where(Vr1 gt -9990.0 and Vr_std(*,ira,ip) ge dvr/2)
                 
;;; If too few values, don't bother even trying all this VAD stuff
                 n=n_elements(where_Vr_is_good)
                 if(n lt 5) then PRINT, 'n lt 5, JUMPING TO 26'
                 if(n lt 5) then goto,jump26
                 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;unfolding stage 2 - absolute value of the radial wind. 
;;; This is ONLY done once, for irpool=3, to derive the grant mean wind
;;; profile used to unfold all subsequent Vr data. In this section of
;;; the code, first guess radial wind Vr_guess, derived from u_guess(p), v_guess(p)
;;; is used as a hypothesis for the true Vr values as a function of azimuth at
;;; this altitude. It is based on a 1-harmonic fit (mean flow +
;;; divergence). We hope the convective turbulence within the ~100km
;;; range of the grand range pool is not so intense that it deviates
;;; from such a hypothetical wind field by more than 2*Vnyq. If so we
;;; are screwed in terms of histogram-based unfolding anyway. For each
;;; hypothesis (trial) Vr_guess, shift the radar Vr_p to within
;;; Vr_guess +- Vnyq and run the VAD curve fitting. If the score
;;; (chisq) is the best yet, keep it. When the trials run out, the
;;; smallest error that came up in those trials WINS and determines
;;; the wind profile. So the big question is, what to use for
;;; u_guess(p), v_guess(p)? Work upward, using the previous (lower)
;;; level and making trial guesses around that. The wind shear from
;;; one level to the next is bounded by a parameter (SHEAR_MAX(p)) set
;;; in the parameters file. 
                 
;;; The strategy is to try all wind SHEAR magnitude and direction 
;;; combinations, relative to the previous level (lower altitude), 
;;; seeking the best VAD fit to the radar data. Not all, but
;;; up to SHEAR_MAX(ip) in 0.25 VNYQ steps, in 30 degree WD sectors. 
;;; Note that all this is Still inside range, pressure loops.
                 
;;; Set up how many discrete shear mags and dirs to try
                 dshearmag=0.25*vnyq
                 nshearmag=SHEAR_MAX(ip)/dshearmag
                 dsheardir=30.0
                 nsheardir=360.0/dsheardir
                 
;;; The strategy is to minimize chi_min, so initialize it with a huge
;;; value and set other stuff to zeros
                 chi_min=1.0E+20
                 coeff_best=[0.0,0.0,0.0]
                 coeff_stdev_best=[0.0,0.0,0.0] ; Arunchandra 
                 Vr_fit_min=fltarr(naz)
                 Vr1_min=fltarr(naz)
                 
;;; Loop over all possible wind shearmag trial values
                 for ishearmag=0,nshearmag-1 do begin
                    shearmag=dshearmag*ishearmag
                    
                    isheardir1=0 ;;; LOOP LIMITS FOR DIRECTIONS
                    isheardir2=nsheardir-1
;;; save time by not trying different directions for 0 magnitude
                    if(ishearmag eq 0) then isheardir2=0
                    
;;; Loop over all possible wind shear directions
                    for isheardir=isheardir1,isheardir2 do begin	
                       sheardir=dsheardir*isheardir
                       
;;; TRIAL VR(az) CURVE IS SET HERE: Previous value (below) plus shear. 
;;; Previous value will be set to 0 for lowest level. 
                       ubelow = u_p[ira,(ip-1)>0] 
                       vbelow = v_p[ira,(ip-1)>0]
                       ubelow = ubelow *(ubelow gt -999)
                       vbelow = vbelow *(vbelow gt -999)
                       Vr_guess1=shearmag*cos((az-sheardir)*PI_180) + $
                                 ubelow*sin(az*PI_180) + $
                                 vbelow*cos(az*PI_180) 
                       
;;; Loop over azimuths, unfolding the Vr toward the shearmag, sheardir hypothesis                     
                       for iaz=0,naz-1 do begin
                          if(Vr1(iaz) lt -9990.0) then print, 'Vr1(iaz) lt -9990.0: jump27'
                          if(Vr1(iaz) lt -9990.0) then goto,jump27 ;; bail out
                          
                          temp=(Vr_guess1(iaz)-Vr1(iaz))/(2.0*vnyq)
                          if(abs(temp) ge 0.5) then begin
                             Vr1(iaz)=Vr1(iaz)+2*vnyq*fix(temp+0.5*sign(temp))
                          endif
                       endfor   ; iaz refolding loop
                       
;;; Now try the VAD wind harmonic fit exercise- Velocity (y) vs Az (x). 
                       x=az(where_Vr_is_good)*PI_180
                       y=Vr1(where_Vr_is_good)
                       input_error_std=Vr_std(where_Vr_is_good,ira,ip)
                       u_guess1=shearmag*sin(sheardir*PI_180)
                       v_guess1=shearmag*cos(sheardir*PI_180)
                       coeff1 = VAD_1harm(x, y, input_error_std, $
                                          u_guess1, v_guess1, fitted_output, coeff_output_stderr, $
                                          chisq_goodness_of_fit)
                                          
                       coeff=coeff1([2,1,0])
                                             
                       
;;; IMPROVEMENT: Quantify the mismatch of the DERIVATIVE of the fitted
;;; output with the DERIVATIVE of the Vr(az) data. Produce another
;;; term to be added to the goodness of fit to make a score that
;;; penalizes folding shocks. This is to fix a tenedency for the
;;; normal goodness of fit to work against some strong localized winds
;;; (like behind a gust front on 2006081103.header from
;;; AMMA-Niamey. The score should be just like the chi-squared_gof
;;; returned from SVDfit: a mean square of the normalized deviations
;;; from the fitted curve. Use periodic shift, a right-handed
;;; derivative if you will, and error bars thereupon. 
;                     vraz_gof =   mean( ((fitted_output-y)/input_error_std)^2 )
;                     vraz_gof =   chisq_goodness_of_fit/n_elements(x) ;; these are equal
                       dvrdaz_fit = shift(fitted_output,-1)-fitted_output
                       dvrdaz_dat = shift(            y,-1)-y
                       dvrdaz_std = sqrt( input_error_std^2 + shift(input_error_std,-1)^2 )
                       dvrdaz_chisq = total( ((dvrdaz_fit-dvrdaz_dat)/dvrdaz_std)^2 )
                       
;;; REDEFINE SCORE: UNEQUAL PARTS FIT AND FIT OF DERIV - more on deriv?
                       chisq_goodness_of_fit = chisq_goodness_of_fit + 5*dvrdaz_chisq
;;; END OF IMPROVEMENT: Quantify the mismatch of the DERIVATIVE...
                       
                       if( debug_unfolding eq 2 ) then begin
                          print, 'SHEARMAG GUESS: ', ishearmag, shearmag, p(ip)
                          print, 'SHEARDIR GUESS: ', isheardir, sheardir
                          print, 'current and best chi_sq so far: ', chisq_goodness_of_fit, chi_min
                          
;;; Diamonds
                          ;plot, x, y, tit=str(ira)+' '+str(p(ip))+' '+ $;commented
                                ;str(chisq_goodness_of_fit), psym=4, yra=[-vnyq,vnyq]*1.5, ystyle=1;commented
;;; Whiskers
                          ;for i=0,n_elements(x)-1 do $;commented
                          ;oplot,x(i)+[0,0],y(i)+[-1,1]*input_error_std(i);commented
                          ;oplot, x, vr_guess1+vnyq, color=150;commented
                          ;oplot, x, vr_guess1-vnyq, color=150;commented
                          ;oplot, x, fitted_output, lines=1, color=200;commented
                          
                       endif
                       
;;; If it's the best fit so far, grab it...
                       if(chisq_goodness_of_fit lt chi_min) then begin
                          
;;; And add red to the debug plot for clarity
                          if( debug_unfolding eq 2 ) then begin
                             oplot, x, fitted_output, color=250 ;;; OVERPLOT RED ON THE ONES THAT ARE BEST
                         ;    stop
                          endif
                          
; Grab the best fit so far                        
                          chi_min=chisq_goodness_of_fit
                          coeff_best=coeff
                          coeff_stdev_best=coeff_output_stderr
                          Vr_fit_min(where_Vr_is_good)=fitted_output
                          Vr1_min=Vr1
                       endif
                       
                    endfor ;;; sheardir loop
                 endfor    ;;; shearmag loop
                 
;;; BAIL OUT POINT to jump27 if some Vr data are bad (-9999)
;;; Why not all the way to jump26? When are Vr data -9999? I forget,
;;; so leave code the way i found it. I think this means the present
;;; level will be set to values from the previous level by lines
;;; below between these 2 jump statements. 
jump27:
                 
;;; Having found the best fit VAD, grab that wind as the profile for
;;; all subsequent unfolds. 
                 
                 Vr_p(*,ira,ip)=Vr1_min
                 u_p(ira,ip)=coeff_best(1)
                 v_p(ira,ip)=coeff_best(0)
;save VAD best wind and corresponding error
                 d01(8,ira,ip,ihr)=coeff_best(1) ;; u
                 d01(9,ira,ip,ihr)=coeff_best(0) ;; v
                 d01(10,ira,ip,ihr)=(coeff_best(0)^2.0+coeff_best(1)^2.0)^0.5
                 d01(26,ira,ip,ihr)=coeff_stdev_best(1)
                 d01(27,ira,ip,ihr)=coeff_stdev_best(0)
                 d01(28,ira,ip,ihr)=chi_min
                 
                 jump26:     ;;; If there aren't at least 5 azimuths with some data, skipped whole thing
              endif else begin ;;; END if VAD7 7 range pool mean wind VAD step
                 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; At this stage, we have a mean wind guess profile (the 7 range pool
;;; is done). We are inside the altitude and range and range
;;; pool loops, and now gonna adjust each azimuth to align with the
;;; mean wind guess. No need to loop, do as an array operation 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                 
; here's the mean wind, uu and vv, from the dataset d01 as built so far
                 uu=d01(8,7,ip,ihr)
                 vv=d01(9,7,ip,ihr)
                 
;;; IF THERE IS AN UNFOLDING OVERRIDE, USE THAT
                 if( n_elements(FORCE_UNFOLD_U) eq np and $
                     n_elements(FORCE_UNFOLD_V) eq np) then begin
                    uu=FORCE_UNFOLD_U(ip)
                    vv=FORCE_UNFOLD_V(ip)
                    d01(8,7,ip,ihr) = FORCE_UNFOLD_U(ip)
                    d01(9,7,ip,ihr) = FORCE_UNFOLD_V(ip)
                 endif
                 
;;; Adjust Vr to be within nyquist of the uu,vv wind                  
                 Vr_guess1 = (uu*costheta(*,iz)+vv*sintheta(*,iz)) $
                             *coszen(ira,iz)                  
;;; If folded, unfold
                 temp=(Vr_guess1 - Vr_p(*,ira,ip))/(2.0*vnyq)
                 isfolded = (abs(temp) ge 0.5 and Vr_p(*,ira,ip) gt -9990.)
                 Vr_p(*,ira,ip) = Vr_p(*,ira,ip) + $
                                  2*vnyq*fix(temp+0.5*sign(temp)) *isfolded
                 
              endelse ;;; if irpool not eaual to 3 (i.e., it's a divergence VAD)
              
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;                        
;;; project Vr_p, now all unfolded, onto horizontal direction
              for iaz=0,naz-1 do begin
                 if(Vr_p(iaz,ira,ip) gt -9990.0) then begin
                    Vr_p(iaz,ira,ip) = Vr_p(iaz,ira,ip)/coszen(ira,iz)
                 endif
              endfor
              
           endfor ;;; ira loop
        endfor    ;;; ip loop
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; OK, the unfolding is done. Now the VAD calculation 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
;;; Dive back into pressure (ip) and range (ira) loops...
        
        for ip=0,np-1 do begin ;; for all pressure levels
           iz=ip2iz(ip)        ;; this is a CLOSEST height, used for geometry and Vterm
           
           for ira=nrpool_half,nr-1-nrpool_half do begin ;;; for all ranges
              
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; grab a temporary copy of Vr vs. azimuth and compute how much of
;;; the azimuth circle has data and how much is missing data. there is
;;; a gapmax (the largest angular gap) and an angle_obs (0 to
;;; 360). Output them as quality flags, and if they are good enough
;;; then go ahead and do the VAD.
              
              gapmax=360.0  ;; initial, worst case
              
              Vr1=Vr_p(*,ira,ip)
              where_Vr_is_good=where(Vr1 gt -9990.0)
              
;;; Surface clutter is characterized by too narrow a histogram of vr,
;;; try killing it this way
;if (ip lt 2) then where_Vr_is_good=where(Vr1 gt -9990.0 and Vr_std(*,ira,ip) ge dvr/2)
              
              n=n_elements(where_Vr_is_good)
              
;;; Stats on missing data
              angle_obs=n*daz
              if (n ge 2) then gapmax = daz * max( deriv( $
                                        [where_Vr_is_good,where_Vr_is_good +naz] ) )
;; save these as qualtiy control flags for later analyses
              if(irpool lt 3) then begin 
                 d01(m+9,ira,ip,ihr)=angle_obs
                 d01(m+12,ira,ip,ihr)=gapmax
              endif
              
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;            
;;; VAD harmonic fitting crashes if there aren't at least 4 datapoints
;;; to fit the 3 parameters (speed, direction, divergence). SO, don't
;;; bother if n < 4.
              
              if(n ge 4) then begin
                 x=az(where_Vr_is_good)*PI_180
                 y=Vr1(where_Vr_is_good)
                 input_error_std=Vr_std(where_Vr_is_good,ira,ip)
                 
;; These "guesses" are starting points for a SVD fitting of a harmonic
;; function (the VAD_1harm). The fit converges, so the guess shouldn't
;; affect the outcome, just the rapidity of convergence. Start with 0 for
;; irpool=3, when we have no wind guess. Afterward, use the VAD wind
;; guess. When that is -9999, there is no data anyway, so no worries. 
                 
                 u_guess1=0.
                 v_guess1=0.
                 if(irpool le 2) then begin
                    u_guess1=d01(8,7,ip,ihr)
                    v_guess1=d01(9,7,ip,ihr)
                 endif
                 
;;; THE VAD STEP
                 coeff1=VAD_1harm (x, y, input_error_std, $
                                   u_guess1, v_guess1, fitted_output, coeff_output_stderr1, $
                                   chisq_goodness_of_fit)
                 coeff=coeff1([2,1,0])
                 coeff_output_stderr=coeff_output_stderr1([2,1,0])
                 
;;; If it worked, save the results
                 
                 if(n_elements(coeff) gt 2) then begin
                    div=2./(range_km(ira)*1000.0) * coeff(2)
                    std_div=2./(range_km(ira)*1000.0) * coeff_output_stderr(2)
                    if(irpool lt 3) then begin
                       d01(m,ira,ip,ihr)=div
                       d01(m+3,ira,ip,ihr)=std_div
                       d01(m+6,ira,ip,ihr)=chisq_goodness_of_fit
                    endif
                    yfit=coeff(0)*cos(az*PI_180)+coeff(1)*sin(az*PI_180)+coeff(2)
                    Vr_fit(*,ira,ip)=yfit
                    div_fit(ira,ip)=coeff(2)
                    div_std(ira,ip)=coeff_output_stderr(2)
                    
                 endif          ; if coeff has >2 elements (a result was returned)
              endif             ; if n>4, a requirement for not crashing the VAD harmfit procedure
           endfor               ; ira loop
        endfor                  ; ip loop
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Done with VAD for all ranges and pressures. Note that we are still
;;; inside the range POOL loop. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Rest of code: plot outputs, and save them to files
        
        print, 'IRPOOL ',irpool
        ;if (debug_unfolding eq 2) then stop ;Arunchandra (removed stop)
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Plotting code: pages and pages. 1. postage stamps. 2. Hourly
;;; summary.
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; 1. IF irpool =1, plot little VAD postage stamps, many to a page
        
        if(plot_output_device gt 0 and irpool eq 1 and onepageonly ne 1) then begin
           if(BW_plots_flag eq 0) then loadct,39
           if(BW_plots_flag eq 1) then begin
              loadct,1,file='~bem/colors1.tbl'
              !p.color = 255
           endif         
           
;;;; Set page size
           !p.position=[1,1,1,1]*0.1
           
;;;;;;;;;;;;;;;;;;;;;; ;;;;;;;;;;;;;;;;;;;;;; ;;;;;;;;;;;;;;;;;;;;;; 
;;; Each postage stamp has an outer box, and other info packed into it
;;; Set gross outer plot ranges         
           xr=[0,360] ;; degrees az range
           dy0=0.499    ;; How much room at top for little sideways histograms
           yr00=20.0    ;;; m/s Vr range
;;; Outer box size
           yr=[-1,1+dy0]*yr00
           
;;; Some character size adjusters
           siz1=0.6 ;;; a character size
           siz2=0.6 ;;; a character size
           
;;; WHICH R and P LEVELS TO PLOT -- set in params file now
           !p.multi=pmulti
           
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; LOOP  over ip and ira for postage-stamp VAD plots
           
           ipanel=-1 ;;; a counter of panels: only 1 gets the page title
           for iip=0, n_elements(plot_p)-1 do begin
              ip =plot_p[iip]
              iz =ip2iz[ip]
              for ira=plot_ra1,plot_ra2,plot_dra do begin
                 
                 ipanel=ipanel+1
                 xticn=strarr(20)
;               xticn(*)=' ' ;; suppress tic marks on most postage stamps
                 xticn(*)=''
;;; On bottom plot, use x ticks
                 if(iip eq n_elements(plot_p)-1) then xticn(*)=''
                 
;;; PLOT outer box (no data yet)
                 
                 plot,xr,yr,/nodata, $
                      xstyle=1,xticks=6,xtickname=xticn,ystyle=1, $
                      charsize=siz1, ytit='m/s'
                 if(p_or_z_coord eq 0) then begin
                    tit='P='+string(p(ip),format='(i4)')+'mb, '+ $
                        'R='+string(r(ira),format='(i2)')+'km'
                 endif
                 
;;; ADD TITLE (range, p) 
                 dyr=yr(1)-yr(0)
                 xyouts,xr(0),yr(1)+dyr*0.03,tit,size=siz2
                 
;;; PAGE TITLE
                 if(ipanel eq 0) then begin
                    tit=exp_name+'  hour '+hour + '  '+$
                        string(nrpool,format='(i1)')+'-range pooled data w/ dr='+ $
                        string(dr/1000,format='(i2)')+'km'
                    xyouts,xr(0),yr(1)+dyr*0.25,tit,size=siz2*1.2
                 endif
                 
;;; Horizontal lines
                 oplot,xr,[0,0];commented
                 oplot,xr,[0,0]+yr00;commented
                 
;;;;;;;;;;  ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;;
;;; Plot little diamonds and error whiskers where obs are numerous enough               
                 temp=vr_p(*,ira,ip)
                 where_Vr_is_good=where(temp gt -9990.0)
                 if(where_Vr_is_good(0) lt 0) then goto,jump11
                 n=n_elements(where_Vr_is_good)
                 x=az(where_Vr_is_good)
                 y=Vr_p(where_Vr_is_good,ira,ip)
                 dy=Vr_std(where_Vr_is_good,ira,ip)
;;; Diamonds
                 oplot,x,y,psym=4,symsize=0.5 ;commented
;;; Whiskers
                 for i=0,n-1 do begin
                    oplot,x(i)+[0,0],y(i)+[-1,1]*dy(i);commented
                 endfor
                 
;;;;;;;;;;  ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;;
;;; Overplot VAD fitted results
                 if(n gt 4) then begin
;;; Cyan sinusoid
                    oplot,az,Vr_fit(*,ira,ip),color=96,line=2 
;;; Red div lines
                    oplot,xr,[0,0]+div_fit(ira,ip),color=250,line=2
                    oplot,xr,[0,0]+div_fit(ira,ip)+div_std(ira,ip),color=250,line=1, thick=0.5
                    oplot,xr,[0,0]+div_fit(ira,ip)-div_std(ira,ip),color=250,line=1, thick=0.5
                 endif
                 
;;;;;;;;;;  ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;;
;;;; Spectral widths are orange boxes across the bottom
                 
                 wid1=wid_p(*,ira,ip)
                 for iaz=0,naz-1 do begin
                    if(wid1(iaz) gt -9990.0) then begin
                       w=daz
                       x=az(iaz)
                       h=wid1(iaz)
                       y=yr(0)+h*0.5
;;; Orange boxes
                       box,x,y,w,h,200
                       
                    endif
                 endfor         ; iaz
                 
;;;;;;;;;;  ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;;
;;;; SDZ indications - may be a distraction on a Vr scale, 
;;; but let's plot them just to see what's there. on bottom plot only
;;; since sdz is not a function of height. 
                 
                 SDZ1=SDZ(*,ira)
                 for iaz=0,naz-1 do begin
                    if(SDZ1(iaz) gt 0 and iip eq n_elements(plot_p)-1) then begin
                       w=daz /2 ;; half the width of the orange width boxes
                       x=az(iaz)
;;; convective spikes only, subtract off a floor value. Since dbz
;;; units are being plotted on a Vr scale anyway, units are
;;; meaningless, it's just a flag of convectiony-ness. 
                       h=(SDZ1(iaz) - SDZ_CONV_THRESH)>0 
                       y=yr(0)+h*0.5
;;; Green boxes
                       box,x,y,w,h,128
                    endif
                 endfor   ;; iaz
                 
;;;;;;;;;;  ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;;
;;; Sideways red histograms (logarithmic)
                 
                 dy1=dy0*yr00
                 for iaz=0,naz-1 do begin
                    dx1=az(iaz)+daz*0.5
                    oplot,[0,0]+dx1,yr00+[0,dy1]
                    if(vr_p(iaz,ira,ip) gt -9990.0) then begin
                       temp=fvrhist_p(1:nVrbins,iaz,ira,ip)
                       temp=alog(temp)/alog(10.0)*5
                       temp1=where(temp lt 0)
                       if(temp1(0) ge 0) then temp(temp1)=0.0
                       oplot,dx1-temp,findgen(nVrbins)/(nVrbins-1.0)*dy1+yr00,psym=10, $
                             color=250
                    endif
                 endfor ;; little loop over iaz
                 
;;;;;;;;;;  ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;; ;;;;;;;;;;
;;; Prevailing VAD wind guess (Green sinusoid through sideways red histograms)
                 
                 uu=d01(8,7,ip,ihr)
                 vv=d01(9,7,ip,ihr)
                 if(uu gt -9990.0) then begin
                    Vr_guess1=uu*sin(az*PI_180)+vv*cos(az*PI_180)
                    Vr_guess1=Vr_guess1*coszen(ira,iz)
                    temp=dy0*yr00/(2.0*vnyq)
;                  oplot,az,yr00*(1+0.5*dy0)+Vr_guess1*temp,color=125
;; better on main plot - green +/- Vnyq swath
                    oplot,az,Vr_guess1+vnyq,color=125
                    oplot,az,Vr_guess1-vnyq,color=125
                 endif
                 
jump11:
              endfor ;; ira loop for plotting
           endfor    ;; ip loop for plotting
        endif        ;; IF plotting to screen or postscript and not just one page
        
     endfor ;;; END RANGE POOL LOOP, FINALLY
     
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;HUGE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;RPOOL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;LOOP
     
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Plot code continues: 
; Making other page of hourly postscript: VAD profiles & echo stats
;;; Plotting code page 2   
     
     if(plot_output_device gt 0) then begin      
        !p.position=[1,1,1,1]*0.1
        !p.multi=[0,2,3,0,0]  ;;; 6 plots to a page
;;; char sizes for plot
        siz1=2.  ;; 2.4
        siz3=1.2 ;;1.2
        
;;;;;;;;;;;;;;;;       ;;;;;;;;;;;;;;;;       ;;;;;;;;;;;;;;;;       
;;;;; Plot echo cover fraction, overlaid with rainrate contours
        
;;; cover at plotalt level. Where is that? 
        
        plotalt = iz2ip[ where(abs(z-2) eq min(abs(z-2))) ]
        covermap_plotalt = total(fdbzhist_p(1:*,*,*,plotalt),1)*100 < 90
;; cap at 90% so color stays red rather than wrapping to white
        
;;; rain rate - map, and averages over circular areas of different sizes
        rain = rrate_GATEZR(2:*,*,*)*3600. *fdbzhist_p(2:*,*,*,plotalt) ;;; mm/hr, by dbz bin
        rainrate = total(rain(*,*,*),1)                                 ;;; a 2D (az, range) map of rainrate mm/hr
        
     
        
;;; arearain is an array of mean rainrates as a function of radius
        arearain =   total( (rainrate*area2d)[*,0] ) / $
                     total( (         area2d)[*,0] )    ;;; mm/hr within radius 1
        for iira= 1, n_elements(r)-1 do $
           arearain = [arearain, $
                       total( (rainrate*area2d)[*,0:iira] ) / $
                       total( (         area2d)[*,0:iira] )]   ;;; mm/hr, avg. over area inside each radius
;;; Save this off in the dataset d01, at the appropriate pressure level 
        d01(5,*,plotalt,ihr)=arearain ;;; mm/h
        
;;;;;;;;;;;;;;;;; BEM 2009 July
;;; stratrain and convrain are arrays of rainrates as a function of
;;; radius. These come from dividing the sdbzhist across the 4:1 line
        
        STD_DBZ = indgen(36,66) mod 36
        ABS_DBZ = indgen(36,66) /36
        strat_mask = ( ABS_DBZ gt STD_DBZ*4 ) ;;; above the 4:1 line
        conv_mask =  ( ABS_DBZ le STD_DBZ*4 ) ;;; below the 4:1 line
        strat_dbzs = total( fsdzhist*strat_mask, 1) ;; result is a 66, 24, 12 array (dbz,az,range)
        conv_dbzs = total(  fsdzhist* conv_mask, 1) ;; result is a 66, 24, 12 array (dbz,az,range)
;;; Convert dbZ to rainrate
        stratrain = rrate_GATEZR(2:*,*,*)*strat_dbzs(2:*,*,*) 
        convrain = rrate_GATEZR(2:*,*,*)*conv_dbzs(2:*,*,*) 
        stratrainrate = total(stratrain(*,*,*),1) ;;; sum over dBZs, a 2D (az, range) map of rainrate        
        convrainrate = total( convrain(*,*,*),1)  ;;; sum over dBZs, a 2D (az, range) map of rainrate 
        
;;; Now average over area within each ccircle of radius iira
        for iira= 0, n_elements(r)-1 do $
           srain =   total( (stratrainrate*area2d)[*,0:iira] ) / $
           total( (         area2d)[*,0:iira] )    ;;; mm/hr within radii
        for iira= 0, n_elements(r)-1 do $
           crain =   total( ( convrainrate*area2d)[*,0:iira] ) / $
           total( (         area2d)[*,0:iira] )    ;;; mm/hr within radii
        
;;; Now ignore the absolute values and just redistribute C-S fraction
        sfrac = srain/( (srain+crain)>1e-10)
;;; Save this off in the dataset d01, at the appropriate pressure level 
        d01(29,*,plotalt,ihr)=sfrac ;;; mm/h
        
        
        temp = where(covermap_plotalt le 0.001)
        if (max(temp) ge 0) then rainrate(temp) = 0.0
        
        if(BW_plots_flag eq 0) then coverlev=findgen(20)*5 ;/5.
        if(BW_plots_flag eq 1) then coverlev=findgen(6)*25
        
;;; size of box (km)
        temp0 = max(range_km)
        xr=[-1,1]*temp0
        yr=xr
        
;;; Contours (filled) color = cover
        contour, [covermap_plotalt,covermap_plotalt(0,*)],xx,yy,$
                 levels = coverlev,/fill, $
                 xrange=xr,yrange=yr,xstyle=1,ystyle=1, $
                 xtit='east (km)',ytit='north (km)', $
                 charsize=siz1
                 
                 
        
;;; Overlay contours of rainrate - too many?? White
;      lev2= (0.015+findgen(10)*0.05) ;; adjusted for looks
;      contour, /overplot, ([rainrate,rainrate(0,*)]),xx,yy,$
;        c_thick=(findgen(10)+1)*0.8,levels=lev2, color=255
        
;;; Overlay contours of SDZ (convectiveness)
        lev2= SDZ_CONV_THRESH +indgen(23)
        contour, /overplot, ([sdz,sdz(0,*)]),xx,yy,$
                 c_thick=(findgen(23)+1)*0.8,levels=lev2, color=255
        
;;; Title      
        dyr=yr(1)-yr(0)
;      tit='cover & R(Z) at '+str(fix(p(plotalt)))
        tit='cover at '+str(fix(p(plotalt)))+', sDBZ'
        xyouts,xr(0),yr(1)+dyr*0.03,tit,size=siz3
        
;;;;; oplot maximum range circle
        temp=(findgen(temp0*2)-temp0)/(temp0+0.0)*!pi
        oplot,temp0*sin(temp),temp0*cos(temp)
        
;print, coverlev
        
;;; color key
      ;;;MAKE_KEY,COLORS=findgen(n_elements(coverlev))*!D.N_colors/n_elements(coverlev),$
        ;;;labels = round_any([min(coverlev),(min(coverlev)+max(coverlev))/2.+mean(deriv(coverlev))/2.],sig=1),$
        ;;;units = ' %',/orientation
        
;;; PAGE TITLE
        tit=exp_name+'  hour '+hour +'; ZR-rain(r<radii):' +string(arearain[radii_to_plot], $
                                                                   format='(4f6.2)')+' mm/h'
        xyouts,xr(0),yr(1)+dyr*0.25,tit,size=0.6*1.2
        
;;;;;;;;;;;;;;    ;;;;;;;;;;;;;;     ;;;;;;;;;;;;;;       ;;;;;;;;;;;;;;              
;;;;;;; CFAD of echo: contours fixed
        if(BW_plots_flag eq 0) then lev=findgen(20)*0.5-8 ;; log10 of cfad
        if(BW_plots_flag eq 0) then lev=findgen(16)*0.5-6 ;; log10 of cfad
        if(BW_plots_flag eq 1) then lev=findgen(6)*2-8
        
;;; set plot ranges
        xr=[0,max(dbz)]
        if(p_or_z_coord eq 0) then begin
           
           yr=[max(p)-dp/2,min(p)+dp/2] ;[1000,100]
           ytit='P (mb)'
        endif
        if(p_or_z_coord eq 1) then begin
           yr=[0,18]
           ytit='Z (km)'
        endif
        
;;; Title set, contour CFAD
        tit='CFAD r<' + string(max(range(CFAD_RANGE_BINS2))/1000,format='(i2)') +'km'
        contour, alog(cfad_p*100.0>1e-30)/alog(10.0), r03, p, $
                 levels=lev, /fill, $
                 xrange=xr,xstyle=1,yrange=yr,ystyle=1, $
                 xtit='nonzero echo (dBZe)', $ ;,ytit=ytit, $
                 charsize=siz1
        
        
;;; TITLE      
        dyr=yr(1)-yr(0)
        xyouts,xr(0),yr(1)+dyr*0.03,tit,size=siz3
        
;;; COLOR KEY
      ;;;lab= $
        ;;;round_any([min(lev),(min(lev)+max(lev))/2.+mean(deriv(lev))/2.],sig=1)
      ;;;lab=['10!U-6!N','10!U-2!N']
      ;;;MAKE_KEY,COLORS=findgen(n_elements(lev))*!D.N_colors/n_elements(lev),$
        ;;;labels=lab,units=' %',/orientation 
        
;;;;;;;;;;;;;;;;;;;;;;;;  ;;;;;;;;;;;;;;;;;;;;;;;;  ;;;;;;;;;;;;;;;;;;;;;;;;  
;;; Divergence profile plots
        if(p_or_z_coord eq 0) then begin
           yr=[max(p)-dp/2,min(p)+dp/2] ;[1000,100]
           ytit='P (mb)'
        endif
        if(p_or_z_coord eq 1) then begin
           yr=[0,18]
           ytit='Z (km)'
           dy1=1.0
           y10=16.0
        endif
        zplot=p
        
;;; legend locations for 4 lines 
        y10= min(p) + (max(p)-min(p))/10
        dy1= (max(p)-min(p))/15
        
        ctit0=['Divergence','Div (3R pools)', $
               'Div (5R pools)']
        xr=[-1,1]*3.0
        x1=xr(0)*0.9
        dx1=-xr(0)*0.25
        
        m1=11 ;;; Dataset codes
        m2=13
        
;;; m is the data array code for range pools 1,3,5. LOOP
        for m=m1,m2 do begin
           plot,zplot,zplot,xtit='Div (10!U-4!N s!U-1!N)',ytit=ytit,/nodata,$
                xra=xr,xstyle=1,yrange=yr,ystyle=1, $
                charsize=siz1
;;; TITLE
           tit=ctit0(m-m1)
           xyouts,xr(0),yr(1)+dyr*0.03,tit,size=siz3
           
;;; legend
           xyouts,x1+dx1*0.7,y10,'Radius (km)',size=0.6
           
;;; Bar plots of good values at 4 ranges         
           if(BW_plots_flag eq 0) then begin
              
;;; The 4 ranges selcted in parameters
;            for ira=3,9,2 do begin
              for iira= 0, n_elements( radii_to_plot )-1 do begin
                 ira = radii_to_plot[iira]
                 col=(ira-3)*250./6 ;;; colors
                 lin=0              ;;; line styles
                 thi=1              ;;; thicknesses
                 
;;; Grab div values
                 div=d01(m,ira,*,ihr)
                 divstd=d01(m+3,ira,*,ihr)
                 index=where(div gt -9990.0 and divstd lt 999)
                 if(index(0) ge 0) then begin
                    n11=n_elements(index)
                    for i11=0,n11-1 do begin
                       j11=index(i11)
                       w11=div(j11)*1.0E+4
                       h11=dp*0.5
                       x11=w11*0.5
                       y11=zplot(j11)
;;; plot the good ones
;                     box2,x11,y11,w11,h11,col
;;; use plots
                       xxx = [0,w11,w11,0] < max(xr) > min(xr)
                       yyy = [y11-h11,y11-h11,y11+h11,y11+h11]
                       plots, xxx,yyy, color = col
;;; Whiskers for div_std, offset a skinch vertically for seeing
                       oplot, [(div(j11)-divstd(j11))*1.0E+4, $
                               (div(j11)+divstd(j11))*1.0E+4], $
                              [y11,y11] + (iira-1.5) *min(deriv(zplot))/10, color=col
                       
                    endfor      ; i11 (loop over good data levels)
                 endif          ; if there are good levels
                 
;;; Legend for 4 lines               
                 y1=y10+dy1*(ira-1)*0.5
                 oplot,[0,dx1]+x1,[0,0]+y1,color=col,line=lin,thick=thi
                 xyouts,x1+dx1*1.2,y1,string(r(ira),format='(f4.1)'),size=0.6
                 
              endfor ;; iira radii to plot
           endif 
;;; Ver, 0
           oplot,[0,0],yr
           
        endfor  ;; m loop (over range pools 1,3,5 for the 3 div profile plots)
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;      
;;; VAD Wind profiles
        
;;; xrange
        xr=[-1,1]*25.0
        
;;; Wind profile plot
        tit='big VAD u & v'
        plot, zplot, zplot, /nodata, $
              xtit='u,v (m s!U-1!N)',ytit=ytit,$
              xra=xr,xstyle=1, $
              yrange=yr,ystyle=1, $
              charsize=siz1
        dyr=yr(1)-yr(0)
        xyouts,xr(0),yr(1)+dyr*0.03,tit,size=siz3
;;; ver,0      
        oplot,[0,0],yr
        
;;; Grab VAD wind profiles and plot them
        u1=d01(8,7,*,ihr)
        v1=d01(9,7,*,ihr)
        stdev_u = d01(26,7,*,ihr)
        stdev_v = d01(27,7,*,ihr)
        index=where(u1 gt -9990.0)
        
        if(index(0) ge 0) then begin
           oplot,u1(index),zplot(index),thick=4,color=30
           oplot,v1(index),zplot(index),thick=4,color=250,lines=2
        endif
        
;;; Whiskers for wind_std
        
        for i=0,n_elements(index)-1 do begin
           oplot, [(u1-stdev_u)[index[i]], $
                   (u1+stdev_u)[index[i]]], $
                  [zplot[index[i]], zplot[index[i]] ], color=30
           
           oplot, [(v1-stdev_v)[index[i]], $
                   (v1+stdev_v)[index[i]]], $
                  [zplot[index[i]], zplot[index[i]] ], color=250
        endfor
        
        if (onepageonly ne 1) then begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; MORE PLOTS -- cross sections W-E and N-S of mean reflecctivity and Vr
           dbz_hist_values = (findgen(ndbz,naz,nr,nz) mod NDBZ)
           meandbz = total( fdbzhist(2:*,*,*,*) *dbz_hist_values(2:*,*,*,*) ,1) / $
                     (total( fdbzhist(2:*,*,*,*) ,1)  >1) ;; >1 avoids 0/0 error
;;; max dbz 
           maxdbz = meandbz*0
           FOR iz=0,nz-1 do $
              for ira=0,nr-1 do $
                 for iaz=0,naz-1 do $
                    maxdbz(iaz,ira,iz) = max( where( dbzhist(*,iaz,ira,iz) gt 0) )
           
;;; Raw doppler velocity - divide Vrsum by nobsVr
           nobsVr=total(vrhist(1:nVrbins,*,*,*) ,1)
           meanVr =Vrsum / (nobsVr > 1)
           
;;; Azes and levels and labels
           raxis = (findgen(nr*2) - nr +0.5) *DR/1000.
           yaxis = z
           tit=exp_name+'  hour '+hour
           
           !p.multi=[0,1,4]
;;; WE section of dbz - meean or max. How about both -naww just max
           WExsec = fltarr(nr*2,nz)
           WExsec(nr:*,*)   =                meandbz(  naz/4,*,*)
           WExsec(0:nr-1,*) = rotate( reform(meandbz(3*naz/4,*,*)) ,5)
;;;CONTOUR dbz
           lev = (indgen(ndbz))
                                ; pfh lev = (indgen(ndbz))(0:59)
;      contour, WExsec, raxis, yaxis, xtit='W-E distance (km)', ytit='z (km)', $
;        /fill, levels=lev, tit='W-E cross section of mean dBZ ', subtitle = tit
;;; color key
;      MAKE_KEY,COLORS=findgen(n_elements(lev))*!D.N_colors/n_elements(lev),$
;        labels = round_any([min(lev),(min(lev)+max(lev))/2.+mean(deriv(lev))/2.],sig=2),$
;        units = 'dBZ',/orientation
           
;;; max dbz
           WExsec(nr:*,*)   =                maxdbz(  naz/4,*,*)
           WExsec(0:nr-1,*) = rotate( reform(maxdbz(3*naz/4,*,*)) ,5)
;;;CONTOUR dbz
           contour, WExsec, raxis, yaxis, xtit='W-E distance (km)', ytit='z (km)', $
                    /fill, levels=lev, tit='W-E section of max dBZ; curve is sDZ index', subtitle = tit
;;; color key
      ;;;MAKE_KEY,COLORS=findgen(n_elements(lev))*!D.N_colors/n_elements(lev),$
        ;;;labels = round_any([min(lev),(min(lev)+max(lev))/2.+mean(deriv(lev))/2.],sig=2),$
        ;;;units = 'dBZ',/orientation
;;; OVERPLOT SDZ>10 for convective indicator
           xsec_series = fltarr(nr*2)
           xsec_series(nr:*,*)   = sdz(  naz/4,*)
           xsec_series(0:nr-1,*) = rotate(sdz(3*naz/4,*), 1)
           oplot, raxis, yaxis(0) +(xsec_series-SDZ_CONV_THRESH >0), $
                  thick=5,lines=2,color=255 ;;; dbZ-1 in km units - should look OK
           oplot, raxis, yaxis(0) +(xsec_series-SDZ_CONV_THRESH >0), $
                  thick=5,lines=1 ;;; dbZ-1 in km units - should look OK
           
;;; WE section of Vr (raw) and dvr/dx
           levvr = -fix(vnyq) + findgen(2*fix(vnyq)+1)
           WExsec(nr:*,*)   =                meanvr(  naz/4,*,*)
           WExsec(0:nr-1,*) = rotate( reform(meanvr(3*naz/4,*,*)) ,5)
;;;CONTOUR Vr
           contour, WExsec, raxis, yaxis, xtit='W-E distance (km)', ytit='z (km)', $
                    /fill, levels=levvr, tit='W-E cross section of mean Vr (NOT UNFOLDEwD)' ; and d/dx'
;      dvrdx = (shift(WExsec,-1)-shift(WExsec,1))
;      contour, dvrdx, raxis, yaxis, /over, levels = indgen(vnyq)+1
;      contour, dvrdx, raxis, yaxis, /over, levels = indgen(vnyq)-vnyq, c_lines=1
;;; color key
      ;;;MAKE_KEY,COLORS=findgen(n_elements(levvr))*!D.N_colors/n_elements(levvr),$
       ;;; labels = round_any([min(levvr),(min(levvr)+max(levvr))/2.],sig=2),$
       ;;; units = 'm/s',/orientation
           
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; NS section
           
;;; mean dbz
           NSxsec = fltarr(nr*2,nz)
           NSxsec(nr:*,*)   =                meandbz(    0,*,*)
           NSxsec(0:nr-1,*) = rotate( reform(meandbz(naz/2,*,*)) ,5)
;;; CONTOUR
;      contour, NSxsec, raxis, yaxis, xtit='S-N distance (km)', ytit='z (km)', $
;        /fill, levels=lev, tit='S-N cross section of mean dBZ '
;;; color key
;      MAKE_KEY,COLORS=findgen(n_elements(lev))*!D.N_colors/n_elements(lev),$
;        labels = round_any([min(lev),(min(lev)+max(lev))/2.+mean(deriv(lev))/2.],sig=2),$
;        units = 'dBZ',/orientation
           
;;; max
           NSxsec = fltarr(nr*2,nz)
           NSxsec(nr:*,*)   =                maxdbz(    0,*,*)
           NSxsec(0:nr-1,*) = rotate( reform(maxdbz(naz/2,*,*)) ,5)
;;; CONTOUR
           contour, NSxsec, raxis, yaxis, xtit='S-N distance (km)', ytit='z (km)', $
                    /fill, levels=lev, tit='S-N section of max dBZ; curve is sDZ index'
;;; color key
      ;;;MAKE_KEY,COLORS=findgen(n_elements(lev))*!D.N_colors/n_elements(lev),$
        ;;;labels = round_any([min(lev),(min(lev)+max(lev))/2.+mean(deriv(lev))/2.],sig=2),$
       ;;; units = 'dBZ',/orientation
;;; OVERPLOT SDZ>10 for convective indicator
           xsec_series = fltarr(nr*2)
           xsec_series(nr:*,*)   = sdz(    0,*)
           xsec_series(0:nr-1,*) = rotate(sdz(naz/2,*),1)
           oplot, raxis, yaxis(0) +(xsec_series-SDZ_CONV_THRESH >0), $
                  thick=5,lines=2,color=255 ;;; dbZ-1 in km units - should look OK
           oplot, raxis, yaxis(0) +(xsec_series-SDZ_CONV_THRESH >0), $
                  thick=5,lines=1 ;;; dbZ-1 in km units - should look OK
           
;;; NS section of Vr (raw)
           NSxsec(nr:*,*)   =                meanvr(    0,*,*)
           NSxsec(0:nr-1,*) = rotate( reform(meanvr(naz/2,*,*)) ,5)
;;;CONTOUR Vr and dvr/dx
           contour, NSxsec, raxis, yaxis, xtit='S-N distance (km)', ytit='z (km)', $
                    /fill, levels=levvr, tit='S-N cross section of mean Vr (NOT UNFOLDED)' ; and d/dy'
;      dvrdx = (shift(NSxsec,-1)-shift(NSxsec,1))
;      contour, dvrdx, raxis, yaxis, /over, levels = indgen(vnyq)+1
;      contour, dvrdx, raxis, yaxis, /over, levels = indgen(vnyq)-vnyq, c_lines=1
;;; color key
      ;;;MAKE_KEY,COLORS=findgen(n_elements(levvr))*!D.N_colors/n_elements(levvr),$
        ;;;labels = round_any([min(levvr),(min(levvr)+max(levvr))/2.],sig=2),$
        ;;;units = 'm/s',/orientation
           
;;; XXX onepageonly should exclude this plot
;;; Try contours of SDZ (convectiveness), and joint histogram SDZ-DBZ
           !p.multi=[0,1,2]
;;; Joint histogram
           jhist = total(total( sdzhist, 4),3)
           contour, jhist(1:*,2:*), sdz_values(1:*), dbz(2:*), $
                    xtit='dBz gradient', ytit='dBZ', /fill, nlev=23, tit='SDZ - dbZ joint histogram, 2-4km layer'
           oplot, [0,20], [0,80] ;;; suggested line for strat-conv
;;; Joint histogram (short ranges only)
           jhist = total(total( sdzhist(*,*,*,0:nr/2), 4),3)
           contour, jhist(1:*,2:*), sdz_values(1:*), dbz(2:*), $
                    xtit='dBz gradient', ytit='dBZ', /fill, nlev=23, tit='joint histogram inside r< rmax/2'
           oplot, [0,20], [0,80] ;;; suggested line for strat-conv
           
        endif   ;; if (onepageonly ne 1)
     endif       ;; If plotting is wanted
     
     
; Digital data output of results!!!
;;; Grab results for file output
     
     rain_plotalt_GATEZR(ihr,*)=d01(5,*,3,ihr)
     stratiform_fraction(ihr,*)=d01(29,*,3,ihr)
     rain_plotalt_localZR(ihr,*)=d01(6,*,3,ihr)
     cfad0(*,*,ihr)=cfad_p
     divergence(ihr,*,*)=d01(12,*,*,ihr) ;;; rpool = 1 results
     div_stdev(ihr,*,*)=d01(15,*,*,ihr)
     chi_div(ihr,*,*)=d01(18,*,*,ihr)
     uwind(ihr,*)=d01(8,7,*,ihr)
     vwind(ihr,*)=d01(9,7,*,ihr)
     u_std(ihr,*)=d01(26,7,*,ihr)
     v_std(ihr,*)=d01(27,7,*,ihr)
     chi_uv(ihr,*)=d01(28,7,*,ihr)
     
   

  endfor ;;; ihr loop over hours in the dataset
  
  
  ;;; IDL .sav file savefile output
  file1=exp_name+'.sav'

  ;;;;;;;; Commenting saving the file until the code debugging ends;;;;;;;;

  save, file=file1, $
    timestrings,cover_0db,cover_15db,cover_30db,rain_plotalt_GATEZR,rain_plotalt_localZR,$
    stratiform_fraction, $
    divergence,div_stdev,chi_div,uwind,vwind,u_std,v_std,chi_uv, cfad0, r, p, arearain, d01, v01
    
    ;;;;;;;; Commenting saving the output variables into netcdf file ;;;;;;;;
  
  output_flag = netcdf_generation (timestrings, cover_0db, cover_15db, cover_30db, rain_plotalt_GATEZR, rain_plotalt_localZR, $
  stratiform_fraction, divergence, div_stdev, chi_div, uwind, vwind, u_std, v_std, chi_uv, cfad0, r, p, arearain, d01, v01)
  
  stop
  
;;;;;;;;;;;;;;;;;;;;;; DONE WITH ALL PROCESSING OF ALL HOURS. 
;close postscript file
  
  if(plot_output_device eq 2) then begin
     device,/close
     set_plot,'x'
  endif
  
  file='radar_proc.out1_'+exp_name+'_'+qc_string(experimental_qc_flag) $
       +cguess(wind_firstguess)
;opt_grid,file,nv0,nr,np,nt,v01,r,p,t,d01
;spawn,'gzip -f '+file
  
  file='radar_proc.out3_'+exp_name+'_'+qc_string(experimental_qc_flag) $
       +cguess(wind_firstguess)
;opt_grid,file,nv03,nr03,nz,nt,v03,r03,z,t,d03
;spawn,'gzip -f '+file
  

        
  ;;;;;;;; Commenting saving the file until the code debugging ends;;;;;;;; 
  
  
;spawn,'gzip -f '+file1
  
  stop
  
end ;; END OF PROGRAM


