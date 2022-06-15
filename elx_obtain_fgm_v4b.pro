pro ELX_obtain_FGM_v4b
; 
; 2022-06-01 VA: This version corrects the rotation to +spinphase in snwsidt0=sin(+spinphase0) etc.
; 2022-06-01 VA: This version corrects the rotation back to 44 deg.
; 2022-03-26 VA: This version rotates by the correct angle (180+44) deg,
;     and can add an offset to the IGRF to represent the "real" field.
;     This offset can be constant, linear or quadratic with time. 
;
  elf_init
  tplot_options, 'xmargin', [17,12]
  tplot_options, 'ymargin', [2,3]
  cwdirname='C:\My Documents\ucla\Elfin\MAG\CAL_12_params' ; your directory here, if other than default IDL dir
  cwd,cwdirname
  pival=!PI
  Re=6378.1 ; Earth equatorial radius in km
;
; pick an event; here from 2021-12-24 11:40-12:10
;
; two known FGM collections (NA+ND) exist there: 11:48-11:55 (NA) and 12:01 - 12:07 (ND)
;
  tstart='2021-12-24/11:40:00' 
  tend='2021-12-24/12:10:00'
;  
  tstart1='2021-12-24/11:48:00'
  tend1='2021-12-24/11:55:00'
;
;
  tstart2='2021-12-24/12:01:00'
  tend2='2021-12-24/12:07:00'
;
  ; next is a 5nT max amplitude wave in X, 3nT in Z and 1nT in Y direction, arbitrary polarization, gaussian shape, 1st scizone
  Twav1=1.25;sec
  fwhm1=20. ;sec
  A1x=5.
  A1y=1.
  A1z=3.
  temic1=make_array(1201,/index,/double)/10. + time_double(tstart1)+3.*60.; 0.1s res., 1min dur., 3min into scizone
  temic1_cntr=(temic1[0]+temic1[-1])/2.
  emic1x=A1x*cos(-(2.*!PI/Twav1)*temic1)*exp(-((temic1-temic1_cntr)/fwhm1)^2)
  emic1y=A1y*cos(-(2.*!PI/Twav1)*temic1)*exp(-((temic1-temic1_cntr)/fwhm1)^2)
  emic1z=A1z*sin(-(2.*!PI/Twav1)*temic1)*exp(-((temic1-temic1_cntr)/fwhm1)^2)
  store_data,'emic1',data={x:temic1,y:[[emic1x],[emic1y],[emic1z]]}
  ;
  ; next is a 4nT max amplitude wave in X, 3.nT in Z and 0.5nT in Y direction, arbitrary polarization, gaussian shape, 2nd scizone
  Twav2=2.2d;sec
  fwhm2=15.d ;sec
  A2x=4.d
  A2y=0.5d
  A2z=3.d
  temic2=make_array(1201,/index,/double)/10. + time_double(tstart2)+2.*60.; 0.1s res., 2min dur., 2min into sci. zone
  temic2_cntr=(temic2[0]+temic2[-1])/2.
  emic2x=A1x*cos(-(2.*!PI/Twav2)*temic2)*exp(-((temic2-temic2_cntr)/fwhm2)^2)
  emic2y=A1y*cos(-(2.*!PI/Twav2)*temic2)*exp(-((temic2-temic2_cntr)/fwhm2)^2)
  emic2z=A1z*sin(-(2.*!PI/Twav2)*temic2)*exp(-((temic2-temic2_cntr)/fwhm2)^2)
  store_data,'emic2',data={x:temic2,y:[[emic2x],[emic2y],[emic2z]]}
  options,'emic?',colors=['b','g','r']
;
 time2plot=[tstart,tend]
 timeduration=time_double(tend)-time_double(tstart)
 timespan,tstart,timeduration,/seconds ; set the analysis time interval
;
sclet='a'
elf_load_state,probe=sclet ; get position and attitude info. (in the future will have spin rate too)
;
; choose the sidereal spin period, omega sidereal, or: wsid to roughly match
; the spin-period inferred from MRM measurements (listed on the orbit plots)
; Note that this cannot be determined exactly, only approximately, since there is no sun-sensor on board
; and the magnetic field changes along the orbit so zero-crossings change faster or slower than spin period.
; We measure the synodic spin period by counting the delay between zero-crossings of the radial B-field component.
; However, if dataset is only 1 or 2 spins, then B doesn't change much and the wsid can match wsyn within tolerance. 
;
Tsid = 2.84d; seconds (siderial spin period assumed same as synodic listed on the orbit plots)
wsid = 2.*!PI/Tsid
;
elf_load_fgm,probe=sclet ; just for information purposes
;
; establish a time array at 0.1s resolution (same as FGM) to associate the fake data with
; then interpolate position to that fake data time array, so you can find the IGRF B at that position
; then time_clip the final data to two intervals, interval 1 and 2, so you can write those in 2 data files
;
tresfakedata=0.1d ; seconds
nfaketimes=long(double(timeduration)/tresfakedata) + 1
tarray=make_array(nfaketimes,/double,/index)/10.+time_double(tstart) ; time in seconds since 1970
tinterpol_mxn,'el'+sclet+'_pos_gei',tarray,newname='elx_pos_gei'
tinterpol_mxn,'el'+sclet+'_att_gei',tarray,newname='elx_att_gei'
tnormalize,'elx_att_gei',newname='elx_att_gei'
;
; Find the B-field in various coordinates at the satellite location.
cotrans,'elx_pos_gei','elx_pos_gse',/GEI2GSE
cotrans,'elx_pos_gse','elx_pos_gsm',/GSE2GSM
tt89,'elx_pos_gsm',/igrf_only,newname='elx_bt89_gsm',period=0.1; gets IGRF field at ELF location
cotrans,'elx_bt89_gsm','elx_bt89_gse',/GSM2GSE ; cast Bfield into GSE->GEI coords for transforming to body coord's
cotrans,'elx_bt89_gse','elx_bt89_gei',/GSE2GEI ; cast Bfield into GSE->GEI coords for transforming to body coord's
;
; Find B-field in Despun Satellite coordinates (Z=along spin axis, Y=along Bperp vector on spin plane, X=YxZ)
; This is consistent with the DMXL coordinate system defined in the COORDSYS-ADCS-01-R-1_20150604_va.docx
; When the FGM radial component (along Y body axis) crosses peak B field,
;    ... then the X-axis in body coord's crosses zero field in the ascending direction.
; So the despun coordinates have Y axis along the GEI direction of the Bperp field.
; All computations in GEI. To find X unit vector you can just take the cross product of B into spin axis, then normalize.
; Then to find the Y unit vector you can just take ZxX and normalize.
; The assumption is that the IGRF field does not change over short periods during the orbit.
; So we can use instantaneous field (GEI) to determine the rotation matrix (it is same as the avg field over spin period).
; (Else we can just do a low-pass filter or running-average of the data over a 3-6 second window before determining rotmat.)
; It is also assumed that the satellite spin axis is along the Z body axis. This is not entirely true but it does not
;    matter because we dont care about the actual satellite coordinates, just the spin-axis, and spin-plane coordinates.
get_data,'elx_att_gei',data=elx_att_gei
Zaxis=double(elx_att_gei.y) ; this is the unit vector of the spin axis in GEI coord's
tnormalize,'elx_bt89_gei',newname='bhat'
get_data,'bhat',data=bhat
Bunit=double(bhat.y) ; this is the Bfield unit vector in GEI coord's 
Xaxis=crossp2(Bunit,Zaxis) ; this is the unit vector across B and the spin axis
mytot=sqrt(total(Xaxis^2,2))
for j=0,2 do Xaxis[*,j]=Xaxis[*,j]/mytot ; normalize vector
Yaxis=crossp2(Zaxis,Xaxis) ; this is the unit vector along Bperp (perp to spin axis).
npoints=n_elements(Xaxis[*,0])
DMXL2GEI=reform([[Xaxis],[Yaxis],[Zaxis]],npoints,3,3)
store_data,'DMXL2GEI',data={x:elx_att_gei.x,y:DMXL2GEI} ; the inverse of this one is GEI2DMXL
;tvector_rotate,'DMXL2GEI','elx_att_gei',/invert,newname = 'Zaxis_DMXL' ; just testing it gives [0,0,1] everywhere
tvector_rotate,'DMXL2GEI','elx_bt89_gei',/invert,newname = 'elx_bt89_dmxl' ; from GEI to DMXL
options,'elx_bt89_dmxl',labels=['X','Y','Z'],YSUBTITLE='elx_bt89_dmxl'
copy_data,'elx_bt89_dmxl','elx_bt89_dmxl_true'
calc,"'elx_bt89_dmxl_offset'='elx_bt89_dmxl_true'-'elx_bt89_dmxl_true'"
get_data,'elx_bt89_dmxl_offset',data=elx_bt89_dmxl_offset,dlim=mydlim,lim=mylim
;stop
elx_bt89_dmxl_offset_x0=elx_bt89_dmxl_offset.x ; these are destined to be relative time [0,...1] since the start of SZ...
iinSZ1=where(elx_bt89_dmxl_offset.x ge time_double(tstart1) and elx_bt89_dmxl_offset.x le time_double(tend1), jinSZ1)
iinSZ2=where(elx_bt89_dmxl_offset.x ge time_double(tstart2) and elx_bt89_dmxl_offset.x le time_double(tend2), jinSZ2)
elx_bt89_dmxl_offset_x0[*]=0.
if jinSZ1 gt 0 then elx_bt89_dmxl_offset_x0[iinSZ1]= $
  (elx_bt89_dmxl_offset.x[iinSZ1]-elx_bt89_dmxl_offset.x[iinSZ1[0]])/(elx_bt89_dmxl_offset.x[iinSZ1[-1]]-elx_bt89_dmxl_offset.x[iinSZ1[0]])
if jinSZ2 gt 0 then elx_bt89_dmxl_offset_x0[iinSZ2]= $
  (elx_bt89_dmxl_offset.x[iinSZ2]-elx_bt89_dmxl_offset.x[iinSZ2[0]])/(elx_bt89_dmxl_offset.x[iinSZ2[-1]]-elx_bt89_dmxl_offset.x[iinSZ2[0]])
;store_data,'elx_bt89_dmxl_offset_x0',data={x:elx_bt89_dmxl_offset.x,y:elx_bt89_dmxl_offset_x0}
elx_bt89_dmxl_offset_y0=0.  +elx_bt89_dmxl_offset_x0*0.           ; = 0. ; no offset in DMXL-X direction
elx_bt89_dmxl_offset_y1=-100.+(elx_bt89_dmxl_offset_x0-0.5)^2*400. ; quad. offset ; -100.+(elx_bt89_dmxl_offset_x0-0.5)^2*400. ; quad. offset ;+25.-elx_bt89_dmxl_offset_x0*100. ; lin. offset ; = -100. ;const. offset
elx_bt89_dmxl_offset_y2= 150.-(elx_bt89_dmxl_offset_x0-0.5)^2*600.  ; quad. offset;  150.-(elx_bt89_dmxl_offset_x0-0.5)^2*600.  ; quad. offset ;-50.+elx_bt89_dmxl_offset_x0*150. ; lin. offset; = +150. ;const. offset
store_data,'elx_bt89_dmxl_offset',data={x:elx_bt89_dmxl_offset.x, $
   y:[[elx_bt89_dmxl_offset_y0],[elx_bt89_dmxl_offset_y1],[elx_bt89_dmxl_offset_y2]]},dlim=mydlim,lim=mylim
;calc," 'elx_bt89_dmxl'='elx_bt89_dmxl'+'elx_bt89_dmxl_offset'"
;
; add the two waves
;
tinterpol_mxn,'emic1','elx_bt89_dmxl',newname='emic1_int',/repeat_extrapolate ; sets the extrapolated values to near-zero
tinterpol_mxn,'emic2','elx_bt89_dmxl',newname='emic2_int',/repeat_extrapolate ; sets the extrapolated values to near-zero
calc," 'emicwave' = 'emic1_int' + 'emic2_int' "
;
; add the waves to the data in despun coordinates (dmxl)
;
calc," 'elx_igrf_dmxl' ='elx_bt89_dmxl' "
calc," 'elx_bt89_dmxl' = 'elx_bt89_dmxl' + 'emicwave' " ; the wave is imperceptible unless the full field is subtracted
;
; Now spin up the data using the sidereal spin period in a fixed coordinate system rotating with satellite as (X0,Y0,Z)
; where X0 and Y0 are the DMXL coordinates at the start of the interval (arbitrary initial phase), and Z is the spin axis.
; The transformation matrix is (approximately) fixed in time, assuming Z does not change (much) in time.
; Z axis can change (attitude may vary) but this is very small.
; As satellite rotates about Z, the Bperp field rotates in satellite spin plane in opposite direction.
; Sidereal spin period remains constant here, but in future can change with time, based on state file, no problem including.
; Sidereal spinphase = wsid*(t-t0) and BX0=Bperp(t)*sin(-wsid*(t-t0)); BY0=Bperp(t)*cos(-wsid*(t-t0)).
; An arbitraty start time spinphase0 can account for the fact that X0 may not point at Bperp zero-crossing initially
; spinphase_init is assumed to be = 0. now, without loss of generality. (Means t0 starts when X0 is along
;    Bperp(t), and can be found from DMXL (Yaxis) discussed in the previous section.
spinphase_init=0.d
spinphase0 = spinphase_init + wsid*(elx_att_gei.x-elx_att_gei.x[0]) ; this is spinphase for fixed coordinate system (X0,Y0,Z0)
;
; Now get correction angle dsp to spinphase (sp) for a moving coordinate system. The correction d(deltaspinphase) is
;   cumulative-summed in time, as it's due to both spin-axis change and Bperp direction change along sat. path.
;   Simply taking the angle of Y(t+dt) relative to Y(t) will give the new Bperp direction momentary change.
;   This d(dsp) is positive if along the spin period (phi is in the positive spin phase direction), i.e., R.H. around Z.
;
Zaxis_last=Zaxis
Zaxis_last[1:npoints-1,*]=Zaxis[0:npoints-2,*] ; this is the last point's Z axis
Yaxis_proj=crossp2(Zaxis_last,crossp2(Yaxis,Zaxis_last))
mytot=sqrt(total(Yaxis_proj^2,2))
for j=0,2 do Yaxis_proj[*,j]=Yaxis_proj[*,j]/mytot ; this is unit vector of Yaxis projected onto last time's XY plane
;
Yaxis_last=Yaxis
Yaxis_last[1:npoints-1,*]=Yaxis[0:npoints-2,*] ; this is the last point's Y axis
ddsp=asin(total(crossp2(Yaxis_last,Yaxis_proj)*Zaxis_last,2)) ; rotation angle (rads) of Bperp the XY plane, >0 in +Z sense. 
;
dsp=total(ddsp,/cumulative) ; dsp=spinphase angle (rads) relative to the original angle (its deriv. is: wsyn(t)-wsid(t))
store_data,'dsp',data={x:elx_att_gei.x,y:dsp*180./!PI} ; in degrees just for plotting purposes!
spinphase = spinphase0 + dsp ; this is spinphase in local coordinate system (X,Y,Z)
;
snwsidt0=sin(+spinphase0)
cswsidt0=cos(+spinphase0)
snwsidt=sin(+spinphase)
cswsidt=cos(+spinphase)
get_data,'elx_bt89_dmxl',data=elx_bt89_dmxl,dlim=dl_elx_bt89_dmxl,lim=l_elx_bt89_dmxl
BX0 =   elx_bt89_dmxl.y[*,0] * cswsidt0 + elx_bt89_dmxl.y[*,1] * snwsidt0
BY0 = - elx_bt89_dmxl.y[*,0] * snwsidt0 + elx_bt89_dmxl.y[*,1] * cswsidt0
store_data,'elx_bt89_smxl0',data={x:elx_att_gei.x,y:[[BX0],[BY0],[elx_bt89_dmxl.y[*,2]]]},dlim=dl_elx_bt89_dmxl,lim=l_elx_bt89_dmxl
BX =   elx_bt89_dmxl.y[*,0] * cswsidt + elx_bt89_dmxl.y[*,1] * snwsidt
BY = - elx_bt89_dmxl.y[*,0] * snwsidt + elx_bt89_dmxl.y[*,1] * cswsidt
store_data,'elx_bt89_smxl',data={x:elx_att_gei.x,y:[[BX],[BY],[elx_bt89_dmxl.y[*,2]]]},dlim=dl_elx_bt89_dmxl,lim=l_elx_bt89_dmxl
options,'elx_bt89_?mxl',colors=['b','g','r']
options,'elx_bt89_smxl',labels=['X','Y','Z'],YSUBTITLE='elx_bt89_smxl'
options,'elx_bt89_smxl0',labels=['X0','Y0','Z'],YSUBTITLE='elx_bt89_smxl'
;
options, 'elx*','databar',0.,linestyle=2
;
options, 'elx*', 'databar', {yval:[0.], linestyle:2}
;tplot,'elx_bt89_gei emicwave elx_bt89_dmxl dsp elx_bt89_smxl'
;tplot_apply_databar
;timebar,[tstart1,tstart2]
;timebar,[tend1,tend2],linestyle=2
;
; get data in sensor coordinates
; SMXL is approximately along the geometric CubeSat coordinate system
; SMXL X,Y,Z are: X ~along long-axis of CubeSat, Y is ~along FGM boom direction, Z is along spin axis (Lvec).
; FGM Z axis is along boom direction, X and Y axes are rotated about FGMZ (=SMXL Y) on their way out.
; Based on lab and on-orbit tests pulsing the coils, the angle of rotation for EL-A is 44deg and the gain ratio Gx/Gy=0.89
; We did not do a pulse test for EL-B. However the angle and gain may be different anyway. The in-orbit phase
; analysis also shows that using FGM-Z (radial) as reference, the FGMX and FGMY sensor axes are a few deg. off from orthogonal.
; Assuming the above the rotation matrix is as follows
fgm2smxl_rotang=double(44.);deg rot. of FGM-X,Y data about FGMZ to X'Y', that will render the Y' axis opposite to SMXL Z axis (L)
csfgmrotang=cos(fgm2smxl_rotang*!PI/180.)
snfgmrotang=sin(fgm2smxl_rotang*!PI/180.)
fgmx=[csfgmrotang,0,snfgmrotang] ; in smxl system
fgmy=[snfgmrotang,0,-csfgmrotang]; in smxl system
fgmz=[0,1,0]                     ; in smxl system
smxl2fgm = transpose([[fgmx],[fgmy],[fgmz]]); takes a point in SMXL system and rotates it to FGM sensor coord's
smxl2fgm_rotmat=make_array(npoints,3,3,/double)
for j=0,npoints-1 do smxl2fgm_rotmat[j,*,*]=smxl2fgm[*,*]
store_data,'smxl2fgm',data={x:elx_att_gei.x,y:smxl2fgm_rotmat}; for tplot var's (use keyword /inverse for fgm2smxl)
; example: tvector_rotate,'smxl2fgm','fgmx',newname = 'fgmx_in_fgmcoords' ; from GEI to DMXL produces [1,0,0] everywhere
tvector_rotate,'smxl2fgm','elx_bt89_smxl',newname = 'elx_bt89_fgmc' ; in fgm coordinates now
;
; insert sensor noise to the datastream (random noise at 50pT)
;
seed=10
noise_amp=0.05 ; nT
fgmx_noise = RANDOMU(seed,npoints)
fgmy_noise = RANDOMU(seed,npoints)
fgmz_noise = RANDOMU(seed,npoints)
store_data,'fgm_noise',data={x:elx_att_gei.x,y:[[fgmx_noise],[fgmy_noise],[fgmz_noise]]*2.*noise_amp-noise_amp}
options,'fgm_noise',colors=['b','g','r']
options,'fgm_noise',labels=['X','Y','Z'],YSUBTITLE='elx_bt89_dmxl'
options,'fgm_noise','databar', {yval:[0.], linestyle:2}
;
calc,"'elx_bt89_fgmc'='elx_bt89_fgmc' +'fgm_noise'"
;
tplot,'elx_bt89_gei emicwave elx_bt89_dmxl elx_bt89_dmxl_true dsp elx_bt89_smxl fgm_noise elx_bt89_fgmc'
tplot_apply_databar
timebar,[tstart1,tstart2]
timebar,[tend1,tend2],linestyle=2
;
print,'*********************************************************************************'
stop
;makepng,'fig01'
;
tplot_ascii,'elx_att_gei elx_bt89_gei emicwave elx_bt89_dmxl_true elx_bt89_dmxl dsp elx_bt89_smxl fgm_noise elx_bt89_fgmc',trange=[tstart1,tend1],/precise
;tplot_ascii,'elx_att_gei elx_bt89_gei emicwave elx_bt89_dmxl_true elx_bt89_dmxl dsp elx_bt89_smxl fgm_noise elx_bt89_fgmc',trange=[tstart2,tend2],/precise
;
stop
;
; get igrf during fgs times
;
tinterpol_mxn,'elx_igrf_dmxl','ela_fgs',newname='elx_igrf_dmxl_int',/no_extrapolate ; sets the extrapolated values to near-zero
tplot,'ela_fgs elx_igrf_dmxl_int elx_bt89_gei emicwave elx_bt89_dmxl dsp elx_bt89_smxl fgm_noise elx_bt89_fgmc'
tplot_apply_databar
timebar,[tstart1,tstart2]
timebar,[tend1,tend2],linestyle=2
;makepng,'fig02'
tplot_ascii,'ela_fgs elx_igrf_dmxl_int',trange=[tstart1,tend1],/precise
;tplot_ascii,'ela_fgs elx_igrf_dmxl_int',trange=[tstart2,tend2],/precise
;
stop
end
