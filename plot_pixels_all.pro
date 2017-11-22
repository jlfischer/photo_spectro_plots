PRO plot_pixels_all

;Purpose:
;Plotting MANGA vdisp & stellar velocity (cubes and radial profiles)

dir='/data3/MANGA/MPL-5/'
type='LOGCUBE'

;userdefined symbol
t = findgen(21)/20.0 * 2.0 * !Pi
usersym, cos(t)*3., sin(t)*3., /fill


;table with galaxy photometry
tab=mrdfits('/home/helenado/MANGA/catalogues/tab_ser_r_manga_dr14.fits', 1, hdr_tab)


restore, '/home/helenado/MANGA/test_vdisp/PA_npix.sav', /verbose
pa_hds=pa[*,1] ;position angle derived by HDS
; RESTORE: Restored variable: GALCOUNT.
; RESTORE: Restored variable: NPIX.
; RESTORE: Restored variable: NPIX_OK.
; RESTORE: Restored variable: PA.

;=====================cilce begins================

stop

for i=0, n_elements(tab)-1 do begin

print, "Galaxy", i
print, "Galaxy", i
print, '================='

name=tab(i).plateifu
z=tab(i).z

if z gt 0 then begin

;Reff
re=tab(i).r_tot
fac=red_xcl_factor(z)    ;kpc/arcsec
re_kpc=re*fac ;Effective radious in kpc
cont_vect=re*[0.5,1.,2.]

; Read Galaxy structure
dir1='/data3/MANGA/MPL-5/tables/all/offset/'
file=+strcompress(name, /remove_all)+'-'+type+''

tmp_exist = FILE_TEST(dir1+file+'-veldisp_corr.fits') 
if tmp_exist eq 1 then begin  ; if file doesn't exist skip cicle
gal=mrdfits(dir1+file+'-veldisp_corr.fits',1, hdr1)


;correct stellar velocity
if gal(0).offset gt -100 and gal(0).offset lt 100 then begin
st_vel_corr=gal.st_vel-gal.offset
endif else begin
st_vel_corr=gal.st_vel
endelse

;--------central pixel ra dec -------

ra0=gal[0].ra0   ;central ra,dec derived in sigma_lick_cube_HDS
dec0=gal[0].dec0

if abs(tab(i).ra-ra0) gt 1.e-3 or abs(tab(i).dec-dec0) gt 1.e-3 then begin
print, 'Something wrong with position!'
stop

endif


;=========== Plot Positions=================

;p =  [xlow, ylow, xup,  yup ]
p1  = [0.07, 0.52, 0.45, 0.9]
pb1 = [0.5, 0.50, 0.52, 0.9]

p2  = [0.07, 0.07, 0.45, 0.45]
pb2 = [0.5, 0.07, 0.52, 0.45]

p3  = [0.6, 0.52, 0.94, 0.90]
p4  = [0.6, 0.07, 0.94, 0.45]

;========================= PLOT ==================================================
;=================================================================================

;Definitions
;----------------

if n_elements(gal) gt 1 then begin
ok=where(gal.veldisp gt mean(gal.veldisp)-3*stdev(gal.veldisp) and gal.veldisp lt mean(gal.veldisp)+3*stdev(gal.veldisp) and gal.veldisp lt 400. and gal.st_vel ne min(gal.st_vel))
endif else begin ok=-1
endelse


if n_elements(ok) gt 50 then begin  ;avoid galaxies with less than 50 ok pixels
if min(gal(ok).st_vel) lt 0. and  max(gal(ok).st_vel) gt 0. then begin ;avoid galaxies with not negative or positive velocities

print, 'Empieza plot'
print, file
print, '============'


set_plot,'ps'
 
device,filename='/data3/MANGA/MPL-5/plots/all/offset/radial_profile_'+file+'.ps'
device,xsize=20.0,ysize=20.0,xoffset=0.,yoffset=4.0,/color


;===================
;vdisp
;===================


print, 'ok Veldisp array'
help, gal.veldisp
help, ok


;variables
;---------------------
y=(gal.dec-dec0)*3600
z=gal.veldisp
x=sqrt(gal.dis^2-y^2)

x(where(gal.ra lt ra0))=-x(where(gal.ra lt ra0)) ;include negative sqrt
nan=(where(finite(x) eq 0))
if nan(0) ne -1 then  x(nan)=0 ;avoid nan (when dis ~ y)

xplot=x(ok)
yplot=y(ok)
zplot=z(ok)


;Plot limits
vv=[abs(xplot), abs(yplot)]
xlim1=max(vv)+max(vv)*0.15
xlim2=-xlim1

ylim1=xlim2
ylim2=xlim1

;color bar
a=round(min(zplot))-1
b=round(max(zplot))+1

m=255./(b-a)
r=-a*m
vcol=zplot*m+r

;pixel size
pix_size=0.4

;PLOT
;------------------------
plot, xplot, yplot, psym=8,xtit='dis RA [arcsec]',ytit='dis Dec [arcsec]',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=9,XTICKFORMAT='(F5.1)',position=p1,charsize=0.65,  /nodata , /noerase

;plot all pixels
cgLoadct,0
usersym, cos(t)*3., sin(t)*3.
oplot,x, y, psym=8, symsize=pix_size

;plot color vdisp
loadct, 20
usersym, cos(t)*3., sin(t)*3., /fill
for k=0, n_elements(gal(ok).ra)-1 do begin  &$
  oplot, [xplot(k), xplot(k)], [yplot(k), yplot(k)], symsize=pix_size, col=vcol(k), psym=8  &$
endfor

cgColorbar,divisions=11,range=[a,b],/vertical,tit='Velocity Dispersion [km/s]',tlocation="right",charsize=0.9,position=pb1

;centre
cgLoadct,0
oplot, [0, 0], [0, 0], psym=7, thick=3, symsize=2, col=0


;Contours
cgLoadct,0
if re gt 0 then begin
cgContour, gal(ok).dis, xplot, yplot, /irregular, Color=cgColor('Dark Green'),  LEVELS=cont_vect ,xtit='',ytit='',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=1,ystyle=1,XTICKFORMAT='(F5.1)',position=p1,charsize=0.65,  /noerase
endif


;Galaxy info
;================================================

cgLoadct,0
xyouts, [xlim1, xlim1],[ylim2*1.1, ylim2*1.1], 'ID_Manga='+strcompress(tab(i).plateifu, /remove_all)+'', charthick=3, charsize=1
xyouts, [xlim1, xlim1],[ylim2*1.2, ylim2*1.2] , 'Galcount='+strcompress(STRING(tab(i).galcount, format='(I)'), /remove_all)+'', charthick=3, charsize=1
xyouts, [0.5*xlim2, 0.5*xlim2],[ylim2*1.1, ylim2*1.1 ], 'R_eff='+strcompress(STRING(tab(i).r_tot, format='(F5.2)'), /remove_all)+' [arcsec]', charthick=3, charsize=1
xyouts, [0.5*xlim2, 0.5*xlim2],[ylim2*1.2, ylim2*1.2], 'z='+strcompress(STRING(tab(i).z, format='(F5.3)'), /remove_all)+'', charthick=3, charsize=1.0


;=================================
;stellar velocity
;=================================

;Definitions
;----------------
bad=where(gal.st_vel eq min(gal.st_vel))  ; min values are too small! Probably associated to error
good=where(gal.st_vel ne min(gal.st_vel))

xplot=x(good)
yplot=y(good)
zplot=st_vel_corr(good)

;color bar
a=-200.
b=200.

m=255./(b-a)
r=-a*m
vcol2=zplot*m+r

;PLOT
;----------------------
plot, xplot, yplot, psym=1, xtit='delta RA (arcsec)',ytit='delta Dec (arcsec)',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=9,XTICKFORMAT='(F5.1)',position=p2,charsize=0.65, /noerase, /nodata

;all pixels
cgLoadct,0
usersym, cos(t)*3., sin(t)*3.
oplot, x, y, psym=8, symsize=pix_size


;color stellar vel
cgloadct, 18, /reverse
usersym, cos(t)*3., sin(t)*3., /fill
for k=0, n_elements(xplot)-1 do begin  &$
  oplot, [xplot(k), xplot(k)], [yplot(k), yplot(k)], symsize=pix_size, col=vcol2(k), psym=8  &$
endfor

cgColorbar,divisions=11,range=[a, b],/vertical,tit='Stellar Velocity [km/s]',tlocation="right",charsize=0.9,position=pb2


;contours
cgLoadct,0
if re gt 0 then begin
cgContour, gal(good).dis,xplot, yplot, /irregular, Color=cgColor('Dark Green'),  LEVELS=cont_vect, thick=2,xtit='',ytit='',xr=[xlim1, xlim2],yr=[ylim1,ylim2],xstyle=1,ystyle=1,XTICKFORMAT='(F5.1)',position=p2,charsize=0.65, /noerase
endif


;Rotation angle
alpha=pa_hds(i)
m_pa=tan((90-alpha)*2*!pi/360) ;slope

x1=min(xplot)
x2=max(xplot)

;Pa line
oplot,[x1,x2], [x1,x2]*m_pa, col=cgcolor('Navy'), thick=6

;centre
cgLoadct,0
oplot, [0, 0], [0, 0], psym=7, thick=3, symsize=2, col=0

;Label PA
xyouts, [xlim1*0.9,xlim1*0.9],[ylim2*0.8, ylim2*0.8], 'PA='+strcompress(STRING(alpha, FORMAT='(I)'), /remove_all)+'', charthick=3, charsize=1.0, col=cgcolor('Navy')


;===================== Profiles =====================================================
;====================================================================================

xfit=gal(ok).dis_kpc
yfit= gal(ok).veldisp

;fits to vdisp radial profiles
;--------------------------------------
;pixels mean
meanbin,xfit, yfit,0,ycut1=1.,ycut2=300,xmin=0.1,xmax=10.,xbin=0.5,xmean=xmean,ymean=ymean,rmsmeant=rmsmean,weight=weight,nelem=nelem,sym=8,nmax=1,zbin=zbin,rms68=yper68,rms32=yper32,rms95=yper95,rms5=yper5,/no2sigmacontour,rmsmeanm=rmsmeanm,rmsmeanp=rmsmeanp, /noshow

;pixels linear
lin_cte=linfit(xfit, yfit)


;=================================
; velocity dispersion
;=================================

xlim1=min(xfit)-0.5
xlim2=max(xfit)+1

ylim1=min(yfit)-10
ylim2=max(yfit)+10


;plot data
;-----------------------------------
plot,xfit, yfit, psym=8,xtit='',ytit='', xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=1,XTICKFORMAT='(F4.1)',position=p3,charsize=1.0, /noerase, /nodata

POLYFILL, [xmean,reverse(xmean)] ,[yper68, reverse(yper32)], color=cgColor('RYB3')
oplot,xfit, yfit, psym=8, symsize=0.1, color=cgColor('Gray')

if n_elements(xmean) gt 1 then begin
oplot, xmean, ymean, col=cgColor('RYB1'), thick=6
oplot, xmean, ymean+rmsmeanp, line=2, col=cgColor('RYB2'), thick=3
oplot, xmean, ymean-rmsmeanm, line=2, col=cgColor('RYB2'), thick=3
endif


;Vdisp integrated & fits
;---------------------------
oplot, [0,max(xfit)], [lin_cte[0], (lin_cte[0]+lin_cte[1]*max(xfit))], thick=7,color=cgColor('Charcoal') ;y=ax+b


;Reff lines
;----------------------
oplot, 0.5*[re_kpc, re_kpc], [ylim1, ylim2], line=2, thick=4, Color=cgColor('Dark Grey')
oplot, [re_kpc, re_kpc], [ylim1, ylim2], line=2, thick=4,Color=cgColor('Dark Grey')
oplot, 2.0*[re_kpc, re_kpc], [ylim1, ylim2], line=2, thick=4,Color=cgColor('Dark Grey')

xyouts, 0.52*[re_kpc, re_kpc], [ylim2*0.98],'0.5 Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0
if  1.02*[re_kpc] lt (max(gal.dis_kpc)+1)  then  xyouts, 1.02*[re_kpc, re_kpc], [ylim2*0.98],'Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0
if  2.02*[re_kpc] lt  (max(gal.dis_kpc)+1) then  xyouts, 2.02*[re_kpc, re_kpc], [ylim2*0.98],'2 Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0



;close axis
;----------------------
trange=[xlim1, xlim2]/fac
axis,xaxis=1,xr=trange,xst=1, xtitle='Distance (arcsec)',charsize=1.0



;=================================
;Stellar Velocity
;=================================


xfit=gal(good).dis_kpc
yfit=st_vel_corr(good)

xmin=min(xfit)-0.5
xmax=max(xfit)+1

ymin=a
ymax=b


;define upper-lower rotating vel
;----------------------
bin=findgen(round((xmax-xmin)/0.5))*0.5

st_vel_up=fltarr(n_elements(bin)-1)
st_vel_low=fltarr(n_elements(bin)-1)

for bb=0, n_elements(bin)-2 do begin &$
   
   kk=where((xfit gt bin(bb)) and (xfit lt bin(bb+1)))  &$

   if n_elements(kk) gt 5 then begin  &$

      st_vel_up(bb)=max(yfit(kk)) &$
      st_vel_low(bb)=min(yfit(kk)) &$

   endif else begin &$
        st_vel_low(bb)=-999. & st_vel_up(bb)=-999. &$
   endelse &$


endfor 



;PLOT
;----------------------
plot, xfit,yfit, psym=8,xtit='Distance (kpc)',ytit='', xr=[xmin, xmax],yr=[ymin,ymax],xstyle=9,ystyle=1,XTICKFORMAT='(F4.1)',position=p4,charsize=1.0, /noerase, /nodata
oplot, xfit, yfit, psym=8, symsize=0.1, color=cgColor('Gray')

;plot upper-lower limits
wu=where(st_vel_up ne -999. and st_vel_up gt 0)
wl=where(st_vel_low ne -999. and st_vel_low le 0)

if wl(0) ne -1 then begin
oplot, bin(wl)+0.5, st_vel_low(wl), psym=8, symsize=0.4, col=cgcolor('Dodger Blue')
oplot, bin(wl)+0.5, st_vel_low(wl), col=cgcolor('Dodger Blue')

oplot, bin(wl)+0.5, -1*st_vel_low(wl), psym=8, symsize=0.4, col=cgcolor('Cyan')
oplot, bin(wl)+0.5, -1*st_vel_low(wl), col=cgcolor('Cyan')
endif

if wu(0) ne -1 then begin
oplot, bin(wu)+0.5, st_vel_up(wu), psym=8, symsize=0.4,col=cgcolor('Deep Pink')
oplot, bin(wu)+0.5, st_vel_up(wu),col=cgcolor('Deep Pink')
endif


;Reff
;----------------------
oplot, 0.5*[re_kpc, re_kpc],[ymin*0.9, ymax*1.1] , line=2, thick=4, Color=cgColor('Dark Grey')
oplot, [re_kpc, re_kpc], [ymin*0.9, ymax*1.1], line=2, thick=4,Color=cgColor('Dark Grey')
oplot, 2.0*[re_kpc, re_kpc], [ymin*0.9, ymax*1.1], line=2, thick=4,Color=cgColor('Dark Grey')

xyouts, 0.52*[re_kpc, re_kpc], [ymax*0.95],'0.5 Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0
if  1.02*[re_kpc] lt (max(gal.dis_kpc)+1)  then  xyouts, 1.02*[re_kpc, re_kpc], [ymax*0.95],'Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0
if  2.02*[re_kpc] lt  (max(gal.dis_kpc)+1) then  xyouts, 2.02*[re_kpc, re_kpc], [ymax*0.95],'2 Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0

;0 axis line
;--------------
oplot, [xmin, xmax] , [0,0], line=2, thick=4, Color=cgColor('Dark Grey')

;close axis
;----------------------
trange=[xmin, xmax]/fac
axis,xaxis=1,xr=trange,xst=1, xtitle='',charsize=1.0


device,/close & set_plot, 'X'

endif  ;z==0

endif ; more than 50 pixels

endif ; pos and neg st_vel

endif ; .fits doesn't exist



endfor

stop
stop


END





