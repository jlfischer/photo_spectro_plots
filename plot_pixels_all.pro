PRO plot_pixels_all

;Test for plotting each pixel separately, problems with contours
;input is the same as for plot_contour_profile

;dir='/home/helenado/MANGA/data/Ell/all/'

dir='/data3/MANGA/MPL-5/'
type='LOGCUBE'

t = findgen(21)/20.0 * 2.0 * !Pi
usersym, cos(t)*3., sin(t)*3., /fill
pix_size=0.5

;table with galaxy info
tab=mrdfits('/home/helenado/MANGA/catalogues/tab_ser_r_manga_dr14.fits', 1, hdr_tab)




;=====================cilce begins================

;for i=0, n_elements(plate)-1 do begin
for i=0, 65 do begin
print, "Galaxy", i

;29

name=tab(i).plateifu
z=tab(i).z
re=tab(i).r_tot


fac=red_xcl_factor(z)    ;kpc/arcsec
re_kpc=re*fac ;Effective radious in kpc

cont_vect=re*[0.5,1.,2.]

; Read structure
dir1='/data3/MANGA/MPL-5/tables/all/offset/'
file=+strcompress(name, /remove_all)+'-'+type+''
gal=mrdfits(dir1+file+'-veldisp_corr.fits',1, hdr1)

;correct stellar velocity
st_vel_corr=gal.st_vel-gal.offset

;--------central pixel ra dec -------

ra0=gal[0].ra0   ;central ra,dec derived in sigma_lick_cube_HDS
dec0=gal[0].dec0

if abs(tab(i).ra-ra0) gt 1.e-3 or abs(tab(i).dec-dec0) gt 1.e-3 then begin
print, 'Something wrong with position!'
stop

endif

;=========== Plot Positions=================

p1=[0.09,0.4,0.4,0.65]
pb1=[0.44,0.4,0.47,0.65]

p2=[0.09,0.07,0.4,0.32]
pb2=[0.44,0.07,0.47,0.32]

p3=[0.59,0.4,0.92,0.65]

p4=[0.59,0.07,0.92,0.32]

;========================= PLOT ==================================================
;=================================================================================



print, 'Empieza plot'
print, file
;stop

set_plot,'ps'
 
device,filename='/data3/MANGA/MPL-5/plots/all/offset/radial_profile_'+file+'.ps'
device,xsize=20.0,ysize=23.0,xoffset=0.,yoffset=4.0,/color


;===================
;vdisp
;===================

;Definitions
;----------------
ok=where(gal.veldisp gt mean(gal.veldisp)-3*stdev(gal.veldisp) and gal.veldisp lt mean(gal.veldisp)+3*stdev(gal.veldisp) and gal.st_vel ne min(gal.st_vel)
;if flag(i) eq 3 then ok=where(gal.veldisp gt mean(gal.veldisp)-3*stdev(gal.veldisp) and gal.veldisp lt mean(gal.veldisp)+3*stdev(gal.veldisp) and gal.dis lt re)


print, 'ok Veldisp array'
help, gal.veldisp
help, ok


;color vdisp
;---------------------


xplot=gal(ok).ra
yplot=gal(ok).dec
zplot=gal(ok).veldisp

;a=100
;b=270

a=round(min(gal(ok).veldisp))-1
b=round(max(gal(ok).veldisp))+1

m=255./(b-a)
r=-a*m
vcol=zplot*m+r

xlim1=min(gal(ok).ra)-0.0003
xlim2=max(gal(ok).ra)+0.0003

ylim1=min(gal(ok).dec)-0.0003
ylim2=max(gal(ok).dec)+0.0003

;------------------------

plot, xplot, yplot, psym=8,xtit='RA',ytit='Dec',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=9,XTICKFORMAT='(F8.4)',position=p1,charsize=0.65,  /nodata , /noerase

cgLoadct,0
usersym, cos(t)*3., sin(t)*3.
oplot, gal.ra, gal.dec, psym=8, symsize=pix_size

usersym, cos(t)*3., sin(t)*3., /fill
loadct, 20
for k=0, n_elements(gal(ok).ra)-1 do begin  &$
  oplot, [xplot(k), xplot(k)], [yplot(k), yplot(k)], symsize=pix_size, col=vcol(k), psym=8  &$
endfor

cgColorbar,divisions=11,range=[a,b],/vertical,tit='Velocity Dispersion [km/s]',tlocation="right",charsize=0.9,position=pb1


;Contours
;----------------------
cgLoadct,0
cgContour, gal.dis, gal.ra, gal.dec, /irregular, Color=cgColor('Dark Green'),  LEVELS=cont_vect ,xtit='RA',ytit='Dec',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=1,ystyle=1,XTICKFORMAT='(F8.4)',position=p1,charsize=0.65,  /noerase

cgLoadct,0
oplot, [ra0, ra0], [dec0, dec0], psym=7, thick=3, symsize=2, col=0


;Galaxy info
;================================================

delta=max(yplot)-min(yplot)
print, "Delta value", delta

cgLoadct,0


cgLoadct,0
xyouts, [min(xplot),min(xplot)],[max(yplot)+1.35*delta,max(yplot)+1.35*delta], 'ID_Manga='+strcompress(tab(i).plateifu, /remove_all)+'', charthick=3, charsize=1
xyouts, [min(xplot),min(xplot)],[max(yplot)+1.20*delta,max(yplot)+1.20*delta] , 'Galcount='+strcompress(tab(i).galcount, /remove_all)+'', charthick=3, charsize=1

xyouts, [min(xplot),min(xplot)],[max(yplot)+0.5*delta,max(yplot)+0.5*delta], 'R_eff='+strcompress(tab(i).r_tot, /remove_all)+'', charthick=3, charsize=1
xyouts, [min(xplot),min(xplot)],[max(yplot),max(yplot)], 'z='+strcompress(tab(i).z, /remove_all)+'', charthick=3, charsize=1.0


;=================================
;stellar velocity
;=================================


;Definitions
;----------------
bad=where(gal.st_vel eq min(gal.st_vel))  ; min values are too small! Probably associated to error
good=where(gal.st_vel ne min(gal.st_vel))

;if flag(i) eq 3 then begin 
;   bad=where(gal.st_vel eq min(gal.st_vel) or gal.dis lt re)
;   good=where(gal.st_vel ne min(gal.st_vel) and gal.dis lt re)
;endif


xplot=gal(good).ra
yplot=gal(good).dec
zplot=st_vel_corr(good)

a=-200.
b=200.

m=255./(b-a)
r=-a*m
vcol2=zplot*m+r

;plot pixels
;----------------------
plot, xplot, yplot, psym=1, xtit='RA',ytit='Dec',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=9,XTICKFORMAT='(F8.4)',position=p2,charsize=0.65, /noerase, /nodata

cgLoadct,0
usersym, cos(t)*3., sin(t)*3.
oplot, gal.ra, gal.dec, psym=8, symsize=pix_size

usersym, cos(t)*3., sin(t)*3., /fill
loadct, 18
;loadct, 33
for k=0, n_elements(xplot)-1 do begin  &$
  oplot, [xplot(k), xplot(k)], [yplot(k), yplot(k)], symsize=pix_size, col=vcol2(k), psym=8  &$
endfor

cgColorbar,divisions=11,range=[a, b],/vertical,tit='Stellar Velocity [km/s]',tlocation="right",charsize=0.9,position=pb2

;contours
;----------------------
cgLoadct,0
cgContour, gal.dis, gal.ra, gal.dec, /irregular, Color=cgColor('Dark Green'),  LEVELS=cont_vect, thick=2,xtit='RA',ytit='Dec',xr=[xlim1, xlim2],yr=[ylim1,ylim2],xstyle=1,ystyle=1,XTICKFORMAT='(F8.4)',position=p2,charsize=0.65, /noerase
oplot, [ra0, ra0], [dec0, dec0], psym=7, thick=3, symsize=2, col=0


;===================== Profiles =====================================================
;====================================================================================

xfit=gal(ok).dis_kpc
yfit= gal(ok).veldisp

;fits to vdisp radial profiles
;--------------------------------------

;pixels mean
meanbin,xfit, yfit,0,ycut1=20,ycut2=300,xmin=0.1,xmax=10.,xbin=0.8,xmean=xmean,ymean=ymean,rmsmeant=rmsmean,weight=weight,nelem=nelem,sym=8,nmax=1,zbin=zbin,rms68=yper68,rms32=yper32,rms95=yper95,rms5=yper5,/no2sigmacontour,rmsmeanm=rmsmeanm,rmsmeanp=rmsmeanp, /noshow

;pixels linear
non_zero=where(xfit ne 0)
lin_cte=linfit(xfit, yfit)
log_lin_cte=linfit(alog10(xfit(non_zero)), alog10(yfit(non_zero)))


;=================================
; velocity dispersion
;=================================

xlim1=min(xfit)-0.5
xlim2=max(xfit)+1

ylim1=min(yfit)-10
ylim2=max(yfit)+10


;plot data
;-----------------------------------
plot,xfit, yfit, psym=8,xtit='',ytit='Velocity Dispersion', xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=1,XTICKFORMAT='(F4.1)',position=p3,charsize=1.0, /noerase, /nodata

POLYFILL, [xmean,reverse(xmean)] ,[yper68, reverse(yper32)], color=cgColor('RYB3')
oplot,xfit, yfit, psym=8, symsize=0.1, color=cgColor('Gray')
oplot, xmean, ymean, col=cgColor('RYB1'), thick=6
oplot, xmean, ymean+rmsmeanp, line=2, col=cgColor('RYB2'), thick=3
oplot, xmean, ymean-rmsmeanm, line=2, col=cgColor('RYB2'), thick=3


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

; Slope values
;-----------------
delta=max(yfit)-min(yfit)
xyouts, [min(xfit),min(xfit)],[max(yfit)+0.45*delta,max(yfit)+0.45*delta], 'Slope pix='+strcompress(STRING(lin_cte(1), FORMAT='(F6.3)'), /remove_all)+'',charthick=3, charsize=1.0

;close axis
;----------------------
trange=[xlim1, xlim2]/fac
axis,xaxis=1,xr=trange,xst=1, xtitle='Distance (arcsec)',charsize=1.0



;=================================
;Stellar Velocity
;=================================

;stop

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

   print, bb &$
   kk=where((xfit gt bin(bb)) and (xfit lt bin(bb+1)))  &$
   print, n_elements(kk) &$

   if n_elements(kk) gt 5 then begin  &$

      st_vel_up(bb)=max(yfit(kk)) &$
      st_vel_low(bb)=min(yfit(kk)) &$

   endif else begin &$
        st_vel_low(bb)=-999. & st_vel_up(bb)=-999. &$
   endelse &$


endfor 

;plot data
;----------------------
plot, xfit,yfit, psym=8,xtit='Distance (kpc)',ytit='Stellar Velocity', xr=[xmin, xmax],yr=[ymin,ymax],xstyle=9,ystyle=9,XTICKFORMAT='(F4.1)',position=p4,charsize=1.0, /noerase, /nodata

oplot, xfit, yfit, psym=8, symsize=0.1, color=cgColor('Gray')

oplot, bin(where(st_vel_low ne -999.))+0.5, st_vel_low(where(st_vel_low ne -999.)), psym=8, symsize=0.4, col=cgcolor('Dodger Blue')
oplot, bin(where(st_vel_low ne -999.))+0.5, st_vel_low(where(st_vel_low ne -999.)), col=cgcolor('Dodger Blue')

oplot, bin(where(st_vel_low ne -999.))+0.5, -1*st_vel_low(where(st_vel_low ne -999.)), psym=8, symsize=0.4, col=cgcolor('Cyan')
oplot, bin(where(st_vel_low ne -999.))+0.5, -1*st_vel_low(where(st_vel_low ne -999.)), col=cgcolor('Cyan')

oplot, bin(where(st_vel_up ne -999.))+0.5, st_vel_up(where(st_vel_up ne -999.)), psym=8, symsize=0.4,col=cgcolor('Deep Pink')
oplot, bin(where(st_vel_up ne -999.))+0.5, st_vel_up(where(st_vel_up ne -999.)),col=cgcolor('Deep Pink')



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

trange=[-2, 2]
axis,yaxis=1,yr=trange,yst=1, ytitle='B/T',charsize=1.0


device,/close & set_plot, 'X'



;endif



endfor

stop




END
