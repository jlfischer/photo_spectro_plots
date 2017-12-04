PRO plot_pixels_sn

;Purpose:
;Plotting MANGA vdisp & stellar velocity (cubes and radial profiles)
;Written by HDS on 11/2017


;userdefined symbol
t = findgen(21)/20.0 * 2.0 * !Pi
usersym, cos(t)*3., sin(t)*3., /fill


dir='/data3/MANGA/MPL-5/'

type='LOGCUBE'
type1='LOGCUBE-VOR10-GAU-MILESHC'
type2='MAPS-SPX-GAU-MILESHC'



;table with galaxy photometry
tab=mrdfits('/home/helenado/MANGA/catalogues/tab_ser_r_manga_dr14.fits', 1, hdr_tab)



;=====================cilce begins================



;for i=0, n_elements(tab)-1 do begin
for i=0, 0 do begin

galcount=tab(i).galcount
name=tab(i).plateifu
z=tab(i).z

print, "Galaxy", i
print, "Galaxy", galcount,'  ',  name
print, '========================'



if z gt 0 then begin

;Reff
re=tab(i).r_tot
fac=red_xcl_factor(z)    ;kpc/arcsec
re_kpc=re*fac ;Effective radious in kpc
cont_vect=re*[0.5,1.,2.]


;Read logcube
;-------------
dir='/data3/MANGA/MPL-5/data_cubes/logcube/'
file='manga-'+strcompress(name, /remove_all)+'-'+type+'.fits.gz'

data1=mrdfits(dir+file, 1, hdr1)    ;flux  + Need  hdr with central ra, dec !!!
data2=mrdfits(dir+file, 2, hdr)     ;ivar
;mask=mrdfits(dir+file, 3, hdr)      ;mask
wave_lin=mrdfits(dir+file, 4, hdr)  ;wavelength

dims=size(data1, /dimens)   ;cube pixels

print, 'Cube size'
print, '----------'
print, dims

;definitions
;---------------------------
wave=alog10(wave_lin)
err=sqrt(1./data2)          ;err=sqrt(1/ivar)

;--------central pixel ra dec -------

;central ra and dec from header
ra0=double(sxpar(hdr1,'OBJRA')) 
dec0=double(sxpar(hdr1,'OBJDEC'))


;Read stellar velocity & positions
;----------------------------------
dir2='/data3/MANGA/MPL-5/data_cubes/maps/'
file2='manga-'+strcompress(name, /remove_all)+'-'+type2+'.fits.gz'

tmp_exist = FILE_TEST(dir2+file2) 

if tmp_exist eq 1 then begin  ; if file doesn't exist skip cicle

spx_coo=mrdfits(dir2+file2, 1, hdr2)  
stellar_vel=mrdfits(dir2+file2, 15, hdr2)
sigma_manga=mrdfits(dir2+file2, 18, hdr2)
sigma_corr=mrdfits(dir2+file2, 21, hdr2) 

sigma_manga_corr=sqrt(sigma_manga^2-sigma_corr^2)
sigma_manga_corr(where(sigma_manga lt sigma_corr))=0

;derive, coord, dis, flux, sn
;------------------------------------
coord=dis_pix(dims, hdr1, data1, err)     ;function to derive dis in arcsec


;Definitions
;----------------

vdisp=reform(sigma_manga_corr, dims[0]*dims[1])
st_vel=reform(stellar_vel, dims[0]*dims[1])

ra=reform(coord[*,*,0], dims[0]*dims[1])
dec=reform(coord[*,*,1], dims[0]*dims[1])
dis=reform(coord[*,*,2], dims[0]*dims[1])
flux=reform(coord[*,*,3], dims[0]*dims[1])
sn=reform(coord[*,*,4], dims[0]*dims[1])



;correct by offset
;---------------
off=offset(dis,st_vel)   ;function to derive st_vel off


if off gt -100 and off lt 100 then begin
st_vel_corr=st_vel-off
endif else begin
st_vel_corr=st_vel
endelse
;vdisp corrected by stellar vel

sig_corr=vdisp*sqrt(1+(st_vel_corr/vdisp)^2)



;pixels to be plotted
sn_lim=5.
sn_lim1=0.

ok1=where(sn gt sn_lim1)  ;all pixels
ok=where(sn gt sn_lim)    ;sn pixels

ok10=where(sn gt 10)    ;sn pixels
ok8=where(sn gt 8)    ;sn pixels
ok6=where(sn gt 6)    ;sn pixels


mvd=median(vdisp(ok))
std=stdev(vdisp(ok))

print, "Number pixels ok "
help, ok1, ok


if n_elements(ok) gt 50 then begin  ;avoid galaxies with less than 50 ok pixels
if min(st_vel_corr) lt 0. and  max(st_vel_corr) gt 0. then begin ;avoid galaxies with not negative or positive velocities

print, 'Empieza plot'
print, file
print, '============'


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



;dir_plot='/home/helenado/MANGA/test_vdisp/panel/'
dir_plot='/data3/MANGA/MPL-5/plots/all/offset/SN5'

set_plot,'ps'
device,filename=dir_plot+'/radial_profile_'+strcompress(name, /remove_all)+'_new.ps'
device,xsize=20.0,ysize=20.0,xoffset=0.,yoffset=4.0,/color

;set_plot, 'X'



;vdisp
;===================


;variables
;---------------------
y=(dec-dec0)*3600
x=double(sqrt(dis^2-y^2))
z=vdisp


x(where(ra lt ra0))=-x(where(ra lt ra0)) ;include negative sqrt
nan=(where(finite(x) eq 0))


if nan(0) ne -1 then  begin

print,"NAN values xpos"
xx=fltarr(n_elements(nan))
for n=0, n_elements(nan)-1 do begin
   tmp=where(y eq y(nan(n)))
   tmp2=where(finite(x(tmp)) eq 0)
   xx(n)=(x(tmp(tmp2+1))+x(tmp(tmp2-1)))/2.
endfor

x(nan)=xx
endif



xplot=x(ok)
yplot=y(ok)
zplot=z(ok)


;Plot limits zoom
vv=[abs(x(ok1)), abs(y(ok1))]
xlim1=max(vv)
xlim2=-xlim1

ylim1=xlim2
ylim2=xlim1


a=round(min(zplot))
b=round(mean(zplot)+3*stdev(zplot))


m=255./(b-a)
r=-a*m
vcol=zplot*m+r



;pixel size
pix_size=0.4
chs=0.9

;PLOT
;------------------------
plot, xplot, yplot, psym=8,xtit='arcsec',ytit='arcsec',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=9,XTICKFORMAT='(F5.1)',position=p1,charsize=chs,  /nodata , /noerase

;plot all pixels
cgLoadct,0
usersym, cos(t)*3., sin(t)*3.
oplot,x(ok1), y(ok1), psym=8, symsize=pix_size*0.5

;plot color vdisp
loadct, 20
usersym, cos(t)*3., sin(t)*3., /fill
for k=0, n_elements(xplot)-1 do begin  &$
  oplot, [xplot(k), xplot(k)], [yplot(k), yplot(k)], symsize=pix_size, col=vcol(k), psym=8  &$
endfor

print, "color bar limits"
print, a, b
cgLoadct,0
cgColorbar,divisions=11,range=[a,b],/vertical,tit='Velocity Dispersion [km/s]',tlocation="right",charsize=0.9,position=pb1
loadct, 20
cgColorbar,divisions=11,range=[a,b],/vertical,tit='Velocity Dispersion [km/s]',tlocation="right",charsize=0.9,position=pb1

;centre
cgLoadct,0
oplot, [0, 0], [0, 0], psym=7, thick=3, symsize=2, col=0


;N-E arrow
width=(xlim1-xlim2)/2
yarr=width*0.7
xarr=-yarr

cgarrow, xarr, yarr, xarr, yarr + 0.2*width, hsize=200, /data , thick=5, hthick=5
cgarrow, xarr, yarr, xarr + 0.2*width, yarr, hsize=250., /data , thick=5
xyouts, [xarr - 0.1*width, xarr - 0.1*width], [yarr + 0.1*width, yarr + 0.1*width], 'N', charsize=chs, charthick=2
xyouts, [xarr + 0.1*width, xarr + 0.1*width], [yarr - 0.1*width, yarr - 0.1*width], 'E', charsize=chs, charthick=2


;Contours
cgLoadct,0
if re gt 0 then begin
cgContour, dis(ok), xplot, yplot, /irregular, Color=cgColor('Dark Green'),  LEVELS=cont_vect ,xtit='',ytit='',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=1,ystyle=1,XTICKFORMAT='(F5.1)',position=p1,charsize=chs,  /noerase
endif



;Galaxy info
;================================================

cgLoadct,0
xyouts, [xlim1, xlim1],[ylim2*1.1, ylim2*1.1], 'ID_Manga='+strcompress(tab(i).plateifu, /remove_all)+'', charthick=3, charsize=1
xyouts, [xlim1, xlim1],[ylim2*1.2, ylim2*1.2] , 'Galcount='+strcompress(STRING(tab(i).galcount, format='(I)'), /remove_all)+'', charthick=3, charsize=1
xyouts, [0.5*xlim2, 0.5*xlim2],[ylim2*1.1, ylim2*1.1 ], 'R_eff='+strcompress(STRING(tab(i).r_tot, format='(F5.2)'), /remove_all)+'', charthick=3, charsize=1
xyouts, [0.5*xlim2, 0.5*xlim2],[ylim2*1.2, ylim2*1.2], 'z='+strcompress(STRING(tab(i).z, format='(F5.3)'), /remove_all)+'', charthick=3, charsize=1.0


;=================================
;stellar velocity
;=================================

;Definitions
;----------------
bad=where(st_vel eq min(st_vel))  ; min values are too small! Probably associated to error
good=where(st_vel ne min(st_vel) and sn gt sn_lim)

xplot=x(good)
yplot=y(good)
zplot=st_vel_corr(good)


;color bar
a=-200.
b=200.

m=255./(b-a)
r=-a*m
vcol2=zplot*m+r



;position angle
if min(zplot) lt 0 and max(zplot) gt 0 then begin

   print,"Derive PA"
   PA=pos_angle(zplot, xplot, yplot, dis(good), re)

   alpha0=pa(0)
   alpha=pa(1)

endif else begin
   alpha0=-999.
   alpha=-999.
endelse


m_pa0=tan((90-alpha0)*2*!pi/360) ;slope
m_pa=tan((90-alpha)*2*!pi/360)


;PA from Manga
dir2='/data3/MANGA/MPL-5/data_cubes/maps/'
file2='manga-'+strcompress(name, /remove_all)+'-MAPS-SPX-GAU-MILESHC.fits.gz'

tmp=mrdfits(dir2+file2, 0, hdr)
str=strsplit(hdr(78), /extract, "   ")
alpha_mang=float(str(2))

 
m_mang=tan((90-alpha_mang)*2*!pi/360)
if pa(2) eq 1  then alpha_mang = alpha_mang + 180 ; if ra1 lt ra2

print, alpha, alpha0, alpha_mang


;PLOT
;----------------------
plot, xplot, yplot, psym=1, xtit='arcsec',ytit='arcsec',xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=9,XTICKFORMAT='(F5.1)',position=p2,charsize=chs, /noerase, /nodata

;all pixels
cgLoadct,0
usersym, cos(t)*3., sin(t)*3.
oplot, x(ok1), y(ok1), psym=8, symsize=pix_size*0.5


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
cgContour, dis(good),xplot, yplot, /irregular, Color=cgColor('Dark Green'),  LEVELS=cont_vect, thick=2,xtit='',ytit='',xr=[xlim1, xlim2],yr=[ylim1,ylim2],xstyle=1,ystyle=1,XTICKFORMAT='(F5.1)',position=p2,charsize=chs, /noerase
endif


;Rotation angle
x1=min(xplot)
x2=max(xplot)

;Pa line
oplot,[x1,x2], [x1,x2]*m_pa0, col=cgcolor('Navy'), thick=6
oplot,[x1,x2], [x1,x2]*m_pa, col=cgcolor('RED7'), thick=6
oplot,[x1,x2], [x1,x2]*m_mang, col=cgcolor('GRN7'), thick=6

;centre
cgLoadct,0
oplot, [0, 0], [0, 0], psym=7, thick=3, symsize=2, col=0

;Label PA
xyouts, [xlim1*0.9,xlim1*0.9],[ylim2*0.8, ylim2*0.8], 'PA0='+strcompress(STRING(alpha0, FORMAT='(F4.0)'), /remove_all)+'', charthick=3, charsize=chs, col=cgcolor('Navy')

xyouts, [xlim1*0.9,xlim1*0.9],[ylim2*0.7, ylim2*0.7], 'PA='+strcompress(STRING(alpha, FORMAT='(F4.0)'), /remove_all)+'', charthick=3, charsize=chs, col=cgcolor('RED7')


xyouts, [xlim1*0.9,xlim1*0.9],[ylim2*0.6, ylim2*0.6], 'PA_mang='+strcompress(STRING(alpha_mang, FORMAT='(F4.0)'), /remove_all)+'', charthick=3, charsize=chs, col=cgcolor('GRN7')

;=====================   Profiles   =================================================
;====================================================================================


;=================================
; velocity dispersion
;=================================

ok_vd=where(vdisp gt 0.5 and vdisp lt mvd+3*std and vdisp lt 400. and st_vel ne min(st_vel) and sn gt sn_lim)

print, "veldips fit"
help, ok, ok_vd
print, "============="

dis_kpc=dis*fac
xfit=dis_kpc(ok_vd)  ; dis in kcp
yfit=vdisp(ok_vd)



;fits to vdisp radial profiles
;--------------------------------------
;pixels mean
meanbin,xfit, yfit,0,ycut1=40.,ycut2=300,xmin=0.1,xmax=10.,xbin=0.5,xmean=xmean,ymean=ymean,rmsmeant=rmsmean,weight=weight,nelem=nelem,sym=8,nmax=1,zbin=zbin,rms68=yper68,rms32=yper32,rms95=yper95,rms5=yper5,/no2sigmacontour,rmsmeanm=rmsmeanm,rmsmeanp=rmsmeanp, /noshow



;compute integrated vdisp
r_vect=findgen(40)/2.+0.5
vdisp_int=fltarr(n_elements(r_vect))
vdisp_int_corr=fltarr(n_elements(r_vect))

for rr=0, n_elements(r_vect)-1 do begin
       
   p=where(dis_kpc le r_vect[rr] and vdisp gt 40. and vdisp lt mvd+3*std and sn gt sn_lim)

   num=0
   den=0
   
   if p(0) ne -1 then begin

   num_corr=total(sig_corr(p)*flux(p))
   num=total(vdisp(p)*flux(p))
   den=total(flux(p))
   vdisp_int(rr)=num/den
   vdisp_int_corr(rr)=num_corr/den
   

   endif else begin   
   vdisp_int(rr)=-999.
   vdisp_int_corr(rr)=-999.
   endelse
endfor



;define plot limits
xlim1=min(dis_kpc(good))-0.5 ;same x lim as stellar vel plot
xlim2=max(dis_kpc(good))+1

ylim1=min(yfit)-10
ylim2=max(yfit)+10


;plot data
;-----------------------------------
plot,xfit, yfit, psym=8,xtit='',ytit='', xr=[xlim1, xlim2],yr=[ylim1, ylim2],xstyle=9,ystyle=1,XTICKFORMAT='(F4.1)',position=p3,charsize=1.0, /noerase, /nodata

;plot mean
POLYFILL, [xmean,reverse(xmean)] ,[yper68, reverse(yper32)], color=cgColor('RYB3')
oplot,dis_kpc(ok), vdisp(ok), psym=8, symsize=0.1, color=cgColor('Light Gray')
oplot,dis_kpc(ok10), vdisp(ok10), psym=8, symsize=0.1, color=cgColor('Orange')


if n_elements(xmean) gt 1 then begin
oplot, xmean, ymean, col=cgColor('RYB1'), thick=6
oplot, xmean, ymean+rmsmeanp, line=2, col=cgColor('RYB2'), thick=3
oplot, xmean, ymean-rmsmeanm, line=2, col=cgColor('RYB2'), thick=3
endif



;Vdisp integrated
wr=where(r_vect lt 2.5*re_kpc and vdisp_int gt 0.)
if wr(0) ne -1  then oplot, r_vect[wr], vdisp_int[wr], psym=-8, symsize=0.4, col=cgcolor('Dark Green')

wr=where(r_vect lt 2.5*re_kpc and vdisp_int_corr gt 0.)
if wr(0) ne -1  then oplot, r_vect[wr], vdisp_int_corr[wr], psym=-8, symsize=0.3, col=cgcolor('Lime Green')


;Reff lines
oplot, 0.5*[re_kpc, re_kpc], [ylim1, ylim2], line=2, thick=4, Color=cgColor('Dark Grey')
oplot, [re_kpc, re_kpc], [ylim1, ylim2], line=2, thick=4,Color=cgColor('Dark Grey')
oplot, 2.0*[re_kpc, re_kpc], [ylim1, ylim2], line=2, thick=4,Color=cgColor('Dark Grey')

xyouts, 0.52*[re_kpc, re_kpc], [ylim2*0.98],'0.5 Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0
if  1.02*[re_kpc] lt (max(dis_kpc)+1)  then  xyouts, 1.02*[re_kpc, re_kpc], [ylim2*0.98],'Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0
if  2.02*[re_kpc] lt  (max(dis_kpc)+1) then  xyouts, 2.02*[re_kpc, re_kpc], [ylim2*0.98],'2 Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0



;close axis
trange=[xlim1, xlim2]/fac
axis,xaxis=1,xr=trange,xst=1, xtitle='Distance (arcsec)',charsize=1.0



;=================================
;Stellar Velocity
;=================================


xfit=dis_kpc(good)
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



;Plot data
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


;Reff lines
oplot, 0.5*[re_kpc, re_kpc],[ymin*0.9, ymax*1.1] , line=2, thick=4, Color=cgColor('Dark Grey')
oplot, [re_kpc, re_kpc], [ymin*0.9, ymax*1.1], line=2, thick=4,Color=cgColor('Dark Grey')
oplot, 2.0*[re_kpc, re_kpc], [ymin*0.9, ymax*1.1], line=2, thick=4,Color=cgColor('Dark Grey')

xyouts, 0.52*[re_kpc, re_kpc], [ymax*0.95],'0.5 Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0
if  1.02*[re_kpc] lt (max(dis_kpc)+1)  then  xyouts, 1.02*[re_kpc, re_kpc], [ymax*0.95],'Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0
if  2.02*[re_kpc] lt  (max(dis_kpc)+1) then  xyouts, 2.02*[re_kpc, re_kpc], [ymax*0.95],'2 Re',  charthick=3, charsize=1,  Color=cgColor('Dark Grey'),  Orientation=-90.0

;0 axis line
oplot, [xmin, xmax] , [0,0], line=2, thick=4, Color=cgColor('Dark Grey')

;close x axis
trange=[xmin, xmax]/fac
axis,xaxis=1,xr=trange,xst=1, xtitle='',charsize=1.0


;close plot
device, /close & set_plot, 'X'


endif  ;z==0

endif ; more than 50 pixels

endif ; pos and neg st_vel

endif ; .fits doesn't exist


endfor

;stop
;stop


END


;   FUNCTIONS
;==================================================================

;derive distance & coord for each pixel
;--------------------------------
Function dis_pix, dims, hdr, flux_cube, error_cube

ra0=double(sxpar(hdr,'OBJRA')) 
dec0=double(sxpar(hdr,'OBJDEC'))

dis_all=dblarr(dims[0], dims[1], 5)

for i=0, dims[0]-1 do begin
   x=i
   ;print, 'i='
   ;print, i

   for j=0, dims[1]-1 do begin
     
      y=j
  
      ;Distance to central pixel
      xyad, hdr, i, j, ra, dec          ;ra, dec of this pixel
      GCIRC, 2, ra0, dec0, ra, Dec, DIS ; dis in arcsec if ra & dec are in degrees


      ;flux, err at each pix
      flux=reform(flux_cube(i,j, *))
      error=reform(error_cube(i,j, *))
      sn=flux/error
      w=1/error^2
                      
      dis_all[i,j,0] = ra
      dis_all[i,j,1] = dec
      dis_all[i,j,2] = dis
      

     ; dis_all[i,j,3] = total(flux)  
     ; dis_all[i,j,4] = median(flux/error) 
      
      dis_all[i,j,3] = total(flux*w)/total(w)
      dis_all[i,j,4] = total(sn*w)/total(w)


  endfor
endfor

return, dis_all

END




;calculate stellar vel offset
;-------------------------------
Function offset, dis, stellar_vel

cc=where(dis lt max(dis)*0.05) ;pixels within 5% max dis

print, "calculating offset within", max(dis)*0.05, "arcsec"

if n_elements(cc) gt 1 then begin

   offset=median(stellar_vel(cc))
   offset_mean=mean(stellar_vel(cc))
   stellar_vel_corr=stellar_vel-offset

endif else begin

   offset=-99.
   offset_mean=-99.
   stellar_vel_corr=stellar_vel

endelse

print, "Offset:", offset, offset_mean

return, offset

END


;Pruebas neigbourghs

;xx=x((ok(sort(x(ok)))))  ;order ra

;tmp=fltarr(n_elements(ok)-1)
;for t=0, n_elements(xx)-2 do tmp(t)= xx(t+1)-xx(t) ;distance between pixels

;av=median(tmp)
