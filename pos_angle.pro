
FUNCTION  pos_angle, vel, xplot, yplot, dis, re

;Derive Rotation Angle
;-----------------------


;blue & red
w1=where(vel gt 0.75*max(vel)) ;redshift
ra1=median(xplot(w1)) 
dec1=median(yplot(w1))

w2=where(vel lt 0.75*min(vel)) ;blueshift
ra2=median(xplot(w2))
dec2=median(yplot(w2))

rot_line=linfit([ra1,0,ra2],[dec1,0,dec2])
m_rot=rot_line(1)
beta=atan(m_rot)*360/(2*!pi) ;in degrees

if ra1 gt ra2 then  alpha=90-beta      ;redshift is upper right part
if ra1 lt ra2 then  alpha=180+(90-beta)



;no rotation
w0=where(abs(vel) lt 10 and dis lt 0.5*re) 
if w0(0) ne -1 then begin
central_line=linfit(xplot(w0), yplot(w0))
m0=-1/central_line(1)
beta0=atan(m0)*360/(2*!pi) ;in degree

if ra1 gt ra2 then  alpha0=90-beta0      ;redshift is upper right part
if ra1 lt ra2 then  alpha0=180+(90-beta0)

endif else begin alpha0=-1000.
endelse

if ra1 lt ra2 then add=1 else add=0

PA=[alpha0, alpha, add ]

return, PA


END
