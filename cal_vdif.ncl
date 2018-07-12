load "./set_sigma.ncl"

begin

  ;print(z)
  ;print(zz)
  ;print(dzz)
  data2 = asciiread("cross-section_elements.dat",-1,"string")
  eles  = stringtoint(str_get_field(data2,1," "))
  dime = dimsizes(eles)
  ;print(eles)

  data = asciiread("info_cell_art.dat",-1,"string")
  art  = stringtodouble(str_get_field(data,1," "))

  ;iele = 162184
  f = addfile("../netcdf/BaratariaPass_0001.nc","r")
  nv = floattoint(f->nv)
  vx = flt2dble(f->x)
  vy = flt2dble(f->y)
  h  = flt2dble(f->h)
  a1u = flt2dble(f->a1u)
  a2u = flt2dble(f->a2u)

  dim=dimsizes(nv)
  n=dim(1)
  m=dimsizes(f->x)

  xc=new((/n/),double)
  yc=new((/n/),double)
  do i=0,n-1
    xc(i)=(vx(nv(0,i)-1)+vx(nv(1,i)-1)+vx(nv(2,i)-1))/3.0d
    yc(i)=(vy(nv(0,i)-1)+vy(nv(1,i)-1)+vy(nv(2,i)-1))/3.0d
  end do

; time = 72 (2008-07-31_12:00)
; in netcdf it is 71*2
; u - 142
; uf - 143
 do it=142,142
  ; u(n)
   u = flt2dble(f->u(it,:,eles-1))
   v = flt2dble(f->v(it,:,eles-1))
  ; u(n+1)
   uf = flt2dble(f->u(it+1,:,eles-1))
   vf = flt2dble(f->v(it+1,:,eles-1))
   
   kmz = flt2dble(f->km(it,:,:))

   el = flt2dble(f->zeta(it,:))
   elT = flt2dble(f->zeta(it,:))
   D = h + el
   DT = h + elT

   kmz2 = new((/19,n/),double)
   DT1 = new((/n/),double)
   D1 = new((/n/),double)
   do i=0,n-1
    kmz2(:,i)=(kmz(:,nv(0,i)-1)+kmz(:,nv(1,i)-1)+kmz(:,nv(2,i)-1))/3.0d
    DT1(i)=(DT(nv(0,i)-1)+DT(nv(1,i)-1)+DT(nv(2,i)-1))/3.0d
    D1(i)=(D(nv(0,i)-1)+D(nv(1,i)-1)+D(nv(2,i)-1))/3.0d
   end do

   kmz3 = kmz2(:,eles-1)
   DT2 = DT1(eles-1)
   D2  = D1(eles-1)

  ; kmzz2 = new((/18,dime/),double)
  ; do k=0,17
  ;    kmzz2(k,:) =0.5*(kmzz(k,:)+kmzz(k+1,:))
  ; end do

;calculate term1 and  term 2
  ; loop from k=2 to kb-2
   termu1 = new((/18,dime/),double)
   termu2 = new((/18,dime/),double)
   termv1 = new((/18,dime/),double)
   termv2 = new((/18,dime/),double)

   do k=1,16
     termu1(k,:) = kmz3(k,:)*(uf(k-1,:)-uf(k,:))/(z(k)-z(k+1)*(z(k-1)-z(k+1)))
     termu2(k,:) = kmz3(k+1,:)*(uf(k,:)-uf(k+1,:))/((z(k)-z(k+1))*(z(k)-z(k+2)))
     termv1(k,:) = kmz3(k,:)*(vf(k-1,:)-vf(k,:))/(z(k)-z(k+1)*(z(k-1)-z(k+1)))
     termv2(k,:) = kmz3(k+1,:)*(vf(k,:)-vf(k+1,:))/((z(k)-z(k+1))*(z(k)-z(k+2)))
   end do

 ; 2*(term1-term2)/DT2
   vdifU = new((/18,dime/),double)
   vdifV = new((/18,dime/),double)
   vdifU = 0.0d
   vdifV = 0.0d

   do k=1,16
     vdifU(k,:) = 2*(termu1(k,:)-termu2(k,:))/DT2   ;/art(eles-1)
     vdifV(k,:) = 2*(termv1(k,:)-termv2(k,:))/DT2   ;/art(eles-1)
   end do

;Surface wind stress
;add later

;Bottom friction
   CBCMIN = 0.0025d
   Z0     = 0.003d
   VK2    = 0.16d

   wubot = new((/dime/),double)
   wvbot = new((/dime/),double)
   do i=0,dime-1
       if(DT2(i) .gt. 3.0d)then
          ZTEMP=(zz(kbm1-1)-z(kb-1))*DT2(i)/Z0
          CBC = max((/CBCMIN,VK2/(log(ZTEMP))^2/))
       else
          ZTEMP=(zz(kbm1-1)-z(kb-1))*3.0/Z0
          CBC = max((/CBCMIN,VK2/(log(ZTEMP))^2/))
       end if

       TPS = CBC*sqrt(uf(kbm1-1,i)^2+vf(kbm1-1,i)^2)
       wubot(i) = -TPS*u(kbm1-1,i)
       wvbot(i) = -TPS*v(kbm1-1,i)
   end do 
;   print(wubot)
;   print(wvbot)

   vdifU(17,:) = 2*(kmz3(kbm1-1,:)*(uf(kbm2-1,:)-uf(kbm1-1,:))/(z(kbm2-1)-z(kb-1))- \
                 D2*wubot/1000.0d)/DT2/(z(kbm1-1)-z(kb-1))
   vdifV(17,:) = 2*(kmz3(kbm1-1,:)*(vf(kbm2-1,:)-vf(kbm1-1,:))/(z(kbm2-1)-z(kb-1))- \
                 D2*wvbot/1000.0d)/DT2/(z(kbm1-1)-z(kb-1))
      
 end do  ; loop it

 write_table("cross-section_vdif_flood_time0072.dat","w",[/vdifU,vdifV/], \
              "%12.6f   %12.6f")

end
