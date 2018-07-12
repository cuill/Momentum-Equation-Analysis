begin

  data2 = asciiread("cross-section_elements.dat",-1,"string")
  eles  = stringtoint(str_get_field(data2,1," "))
  dime = dimsizes(eles)
  ;print(eles)

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

; time 72 (it=71*2, 2008-07-31_12:00)
 do it=142,142
   u = flt2dble(f->u(it,:,eles-1))
   v = flt2dble(f->v(it,:,eles-1))
   uf = flt2dble(f->u(it+1,:,eles-1))
   vf = flt2dble(f->v(it+1,:,eles-1))
   el = flt2dble(f->zeta(it,:))
   D = h + el
   elf = flt2dble(f->zeta(it+1,:))
   DT = h+ elf

   D1=new((/n/),double)
   DT1=new((/n/),double)
   do i=0,n-1
        D1(i)=(D(nv(0,i)-1)+D(nv(1,i)-1)+D(nv(2,i)-1))/3.0d
        DT1(i)=(DT(nv(0,i)-1)+DT(nv(1,i)-1)+DT(nv(2,i)-1))/3.0d
   end do
 end do  ; loop it
  
  D11 = D1(eles-1)
  DT11 = DT1(eles-1)

; Vertical diffusion
  ddtX = new((/18,dime/),double)
  ddtY = new((/18,dime/),double)
  ddtX = 0.0d
  ddtY = 0.0d

  do k=0,17
     ddtX(k,:) = (DT11*uf(k,:)-D11*u(k,:))/0.2d
     ddtY(k,:) = (DT11*vf(k,:)-D11*v(k,:))/0.2d
  end do

  write_table("cross-section_ddt_flood_time0072.dat","w",[/ddtX,ddtY/], \
              "%12.6f   %12.6f")
end
