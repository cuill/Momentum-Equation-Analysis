load "./dens2.ncl"
load "./set_sigma.ncl"
external BAROC "./baropg90.so"

begin
 
  data = asciiread("info_cell_art.dat",-1,"string")
  art  = stringtodouble(str_get_field(data,1," "))
  
  data2 = asciiread("cross-section_elements.dat",-1,"string")
  eles  = stringtoint(str_get_field(data2,1," "))
  dime = dimsizes(eles)
  print(eles)

;read coriolis parameters
  data3 = asciiread("BaratariaPass_cor.dat",-1,"string")
  cor = stringtodouble(str_get_field(data3,3," "))

  iele = 162184
  f = addfile("../netcdf/BaratariaPass_0001.nc","r")
  nv = floattoint(f->nv)
  vx = flt2dble(f->x)
  vy = flt2dble(f->y)
  h  = flt2dble(f->h)

  dim=dimsizes(nv)
  n=dim(1)
  m=dimsizes(f->x)

  xc=new((/n/),double)
  yc=new((/n/),double)
  corr=new((/n/),double)
  do i=0,n-1
    xc(i)=(vx(nv(0,i)-1)+vx(nv(1,i)-1)+vx(nv(2,i)-1))/3.0d
    yc(i)=(vy(nv(0,i)-1)+vy(nv(1,i)-1)+vy(nv(2,i)-1))/3.0d
    corr(i)=(cor(nv(0,i)-1)+cor(nv(1,i)-1)+cor(nv(2,i)-1))/3.0d
  end do
  corr = 2.*7.292e-5*sin(corr*2.0*3.14159/360.0)

  T = new((/18,m/),double)
  T = 20.0d

; time 72 (it=71, 2008-07-31_12:00)
  it = 71
  S1 = flt2dble(f->salinity(it,:,:))
; density at nodes
  rho1 = dens2(T,S1)
 ; print(rho1(:,0))

; density at elements
  rho = new((/18,n/),double)
  do i=0,n-1
   rho(:,i)=(rho1(:,nv(0,i)-1)+rho1(:,nv(1,i)-1)+rho1(:,nv(2,i)-1))/3.0d
  end do
 ; print(rho(:,0))

  el = flt2dble(f->zeta(it,:))
; depth at nodes
  DT = el + h
; depth at elements
  DT1=new((/n/),double)
  do i=0,n-1
    DT1(i)=(DT(nv(0,i)-1)+DT(nv(1,i)-1)+DT(nv(2,i)-1))/3.0d
  end do

  dpbcx = new((/kbm1,n/),double)
  dpbcy = new((/kbm1,n/),double)
  dpbcx = 0.0d
  dpbcy = 0.0d

  print(kbm1)
  print(zz)
  BAROC::BAROPG(rho1,rho,vx,vy,zz(0:17),nv,n,m,kbm1,art,DT,DT1,dpbcx,dpbcy)
  dpbcx = dpbcx/1000.0d
  dpbcy = dpbcy/1000.0d
  ;print(dpbcx(:,ele-1))
  ;print(dpbcy(:,ele-1))

;calculate coriolis forcing
  corx = new((/kbm1,dime/),double)
  cory = new((/kbm1,dime/),double)

  u = flt2dble(f->u(it,:,eles-1))
  v = flt2dble(f->v(it,:,eles-1))
  do k=0,kbm1-1
    corx(k,:) =  corr(eles-1)*v(k,:)*DT1(eles-1)
    cory(k,:) = -corr(eles-1)*u(k,:)*DT1(eles-1)
  end do

;calculate vertical advection
  VadvX = new((/kbm1,dime/),double)
  VadvY = new((/kbm1,dime/),double)
  w = flt2dble(f->w(it,:,eles-1))
; k=0
  VadvX(0,:) = -w(1,:)*(u(0,:)*dz(1)+u(1,:)*dz(0))/(dz(0)+dz(1))
  VadvY(0,:) = -w(1,:)*(v(0,:)*dz(1)+v(1,:)*dz(0))/(dz(0)+dz(1))
; k=17
  VadvX(17,:) = w(17,:)*(u(17,:)*dz(16)+u(16,:)*dz(17))/(dz(17)+dz(16))
  VadvY(17,:) = w(17,:)*(v(17,:)*dz(16)+v(16,:)*dz(17))/(dz(17)+dz(16))
; k=1,16
  do k=1,kbm1-2
     VadvX(k,:) = w(k,:)*(u(k,:)*dz(k-1)+u(k-1,:)*dz(k))/(dz(k)+dz(k-1)) - \
                  w(k+1,:)*(u(k,:)*dz(k+1)+u(k+1,:)*dz(k))/(dz(k)+dz(k+1))
     VadvY(k,:) = w(k,:)*(v(k,:)*dz(k-1)+v(k-1,:)*dz(k))/(dz(k)+dz(k-1)) - \
                  w(k+1,:)*(v(k,:)*dz(k+1)+v(k+1,:)*dz(k))/(dz(k)+dz(k+1))
  end do
  
  
  write_table("cross-section_barocp_flood_time0072.dat","w",[/dpbcx(:,eles-1),dpbcy(:,eles-1),corx,cory,VadvX,VadvY/], \
              "%12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f")
end
