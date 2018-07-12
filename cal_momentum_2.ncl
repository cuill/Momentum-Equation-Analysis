external TGE "./tge90.so"
external ADVUV "./adv_uv_edg_gcn90.so"

begin

  data = asciiread("info_cell_art.dat",-1,"string")
  art  = stringtodouble(str_get_field(data,1," "))

  data2 = asciiread("cross-section_elements.dat",-1,"string")
  eles  = stringtoint(str_get_field(data2,1," "))
  dime = dimsizes(eles)
  print(eles)

  iele = 162184
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

  ne = 430057
  nct = 3*n
  ntsn = new((/m/),integer)
  nbsn = new((/11,m/),integer)
  ntve = new((/m/),integer)
  nbve = new((/9,m/),integer)
  nbvt = new((/9,m/),integer)
     
  ntsn = 0
  nbsn = 0
  ntve = 0
  nbve = 0
  nbvt = 0

  iec = new((/2,ne/),integer)
  ienode = new((/2,ne/),integer)
  nbe = new((/3,n/),integer)
  xijc = new((/ne/),double)
  yijc = new((/ne/),double)
  dltxc = new((/ne/),double)
  dltyc = new((/ne/),double)
  iec = 0
  ienode = 0
  nbe = 0
  xijc = 0.0d
  yijc = 0.0d
  dltxc = 0.0d
  dltyc = 0.0d
   
  TGE::TRIANGLE_GRID_EDGE(vx,vy,xc,yc,nv,n,nct,m,ne,ntsn,nbsn,ntve,nbve,nbvt,iec,ienode,nbe,xijc,yijc,dltxc,dltyc)

; Momentum equation.
;time 72 (it=71,2008-07-31_12:00)
  do it=71,71

    kbm1 = 18
    advXU = new((/kbm1,dime/),double)
    advXV = new((/kbm1,dime/),double)
    advYU = new((/kbm1,dime/),double)
    advYV = new((/kbm1,dime/),double)
    viscX = new((/kbm1,dime/),double)
    viscY = new((/kbm1,dime/),double)
    dpbpX = new((/kbm1,dime/),double)
    dpbpY = new((/kbm1,dime/),double)
    advXU = 0.0d
    advXV = 0.0d
    advYU = 0.0d
    advYV = 0.0d
    viscX = 0.0d
    viscY = 0.0d
    dpbpX = 0.0d
    dpbpY = 0.0d

    do k=0,kbm1-1

      u = flt2dble(f->u(it,k,:))
      v = flt2dble(f->v(it,k,:))
      el = flt2dble(f->zeta(it,:))
      DT = el + h
      DT1=new((/n/),double)
      do i=0,n-1
        DT1(i)=(DT(nv(0,i)-1)+DT(nv(1,i)-1)+DT(nv(2,i)-1))/3.0d
      end do

    advxu = new((/n/),double)
    advxv = new((/n/),double)
    advyu = new((/n/),double)
    advyv = new((/n/),double)
    viscx = new((/n/),double)
    viscy = new((/n/),double)
    dpbpx = new((/n/),double)
    dpbpy = new((/n/),double)
    advxu = 0.0d
    advxv = 0.0d
    advyu = 0.0d
    advyv = 0.0d
    viscx = 0.0d
    viscy = 0.0d
    dpbpx = 0.0d
    dpbpy = 0.0d
    ADVUV::ADV_UV_EDGE_GCN(u,v,xc,yc,n,m,art,ne,iec,ienode,nbe,DT,DT1,el,xijc,yijc,a1u,a2u, \
                         dltyc,dltxc,advxu,advxv,advyu,advyv,viscx,viscy,dpbpx,dpbpy)
    advXU(k,:) = advxu(eles-1)
    advXV(k,:) = advxv(eles-1)
    advYU(k,:) = advyu(eles-1)
    advYV(k,:) = advyv(eles-1)
    viscX(k,:) = viscx(eles-1)
    viscY(k,:) = viscy(eles-1)
    dpbpX(k,:) = dpbpx(eles-1)
    dpbpY(k,:) = dpbpy(eles-1)
   end do  ;loop k
 end do  ;loop it

 write_table("cross-section_momentum_flood_time0072.dat","w",[/advXU,advXV,advYU,advYV, \
             viscX,viscY,dpbpX,dpbpY/], \
             "%12.6f   %12.6f  %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f")
end
