external TGE "./tge90.so"

begin

  data = asciiread("info_cell_art.dat",-1,"string")
  art  = stringtofloat(str_get_field(data,1," "))

  f = addfile("../netcdf/BaratariaPass_0001.nc","r")
     nv = floattoint(f->nv)
     vx  = flt2dble(f->x)
     vy  = flt2dble(f->y)

     dim=dimsizes(nv)
  ;   printVarSummary(nv)
     n=dim(1)
     m=dimsizes(f->x)
     three=3

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
 
   
; NCL array starting from 0
; FORTRAN array starting from 1
     TGE::TRIANGLE_GRID_EDGE(vx,vy,xc,yc,nv,n,nct,m,ne,ntsn,nbsn,ntve,nbve,nbvt,iec,ienode,nbe,xijc,yijc,dltxc,dltyc)

end
