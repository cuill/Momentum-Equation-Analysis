begin
 kb = 19
 kbm1 = kb-1
 kbm2 = kb-2
 p_sigma = 1.0
; size(z) = 19
; size(zz) = 18 (negative)
; size(dz) = 18
; size(dzz) = 17 (positive)

 z = new((/kb/),double)
 dz = new((/kb/),double)
 zz = new((/kb/),double)
 dzz = new((/kb/),double)
 
 do k=0,kb-1
    z(k) = -(k/(kb-1.0))^p_sigma
 end do

 do k=0,kbm1-1
    dz(k) = z(k)-z(k+1)
    zz(k) = 0.5*(z(k)+z(k+1))
 end do
 zz(kb-1) = 2.0*zz(kb-2)-zz(kb-3)
 ;print(z)
 ;print(dz)
 ;print(zz)

 do k=0,kbm2-1
    dzz(k) = zz(k)-zz(k+1)
 end do
;print(dzz)

 dzz(kbm1-1) = 0.0
 dzz(kb-1) = 0.0
 dz(kb-1)    = 0.0
;WRITE(11,80) K,Z(K),ZZ(K),DZ(K),DZZ(K) 

end
