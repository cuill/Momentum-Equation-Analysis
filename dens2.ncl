;****************
;dens2
;****************
;calculate density using salinity and temperature
;************************************************

undef("dens2")
function dens2(T:numeric, S:numeric)
;T - temperature (Celsius)
;s - salinity (psu) 

begin

  rhof =   S^3*6.76786136e-6 - S*S*4.8249614e-4 +  \
           S*8.14876577e-1 - 0.22584586e0

  rhof = rhof* (T^3*1.667e-8-T*T*8.164e-7+ \
           T*1.803e-5)

  rhof = rhof + 1. - T^3*1.0843e-6 + T*T*9.8185e-5 - T*4.786e-3

  rhof = rhof*(S^3*6.76786136e-6-S*S*4.8249614e-4+S*8.14876577e-1+3.895414e-2)

  rhof = rhof - (T-3.98)^2*(T+283.0)/(503.57*(T+67.26))
 ; rhof1 = (4.2725-rhof)/(4.2725+2.8931)
 ; rhof1 = (9.54-S)/9.54
  
  return (rhof)

end
