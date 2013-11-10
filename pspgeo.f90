! Pavan (pavan_nandan_racherla@alumni.cmu.edu)
! Camp: University of The Aegean; June 23-27 2008
! Source: WRF V3 -> $WRFV3/share/module_llxy.F and $WRFV3/wrf_timeseries.F
! Updated October 13 2008

subroutine pspgeo(i, j, gs, psp_sl, psp_tl, lon_sw, lat_sw, lat, lon)

  use ProjUtils

  implicit none
  type(proj_info) :: proj_struct
  real, intent(in) :: i, j, gs, psp_sl, psp_tl, lon_sw, lat_sw
  real, intent(out) :: lat, lon

  call map_init(proj_struct) ! set up the map transformation structure.

  call map_set(PROJ_PS, &
       &proj_struct, &
       &stdlon = psp_sl, &
       &truelat1 = psp_tl, &
       &lon1 = lon_sw, &
       &lat1 = lat_sw, &
       &knowni = 1., &
       &knownj = 1., &
       &dx = gs)
  
  call ij_to_latlon( proj_struct, i, j, lat, lon )

  return
end subroutine pspgeo
